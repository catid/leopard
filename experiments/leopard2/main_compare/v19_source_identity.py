#!/usr/bin/python3
"""Authenticate held v19 worktrees; leopard-79h.38.5.4.8.2.2.2.1.

Reconstruct trees from actual file bytes, not Git status or index assertions.
The caller owns the enclosing source/preflight/host/lock lifetimes and later
creation, build and runtime gates. No clone, compiler or workload runs here.
"""
from __future__ import annotations

from contextlib import ExitStack
import base64
import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import stat

HERE = Path(__file__).resolve().parent
_dependency = HERE / "v19_streamed_sources.py"
if _dependency.resolve(strict=True) != _dependency:
    raise RuntimeError("source identity dependency is not canonical")
_spec = importlib.util.spec_from_file_location("v19_identity_streamed", _dependency)
streamed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(streamed)
preflight, owners = streamed.preflight, streamed.owners
host, provenance, require = streamed.host, streamed.provenance, streamed.require

SUBMODULE_COMMIT = "cad518a93b326f0f644b7972d488d04eaa2b0475"
REPOSITORIES = ("candidate-source", "candidate-source/sse2neon", "leopard1-source")
MAX_METADATA_BYTES = 1 << 20
MAX_PATH_TOTAL_BYTES = 1 << 20
HEX40 = re.compile(r"[0-9a-f]{40}")


def _metadata_path(path: str) -> bool:
    return any(path == role + "/.git" or path.startswith(role + "/.git/") for role in REPOSITORIES)


def object_id(kind: str, data: bytes) -> str:
    return hashlib.sha1(f"{kind} {len(data)}\0".encode() + data,
                        usedforsecurity=False).hexdigest()


def tree_id(entries: list[dict]) -> str:
    """Serialize one bounded ordinary-file/gitlink tree with Git byte ordering."""
    require(type(entries) is list and 0 < len(entries) <= streamed.MAX_FILES,
            "source tree inventory is empty or oversized")
    root = {}
    total = 0
    for entry in entries:
        require(type(entry) is dict and set(entry) == {"path", "mode", "object_id"},
                "source tree entry has unexpected fields")
        path, mode, oid = entry["path"], entry["mode"], entry["object_id"]
        require(type(path) is str and type(mode) is str and type(oid) is str and
                mode in ("100644", "100755", "160000") and HEX40.fullmatch(oid),
                "source tree entry has invalid type, mode or object ID")
        provenance._require_safe_unicode(path, "source tree path")
        components = path.split("/")
        total += len(path.encode("utf-8"))
        require(total <= MAX_PATH_TOTAL_BYTES and len(components) <= streamed.MAX_DEPTH and
                all(part not in ("", ".", "..", ".git") and "\\" not in part and
                    len(part.encode("utf-8")) <= 255 for part in components),
                "source tree path is unsafe or oversized")
        current = root
        for part in components[:-1]:
            current = current.setdefault(part, {})
            require(type(current) is dict, "source tree file/directory paths overlap")
        require(components[-1] not in current, "source tree path is duplicated or overlaps")
        current[components[-1]] = mode, oid

    def digest(directory):
        body = bytearray()
        for name, value in sorted(directory.items(), key=lambda row:
                row[0].encode("utf-8") + (b"/" if type(row[1]) is dict else b"")):
            mode, oid = ("40000", digest(value)) if type(value) is dict else value
            body.extend(mode.encode("ascii") + b" " + name.encode("utf-8") + b"\0" + bytes.fromhex(oid))
        require(len(body) <= MAX_METADATA_BYTES, "source tree object exceeds byte bound")
        return object_id("tree", body)
    return digest(root)


def clone_config(role: str, canonical_root: str) -> bytes:
    """Exact config dialect of the pinned native fresh-clone procedure."""
    require(role in REPOSITORIES, "unexpected source repository role")
    host.canonical_path(canonical_root, nonroot=True)
    config = ("[core]\n\trepositoryformatversion = 0\n\tfilemode = true\n"
              "\tbare = false\n\tlogallrefupdates = true\n[remote \"origin\"]\n"
              "\turl = " + canonical_root + "/" + role + "\n")
    if role == "leopard1-source":
        config += ("\tfetch = +refs/heads/master:refs/remotes/origin/master\n"
                   "[branch \"master\"]\n\tremote = origin\n\tmerge = refs/heads/master\n")
    return config.encode("utf-8")


class PinnedSourceIdentity:
    """Retain the inspected Git tool and borrow live source/preflight owners.

    Initial bytes must reconstruct the pinned commit trees; the candidate also
    matches the pinned SHA-256 manifest including permission modes. Git queries
    use the sealed executable and retained root/gitdir fds. This is a narrow
    ordinary detached-clone dialect, not general repository compatibility or
    proof of independent creation, trusted runtime libraries, or atomicity.
    """

    def __init__(self, source_owner, retained_preflight):
        require(type(source_owner) is streamed.StreamingSourceOwner and
                type(retained_preflight) is preflight.PinnedPreflight,
                "source authentication requires live owners, not records")
        self.source = source_owner
        self.retained = retained_preflight
        self._stack = ExitStack()
        self._state = "new"
        self._git = None
        self._record = None

    def _read(self, relative):
        require(relative in self.source._files, "required source metadata is missing: " + relative)
        entry = self.source._files[relative]
        self.source._metadata(entry)
        require(entry[3][6] <= MAX_METADATA_BYTES, "source metadata exceeds byte bound")
        data = os.pread(entry[0], entry[3][6] + 1, 0)
        require(len(data) == entry[3][6] and hashlib.sha256(data).hexdigest() == entry[4],
                "source metadata bytes changed")
        self.source._metadata(entry)
        return data

    def _validate_inputs(self, *, evict_cache=False):
        self.source.validate_current(evict_cache=evict_cache)
        self.retained.validate_current()
        if self._git is not None:
            owners.verify_current_bytes(self._git, self._git.identity["sha256"])

    def _metadata_policy(self, role, commit, canonical_root):
        prefix = role + "/.git/"
        require(self._read(prefix + "HEAD") == (commit + "\n").encode("ascii"),
                "source HEAD is not the pinned detached commit")
        require(self._read(prefix + "config") == clone_config(role, canonical_root),
                "source Git config differs from the exact clone dialect")
        forbidden = ("commondir", "gitdir", "config.worktree", "info/grafts",
                     "objects/info/alternates", "objects/info/http-alternates", "refs/replace")
        for path in (*self.source._files, *self.source._directories):
            if not path.startswith(prefix):
                continue
            name = path[len(prefix):]
            require(not any(name == value or name.startswith(value + "/") for value in forbidden),
                    "source Git metadata redirects or replaces object identity")
            if path in self.source._files and name.startswith("hooks/"):
                require(name.endswith(".sample"), "source Git metadata contains an active hook")
        if prefix + "shallow" in self.source._files:
            require(self._read(prefix + "shallow") == (commit + "\n").encode("ascii"),
                    "source shallow boundary differs from the pinned commit")
        if prefix + "packed-refs" in self.source._files:
            require(b" refs/replace/" not in self._read(prefix + "packed-refs"),
                    "source packed refs contain replacement objects")

    def _query(self, role, arguments):
        self._validate_inputs(evict_cache=True)
        root_fd = self.source._directories[role][0]
        git_fd = self.source._directories[role + "/.git"][0]
        command = ("/usr/bin/git", "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
                   "-c", "core.hooksPath=/dev/null", "-c", "core.excludesFile=/dev/null",
                   f"--git-dir=/proc/self/fd/{git_fd}", f"--work-tree=/proc/self/fd/{root_fd}",
                   *arguments)
        data = provenance._run(command, "v19 source " + role + " " + arguments[0],
            maximum_bytes=MAX_METADATA_BYTES, timeout=30,
            inherited_descriptors=(root_fd, git_fd), executable_descriptor=self._git.executable_descriptor)
        self._validate_inputs(evict_cache=True)
        return data

    def _worktree(self, role):
        entries, files = [], []
        for path, entry in sorted(self.source._files.items()):
            if not path.startswith(role + "/"):
                continue
            relative = path[len(role) + 1:]
            if _metadata_path(path):
                continue
            files.append({"path": relative, "size": entry[3][6],
                          "mode": stat.S_IMODE(entry[3][2]), "sha256": entry[4]})
            if role == "candidate-source" and relative.startswith("sse2neon/"):
                continue
            self.source._metadata(entry)
            size, offset = entry[3][6], 0
            digest = hashlib.sha1(f"blob {size}\0".encode("ascii"), usedforsecurity=False)
            actual = hashlib.sha256()
            while offset < size:
                part = os.pread(entry[0], min(streamed.HASH_BLOCK_BYTES, size - offset), offset)
                require(part, "source blob ended during authentication")
                digest.update(part)
                actual.update(part)
                offset += len(part)
            require(not os.pread(entry[0], 1, size) and actual.hexdigest() == entry[4],
                    "source blob bytes changed during authentication")
            self.source._metadata(entry)
            entries.append({"path": relative, "mode": "100755" if entry[3][2] & stat.S_IXUSR else "100644",
                            "object_id": digest.hexdigest()})
        if role in ("candidate-source", "leopard1-source"):
            entries.append({"path": "sse2neon", "mode": "160000", "object_id": SUBMODULE_COMMIT})
        entries.sort(key=lambda entry: entry["path"].encode("utf-8"))
        return entries, files

    def __enter__(self):
        require(self._state == "new", "source identity owner cannot be reused")
        self._state = "entering"
        try:
            self._validate_inputs()
            historical = self.retained.record()
            pins = historical["source_pins"]
            require(set(pins) == {"source_commit", "source_tree", "baseline_commit", "baseline_tree"} and
                    all(type(value) is str and HEX40.fullmatch(value) for value in pins.values()),
                    "retained source pins are malformed")
            manifest = preflight._json(self.retained._bytes("candidate-build-provenance.json"))["tracked_source_manifest"]
            git_roots = {path for path in self.source._directories if path.split("/")[-1] == ".git"}
            require(git_roots == {role + "/.git" for role in REPOSITORIES},
                    "source Git root inventory differs")
            # The exact-main native build leaves its committed sse2neon gitlink
            # uninitialized. Bind precisely that empty directory and its object
            # ID in the reconstructed baseline tree; do not allow other empties.
            baseline_gitlink = "leopard1-source/sse2neon"
            require(baseline_gitlink in self.source._directories and
                    self.source._directories[baseline_gitlink][4] == set(),
                    "baseline gitlink is not the exact uninitialized directory")
            expected_directories = {*REPOSITORIES, baseline_gitlink}
            for path in self.source._files:
                if not _metadata_path(path):
                    components = path.split("/")
                    expected_directories.update("/".join(components[:index]) for index in range(1, len(components)))
            require({path for path in self.source._directories if not _metadata_path(path)} == expected_directories,
                    "source directory topology contains untracked empty directories")
            commits = {"candidate-source": pins["source_commit"], "leopard1-source": pins["baseline_commit"],
                       "candidate-source/sse2neon": SUBMODULE_COMMIT}
            for role in REPOSITORIES:
                self._metadata_policy(role, commits[role], historical["canonical_run_root"])
            self._git = self._stack.enter_context(provenance._RetainedFileSnapshot(
                Path("/usr/bin/git"), "v19 authentication Git", maximum_bytes=owners.MAX_ARTIFACT_BYTES))
            self._git.executable_descriptor  # Seal before the first query.
            repositories = {}
            for role in REPOSITORIES:
                entries, files = self._worktree(role)
                actual_tree = tree_id(entries)
                if role == "candidate-source":
                    require(actual_tree == pins["source_tree"], "candidate worktree differs from pinned tree")
                    require(manifest.get("schema") == "leopard2-tracked-source-tree/v1" and
                            files == manifest.get("files") and
                            sum(row["size"] for row in files) == manifest.get("total_bytes"),
                            "candidate worktree differs from pinned SHA-256 manifest")
                elif role == "leopard1-source":
                    require(actual_tree == pins["baseline_tree"], "baseline worktree differs from pinned tree")
                raw_commit = self._query(role, ("cat-file", "commit", commits[role]))
                require(object_id("commit", raw_commit) == commits[role] and b"\n\n" in raw_commit,
                        "source raw commit does not authenticate its pinned object ID")
                headers = raw_commit.split(b"\n\n", 1)[0].split(b"\n")
                tree_header = b"tree " + actual_tree.encode("ascii")
                require(headers[0] == tree_header and
                        [line for line in headers if line.startswith(b"tree ")] == [tree_header],
                        "source raw commit names a different content tree")
                index = b"".join((row["mode"] + " " + row["object_id"] + " 0\t" + row["path"]).encode("utf-8") +
                                 b"\0" for row in entries)
                flags = b"".join(b"H " + row["path"].encode("utf-8") + b"\0" for row in entries)
                require(len(index) <= MAX_METADATA_BYTES and len(flags) <= MAX_METADATA_BYTES,
                        "source index transcript exceeds byte bound")
                require(self._query(role, ("ls-files", "--stage", "-z")) == index,
                        "source index differs from authenticated worktree")
                for option in ("-v", "-f"):
                    require(self._query(role, ("ls-files", option, "-z")) == flags,
                            "source index contains non-default flags")
                repositories[role] = {"commit": commits[role], "tree": actual_tree,
                    "raw_commit_sha256": hashlib.sha256(raw_commit).hexdigest(),
                    "raw_commit_base64": base64.b64encode(raw_commit).decode("ascii"),
                    "index_entries_sha256": hashlib.sha256(index).hexdigest(),
                    "entries": entries, "entry_count": len(entries),
                    "uninitialized_gitlinks": {"sse2neon": SUBMODULE_COMMIT} if role == "leopard1-source" else {}}
            self._record = {"schema": "leopard2-v19-pinned-source-identity/v1",
                "workspace": str(self.source.workspace), "repositories": repositories,
                "candidate_manifest_sha256": hashlib.sha256(host.canonical_bytes(manifest["files"])).hexdigest(),
                "git_executable": copy.deepcopy(self._git.executable_identity),
                "source_content_verified": True, "commit_and_index_verified": True,
                "source_creation_verified": False, "runtime_closure_verified": False,
                "atomic_snapshot": False, "live_acquisition_armed": False}
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "closed"
            self._stack.close()
            raise

    def validate_current(self, *, evict_cache=False):
        require(self._state == "held", "source identity owner is not usable")
        try:
            self._validate_inputs(evict_cache=evict_cache)
        except BaseException:
            self._state = "failed"
            raise

    def record(self, *, evict_cache=False):
        self.validate_current(evict_cache=evict_cache)
        return copy.deepcopy(self._record)

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held":
                self.validate_current()
            elif value is None:
                raise host.PreflightError("failed source identity owner cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
