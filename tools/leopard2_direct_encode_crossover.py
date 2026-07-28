#!/usr/bin/env python3
"""Reproduce Leopard2 tiny direct-encoder crossover measurements.

Screening and generic pinned modes consume already-built
``bench_leopard2_direct_encode`` binaries.  ``historical-avx2`` instead creates
or strictly resumes a retained fresh Release/explicit-AVX2 build while holding
the canonical authoritative lock, freezes that executable, and only then begins
isolated ABBA measurement.  Every cell is a resumable JSON job with
deterministic input data and retained logs.

Typical use::

    tools/leopard2_direct_encode_crossover.py screen --build-root build
    tools/leopard2_direct_encode_crossover.py pinned --build-root build \
        --cpu 16 --sibling 80
    tools/leopard2_direct_encode_crossover.py analyze \
        --result-dir results/leopard2/direct-encode-crossover/pinned

The default build-root lookup accepts the repository's conventional trees
``build/direct-encode-SCALAR`` and ``build/SCALAR``.  ``--executable-root`` may
instead name a different root or contain a literal ``{backend}`` placeholder.
The forced-path benchmark is intentionally bounded to ``K,R <= 16``; the
adjacent 17 fallback is a correctness/dispatch test, not a meaningful
forced-direct crossover cell.
"""

from __future__ import print_function

import argparse
import concurrent.futures
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import shutil
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


SCHEMA = "leopard2-direct-encode-crossover/v2"
JOB_SCHEMA = "leopard2-direct-encode-crossover-job/v2"
ANALYSIS_SCHEMA = "leopard2-direct-encode-crossover-analysis/v2"
BENCHMARK_SCHEMA = "leopard2-direct-encode-benchmark-v2"
KNOWN_BACKENDS = ("scalar", "ssse3", "avx2", "avx512")
FROZEN_EXECUTABLE_SCHEMA = "leopard2-frozen-executable/v1"
AUTHORITATIVE_COMMANDS = ("pinned", "historical-avx2")
AUTHORITATIVE_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
CANONICAL_GIT = Path("/usr/bin/git")
CONTROLLED_BUILD_SCHEMA = "leopard2-direct-controlled-build/v2"
BENCHMARK_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
GIT_ENVIRONMENT = dict(BENCHMARK_ENVIRONMENT)
GIT_ENVIRONMENT.update({
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_SYSTEM": "/dev/null",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
})


class CrossoverError(Exception):
    """An actionable configuration, execution, or result error."""


class AuthoritativeGuardError(CrossoverError):
    """A held canonical or physical-pair guard was lost or replaced."""


def load_isolation_support(source):
    path = (
        Path(source) / "experiments" / "leopard2" /
        "main_compare" / "run_abba.py"
    ).resolve()
    if not path.is_file():
        raise CrossoverError(
            "canonical isolation support is missing: {}".format(path)
        )
    name = "leopard2_direct_crossover_isolation_support"
    specification = importlib.util.spec_from_file_location(name, str(path))
    if specification is None or specification.loader is None:
        raise CrossoverError(
            "cannot load canonical isolation support: {}".format(path)
        )
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    try:
        specification.loader.exec_module(module)
    except Exception as error:
        raise CrossoverError(
            "cannot initialize canonical isolation support: {}".format(error)
        )
    for function in (
            "cpu_stat_snapshot", "isolation_record", "PairLease",
            "validate_isolation", "validate_pair_lease_identity",
            "validate_topology"):
        if not callable(getattr(module, function, None)):
            raise CrossoverError(
                "canonical isolation support omits {}".format(function)
            )
    if getattr(module, "MAX_SIBLING_NONIDLE_JIFFIES", None) != 0:
        raise CrossoverError(
            "canonical isolation support no longer requires an idle sibling"
        )
    return module


class CanonicalAuthoritativeLock:
    def __init__(self, path=AUTHORITATIVE_LOCK):
        self.path = Path(path)
        self.descriptor = None
        self.identity = None

    def validate_current(self):
        if self.descriptor is None or self.identity is None:
            raise AuthoritativeGuardError(
                "canonical authoritative lock is not held"
            )
        try:
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(str(self.path))
        except OSError as error:
            raise AuthoritativeGuardError(
                "canonical authoritative lock revalidation failed: {}".format(
                    error
                )
            )
        current = {
            "device": descriptor.st_dev,
            "inode": descriptor.st_ino,
            "lock": "exclusive",
            "mode": stat.S_IMODE(descriptor.st_mode),
            "path": str(self.path),
            "uid": descriptor.st_uid,
        }
        if (not stat.S_ISREG(descriptor.st_mode) or
                descriptor.st_nlink != 1 or
                (descriptor.st_dev, descriptor.st_ino) !=
                (path.st_dev, path.st_ino) or current != self.identity):
            raise AuthoritativeGuardError(
                "canonical authoritative lock path was replaced or changed"
            )
        return self.identity

    def __enter__(self):
        if self.path != AUTHORITATIVE_LOCK:
            raise CrossoverError(
                "authoritative runs require canonical lock {}".format(
                    AUTHORITATIVE_LOCK
                )
            )
        flags = os.O_RDWR | os.O_CREAT
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            self.descriptor = os.open(str(self.path), flags, 0o600)
            metadata = os.fstat(self.descriptor)
            path_metadata = os.lstat(str(self.path))
            if (not stat.S_ISREG(metadata.st_mode) or
                    metadata.st_uid != os.getuid() or metadata.st_nlink != 1 or
                    stat.S_IMODE(metadata.st_mode) != 0o600 or
                    (metadata.st_dev, metadata.st_ino) !=
                    (path_metadata.st_dev, path_metadata.st_ino)):
                raise CrossoverError(
                    "canonical authoritative lock has unsafe metadata"
                )
            fcntl.flock(self.descriptor, fcntl.LOCK_EX)
            locked = os.fstat(self.descriptor)
            self.identity = {
                "device": locked.st_dev,
                "inode": locked.st_ino,
                "lock": "exclusive",
                "mode": stat.S_IMODE(locked.st_mode),
                "path": str(self.path),
                "uid": locked.st_uid,
            }
            self.validate_current()
            return self
        except Exception:
            if self.descriptor is not None:
                try:
                    fcntl.flock(self.descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(self.descriptor)
                    self.descriptor = None
                    self.identity = None
            raise

    def __exit__(self, exc_type, exc, traceback):
        validation_error = None
        try:
            self.validate_current()
        except Exception as error:
            validation_error = error
        descriptor = self.descriptor
        self.descriptor = None
        self.identity = None
        try:
            if descriptor is not None:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            if descriptor is not None:
                os.close(descriptor)
        if validation_error is not None:
            raise validation_error


def canonical_authoritative_lock(path=AUTHORITATIVE_LOCK):
    return CanonicalAuthoritativeLock(path)


class AuthoritativePairGuard:
    """Add fail-closed exit validation to the shared pair-lease helper."""

    def __init__(self, guard):
        self.guard = guard
        self.identity = None

    def validate_current(self):
        try:
            self.guard.validate_current()
        except Exception as error:
            raise AuthoritativeGuardError(
                "physical-pair guard revalidation failed: {}".format(error)
            )
        return self.identity

    def __enter__(self):
        self.identity = self.guard.__enter__()
        try:
            self.validate_current()
        except Exception:
            try:
                self.guard.__exit__(None, None, None)
            finally:
                self.identity = None
            raise
        return self.identity

    def __exit__(self, exc_type, exc, traceback):
        validation_error = None
        try:
            self.validate_current()
        except Exception as error:
            validation_error = error
        try:
            result = self.guard.__exit__(exc_type, exc, traceback)
        finally:
            self.identity = None
        if validation_error is not None:
            raise validation_error
        return result


def validate_authoritative_guards(lock_guard, pair_guard):
    try:
        lock_guard.validate_current()
        pair_guard.validate_current()
    except AuthoritativeGuardError:
        raise
    except Exception as error:
        raise AuthoritativeGuardError(
            "authoritative guard revalidation failed: {}".format(error)
        )


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False
    ).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical_bytes(value))


def evidence_contract(frozen_executable_required):
    return {
        "benchmark_schema": BENCHMARK_SCHEMA,
        "execution_mad_path": (
            "metrics.encode_execution.mad_us_per_batch_call"
        ),
        "execution_median_path": (
            "metrics.encode_execution.median_us_per_batch_call"
        ),
        "frozen_executable_required": bool(frozen_executable_required),
        "numeric_policy": "finite JSON numbers; median > 0 and MAD >= 0",
        "parity_identity_fields": [
            "algorithm", "digest", "hashed_bytes",
            "requested_parity_indices",
        ],
        "parity_oracles": [
            "selected_transform_reference_parity_match",
            "direct_transform_parity_match on force-direct invocation",
            "same non-vacuous parity identity across direct/transform ABBA",
            "unrequested_outputs_untouched",
        ],
    }


def decode_json_bytes(value, description):
    def reject_duplicate_keys(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise CrossoverError(
                    "{} contains duplicate JSON key {!r}".format(description, key)
                )
            result[key] = item
        return result

    def reject_nonstandard_constant(constant):
        raise CrossoverError(
            "{} contains non-standard JSON number {}".format(
                description, constant
            )
        )

    def parse_finite_float(text):
        try:
            number = float(text)
        except (OverflowError, ValueError):
            number = float("nan")
        if not math.isfinite(number):
            raise CrossoverError(
                "{} contains a non-finite JSON number {}".format(
                    description, text
                )
            )
        return number

    def parse_bounded_integer(text):
        try:
            number = int(text)
        except ValueError:
            raise CrossoverError(
                "{} contains an invalid JSON integer".format(description)
            )
        if abs(number) > (1 << 64) - 1:
            raise CrossoverError(
                "{} contains a JSON integer outside the 64-bit evidence "
                "domain".format(description)
            )
        return number

    try:
        text = value.decode("utf-8")
    except (AttributeError, UnicodeError) as error:
        raise CrossoverError("cannot decode {}: {}".format(description, error))
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonstandard_constant,
            parse_float=parse_finite_float,
            parse_int=parse_bounded_integer,
        )
    except CrossoverError:
        raise
    except ValueError as error:
        raise CrossoverError("cannot parse {}: {}".format(description, error))


def normalized_output(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(
                value, output, indent=2, sort_keys=True, ensure_ascii=True,
                allow_nan=False
            )
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, str(path))
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def compact_cpu_list(cpus):
    values = sorted(set(int(cpu) for cpu in cpus))
    if not values:
        return ""
    result = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        result.append(str(first) if first == previous else "{}-{}".format(first, previous))
        first = previous = value
    result.append(str(first) if first == previous else "{}-{}".format(first, previous))
    return ",".join(result)


def allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            result = sorted(os.sched_getaffinity(0))
            if result:
                return result
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def read_optional(path):
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def cpu_description(cpus):
    model = None
    flags = []
    cpuinfo = read_optional("/proc/cpuinfo")
    if cpuinfo:
        wanted = set(cpus)
        fallback = None
        for block in cpuinfo.split("\n\n"):
            fields = {}
            for line in block.splitlines():
                if ":" in line:
                    key, value = line.split(":", 1)
                    fields[key.strip().lower()] = value.strip()
            if fallback is None:
                fallback = fields
            try:
                processor = int(fields.get("processor", "-1"))
            except ValueError:
                processor = -1
            if processor in wanted:
                model = fields.get("model name", fields.get("hardware"))
                flags = fields.get("flags", fields.get("features", "")).split()
                break
        if model is None and fallback:
            model = fallback.get("model name", fallback.get("hardware"))
            flags = fallback.get("flags", fallback.get("features", "")).split()
    return model, sorted(set(flags))


def cpu_topology(cpu):
    root = Path("/sys/devices/system/cpu/cpu{}".format(cpu))
    result = {"cpu": cpu}
    fields = {
        "core_id": root / "topology/core_id",
        "package_id": root / "topology/physical_package_id",
        "thread_siblings": root / "topology/thread_siblings_list",
        "scaling_governor": root / "cpufreq/scaling_governor",
    }
    for key, path in fields.items():
        value = read_optional(path)
        if value is not None:
            result[key] = value
    return result


def machine_identity(cpus):
    model, flags = cpu_description(cpus)
    uname = platform.uname()
    result = {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_flags": flags,
        "cpu_model": model,
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "uname": {
            "machine": uname.machine,
            "node": uname.node,
            "release": uname.release,
            "system": uname.system,
            "version": uname.version,
        },
    }
    memory = read_optional("/proc/meminfo")
    if memory:
        for line in memory.splitlines():
            if line.startswith("MemTotal:"):
                result["memory_total"] = line.split(":", 1)[1].strip()
                break
    no_turbo = read_optional("/sys/devices/system/cpu/intel_pstate/no_turbo")
    boost = read_optional("/sys/devices/system/cpu/cpufreq/boost")
    if no_turbo is not None:
        result["intel_pstate_no_turbo"] = no_turbo
    if boost is not None:
        result["cpufreq_boost"] = boost
    return result


def git_command_bytes(source, arguments, description):
    try:
        completed = subprocess.run(
            [str(CANONICAL_GIT)] + list(arguments), cwd=str(source),
            env=GIT_ENVIRONMENT,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=30, check=False
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise CrossoverError("cannot {}: {}".format(description, error))
    if completed.returncode != 0:
        raise CrossoverError(
            "cannot {}: {}".format(
                description,
                normalized_output(completed.stderr).decode(
                    "utf-8", errors="replace"
                ).strip()
            )
        )
    return completed.stdout


def canonical_git_identity():
    try:
        before = CANONICAL_GIT.stat()
        value = CANONICAL_GIT.read_bytes()
        after = CANONICAL_GIT.stat()
    except OSError as error:
        raise CrossoverError(
            "canonical Git executable is unavailable: {}".format(error)
        )
    if (CANONICAL_GIT.resolve() != CANONICAL_GIT or
            not stat.S_ISREG(before.st_mode) or
            not os.access(str(CANONICAL_GIT), os.X_OK) or
            (before.st_dev, before.st_ino, before.st_size,
             before.st_mtime_ns, before.st_ctime_ns) !=
            (after.st_dev, after.st_ino, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns)):
        raise CrossoverError(
            "canonical Git executable is not a stable regular executable"
        )
    return {
        "environment": dict(GIT_ENVIRONMENT),
        "path": str(CANONICAL_GIT),
        "sha256": digest_bytes(value),
    }


def parse_stage_zero_entries(value, description):
    entries = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            header, encoded_path = encoded.split(b"\t", 1)
            mode, object_id, stage = header.decode("ascii").split(" ")
            relative = encoded_path.decode("utf-8")
        except (UnicodeError, ValueError) as error:
            raise CrossoverError(
                "{} contains an invalid entry: {}".format(description, error)
            )
        if (stage != "0" or not re.fullmatch(r"[0-9a-f]{40}", object_id) or
                mode not in ("100644", "100755", "120000", "160000") or
                not relative or relative in entries):
            raise CrossoverError(
                "{} contains a conflicted, unsupported, or duplicate entry "
                "{!r}".format(description, relative)
            )
        entries[relative] = (mode, object_id)
    return entries


def parse_head_tree_entries(value, description):
    entries = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            header, encoded_path = encoded.split(b"\t", 1)
            mode, object_type, object_id = header.decode("ascii").split(" ")
            relative = encoded_path.decode("utf-8")
        except (UnicodeError, ValueError) as error:
            raise CrossoverError(
                "{} contains an invalid entry: {}".format(description, error)
            )
        if (object_type not in ("blob", "commit") or
                not re.fullmatch(r"[0-9a-f]{40}", object_id) or
                mode not in ("100644", "100755", "120000", "160000") or
                not relative or relative in entries or
                (mode == "160000") != (object_type == "commit")):
            raise CrossoverError(
                "{} contains an unsupported or duplicate entry {!r}".format(
                    description, relative
                )
            )
        entries[relative] = (mode, object_id)
    return entries


def parse_index_tags(value, description):
    tags = {}
    for encoded in value.split(b"\0"):
        if not encoded:
            continue
        try:
            tag = encoded[:1].decode("ascii")
            separator = encoded[1:2]
            relative = encoded[2:].decode("utf-8")
        except UnicodeError as error:
            raise CrossoverError(
                "{} contains an invalid tag: {}".format(description, error)
            )
        if (separator != b" " or not relative or relative in tags or
                tag != "H"):
            raise CrossoverError(
                "{} found non-default index flag tag {!r} for {!r}".format(
                    description, tag, relative
                )
            )
        tags[relative] = tag
    return tags


def repository_source_controls(source, stage_entries):
    replacements = normalized_output(git_command_bytes(
        source, (
            "for-each-ref", "--format=%(refname)", "refs/replace/",
        ), "enumerate Git replacement refs"
    )).decode("utf-8", errors="strict").rstrip("\n")
    replacement_refs = replacements.splitlines() if replacements else []
    if replacement_refs:
        raise CrossoverError(
            "Git replacement refs are forbidden: {}".format(
                ", ".join(replacement_refs)
            )
        )

    graft_name = normalized_output(git_command_bytes(
        source, ("rev-parse", "--git-path", "info/grafts"),
        "locate legacy Git grafts"
    )).decode("utf-8", errors="strict").strip()
    graft_path = Path(graft_name)
    if not graft_path.is_absolute():
        graft_path = (source / graft_path).resolve()
    if graft_path.exists():
        raise CrossoverError(
            "legacy Git grafts are forbidden: {}".format(graft_path)
        )

    visible_tags = parse_index_tags(git_command_bytes(
        source, ("ls-files", "-v", "-z"),
        "enumerate assume-unchanged/skip-worktree index flags"
    ), "git ls-files -v")
    fsmonitor_tags = parse_index_tags(git_command_bytes(
        source, ("ls-files", "-f", "-z"),
        "enumerate fsmonitor-valid index flags"
    ), "git ls-files -f")
    if (set(visible_tags) != set(stage_entries) or
            set(fsmonitor_tags) != set(stage_entries)):
        raise CrossoverError(
            "Git index flag enumeration differs from the stage-0 closure"
        )

    head_entries = parse_head_tree_entries(git_command_bytes(
        source, ("ls-tree", "-r", "-z", "HEAD"),
        "enumerate the exact HEAD tree"
    ), "git ls-tree HEAD")
    return {
        "head_index_match": head_entries == stage_entries,
        "head_tree_sha256": digest_value({
            name: list(head_entries[name]) for name in sorted(head_entries)
        }),
        "index_entry_count": len(stage_entries),
        "index_flags": "all ls-files -v/-f tags are canonical H",
        "index_sha256": digest_value({
            name: list(stage_entries[name]) for name in sorted(stage_entries)
        }),
        "legacy_grafts_absent": True,
        "replace_refs": [],
    }


def git_blob_object_id(value):
    header = b"blob " + str(len(value)).encode("ascii") + b"\0"
    return hashlib.sha1(header + value).hexdigest()


def tracked_git_closure(source):
    source = Path(source).resolve()
    output = git_command_bytes(
        source, ("ls-files", "-s", "-z"), "enumerate tracked source closure"
    )
    stage_entries = parse_stage_zero_entries(
        output, "Git stage-0 source closure"
    )
    controls = repository_source_controls(source, stage_entries)
    files = {}
    for relative, (mode, object_id) in stage_entries.items():
        path = source / relative
        try:
            if mode == "160000":
                resolved = path.resolve(strict=True)
                common = os.path.commonpath((str(source), str(resolved)))
                if common != str(source) or not resolved.is_dir():
                    raise CrossoverError(
                        "tracked gitlink escapes or is not initialized: {}".format(
                            relative
                        )
                    )
                head = normalized_output(git_command_bytes(
                    resolved, ("rev-parse", "--verify", "HEAD"),
                    "read submodule HEAD {}".format(relative)
                )).decode("ascii").strip()
                status_text = normalized_output(git_command_bytes(
                    resolved, (
                        "status", "--porcelain=v1", "--untracked-files=all",
                        "--ignore-submodules=none",
                    ), "read submodule status {}".format(relative)
                )).decode("utf-8").rstrip("\n")
                if head != object_id:
                    raise CrossoverError(
                        "initialized submodule {} HEAD {} differs from index "
                        "{}".format(relative, head, object_id)
                    )
                closure = tracked_git_closure(resolved)
                record = {
                    "index_mode": mode,
                    "index_object": object_id,
                    "submodule": {
                        "digest": closure["digest"],
                        "files": closure["files"],
                        "head": head,
                        "repository_controls":
                            closure["repository_controls"],
                        "status": (
                            status_text.splitlines() if status_text else []
                        ),
                        "worktree_clean": not bool(status_text),
                    },
                }
            elif mode in ("100644", "100755", "120000"):
                resolved = path.resolve(strict=True)
                common = os.path.commonpath((str(source), str(resolved)))
                if common != str(source):
                    raise CrossoverError(
                        "tracked path escapes the source root: {}".format(relative)
                    )
                if mode == "120000":
                    if not path.is_symlink():
                        raise CrossoverError(
                            "tracked symlink changed type: {}".format(relative)
                        )
                    before = path.lstat()
                    if not stat.S_ISLNK(before.st_mode):
                        raise CrossoverError(
                            "tracked symlink changed type while hashing: "
                            "{}".format(relative)
                        )
                    value = os.readlink(str(path)).encode("utf-8")
                    after = path.lstat()
                    worktree_mode = "120000"
                else:
                    if not path.is_file() or path.is_symlink():
                        raise CrossoverError(
                            "tracked file changed type: {}".format(relative)
                        )
                    before = path.lstat()
                    if not stat.S_ISREG(before.st_mode):
                        raise CrossoverError(
                            "tracked file changed type while hashing: "
                            "{}".format(relative)
                        )
                    value = path.read_bytes()
                    after = path.lstat()
                    worktree_mode = (
                        "100755" if before.st_mode & 0o111 else "100644"
                    )
                if (
                    before.st_dev, before.st_ino, before.st_mode,
                    before.st_size, before.st_mtime_ns, before.st_ctime_ns
                ) != (
                    after.st_dev, after.st_ino, after.st_mode,
                    after.st_size, after.st_mtime_ns, after.st_ctime_ns
                ):
                    raise CrossoverError(
                        "tracked source changed while hashing: {}".format(
                            relative
                        )
                    )
                record = {
                    "index_mode": mode,
                    "index_object": object_id,
                    "worktree_git_object": git_blob_object_id(value),
                    "worktree_mode": worktree_mode,
                    "worktree_sha256": digest_bytes(value),
                }
            else:
                raise CrossoverError(
                    "unsupported Git index mode {} for {}".format(mode, relative)
                )
        except CrossoverError:
            raise
        except (OSError, UnicodeError) as error:
            raise CrossoverError(
                "cannot hash tracked source {}: {}".format(relative, error)
            )
        files[relative] = record
    if not files:
        raise CrossoverError("Git tracked source closure is empty")
    ordered = {name: files[name] for name in sorted(files)}
    material = {
        "files": ordered,
        "repository_controls": controls,
    }
    return {
        "digest": digest_value(material),
        "files": ordered,
        "repository_controls": controls,
    }


def source_fingerprint(source):
    return tracked_git_closure(Path(source).resolve())


def git_source_identity(source, require_clean):
    def git(arguments):
        return normalized_output(git_command_bytes(
            source, arguments, "inspect Git source identity"
        )).decode(
            "utf-8", errors="strict"
        ).strip()

    head = git(("rev-parse", "--verify", "HEAD"))
    tree = git(("rev-parse", "--verify", "HEAD^{tree}"))
    branch = git(("branch", "--show-current"))
    status = normalized_output(git_command_bytes(
        source, (
            "status", "--porcelain=v1", "--untracked-files=all",
            "--ignore-submodules=none",
        ), "read Git worktree status"
    )).decode("utf-8", errors="strict").rstrip("\n")
    for label, value in (("HEAD", head), ("tree", tree)):
        if not re.fullmatch(r"[0-9a-f]{40}", value):
            raise CrossoverError("Git {} is not a full lowercase SHA-1".format(label))
    if require_clean and status:
        raise CrossoverError(
            "authoritative source tree must be clean in Git: {}".format(
                status.replace("\n", "; ")
            )
        )
    return {
        "branch": branch or None,
        "head": head,
        "status": status.splitlines() if status else [],
        "worktree_clean": not bool(status),
        "tree": tree,
    }


def source_identity(source, require_clean):
    result = source_fingerprint(source)
    result["git"] = git_source_identity(source, require_clean)
    result["git_tool"] = canonical_git_identity()
    return validate_source_state(
        result, "captured source identity", require_clean
    )


def validate_source_state(value, description, require_clean):
    if (not isinstance(value, dict) or
            set(value) != {
                "digest", "files", "git", "git_tool",
                "repository_controls"}):
        raise CrossoverError("{} has an invalid source identity".format(description))
    files = value.get("files")
    git = value.get("git")
    git_tool = value.get("git_tool")
    controls = value.get("repository_controls")
    if (not isinstance(files, dict) or not files or
            value.get("digest") != digest_value({
                "files": files,
                "repository_controls": controls,
            })):
        raise CrossoverError(
            "{} has an invalid tracked-source closure".format(description)
        )
    submodules_clean = True
    tracked_files_clean = True

    def validate_controls(control_value, entries, label):
        expected_keys = {
            "head_index_match", "head_tree_sha256", "index_entry_count",
            "index_flags", "index_sha256", "legacy_grafts_absent",
            "replace_refs",
        }
        index_material = {
            name: [entries[name]["index_mode"], entries[name]["index_object"]]
            for name in sorted(entries)
        }
        if (not isinstance(control_value, dict) or
                set(control_value) != expected_keys or
                type(control_value.get("head_index_match")) is not bool or
                type(control_value.get("legacy_grafts_absent")) is not bool or
                control_value["legacy_grafts_absent"] is not True or
                control_value.get("replace_refs") != [] or
                control_value.get("index_flags") !=
                "all ls-files -v/-f tags are canonical H" or
                type(control_value.get("index_entry_count")) is not int or
                control_value["index_entry_count"] != len(entries) or
                control_value.get("index_sha256") !=
                digest_value(index_material) or
                not isinstance(control_value.get("head_tree_sha256"), str) or
                not re.fullmatch(
                    r"[0-9a-f]{64}", control_value["head_tree_sha256"]) or
                (control_value["head_index_match"] and
                 control_value["head_tree_sha256"] !=
                 control_value["index_sha256"])):
            raise CrossoverError(
                "{} has invalid repository controls for {}".format(
                    description, label
                )
            )
        return control_value["head_index_match"]

    def validate_files(entries, prefix):
        nonlocal submodules_clean, tracked_files_clean
        if not isinstance(entries, dict) or not entries:
            raise CrossoverError(
                "{} {} closure is empty".format(description, prefix)
            )
        if list(entries) != sorted(entries):
            raise CrossoverError(
                "{} {} closure is not canonically ordered".format(
                    description, prefix
                )
            )
        for name, record in entries.items():
            label = "{}{}".format(prefix, name)
            if (not isinstance(name, str) or not name or
                    not isinstance(record, dict) or
                    record.get("index_mode") not in (
                        "100644", "100755", "120000", "160000") or
                    not isinstance(record.get("index_object"), str) or
                    not re.fullmatch(
                        r"[0-9a-f]{40}", record["index_object"])):
                raise CrossoverError(
                    "{} has invalid tracked entry {}".format(description, label)
                )
            if record["index_mode"] == "160000":
                if set(record) != {
                        "index_mode", "index_object", "submodule"}:
                    raise CrossoverError(
                        "{} has invalid gitlink {}".format(description, label)
                    )
                submodule = record.get("submodule")
                if (not isinstance(submodule, dict) or set(submodule) != {
                        "digest", "files", "head", "repository_controls",
                        "status", "worktree_clean"} or
                        submodule.get("head") != record["index_object"] or
                        not isinstance(submodule.get("status"), list) or
                        any(not isinstance(line, str) or not line
                            for line in submodule["status"]) or
                        type(submodule.get("worktree_clean")) is not bool or
                        submodule["worktree_clean"] !=
                        (not bool(submodule["status"])) or
                        submodule.get("digest") != digest_value({
                            "files": submodule.get("files"),
                            "repository_controls":
                                submodule.get("repository_controls"),
                        })):
                    raise CrossoverError(
                        "{} has invalid submodule state {}".format(
                            description, label
                        )
                    )
                submodule_controls_clean = validate_controls(
                    submodule["repository_controls"],
                    submodule["files"], label
                )
                submodules_clean = (
                    submodules_clean and submodule["worktree_clean"] and
                    submodule_controls_clean
                )
                validate_files(submodule["files"], label + "/")
            else:
                if (set(record) != {
                        "index_mode", "index_object", "worktree_git_object",
                        "worktree_mode", "worktree_sha256"} or
                        not isinstance(record.get("worktree_sha256"), str) or
                        not re.fullmatch(
                            r"[0-9a-f]{64}", record["worktree_sha256"]) or
                        not isinstance(
                            record.get("worktree_git_object"), str) or
                        not re.fullmatch(
                            r"[0-9a-f]{40}",
                            record["worktree_git_object"]) or
                        record.get("worktree_mode") not in (
                            "100644", "100755", "120000")):
                    raise CrossoverError(
                        "{} has invalid file state {}".format(
                            description, label
                        )
                    )
                tracked_files_clean = (
                    tracked_files_clean and
                    record["worktree_git_object"] ==
                        record["index_object"] and
                    record["worktree_mode"] == record["index_mode"]
                )

    controls_clean = validate_controls(controls, files, "top-level source")
    validate_files(files, "")
    if (not isinstance(git, dict) or set(git) != {
            "branch", "head", "status", "tree", "worktree_clean"} or
            not isinstance(git.get("status"), list) or
            any(not isinstance(line, str) or not line for line in git["status"]) or
            type(git.get("worktree_clean")) is not bool or
            git["worktree_clean"] != (not bool(git["status"])) or
            not isinstance(git.get("branch"), (str, type(None))) or
            not all(isinstance(git.get(name), str) and
                    re.fullmatch(r"[0-9a-f]{40}", git[name])
                    for name in ("head", "tree"))):
        raise CrossoverError(
            "{} has an invalid exact Git identity".format(description)
        )
    if require_clean and (
            git["worktree_clean"] is not True or not submodules_clean or
            not tracked_files_clean or not controls_clean):
        raise CrossoverError("{} is not a clean Git source".format(description))
    if canonical_bytes(git_tool) != canonical_bytes(canonical_git_identity()):
        raise CrossoverError(
            "{} used a different canonical Git executable".format(description)
        )
    return value


def parse_backends(text):
    result = []
    for item in text.split(","):
        backend = item.strip().lower()
        if backend and backend not in result:
            result.append(backend)
    invalid = [item for item in result if item not in KNOWN_BACKENDS]
    if not result or invalid:
        raise CrossoverError(
            "backends must be a comma-separated subset of {}".format(
                ",".join(KNOWN_BACKENDS)
            )
        )
    return result


def executable_candidates(root, backend):
    rendered = Path(str(root).format(backend=backend))
    roots = [
        rendered,
        rendered / backend,
        rendered / ("direct-encode-" + backend),
    ]
    names = ("bench_leopard2_direct_encode", "bench_leopard2_direct_encode.exe")
    candidates = []
    for directory in roots:
        for name in names:
            candidates.extend((directory / name, directory / "Release" / name))
    unique = []
    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def find_executable(root, backend):
    candidates = executable_candidates(root, backend)
    matches = []
    seen = set()
    for candidate in candidates:
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            resolved = candidate.resolve()
            if str(resolved) not in seen:
                seen.add(str(resolved))
                matches.append(resolved)
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise CrossoverError(
            "benchmark executable lookup for {} is ambiguous: {}".format(
                backend, ", ".join(str(item) for item in matches)
            )
        )
    raise CrossoverError(
        "benchmark executable for {} was not found; checked {}".format(
            backend, ", ".join(str(item) for item in candidates)
        )
    )


def cmake_build_metadata(executable):
    """Fingerprint the CMake inputs nearest an already-built benchmark."""
    cache = None
    directory = executable.parent
    for candidate_root in (directory, directory.parent, directory.parent.parent):
        candidate = candidate_root / "CMakeCache.txt"
        if candidate.is_file():
            cache = candidate
            break
    if cache is None:
        raise CrossoverError(
            "cannot find CMakeCache.txt near benchmark executable {}".format(
                executable
            )
        )
    try:
        cache_bytes = cache.read_bytes()
        cache_lines = cache_bytes.decode("utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise CrossoverError("cannot read {}: {}".format(cache, error))
    prefixes = (
        "CMAKE_BINARY_DIR", "CMAKE_BUILD_TYPE", "CMAKE_C_COMPILER", "CMAKE_C_FLAGS",
        "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS", "CMAKE_GENERATOR",
        "CMAKE_HOME_DIRECTORY",
        "ENABLE_OPENMP", "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
        "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS",
        "LEO2_BENCHMARK_GIT_EXECUTABLE",
        "LEO2_BUILD_FUZZERS", "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    )
    entries = {}
    for line in cache_lines:
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        typed_key, value = line.split("=", 1)
        key = typed_key.split(":", 1)[0]
        if any(key == prefix or key.startswith(prefix + "_") for prefix in prefixes):
            entries[key] = value
    build_root = cache.parent
    extra_files = {}
    for relative in (
            "compile_commands.json", "CMakeFiles/CMakeConfigureLog.yaml",
            "CMakeFiles/CMakeOutput.log", "build.ninja",
            "generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h",
            "CMakeFiles/leopard.dir/flags.make",
            "CMakeFiles/bench_leopard2_direct_encode.dir/DependInfo.cmake",
            "CMakeFiles/bench_leopard2_direct_encode.dir/flags.make",
            "CMakeFiles/bench_leopard2_direct_encode.dir/link.txt",
            "CMakeFiles/bench_leopard2_direct_encode.dir/"
            "bench/leopard2/direct_encode_benchmark.cpp.o",
            "CMakeFiles/leopard_test_hooks.dir/flags.make",
            "CMakeFiles/leopard_test_hooks.dir/link.txt",
            "CMakeFiles/leopard2_backend_avx2_test_hooks.dir/flags.make",
            "libleopard.a", "libleopard_test_hooks.a"):
        path = build_root / relative
        if path.is_file():
            extra_files[relative] = digest_bytes(path.read_bytes())
    for directory in (
            "CMakeFiles/leopard_test_hooks.dir",
            "CMakeFiles/leopard2_backend_ssse3_test_hooks.dir",
            "CMakeFiles/leopard2_backend_avx2_test_hooks.dir",
            "CMakeFiles/leopard2_backend_avx512_test_hooks.dir"):
        root = build_root / directory
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.o")):
            relative = str(path.relative_to(build_root))
            extra_files[relative] = digest_bytes(path.read_bytes())
    executable_stat = executable.stat()
    return {
        "binding_scope": (
            "exact executable, CMake cache, compile database/generator graph, "
            "direct benchmark object/link recipe, and present test-hook "
            "archive/object hashes; embedded clean Git commit/tree attestation "
            "and full tracked-worktree hashes are validated separately"
        ),
        "build_root": str(build_root.resolve()),
        "cmake_cache": str(cache.resolve()),
        "cmake_cache_sha256": digest_bytes(cache_bytes),
        "entries": entries,
        "executable": {
            "mtime_ns": executable_stat.st_mtime_ns,
            "path": str(executable.resolve()),
            "sha256": digest_bytes(executable.read_bytes()),
            "size": executable_stat.st_size,
        },
        "extra_file_sha256": extra_files,
    }


def validate_build_source_binding(
        metadata, source, source_state, backend, require_fresh):
    if not isinstance(metadata, dict):
        raise CrossoverError("CMake build provenance is not an object")
    entries = metadata.get("entries")
    extra = metadata.get("extra_file_sha256")
    executable = metadata.get("executable")
    if not all(isinstance(value, dict) for value in (entries, extra, executable)):
        raise CrossoverError("CMake build provenance is incomplete")
    home_directory = entries.get("CMAKE_HOME_DIRECTORY")
    if (not isinstance(home_directory, str) or not home_directory or
            not Path(home_directory).is_absolute() or
            Path(home_directory).resolve() != source.resolve()):
        raise CrossoverError("CMake build source root differs from --source")
    if entries.get("LEO2_BACKEND_VARIANT", "").lower() != backend:
        raise CrossoverError(
            "CMake LEO2_BACKEND_VARIANT is {!r}, expected {!r}".format(
                entries.get("LEO2_BACKEND_VARIANT"), backend
            )
        )
    if entries.get("LEO2_BUILD_BENCHMARKS", "").upper() not in ("1", "ON", "TRUE"):
        raise CrossoverError("CMake build did not explicitly enable benchmarks")
    if not require_fresh:
        return
    git = source_state.get("git", {})
    if git.get("worktree_clean") is not True:
        raise CrossoverError("authoritative build requires a clean Git source tree")
    if entries.get("CMAKE_BUILD_TYPE", "").lower() != "release":
        raise CrossoverError("authoritative build must use CMAKE_BUILD_TYPE=Release")
    if entries.get("LEO2_BUILD_TESTS", "").upper() not in ("1", "ON", "TRUE"):
        raise CrossoverError("authoritative build did not explicitly enable tests")
    embedded_git = entries.get("LEO2_BENCHMARK_GIT_EXECUTABLE")
    if (not isinstance(embedded_git, str) or not embedded_git or
            Path(embedded_git).resolve() != CANONICAL_GIT):
        raise CrossoverError(
            "authoritative benchmark attestation did not use canonical "
            "{}".format(CANONICAL_GIT)
        )
    required = {
        "compile_commands.json",
        "CMakeFiles/bench_leopard2_direct_encode.dir/"
        "bench/leopard2/direct_encode_benchmark.cpp.o",
        "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h",
    }
    missing = sorted(required - set(extra))
    if missing:
        raise CrossoverError(
            "authoritative build provenance omits {}".format(", ".join(missing))
        )
    if not any(name in extra for name in (
            "build.ninja",
            "CMakeFiles/bench_leopard2_direct_encode.dir/link.txt")):
        raise CrossoverError("authoritative build provenance omits the link graph")
    if not any(name in extra for name in (
            "libleopard.a", "libleopard_test_hooks.a")):
        raise CrossoverError("authoritative build provenance omits the linked archive")
    def source_mtimes(root, entries):
        for relative, record in entries.items():
            path = root / relative
            yield path.lstat().st_mtime_ns
            if record["index_mode"] == "160000":
                yield from source_mtimes(
                    path.resolve(), record["submodule"]["files"]
                )

    try:
        newest_source_mtime = max(
            source_mtimes(source, source_state["files"])
        )
    except OSError as error:
        raise CrossoverError(
            "cannot establish authoritative source freshness: {}".format(error)
        )
    if executable.get("mtime_ns", 0) < newest_source_mtime:
        raise CrossoverError(
            "authoritative executable predates the captured source closure; "
            "use a fresh dedicated build directory"
        )


def frozen_executable_identity(
        backend, origin_executable, executable_sha256, build_metadata,
        source_state):
    return {
        "backend": backend,
        "build_metadata": build_metadata,
        "executable_sha256": executable_sha256,
        "origin_executable": str(Path(origin_executable).resolve()),
        "source_fingerprint": source_state,
    }


def validate_frozen_executable(artifact, build_metadata=None, source_state=None):
    expected_keys = {
        "artifact_id", "backend", "directory", "executable",
        "executable_sha256", "provenance", "provenance_sha256", "schema",
    }
    if (not isinstance(artifact, dict) or set(artifact) != expected_keys or
            artifact.get("schema") != FROZEN_EXECUTABLE_SCHEMA):
        raise CrossoverError("frozen executable artifact has an invalid schema")
    artifact_id = artifact.get("artifact_id")
    if (not isinstance(artifact_id, str) or
            not re.fullmatch(r"[0-9a-f]{64}", artifact_id)):
        raise CrossoverError("frozen executable artifact ID is invalid")
    directory = Path(artifact["directory"]).resolve()
    executable = Path(artifact["executable"]).resolve()
    provenance_path = Path(artifact["provenance"]).resolve()
    if executable.parent != directory or provenance_path.parent != directory:
        raise CrossoverError("frozen executable paths escape the artifact directory")
    try:
        executable_bytes = executable.read_bytes()
        provenance_bytes = provenance_path.read_bytes()
        directory_mode = directory.stat().st_mode
        executable_mode = executable.stat().st_mode
        provenance_mode = provenance_path.stat().st_mode
    except OSError as error:
        raise CrossoverError("cannot read frozen executable artifact: {}".format(error))
    if any(mode & 0o222 for mode in (
            directory_mode, executable_mode, provenance_mode)):
        raise CrossoverError("frozen executable artifact remains writable")
    executable_sha256 = digest_bytes(executable_bytes)
    if (executable_sha256 != artifact.get("executable_sha256") or
            not os.access(str(executable), os.X_OK)):
        raise CrossoverError("frozen executable hash or mode is invalid")
    if digest_bytes(provenance_bytes) != artifact.get("provenance_sha256"):
        raise CrossoverError("frozen executable provenance hash is invalid")
    provenance = decode_json_bytes(provenance_bytes, str(provenance_path))
    provenance_keys = {
        "artifact_id", "backend", "build_metadata", "executable",
        "origin_executable", "schema", "source_fingerprint",
    }
    if (not isinstance(provenance, dict) or set(provenance) != provenance_keys or
            provenance.get("schema") != FROZEN_EXECUTABLE_SCHEMA or
            provenance.get("artifact_id") != artifact_id or
            provenance.get("backend") != artifact.get("backend") or
            provenance.get("executable") != {
                "name": executable.name, "sha256": executable_sha256,
            }):
        raise CrossoverError("frozen executable provenance is inconsistent")
    origin = provenance.get("origin_executable")
    if (not isinstance(origin, dict) or set(origin) != {"path", "sha256"} or
            origin.get("sha256") != executable_sha256):
        raise CrossoverError("frozen executable origin provenance is invalid")
    identity = frozen_executable_identity(
        artifact["backend"], origin["path"], executable_sha256,
        provenance["build_metadata"], provenance["source_fingerprint"]
    )
    if digest_value(identity) != artifact_id:
        raise CrossoverError("frozen executable artifact ID does not match provenance")
    if (build_metadata is not None and
            canonical_bytes(provenance["build_metadata"]) !=
            canonical_bytes(build_metadata)):
        raise CrossoverError("frozen executable build provenance is stale")
    if (source_state is not None and
            canonical_bytes(provenance["source_fingerprint"]) !=
            canonical_bytes(source_state)):
        raise CrossoverError("frozen executable source provenance is stale")
    return artifact


def freeze_executable(
        result_dir, backend, executable, build_metadata, source_state):
    executable = Path(executable).resolve()
    try:
        executable_bytes = executable.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot read executable to freeze {}: {}".format(executable, error)
        )
    current_build_metadata = cmake_build_metadata(executable)
    if canonical_bytes(current_build_metadata) != canonical_bytes(build_metadata):
        raise CrossoverError(
            "origin executable or CMake metadata changed before freezing"
        )
    executable_sha256 = digest_bytes(executable_bytes)
    identity = frozen_executable_identity(
        backend, executable, executable_sha256, build_metadata, source_state
    )
    artifact_id = digest_value(identity)
    frozen_root = Path(result_dir).resolve() / "frozen-executables"
    artifact_dir = frozen_root / "{}-{}".format(backend, artifact_id[:24])
    frozen_name = "bench_leopard2_direct_encode"
    frozen_path = artifact_dir / frozen_name
    provenance_path = artifact_dir / "provenance.json"
    provenance = {
        "artifact_id": artifact_id,
        "backend": backend,
        "build_metadata": build_metadata,
        "executable": {
            "name": frozen_name,
            "sha256": executable_sha256,
        },
        "origin_executable": {
            "path": str(executable),
            "sha256": executable_sha256,
        },
        "schema": FROZEN_EXECUTABLE_SCHEMA,
        "source_fingerprint": source_state,
    }

    if not artifact_dir.exists():
        frozen_root.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(
            prefix=".{}-".format(backend), dir=str(frozen_root)
        ))
        try:
            temporary_executable = temporary / frozen_name
            temporary_executable.write_bytes(executable_bytes)
            atomic_write_json(temporary / "provenance.json", provenance)
            temporary_executable.chmod(0o555)
            (temporary / "provenance.json").chmod(0o444)
            temporary.chmod(0o555)
            try:
                os.replace(str(temporary), str(artifact_dir))
            except OSError:
                if not artifact_dir.is_dir():
                    raise
                temporary.chmod(0o755)
                shutil.rmtree(str(temporary))
        except BaseException:
            if temporary.exists():
                try:
                    temporary.chmod(0o755)
                except OSError:
                    pass
                shutil.rmtree(str(temporary), ignore_errors=True)
            raise

    try:
        origin_bytes_after = executable.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot re-read origin executable after freezing: {}".format(error)
        )
    if (origin_bytes_after != executable_bytes or
            canonical_bytes(cmake_build_metadata(executable)) !=
            canonical_bytes(build_metadata)):
        raise CrossoverError(
            "origin executable or CMake metadata changed while freezing"
        )

    try:
        provenance_sha256 = digest_bytes(provenance_path.read_bytes())
    except OSError as error:
        raise CrossoverError(
            "cannot read frozen provenance {}: {}".format(provenance_path, error)
        )
    artifact = {
        "artifact_id": artifact_id,
        "backend": backend,
        "directory": str(artifact_dir),
        "executable": str(frozen_path),
        "executable_sha256": executable_sha256,
        "provenance": str(provenance_path),
        "provenance_sha256": provenance_sha256,
        "schema": FROZEN_EXECUTABLE_SCHEMA,
    }
    return validate_frozen_executable(artifact, build_metadata, source_state)


def cell(region, backend, k, r, profile, field, shard_bytes, q, mask):
    return {
        "backend": backend,
        "field": field,
        "k": k,
        "mask": mask,
        "profile": profile,
        "q": q,
        "r": r,
        "region": region,
        "shard_bytes": shard_bytes,
    }


def threshold_k(backend):
    return 3 if backend == "scalar" else 2


def compact_grid(backends, r):
    """Cells at each promoted condition and its deliberately excluded neighbors."""
    result = []
    for backend in backends:
        minimum_k = threshold_k(backend)
        # Threshold, large-shard behavior, and the bounded edge in both fields.
        for field in ("gf8", "gf16"):
            for k, shard_bytes in (
                    (minimum_k, 1024), (minimum_k, 65536), (16, 4096)):
                result.append(cell(
                    "candidate", backend, k, r, "low", field,
                    shard_bytes, 1, "0"))
        # One adjacent complete tile and an equal-Q sparse mask expose the two
        # AUTO inputs (byte shape and actual output count) independently.
        result.append(cell(
            "candidate", backend, minimum_k, r, "low", "gf8",
            1088, 1, "0"))
        if r > 1:
            result.append(cell(
                "candidate_sparse_output", backend, minimum_k, r,
                "low", "gf8", 4096, 1, str(r - 1)))
        # These neighbors are intentionally outside AUTO and guard each clause.
        if backend == "scalar":
            for field, shard_bytes in (
                    ("gf8", 1024), ("gf8", 65536), ("gf16", 4096)):
                result.append(cell(
                    "excluded_scalar_k2", backend, 2, r, "low", field,
                    shard_bytes, 1, "0"))
        if r >= 2:
            result.append(cell(
                "excluded_q2", backend, minimum_k, r, "low", "gf8",
                4096, 2, "0-1"))
        if r >= 3:
            result.append(cell(
                "excluded_q2_holey", backend, minimum_k, r, "low", "gf8",
                4096, 2, "0,{}".format(r - 1)))
        result.append(cell(
            "excluded_high_profile", backend, minimum_k, r, "high", "gf8",
            4096, 1, "0"))
        for shard_bytes in (63, 1023):
            result.append(cell(
                "excluded_ragged_tail", backend, minimum_k, r, "low", "gf8",
                shard_bytes, 1, "0"))
    return result


def historical_avx2_grid():
    """Corrected high-profile AVX2 cells, including the omitted 17.88x cell."""
    values = [
        # Exact historical cell omitted by the invalid 2026-07-28 screen.
        cell(
            "candidate_historical_exact", "avx2", 2, 16, "high",
            "gf8", 4096, 1, "0"
        ),
        # The nine retained screen cells, now interpreted only as campaign shape.
        cell("excluded_historical_single_side_control", "avx2", 2, 1, "high", "gf8", 1024, 1, "0"),
        cell("excluded_historical_single_side_control", "avx2", 3, 1, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 4, 4, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 8, 8, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 1024, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 8, 8, "high", "gf8", 65536, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 4, "high", "gf8", 4096, 1, "0"),
        cell("candidate_historical_prior_screen", "avx2", 16, 16, "high", "gf8", 65536, 1, "0"),
        # One-coordinate neighbors make K, R, byte size, mask position, and Q
        # independently visible around the historical cell.
        cell("historical_neighbor_k", "avx2", 3, 16, "high", "gf8", 4096, 1, "0"),
        cell("historical_neighbor_r", "avx2", 2, 15, "high", "gf8", 4096, 1, "0"),
        cell(
            "historical_neighbor_bytes_lower", "avx2", 2, 16, "high", "gf8",
            4032, 1, "0"
        ),
        cell(
            "historical_neighbor_bytes_upper", "avx2", 2, 16, "high", "gf8",
            4160, 1, "0"
        ),
        cell(
            "historical_neighbor_sparse_output", "avx2", 2, 16, "high",
            "gf8", 4096, 1, "15"
        ),
        cell(
            "historical_neighbor_q2", "avx2", 2, 16, "high", "gf8",
            4096, 2, "0-1"
        ),
    ]
    identities = {digest_value(value) for value in values}
    if len(identities) != len(values):
        raise CrossoverError("historical AVX2 campaign contains duplicate cells")
    return values


def full_grid(backends, r):
    result = compact_grid(backends, r)
    seen = {digest_value(item) for item in result}

    def add(item):
        identity = digest_value(item)
        if identity not in seen:
            seen.add(identity)
            result.append(item)

    for backend in backends:
        minimum_k = threshold_k(backend)
        candidate_ks = list(range(minimum_k, 17))
        for field in ("gf8", "gf16"):
            for k in candidate_ks:
                for shard_bytes in (1024, 1088, 2048, 4096, 16384, 65536, 1048576):
                    add(cell(
                        "candidate", backend, k, r, "low", field,
                        shard_bytes, 1, "0"))
            # A one-output sparse request exercises a non-prefix parity row.
            for k in (minimum_k, 4, 8, 16):
                if k < minimum_k:
                    continue
                add(cell(
                    "candidate_sparse_output", backend, k, r, "low", field,
                    4096, 1, str(r - 1)))
        for k in (minimum_k, 4, 8, 16):
            if k < minimum_k:
                continue
            for q in (2, 4, 8):
                if q > r:
                    continue
                add(cell(
                    "excluded_q_gt_1", backend, k, r, "low", "gf8",
                    4096, q, "0-{}".format(q - 1)))
            for profile in ("high",):
                add(cell(
                    "excluded_high_profile", backend, k, r, profile, "gf8",
                    4096, 1, "0"))
        for shard_bytes in (1, 63, 960, 1023, 1025, 1087, 1089, 4095, 4097):
            add(cell(
                "excluded_byte_boundary", backend, minimum_k, r, "low", "gf8",
                shard_bytes, 1, "0"))
        if backend == "scalar":
            for field in ("gf8", "gf16"):
                for shard_bytes in (1024, 4096, 65536, 1048576):
                    add(cell(
                        "excluded_scalar_k2", backend, 2, r, "low", field,
                        shard_bytes, 1, "0"))
        for k in range(1, minimum_k):
            for field in ("gf8", "gf16"):
                add(cell(
                    "excluded_k_below_auto", backend, k, r, "low", field,
                    4096, 1, "0"))
    return result


def parse_r_values(text):
    if text.strip().lower() == "all":
        return list(range(1, 17))
    result = []
    for item in text.split(","):
        try:
            value = int(item.strip())
        except ValueError:
            raise CrossoverError("R must be 'all' or a comma-separated integer list")
        if value < 1 or value > 16:
            raise CrossoverError("every R must be in [1,16]")
        if value not in result:
            result.append(value)
    if not result:
        raise CrossoverError("at least one R value is required")
    return sorted(result)


def sorted_grid(backends, r_values, full):
    values = []
    for r in r_values:
        values.extend(full_grid(backends, r) if full else compact_grid(backends, r))
    return sorted(values, key=lambda item: (
        item["backend"], item["region"], item["profile"], item["field"],
        item["r"], item["k"], item["q"], item["shard_bytes"], item["mask"]
    ))


def stable_seed(cell_value):
    # Keep zero reserved and remain within the benchmark's uint64 parser.
    return int(digest_value({"seed_for": cell_value})[:16], 16) or 1


def invocation_order(mode, job_id, abba_rounds):
    if mode in AUTHORITATIVE_COMMANDS:
        return ("direct", "transform", "transform", "direct") * abba_rounds
    # Alternating the two-cell screening order prevents one path from always
    # receiving the colder process launch while keeping the order deterministic.
    return (("direct", "transform") if int(job_id[:8], 16) % 2 == 0
            else ("transform", "direct"))


def job_identity(
        cell_value, executable, executable_artifact, executable_sha256,
        build_metadata, source_state, machine, settings):
    return {
        "cell": cell_value,
        "build_metadata": build_metadata,
        "executable": str(executable),
        "executable_artifact": executable_artifact,
        "executable_sha256": executable_sha256,
        "machine": machine,
        "settings": settings,
        "source_identity": source_state,
    }


def make_jobs(
        cells, executables, build_metadata, source_state, machine, settings,
        executable_artifacts=None):
    executable_artifacts = executable_artifacts or {}
    jobs = []
    for cell_value in cells:
        executable = executables[cell_value["backend"]]
        executable_sha256 = digest_bytes(executable.read_bytes())
        executable_artifact = executable_artifacts.get(cell_value["backend"])
        identity = job_identity(
            cell_value, executable, executable_artifact, executable_sha256,
            build_metadata[cell_value["backend"]], source_state, machine,
            settings
        )
        configuration_id = digest_value(identity)
        job_id = configuration_id[:24]
        jobs.append({
            "cell": cell_value,
            "build_metadata": build_metadata[cell_value["backend"]],
            "configuration_id": configuration_id,
            "executable": str(executable),
            "executable_artifact": executable_artifact,
            "executable_sha256": executable_sha256,
            "invocation_order": list(invocation_order(
                settings["mode"], job_id, settings["abba_rounds"])),
            "job_id": job_id,
            "seed": stable_seed(cell_value),
            "source_identity": source_state,
        })
    return jobs


def benchmark_argv(job, timed_mode, raw_path, settings):
    item = job["cell"]
    benchmark = settings["benchmark"]
    argv = [
        job["executable"],
        "--k", str(item["k"]),
        "--r", str(item["r"]),
        "--profile", item["profile"],
        "--field", item["field"],
        "--bytes", str(item["shard_bytes"]),
        "--q", str(item["q"]),
        "--requested-parity", item["mask"],
        "--batch", str(benchmark["batch"]),
        "--reuse", str(benchmark["reuse"]),
        "--iterations", str(benchmark["iterations"]),
        "--warmups", str(benchmark["warmups"]),
        "--threads", "1",
        "--seed", str(job["seed"]),
        "--mode", timed_mode,
        "--json", str(raw_path),
    ]
    if settings["mode"] in AUTHORITATIVE_COMMANDS:
        argv = [settings["taskset"], "-c", str(settings["pin_cpu"])] + argv
    return argv


def run_command(argv, cwd, stdout_path, stderr_path, timeout, environment):
    timed_out = False
    try:
        completed = subprocess.run(
            [str(item) for item in argv], cwd=str(cwd), env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout
        )
        returncode = completed.returncode
        stdout = normalized_output(completed.stdout)
        stderr = normalized_output(completed.stderr)
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = 124
        stdout = normalized_output(error.stdout or b"")
        stderr = normalized_output(error.stderr or b"")
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    return {
        "argv": [str(item) for item in argv],
        "cwd": str(cwd),
        "returncode": returncode,
        "stderr_log": stderr_path.name,
        "stderr_sha256": digest_bytes(stderr),
        "stdout_log": stdout_path.name,
        "stdout_sha256": digest_bytes(stdout),
        "timed_out": timed_out,
    }


def controlled_avx2_build(
        source, result_dir, source_state, parallel, validate_guard,
        guard_identity):
    if type(parallel) is not int or parallel <= 0 or parallel > 128:
        raise CrossoverError("controlled build parallelism is outside [1,128]")
    if not callable(validate_guard):
        raise CrossoverError("controlled build requires a live guard validator")
    if (not isinstance(guard_identity, dict) or set(guard_identity) != {
            "canonical_lock", "pair_lease"}):
        raise CrossoverError(
            "controlled build requires exact campaign guard identities"
        )
    validate_guard()
    build_root = result_dir / "controlled-build-avx2"
    record_path = result_dir / "controlled-build.json"
    if build_root.exists() or record_path.exists():
        if not build_root.is_dir() or not record_path.is_file():
            raise CrossoverError(
                "controlled build is only partially retained; use a new "
                "--result-dir"
            )
        try:
            record_bytes = record_path.read_bytes()
            record = decode_json_bytes(record_bytes, str(record_path))
            executable = Path(record["executable"]).resolve()
            metadata = cmake_build_metadata(executable)
        except (CrossoverError, KeyError, OSError, TypeError) as error:
            raise CrossoverError(
                "cannot resume controlled AVX2 build: {}".format(error)
            )
        if (not isinstance(record, dict) or set(record) != {
                "backend", "build_metadata", "build_root", "commands",
                "executable", "executable_sha256", "guard_identity",
                "schema", "source_identity"} or
                record.get("schema") != CONTROLLED_BUILD_SCHEMA or
                record.get("backend") != "avx2" or
                record.get("build_root") != str(build_root.resolve()) or
                canonical_bytes(record.get("source_identity")) !=
                canonical_bytes(source_state) or
                canonical_bytes(record.get("guard_identity")) !=
                canonical_bytes(guard_identity) or
                canonical_bytes(record.get("build_metadata")) !=
                canonical_bytes(metadata) or
                record.get("executable_sha256") !=
                digest_bytes(executable.read_bytes())):
            raise CrossoverError(
                "retained controlled AVX2 build differs from this campaign"
            )
        validate_build_source_binding(
            metadata, source, source_state, "avx2", require_fresh=True
        )
        validate_guard()
        return executable, {
            "path": str(record_path.relative_to(result_dir)),
            "schema": CONTROLLED_BUILD_SCHEMA,
            "sha256": digest_bytes(record_bytes),
        }
    cmake = shutil.which("cmake", path=BENCHMARK_ENVIRONMENT["PATH"])
    ninja = shutil.which("ninja", path=BENCHMARK_ENVIRONMENT["PATH"])
    if not cmake or not ninja:
        raise CrossoverError(
            "controlled authoritative build requires /usr/bin cmake and ninja"
        )
    log_dir = result_dir / "build-logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    configure = [
        str(Path(cmake).resolve()),
        "-S", str(source),
        "-B", str(build_root),
        "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DLEO2_BACKEND_VARIANT=avx2",
        "-DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git",
        "-DLEO2_BUILD_BENCHMARKS=ON",
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
        "-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF",
        "-DLEOPARD_ENABLE_GF8=ON",
        "-DLEOPARD_ENABLE_GF16=ON",
    ]
    build = [
        str(Path(cmake).resolve()), "--build", str(build_root),
        "--target", "bench_leopard2_direct_encode",
        "--parallel", str(parallel),
    ]
    commands = []
    for label, argv in (("configure", configure), ("build", build)):
        validate_guard()
        stdout_path = log_dir / (label + ".stdout.log")
        stderr_path = log_dir / (label + ".stderr.log")
        command = run_command(
            argv, source, stdout_path, stderr_path, 1800,
            dict(GIT_ENVIRONMENT)
        )
        command["environment"] = dict(GIT_ENVIRONMENT)
        command["stdout_log"] = str(stdout_path.relative_to(result_dir))
        command["stderr_log"] = str(stderr_path.relative_to(result_dir))
        commands.append(command)
        validate_guard()
        if command["returncode"] != 0 or command["timed_out"]:
            raise CrossoverError(
                "controlled AVX2 {} failed; see {} and {}".format(
                    label, stdout_path, stderr_path
                )
            )
    executable = (build_root / "bench_leopard2_direct_encode").resolve()
    if not executable.is_file() or not os.access(str(executable), os.X_OK):
        raise CrossoverError(
            "controlled AVX2 build omitted {}".format(executable)
        )
    metadata = cmake_build_metadata(executable)
    validate_build_source_binding(
        metadata, source, source_state, "avx2", require_fresh=True
    )
    validate_guard()
    record = {
        "backend": "avx2",
        "build_metadata": metadata,
        "build_root": str(build_root.resolve()),
        "commands": commands,
        "executable": str(executable),
        "executable_sha256": digest_bytes(executable.read_bytes()),
        "guard_identity": guard_identity,
        "schema": CONTROLLED_BUILD_SCHEMA,
        "source_identity": source_state,
    }
    atomic_write_json(record_path, record)
    try:
        validate_guard()
    except Exception:
        try:
            record_path.unlink()
        except FileNotFoundError:
            pass
        raise
    record_bytes = record_path.read_bytes()
    return executable, {
        "path": str(record_path.relative_to(result_dir)),
        "schema": CONTROLLED_BUILD_SCHEMA,
        "sha256": digest_bytes(record_bytes),
    }


def requested_indices(mask, r):
    indices = set()
    for part in mask.split(","):
        bounds = part.split("-", 1)
        try:
            first = int(bounds[0])
            last = int(bounds[1]) if len(bounds) == 2 else first
        except ValueError:
            raise CrossoverError("invalid generated parity mask {!r}".format(mask))
        if first < 0 or last < first or last >= r:
            raise CrossoverError("generated parity mask is out of range: {!r}".format(mask))
        indices.update(range(first, last + 1))
    return sorted(indices)


def required_mapping(value, path):
    if not isinstance(value, dict):
        raise CrossoverError("benchmark {} must be a JSON object".format(path))
    return value


def build_configuration_digest(entries):
    if not isinstance(entries, dict):
        raise CrossoverError("CMake entries for build digest are not an object")
    variables = (
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_FLAGS",
        "CMAKE_CXX_FLAGS_DEBUG",
        "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO",
        "CMAKE_CXX_FLAGS_MINSIZEREL",
        "ENABLE_OPENMP",
        "LEOPARD_ENABLE_GF8",
        "LEOPARD_ENABLE_GF16",
        "LEO2_BACKEND_VARIANT",
        "LEO2_BENCHMARK_GIT_EXECUTABLE",
        "LEO2_BUILD_BENCHMARKS",
        "LEO2_BUILD_TESTS",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    )
    material = "".join(
        "{}={}\n".format(variable, entries.get(variable, ""))
        for variable in variables
    )
    return digest_bytes(material.encode("utf-8"))


def require_exact_keys(value, expected, path):
    mapping = required_mapping(value, path)
    actual = set(mapping)
    expected = set(expected)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise CrossoverError(
            "benchmark {} has invalid v2 fields (missing={}, extra={})".format(
                path, missing, extra
            )
        )
    return mapping


def required_finite_number(value, path):
    number = None
    if not isinstance(value, bool) and isinstance(value, (int, float)):
        try:
            number = float(value)
        except (OverflowError, TypeError, ValueError):
            number = None
    if number is None or not math.isfinite(number):
        raise CrossoverError(
            "benchmark {} must be a finite JSON number".format(path)
        )
    return number


def required_finite_metric(value, path, allow_zero):
    number = required_finite_number(value, path)
    if (number < 0) if allow_zero else (number <= 0):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise CrossoverError(
            "benchmark {} must be a finite {} JSON number".format(path, qualifier)
        )
    return number


def validate_raw(raw, job, timed_mode, settings):
    raw = require_exact_keys(raw, {
        "build", "correctness", "memory", "methodology", "metrics",
        "operation_model", "parameters", "resolved", "schema",
    }, "document")
    if raw.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    build = require_exact_keys(raw.get("build"), {
        "backend_variant", "build_configuration_sha256", "build_type",
        "compiler", "compiler_version", "cplusplus", "source_commit",
        "source_tracked_dirty", "source_tree", "test_hooks",
    }, "build")
    parameters = require_exact_keys(raw.get("parameters"), {
        "K", "Q", "R", "batch", "forced_mode", "iterations",
        "requested_field", "requested_parity_indices", "requested_profile",
        "reuse", "seed", "shard_bytes", "thread_count", "warmups",
    }, "parameters")
    resolved = require_exact_keys(raw.get("resolved"), {
        "backend", "direct_capable", "field", "padded_side", "parent_count",
        "profile", "thread_count", "timed_path_is_direct",
    }, "resolved")
    correctness = require_exact_keys(raw.get("correctness"), {
        "direct_transform_parity_match", "parity_checksum_fnv1a64",
        "selected_transform_reference_parity_match",
        "unrequested_outputs_untouched",
    }, "correctness")
    operation_model = require_exact_keys(raw.get("operation_model"), {
        "direct_model_applies_to_timed_path",
        "direct_output_accumulations", "direct_output_initializations",
        "direct_row_terms",
        "fixed_coefficient_symbol_terms_before_unit_specialization",
        "hardware_counters", "model_scope", "modeled_output_bytes_read",
        "modeled_output_bytes_written", "modeled_source_bytes_read",
        "transform_operation_counts", "xor_accumulation_symbols",
    }, "operation_model")
    required_mapping(raw.get("memory"), "memory")
    required_mapping(raw.get("methodology"), "methodology")
    metrics = required_mapping(raw.get("metrics"), "metrics")
    execution = require_exact_keys(metrics.get("encode_execution"), {
        "logical_input_GB_per_s", "logical_input_plus_output_GB_per_s",
        "mad_us_per_batch_call", "maximum_us_per_batch_call",
        "median_us_per_batch_call", "minimum_us_per_batch_call",
        "requested_parity_output_GB_per_s",
    }, "metrics.encode_execution")
    item = job["cell"]
    try:
        git_identity = job["source_identity"]["git"]
        expected_build_identity = {
            "source_commit": git_identity["head"],
            "source_tree": git_identity["tree"],
            "source_tracked_dirty": not git_identity["worktree_clean"],
        }
    except (KeyError, TypeError):
        raise CrossoverError("job omits its exact Git source identity")
    for key, value in expected_build_identity.items():
        if canonical_bytes(build.get(key)) != canonical_bytes(value):
            raise CrossoverError(
                "benchmark build {} is {!r}, expected {!r}".format(
                    key, build.get(key), value
                )
            )
    if build.get("test_hooks") is not True:
        raise CrossoverError("benchmark build omits required test hooks")
    if build.get("backend_variant") != item["backend"]:
        raise CrossoverError(
            "benchmark executable variant {!r} differs from cell {!r}".format(
                build.get("backend_variant"), item["backend"]
            )
        )
    try:
        entries = job["build_metadata"]["entries"]
    except (KeyError, TypeError):
        raise CrossoverError("job omits CMake entries for build attestation")
    expected_build_type = entries.get("CMAKE_BUILD_TYPE", "")
    if (build.get("build_type") != expected_build_type or
            build.get("build_configuration_sha256") !=
            build_configuration_digest(entries)):
        raise CrossoverError(
            "benchmark executable build configuration differs from CMake"
        )
    if (not isinstance(build.get("compiler"), str) or
            not build["compiler"] or
            not isinstance(build.get("compiler_version"), str) or
            not build["compiler_version"] or
            type(build.get("cplusplus")) is not int or
            build["cplusplus"] < 201103):
        raise CrossoverError("benchmark build metadata has invalid exact types")
    parity_indices = requested_indices(item["mask"], item["r"])
    if not parity_indices or len(parity_indices) != item["q"]:
        raise CrossoverError(
            "generated parity mask must select exactly Q nonempty outputs"
        )
    expected = {
        "K": item["k"], "R": item["r"], "Q": item["q"],
        "forced_mode": "force_" + timed_mode,
        "shard_bytes": item["shard_bytes"],
        "seed": job["seed"],
        "requested_parity_indices": parity_indices,
        "requested_profile": ("low_v1" if item["profile"] == "low"
                              else "legacy_high_v1"),
        "requested_field": item["field"],
        "batch": settings["benchmark"]["batch"],
        "reuse": settings["benchmark"]["reuse"],
        "iterations": settings["benchmark"]["iterations"],
        "warmups": settings["benchmark"]["warmups"],
        "thread_count": 1,
    }
    for key, value in expected.items():
        if canonical_bytes(parameters.get(key)) != canonical_bytes(value):
            raise CrossoverError(
                "benchmark parameter {} is {!r}, expected {!r}".format(
                    key, parameters.get(key), value
                )
            )
    if resolved.get("profile") != ("low_v1" if item["profile"] == "low"
                                    else "legacy_high_v1"):
        raise CrossoverError("benchmark resolved an unexpected profile")
    if resolved.get("field") != item["field"]:
        raise CrossoverError("benchmark resolved an unexpected field")
    if item["backend"] != "auto" and resolved.get("backend") != item["backend"]:
        raise CrossoverError(
            "benchmark resolved backend {!r}, expected {!r}".format(
                resolved.get("backend"), item["backend"]
            )
        )
    if correctness.get("selected_transform_reference_parity_match") is not True:
        raise CrossoverError("benchmark selected/reference parity check failed")
    if (timed_mode == "direct" and
            correctness.get("direct_transform_parity_match") is not True):
        raise CrossoverError("benchmark direct/transform correctness check failed")
    if (timed_mode == "transform" and
            correctness.get("direct_transform_parity_match") is not None):
        raise CrossoverError(
            "force-transform must not claim a direct/transform comparison"
        )
    if correctness.get("unrequested_outputs_untouched") is not True:
        raise CrossoverError("benchmark unrequested-output guard check failed")
    checksum = correctness.get("parity_checksum_fnv1a64")
    if not isinstance(checksum, str) or not re.fullmatch(r"0x[0-9a-fA-F]{16}", checksum):
        raise CrossoverError("benchmark omitted a well-formed parity checksum")
    checksum = checksum.lower()
    if (resolved.get("direct_capable") is not True or
            type(resolved.get("thread_count")) is not int or
            resolved.get("thread_count") != 1 or
            type(resolved.get("parent_count")) is not int or
            resolved.get("parent_count") <= 0 or
            type(resolved.get("padded_side")) is not int or
            resolved.get("padded_side") <= 0):
        raise CrossoverError("benchmark resolved unexpected capability or thread count")
    expected_direct = timed_mode == "direct"
    if resolved.get("timed_path_is_direct") is not expected_direct:
        raise CrossoverError("benchmark timed a path other than the forced path")

    batch = settings["benchmark"]["batch"]
    direct_terms = item["k"] * item["q"] * batch
    direct_initializations = item["q"] * batch
    direct_accumulations = (item["k"] - 1) * item["q"] * batch
    symbol_bytes = 2 if item["field"] == "gf16" else 1
    if item["shard_bytes"] % symbol_bytes != 0:
        raise CrossoverError("benchmark cell contains a partial field symbol")
    symbols_per_shard = item["shard_bytes"] // symbol_bytes
    expected_operation_model = {
        "direct_model_applies_to_timed_path": expected_direct,
        "direct_output_accumulations": direct_accumulations,
        "direct_output_initializations": direct_initializations,
        "direct_row_terms": direct_terms,
        "fixed_coefficient_symbol_terms_before_unit_specialization":
            direct_terms * symbols_per_shard,
        "hardware_counters": None,
        "model_scope": (
            "direct streaming kernels before cache effects; unit "
            "coefficients specialize to copy/XOR"
        ),
        "modeled_output_bytes_read":
            direct_accumulations * item["shard_bytes"],
        "modeled_output_bytes_written":
            direct_terms * item["shard_bytes"],
        "modeled_source_bytes_read":
            direct_terms * item["shard_bytes"],
        "transform_operation_counts": None,
        "xor_accumulation_symbols":
            direct_accumulations * symbols_per_shard,
    }
    if canonical_bytes(operation_model) != canonical_bytes(
            expected_operation_model):
        raise CrossoverError(
            "benchmark operation_model differs from the exact v2 model"
        )

    median_us = required_finite_metric(
        execution.get("median_us_per_batch_call"),
        "metrics.encode_execution.median_us_per_batch_call", False
    )
    mad_us = required_finite_metric(
        execution.get("mad_us_per_batch_call"),
        "metrics.encode_execution.mad_us_per_batch_call", True
    )
    minimum_us = required_finite_metric(
        execution.get("minimum_us_per_batch_call"),
        "metrics.encode_execution.minimum_us_per_batch_call", False
    )
    maximum_us = required_finite_metric(
        execution.get("maximum_us_per_batch_call"),
        "metrics.encode_execution.maximum_us_per_batch_call", False
    )
    if minimum_us > median_us or median_us > maximum_us:
        raise CrossoverError(
            "benchmark encode-execution minimum/median/maximum are inconsistent"
        )
    for key in (
            "logical_input_GB_per_s",
            "logical_input_plus_output_GB_per_s",
            "requested_parity_output_GB_per_s"):
        required_finite_metric(
            execution.get(key), "metrics.encode_execution." + key, True
        )
    hashed_bytes = (
        settings["benchmark"]["batch"] * item["q"] * item["shard_bytes"]
    )
    if hashed_bytes <= 0:
        raise CrossoverError("benchmark parity identity covers no output bytes")
    parity_identity = {
        "algorithm": "fnv1a64",
        "digest": checksum,
        "hashed_bytes": hashed_bytes,
        "requested_parity_indices": parity_indices,
    }
    return median_us, mad_us, parity_identity


def validate_execution_inputs(job):
    executable = Path(job["executable"])
    executable_artifact = job.get("executable_artifact")
    if executable_artifact is not None:
        validate_frozen_executable(
            executable_artifact, build_metadata=job["build_metadata"],
            source_state=job["source_identity"]
        )
        if (executable.resolve() !=
                Path(executable_artifact["executable"]).resolve() or
                job.get("executable_sha256") !=
                executable_artifact["executable_sha256"]):
            raise CrossoverError(
                "job executable differs from its frozen artifact"
            )
        return
    try:
        executable_sha256 = digest_bytes(executable.read_bytes())
    except OSError as error:
        raise CrossoverError("cannot re-read executable {}: {}".format(
            executable, error
        ))
    if executable_sha256 != job["executable_sha256"]:
        raise CrossoverError("benchmark executable changed during the run")
    if canonical_bytes(cmake_build_metadata(executable)) != canonical_bytes(
            job["build_metadata"]):
        raise CrossoverError("benchmark CMake metadata changed during the run")


def artifact_path(root, relative, description):
    if not isinstance(relative, str) or not relative:
        raise CrossoverError("{} path is missing".format(description))
    root = root.resolve()
    candidate = (root / relative).resolve()
    try:
        common = os.path.commonpath((str(root), str(candidate)))
    except ValueError:
        common = ""
    if common != str(root):
        raise CrossoverError("{} escapes the result directory".format(description))
    return candidate


def validate_job_artifacts(result_dir, result, expected_job, settings):
    if not isinstance(result, dict) or result.get("schema") != JOB_SCHEMA:
        raise CrossoverError("job result has a legacy or unknown schema")
    common_keys = {
        "build_metadata", "cell", "commands", "configuration_id",
        "executable", "executable_artifact", "executable_sha256", "job_id",
        "isolation", "measurements", "reason", "resumed", "schema", "seed",
        "source_identity", "status",
    }
    if result.get("status") != "passed":
        if set(result) != common_keys or result.get("status") != "failed":
            raise CrossoverError("failed job result has an invalid v2 shape")
        return
    if set(result) != common_keys | {"parity_identity", "summary"}:
        raise CrossoverError("passed job result has an invalid v2 shape")
    if result.get("reason") != "" or not isinstance(result.get("resumed"), bool):
        raise CrossoverError("passed job contains inconsistent status metadata")
    if settings.get("mode") in AUTHORITATIVE_COMMANDS:
        isolation_settings = settings["isolation"]
        try:
            support = load_isolation_support(
                expected_job["build_metadata"]["entries"][
                    "CMAKE_HOME_DIRECTORY"
                ]
            )
            support.validate_isolation(
                result.get("isolation"),
                isolation_settings["benchmark_cpu"],
                isolation_settings["reserved_sibling"],
                require_accepted=True,
            )
            support.validate_pair_lease_identity(
                isolation_settings["pair_lease"],
                isolation_settings["benchmark_cpu"],
                isolation_settings["reserved_sibling"],
            )
            if canonical_bytes(
                    result["isolation"]["pair_lease"]) != canonical_bytes(
                    isolation_settings["pair_lease"]):
                raise CrossoverError(
                    "job isolation pair lease differs from campaign settings"
                )
        except Exception as error:
            raise CrossoverError(
                "passed job isolation evidence is invalid: {}".format(error)
            )
    elif result.get("isolation") is not None:
        raise CrossoverError("screening job unexpectedly claims isolation")
    for key in (
            "build_metadata", "cell", "configuration_id", "executable",
            "executable_artifact", "executable_sha256", "job_id", "seed",
            "source_identity"):
        if canonical_bytes(result.get(key)) != canonical_bytes(expected_job.get(key)):
            raise CrossoverError(
                "passed job {} has stale {} metadata".format(
                    expected_job["job_id"], key
                )
            )
    validate_execution_inputs(expected_job)
    order = expected_job.get("invocation_order", [])
    commands = result.get("commands")
    measurements = result.get("measurements")
    if (not order or not isinstance(commands, list) or
            not isinstance(measurements, list) or
            len(commands) != len(order) or len(measurements) != len(order)):
        raise CrossoverError(
            "passed job {} has an incomplete command/measurement sequence".format(
                expected_job["job_id"]
            )
        )
    parity_identities = set()
    log_root = result_dir / "logs" / expected_job["job_id"]
    try:
        source = Path(
            expected_job["build_metadata"]["entries"]["CMAKE_HOME_DIRECTORY"]
        ).resolve()
    except (KeyError, TypeError, ValueError):
        raise CrossoverError("job build provenance omits the source root")
    command_keys = {
        "argv", "cwd", "environment", "returncode", "stderr_log",
        "stderr_sha256", "stdout_log", "stdout_sha256", "timed_out",
    }
    measurement_keys = {
        "benchmark_json", "benchmark_json_sha256", "mad_us", "median_us",
        "parity_identity", "sequence_index", "timed_mode",
    }
    for index, timed_mode in enumerate(order):
        command = commands[index]
        measurement = measurements[index]
        if (not isinstance(command, dict) or set(command) != command_keys or
                not isinstance(measurement, dict) or
                set(measurement) != measurement_keys):
            raise CrossoverError(
                "passed job {} has a malformed sequence entry {}".format(
                    expected_job["job_id"], index
                )
            )
        label = "{:02d}-{}".format(index, timed_mode)
        raw_path = (
            result_dir / "raw" / expected_job["job_id"] / (label + ".json")
        ).resolve()
        expected_relative_raw = str(raw_path.relative_to(result_dir.resolve()))
        expected_environment = dict(BENCHMARK_ENVIRONMENT)
        expected_environment["LEO2_EXPECT_BACKEND"] = (
            expected_job["cell"]["backend"]
        )
        if (type(command.get("returncode")) is not int or
                command.get("returncode") != 0 or
                command.get("timed_out") is not False or
                type(measurement.get("sequence_index")) is not int or
                measurement.get("sequence_index") != index or
                measurement.get("timed_mode") != timed_mode or
                command.get("argv") != benchmark_argv(
                    expected_job, timed_mode, raw_path, settings
                ) or command.get("cwd") != str(source) or
                command.get("environment") != expected_environment or
                command.get("stdout_log") != label + ".stdout.log" or
                command.get("stderr_log") != label + ".stderr.log" or
                measurement.get("benchmark_json") != expected_relative_raw):
            raise CrossoverError(
                "passed job {} has an invalid sequence entry {}".format(
                    expected_job["job_id"], index
                )
            )
        for stream in ("stdout", "stderr"):
            name = command.get(stream + "_log")
            if not isinstance(name, str) or Path(name).name != name:
                raise CrossoverError("job log name is missing or unsafe")
            path = artifact_path(log_root, name, stream + " log")
            try:
                value = path.read_bytes()
            except OSError as error:
                raise CrossoverError("cannot read {}: {}".format(path, error))
            if digest_bytes(value) != command.get(stream + "_sha256"):
                raise CrossoverError("{} hash does not match the job record".format(path))
        raw_path = artifact_path(
            result_dir, measurement["benchmark_json"], "benchmark JSON"
        )
        try:
            raw_bytes = raw_path.read_bytes()
        except OSError as error:
            raise CrossoverError("cannot read {}: {}".format(raw_path, error))
        raw = decode_json_bytes(raw_bytes, str(raw_path))
        if digest_bytes(raw_bytes) != measurement.get("benchmark_json_sha256"):
            raise CrossoverError("{} hash does not match the job record".format(raw_path))
        median_us, mad_us, parity_identity = validate_raw(
            raw, expected_job, timed_mode, settings
        )
        recorded_median = required_finite_metric(
            measurement.get("median_us"),
            "job.measurements[{}].median_us".format(index), False
        )
        recorded_mad = required_finite_metric(
            measurement.get("mad_us"),
            "job.measurements[{}].mad_us".format(index), True
        )
        if recorded_median != median_us or recorded_mad != mad_us:
            raise CrossoverError("raw metrics differ from the passed job summary")
        if canonical_bytes(measurement.get("parity_identity")) != canonical_bytes(
                parity_identity):
            raise CrossoverError(
                "raw parity identity differs from the passed job summary"
            )
        parity_identities.add(canonical_bytes(parity_identity))
    if len(parity_identities) != 1:
        raise CrossoverError(
            "passed job raw artifacts contain different parity identities"
        )
    parity_identity = json.loads(next(iter(parity_identities)).decode("utf-8"))
    if canonical_bytes(result.get("parity_identity")) != canonical_bytes(
            parity_identity):
        raise CrossoverError("passed job omits its non-vacuous parity identity")
    recomputed = summarize_measurements(measurements)
    if canonical_bytes(recomputed) != canonical_bytes(result.get("summary")):
        raise CrossoverError("passed job aggregate does not match its raw measurements")


def median(values):
    try:
        result = float(statistics.median(values))
    except (OverflowError, TypeError, ValueError, statistics.StatisticsError) as error:
        raise CrossoverError("cannot compute finite median: {}".format(error))
    if not math.isfinite(result):
        raise CrossoverError("cannot compute finite median")
    return result


def relative_gain_percent(candidate, reference):
    try:
        result = (reference / candidate - 1.0) * 100.0
    except (OverflowError, ZeroDivisionError):
        result = float("nan")
    if not math.isfinite(result):
        raise CrossoverError("cannot compute a finite encoder gain")
    return result


def paired_log_inference(rounds):
    if len(rounds) != 3:
        raise CrossoverError(
            "authoritative ABBA inference requires exactly three rounds"
        )
    contrasts = []
    for index, round_value in enumerate(rounds):
        round_value = require_exact_keys(round_value, {
            "direct_geometric_mean_us", "gain_percent", "log_contrast",
            "transform_geometric_mean_us",
        }, "rounds[{}]".format(index))
        direct = required_finite_metric(
            round_value.get("direct_geometric_mean_us"),
            "rounds[{}].direct_geometric_mean_us".format(index), False
        )
        transform = required_finite_metric(
            round_value.get("transform_geometric_mean_us"),
            "rounds[{}].transform_geometric_mean_us".format(index), False
        )
        contrast = required_finite_number(
            round_value.get("log_contrast"),
            "rounds[{}].log_contrast".format(index)
        )
        expected_contrast = math.log(transform) - math.log(direct)
        recorded_gain = required_finite_number(
            round_value.get("gain_percent"),
            "rounds[{}].gain_percent".format(index)
        )
        expected_gain = relative_gain_percent(direct, transform)
        if not math.isclose(
                contrast, expected_contrast, rel_tol=1e-15, abs_tol=1e-15):
            raise CrossoverError(
                "round {} log contrast does not match its geometric means".format(
                    index
                )
            )
        if not math.isclose(
                recorded_gain, expected_gain, rel_tol=1e-12, abs_tol=1e-12):
            raise CrossoverError(
                "round {} gain does not match its geometric means".format(index)
            )
        contrasts.append(contrast)
    mean = statistics.mean(contrasts)
    standard_error = statistics.stdev(contrasts) / math.sqrt(3.0)
    # Two-sided 95% Student-t critical value at df=2.
    margin = 4.302652729911275 * standard_error
    try:
        lower = math.exp(mean - margin)
        upper = math.exp(mean + margin)
        speedup = math.exp(mean)
    except OverflowError:
        raise CrossoverError(
            "paired-log confidence interval overflowed finite arithmetic"
        )
    if not all(math.isfinite(value) and value > 0 for value in (
            lower, speedup, upper)):
        raise CrossoverError("paired-log confidence interval is not finite")
    return {
        "confidence": 0.95,
        "degrees_of_freedom": 2,
        "gain_percent": (speedup - 1.0) * 100.0,
        "gain_percent_student_t_interval": [
            (lower - 1.0) * 100.0,
            (upper - 1.0) * 100.0,
        ],
        "log_contrasts": contrasts,
        "speedup_geometric_mean": speedup,
        "speedup_student_t_interval": [lower, upper],
    }


def summarize_measurements(measurements):
    by_mode = {"direct": [], "transform": []}
    if not isinstance(measurements, list):
        raise CrossoverError("job measurements must be a JSON array")
    for index, measurement in enumerate(measurements):
        measurement = required_mapping(
            measurement, "job.measurements[{}]".format(index)
        )
        timed_mode = measurement.get("timed_mode")
        if timed_mode not in by_mode:
            raise CrossoverError("job measurement has an unknown timed mode")
        by_mode[timed_mode].append(required_finite_metric(
            measurement.get("median_us"),
            "job.measurements[{}].median_us".format(index), False
        ))
    if not by_mode["direct"] or not by_mode["transform"]:
        raise CrossoverError("job did not measure both encoder paths")
    direct = median(by_mode["direct"])
    transform = median(by_mode["transform"])
    gain = relative_gain_percent(direct, transform)
    rounds = []
    if len(measurements) % 4 == 0 and len(measurements) >= 4:
        for offset in range(0, len(measurements), 4):
            group = measurements[offset:offset + 4]
            if [item["timed_mode"] for item in group] != [
                    "direct", "transform", "transform", "direct"]:
                rounds = []
                break
            direct_logs = [
                math.log(required_finite_metric(
                    group[position]["median_us"],
                    "ABBA direct invocation", False
                )) for position in (0, 3)
            ]
            transform_logs = [
                math.log(required_finite_metric(
                    group[position]["median_us"],
                    "ABBA transform invocation", False
                )) for position in (1, 2)
            ]
            round_direct = math.exp(statistics.mean(direct_logs))
            round_transform = math.exp(statistics.mean(transform_logs))
            log_contrast = (
                statistics.mean(transform_logs) -
                statistics.mean(direct_logs)
            )
            rounds.append({
                "direct_geometric_mean_us": round_direct,
                "gain_percent": relative_gain_percent(
                    round_direct, round_transform
                ),
                "log_contrast": log_contrast,
                "transform_geometric_mean_us": round_transform,
            })
    inference = paired_log_inference(rounds) if len(rounds) == 3 else None
    return {
        "direct_invocation_medians_us": by_mode["direct"],
        "direct_median_us": direct,
        "gain_percent": gain,
        "paired_log_inference": inference,
        "rounds": rounds,
        "transform_invocation_medians_us": by_mode["transform"],
        "transform_median_us": transform,
    }


def run_job(job, context):
    isolation_context = context.get("isolation")

    def validate_guards():
        if isolation_context is not None:
            validate_authoritative_guards(
                isolation_context["canonical_guard"],
                isolation_context["pair_guard"],
            )

    result_path = context["result_dir"] / "jobs" / (job["job_id"] + ".json")
    if context["resume"] and result_path.is_file():
        validate_guards()
        try:
            previous = decode_json_bytes(
                result_path.read_bytes(), str(result_path)
            )
            if (not isinstance(previous, dict) or
                    previous.get("schema") != JOB_SCHEMA):
                raise CrossoverError("resume job has a legacy or unknown schema")
            if (previous.get("configuration_id") == job["configuration_id"] and
                    previous.get("status") == "passed"):
                validate_job_artifacts(
                    context["result_dir"], previous, job, context["settings"]
                )
                validate_guards()
                return previous
        except (CrossoverError, OSError):
            validate_guards()

    log_dir = context["result_dir"] / "logs" / job["job_id"]
    raw_dir = context["result_dir"] / "raw" / job["job_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "build_metadata": job["build_metadata"],
        "cell": job["cell"],
        "commands": [],
        "configuration_id": job["configuration_id"],
        "executable": job["executable"],
        "executable_artifact": job["executable_artifact"],
        "executable_sha256": job["executable_sha256"],
        "job_id": job["job_id"],
        "isolation": None,
        "measurements": [],
        "reason": "",
        "resumed": False,
        "schema": JOB_SCHEMA,
        "seed": job["seed"],
        "source_identity": job["source_identity"],
        "status": "failed",
    }
    environment = dict(BENCHMARK_ENVIRONMENT)
    environment["LEO2_EXPECT_BACKEND"] = job["cell"]["backend"]
    parity_identities = set()
    isolation_before = None
    if isolation_context is not None:
        support = isolation_context["support"]
        validate_guards()
        isolation_before = {
            "monotonic_ns": time.monotonic_ns(),
            "benchmark_cpu": support.cpu_stat_snapshot(
                isolation_context["cpu"]
            ),
            "reserved_sibling": support.cpu_stat_snapshot(
                isolation_context["sibling"]
            ),
        }
    try:
        for index, timed_mode in enumerate(job["invocation_order"]):
            validate_guards()
            validate_execution_inputs(job)
            label = "{:02d}-{}".format(index, timed_mode)
            raw_path = raw_dir / (label + ".json")
            stdout_path = log_dir / (label + ".stdout.log")
            stderr_path = log_dir / (label + ".stderr.log")
            try:
                raw_path.unlink()
            except FileNotFoundError:
                pass
            argv = benchmark_argv(job, timed_mode, raw_path, context["settings"])
            command = run_command(
                argv, context["source"], stdout_path, stderr_path,
                context["timeout"], environment
            )
            command["environment"] = dict(environment)
            result["commands"].append(command)
            validate_guards()
            validate_execution_inputs(job)
            if command["returncode"] != 0:
                raise CrossoverError(
                    "{} exited with status {}".format(label, command["returncode"])
                )
            try:
                raw_bytes = raw_path.read_bytes()
            except OSError as error:
                raise CrossoverError("cannot read {}: {}".format(raw_path, error))
            raw = decode_json_bytes(raw_bytes, str(raw_path))
            median_us, mad_us, parity_identity = validate_raw(
                raw, job, timed_mode, context["settings"]
            )
            parity_identities.add(canonical_bytes(parity_identity))
            result["measurements"].append({
                "benchmark_json": str(raw_path.relative_to(context["result_dir"])),
                "benchmark_json_sha256": digest_bytes(raw_bytes),
                "mad_us": mad_us,
                "median_us": median_us,
                "parity_identity": parity_identity,
                "sequence_index": index,
                "timed_mode": timed_mode,
            })
        if len(parity_identities) != 1:
            raise CrossoverError("forced paths emitted different parity identities")
        result["parity_identity"] = json.loads(
            next(iter(parity_identities)).decode("utf-8")
        )
        result["summary"] = summarize_measurements(result["measurements"])
        result["status"] = "passed"
    except (CrossoverError, OSError, ValueError) as error:
        result["reason"] = str(error)
    finally:
        if isolation_context is not None and isolation_before is not None:
            try:
                support = isolation_context["support"]
                after_cpu = support.cpu_stat_snapshot(
                    isolation_context["cpu"]
                )
                after_sibling = support.cpu_stat_snapshot(
                    isolation_context["sibling"]
                )
                after_monotonic = time.monotonic_ns()
                validate_guards()
                result["isolation"] = support.isolation_record(
                    isolation_context["cpu"],
                    isolation_context["sibling"],
                    isolation_context["pair_lease"],
                    isolation_before["monotonic_ns"],
                    after_monotonic,
                    isolation_before["benchmark_cpu"],
                    after_cpu,
                    isolation_before["reserved_sibling"],
                    after_sibling,
                )
                support.validate_isolation(
                    result["isolation"],
                    isolation_context["cpu"],
                    isolation_context["sibling"],
                    require_accepted=True,
                )
            except Exception as error:
                result["status"] = "failed"
                result["reason"] = (
                    result["reason"] + "; " if result["reason"] else ""
                ) + "isolation validation failed: {}".format(error)
                result.pop("parity_identity", None)
                result.pop("summary", None)
    atomic_write_json(result_path, result)
    return result


def summarize_region(results, promotion_percent):
    passed = [item for item in results if item.get("status") == "passed"]
    gains = []
    round_gains = []
    confidence_lower_gains = []
    confidence_missing_count = 0
    for index, item in enumerate(passed):
        summary = required_mapping(
            item.get("summary"), "jobs[{}].summary".format(index)
        )
        gains.append(required_finite_number(
            summary.get("gain_percent"),
            "jobs[{}].summary.gain_percent".format(index)
        ))
        rounds = summary.get("rounds")
        if not isinstance(rounds, list):
            raise CrossoverError("passed job summary rounds must be a JSON array")
        for round_index, round_value in enumerate(rounds):
            round_value = required_mapping(
                round_value,
                "jobs[{}].summary.rounds[{}]".format(index, round_index)
            )
            round_gains.append(required_finite_number(
                round_value.get("gain_percent"),
                "jobs[{}].summary.rounds[{}].gain_percent".format(
                    index, round_index
                )
            ))
        inference = summary.get("paired_log_inference")
        if inference is None:
            confidence_missing_count += 1
        else:
            inference = required_mapping(
                inference,
                "jobs[{}].summary.paired_log_inference".format(index)
            )
            interval = inference.get("gain_percent_student_t_interval")
            if (not isinstance(interval, list) or len(interval) != 2):
                raise CrossoverError("passed job has an invalid confidence interval")
            lower = required_finite_number(
                interval[0],
                "jobs[{}].summary.paired_log_inference.lower".format(index)
            )
            upper = required_finite_number(
                interval[1],
                "jobs[{}].summary.paired_log_inference.upper".format(index)
            )
            if lower > upper:
                raise CrossoverError("passed job confidence interval is reversed")
            confidence_lower_gains.append(lower)
    gains.sort()
    round_gains.sort()
    return {
        "cell_count": len(results),
        "confidence_interval_missing_count": confidence_missing_count,
        "confidence_promotion_count": sum(
            value >= promotion_percent for value in confidence_lower_gains
        ),
        "confidence_worst_lower_gain_percent": (
            min(confidence_lower_gains) if confidence_lower_gains else None
        ),
        "failed_count": len(results) - len(passed),
        "gain_max_percent": max(gains) if gains else None,
        "gain_median_percent": median(gains) if gains else None,
        "gain_min_percent": min(gains) if gains else None,
        "improvement_count": sum(value > 0 for value in gains),
        "promotion_count": sum(value >= promotion_percent for value in gains),
        "regression_count": sum(value < 0 for value in gains),
        "round_gain_median_percent": median(round_gains) if round_gains else None,
        "round_gain_min_percent": min(round_gains) if round_gains else None,
        "round_regression_count": sum(value < 0 for value in round_gains),
        "severe_regression_count": sum(value < -2.0 for value in gains),
    }


def analyze_results(
        results, promotion_percent, manifest_configuration_id=None,
        run_status=None, source_changed_during_run=None,
        execution_input_error=None):
    ordered = sorted(results, key=lambda item: item.get("job_id", ""))
    regions = {}
    for item in ordered:
        region = item.get("cell", {}).get("region", "unknown")
        regions.setdefault(region, []).append(item)
    candidate = []
    excluded = []
    for region, values in regions.items():
        if region.startswith("candidate"):
            candidate.extend(values)
        else:
            excluded.extend(values)
    candidate_summary = summarize_region(candidate, promotion_percent)
    analysis = {
        "candidate": candidate_summary,
        "excluded_neighbors": summarize_region(excluded, promotion_percent),
        "execution_input_error": execution_input_error,
        "jobs_failed": sum(item.get("status") != "passed" for item in ordered),
        "jobs_passed": sum(item.get("status") == "passed" for item in ordered),
        "jobs_total": len(ordered),
        "manifest_configuration_id": manifest_configuration_id,
        "promotion_percent": promotion_percent,
        "regions": {
            region: summarize_region(values, promotion_percent)
            for region, values in sorted(regions.items())
        },
        "run_status": run_status,
        "schema": ANALYSIS_SCHEMA,
        "source_changed_during_run": source_changed_during_run,
    }
    analysis["candidate"]["all_cells_meet_promotion_threshold"] = bool(candidate) and (
        candidate_summary["failed_count"] == 0 and
        candidate_summary["promotion_count"] == candidate_summary["cell_count"]
    )
    analysis["candidate"][
        "all_cells_confidently_meet_promotion_threshold"
    ] = bool(candidate) and (
        candidate_summary["failed_count"] == 0 and
        candidate_summary["confidence_interval_missing_count"] == 0 and
        candidate_summary["confidence_promotion_count"] ==
        candidate_summary["cell_count"]
    )
    return analysis


def manifest_identity(manifest):
    return {
        "evidence_contract": manifest["evidence_contract"],
        "executables": manifest["executables"],
        "jobs": manifest["jobs"],
        "machine": manifest["machine"],
        "settings": manifest["settings"],
        "source_fingerprint": manifest["source_fingerprint"],
    }


def validate_cell_value(value, description):
    expected_keys = {
        "backend", "field", "k", "mask", "profile", "q", "r", "region",
        "shard_bytes",
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        raise CrossoverError("{} has an invalid cell shape".format(description))
    if (value["backend"] not in KNOWN_BACKENDS or
            value["field"] not in ("gf8", "gf16") or
            value["profile"] not in ("low", "high") or
            not isinstance(value["region"], str) or not value["region"] or
            not isinstance(value["mask"], str)):
        raise CrossoverError("{} has invalid cell labels".format(description))
    for key in ("k", "r", "q", "shard_bytes"):
        if type(value[key]) is not int or value[key] <= 0:
            raise CrossoverError(
                "{} cell {} must be a positive integer".format(description, key)
            )
    if value["k"] > 16 or value["r"] > 16 or value["q"] > value["r"]:
        raise CrossoverError("{} is outside the bounded direct domain".format(
            description
        ))
    if len(requested_indices(value["mask"], value["r"])) != value["q"]:
        raise CrossoverError("{} mask cardinality differs from Q".format(description))


def validate_manifest_settings(settings, jobs, path):
    core_keys = {
        "abba_rounds", "benchmark", "frozen_executable_required", "mode",
        "pin_cpu", "placement_policy", "taskset", "taskset_sha256",
        "timeout_seconds_per_invocation", "workers",
    }
    mode = settings.get("mode")
    expected_keys = set(core_keys)
    if mode in AUTHORITATIVE_COMMANDS:
        expected_keys.add("isolation")
    if mode == "historical-avx2":
        expected_keys.update(("campaign", "controlled_build"))
    if set(settings) != expected_keys:
        raise CrossoverError("{} has invalid v2 settings fields".format(path))
    benchmark = settings.get("benchmark")
    if (not isinstance(benchmark, dict) or
            set(benchmark) != {"batch", "iterations", "reuse", "warmups"}):
        raise CrossoverError("{} has invalid benchmark settings".format(path))
    for key in ("batch", "iterations", "reuse", "warmups"):
        if type(benchmark[key]) is not int or benchmark[key] <= 0:
            raise CrossoverError("{} benchmark {} is invalid".format(path, key))
    for key in ("abba_rounds", "timeout_seconds_per_invocation", "workers"):
        if type(settings.get(key)) is not int or settings[key] < 0:
            raise CrossoverError("{} setting {} is invalid".format(path, key))
    if settings["timeout_seconds_per_invocation"] == 0 or settings["workers"] == 0:
        raise CrossoverError("{} has a zero execution setting".format(path))
    if not isinstance(settings.get("placement_policy"), str):
        raise CrossoverError("{} has invalid placement policy".format(path))
    if mode not in ("screen",) + AUTHORITATIVE_COMMANDS:
        raise CrossoverError("{} has an unknown runner mode".format(path))
    frozen = settings.get("frozen_executable_required")
    if type(frozen) is not bool or frozen != (mode in AUTHORITATIVE_COMMANDS):
        raise CrossoverError("{} has inconsistent frozen policy".format(path))
    if mode in AUTHORITATIVE_COMMANDS:
        if (settings["abba_rounds"] != 3 or settings["workers"] != 1 or
                type(settings.get("pin_cpu")) is not int or
                not isinstance(settings.get("taskset"), str) or
                Path(settings["taskset"]) != Path("/usr/bin/taskset") or
                not isinstance(settings.get("taskset_sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", settings["taskset_sha256"])):
            raise CrossoverError("{} has invalid authoritative ABBA settings".format(
                path
            ))
        try:
            taskset_sha256 = digest_bytes(Path("/usr/bin/taskset").read_bytes())
        except OSError as error:
            raise CrossoverError(
                "{} cannot revalidate /usr/bin/taskset: {}".format(path, error)
            )
        if taskset_sha256 != settings["taskset_sha256"]:
            raise CrossoverError("{} taskset executable changed".format(path))
        isolation = settings.get("isolation")
        isolation_keys = {
            "allowed_cpu_set_at_launch", "benchmark_cpu", "canonical_lock",
            "child_environment", "housekeeping_cpu_set", "pair_lease",
            "reserved_sibling",
        }
        if not isinstance(isolation, dict) or set(isolation) != isolation_keys:
            raise CrossoverError(
                "{} has incomplete authoritative isolation settings".format(path)
            )
        cpu = isolation.get("benchmark_cpu")
        sibling = isolation.get("reserved_sibling")
        allowed = isolation.get("allowed_cpu_set_at_launch")
        housekeeping = isolation.get("housekeeping_cpu_set")
        if (type(cpu) is not int or type(sibling) is not int or cpu == sibling or
                cpu != settings["pin_cpu"] or
                not isinstance(allowed, list) or
                not isinstance(housekeeping, list) or
                any(type(value) is not int or value < 0
                    for value in allowed + housekeeping) or
                allowed != sorted(set(allowed)) or
                housekeeping != sorted(set(housekeeping)) or
                set(housekeeping) != set(allowed) - {cpu, sibling} or
                not housekeeping):
            raise CrossoverError(
                "{} has invalid authoritative CPU-pair settings".format(path)
            )
        if isolation.get("child_environment") != BENCHMARK_ENVIRONMENT:
            raise CrossoverError(
                "{} has an unsanitized benchmark environment".format(path)
            )
        lock = isolation.get("canonical_lock")
        if (not isinstance(lock, dict) or set(lock) != {
                "device", "inode", "lock", "mode", "path", "uid"} or
                lock.get("path") != str(AUTHORITATIVE_LOCK) or
                lock.get("lock") != "exclusive" or lock.get("mode") != 0o600 or
                lock.get("uid") != os.getuid() or
                any(type(lock.get(key)) is not int or lock[key] < 0
                    for key in ("device", "inode"))):
            raise CrossoverError(
                "{} has invalid canonical-lock identity".format(path)
            )
        pair = isolation.get("pair_lease")
        payload = pair.get("payload") if isinstance(pair, dict) else None
        if (not isinstance(pair, dict) or set(pair) != {
                "device", "directory_device", "directory_inode", "inode",
                "lock", "path", "payload", "sha256"} or
                not isinstance(payload, dict) or
                payload.get("cpus") != sorted((cpu, sibling)) or
                payload.get("schema") != "leopard2-cpu-pair-lease/v1"):
            raise CrossoverError(
                "{} has invalid physical-pair lease identity".format(path)
            )
        try:
            source_root = jobs[0]["build_metadata"]["entries"][
                "CMAKE_HOME_DIRECTORY"
            ]
            support = load_isolation_support(source_root)
            support.validate_pair_lease_identity(pair, cpu, sibling)
        except Exception as error:
            raise CrossoverError(
                "{} physical-pair lease failed canonical validation: {}".format(
                    path, error
                )
            )
    elif (settings["abba_rounds"] != 0 or settings.get("pin_cpu") is not None or
          settings.get("taskset") is not None or
          settings.get("taskset_sha256") is not None):
        raise CrossoverError("{} has invalid screening settings".format(path))
    if mode == "historical-avx2":
        cells = [job.get("cell") for job in jobs]
        expected_cells = historical_avx2_grid()
        expected_campaign = {
            "cell_count": len(expected_cells),
            "cells_sha256": digest_value(expected_cells),
            "name": "corrected-high-profile-explicit-avx2",
        }
        if (canonical_bytes(cells) != canonical_bytes(expected_cells) or
                canonical_bytes(settings.get("campaign")) !=
                canonical_bytes(expected_campaign) or benchmark != {
                    "batch": 1, "iterations": 15, "reuse": 64, "warmups": 4,
                }):
            raise CrossoverError(
                "{} is not the exact historical AVX2 campaign".format(path)
            )
        controlled = settings.get("controlled_build")
        if (not isinstance(controlled, dict) or set(controlled) != {
                "path", "schema", "sha256"} or
                controlled.get("schema") != CONTROLLED_BUILD_SCHEMA or
                controlled.get("path") != "controlled-build.json" or
                not isinstance(controlled.get("sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", controlled["sha256"])):
            raise CrossoverError(
                "{} has invalid controlled-build descriptor".format(path)
            )


def validate_controlled_build(settings, jobs, manifest_path):
    if settings.get("mode") != "historical-avx2":
        return
    descriptor = settings["controlled_build"]
    result_dir = Path(manifest_path).resolve().parent
    record_path = artifact_path(
        result_dir, descriptor["path"], "controlled build record"
    )
    try:
        record_bytes = record_path.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot read controlled build record {}: {}".format(
                record_path, error
            )
        )
    if digest_bytes(record_bytes) != descriptor["sha256"]:
        raise CrossoverError("controlled build record hash is stale")
    record = decode_json_bytes(record_bytes, str(record_path))
    expected_keys = {
        "backend", "build_metadata", "build_root", "commands", "executable",
        "executable_sha256", "guard_identity", "schema", "source_identity",
    }
    if (not isinstance(record, dict) or set(record) != expected_keys or
            record.get("schema") != CONTROLLED_BUILD_SCHEMA or
            record.get("backend") != "avx2"):
        raise CrossoverError("controlled build record has an invalid schema")
    build_root = (result_dir / "controlled-build-avx2").resolve()
    executable = (build_root / "bench_leopard2_direct_encode").resolve()
    if (record.get("build_root") != str(build_root) or
            record.get("executable") != str(executable)):
        raise CrossoverError("controlled build paths differ from the campaign")
    try:
        executable_bytes = executable.read_bytes()
    except OSError as error:
        raise CrossoverError(
            "cannot read controlled build executable: {}".format(error)
        )
    if (digest_bytes(executable_bytes) != record.get("executable_sha256") or
            not os.access(str(executable), os.X_OK)):
        raise CrossoverError("controlled build executable changed")
    representative = jobs[0]
    if (canonical_bytes(record.get("build_metadata")) != canonical_bytes(
            representative.get("build_metadata")) or
            canonical_bytes(record.get("source_identity")) != canonical_bytes(
                representative.get("source_identity"))):
        raise CrossoverError(
            "controlled build identity differs from manifest jobs"
        )
    expected_guard_identity = {
        "canonical_lock": settings["isolation"]["canonical_lock"],
        "pair_lease": settings["isolation"]["pair_lease"],
    }
    if canonical_bytes(record.get("guard_identity")) != canonical_bytes(
            expected_guard_identity):
        raise CrossoverError(
            "controlled build guard identity differs from the campaign"
        )
    commands = record.get("commands")
    if not isinstance(commands, list) or len(commands) != 2:
        raise CrossoverError("controlled build command sequence is incomplete")
    command_keys = {
        "argv", "cwd", "environment", "returncode", "stderr_log",
        "stderr_sha256", "stdout_log", "stdout_sha256", "timed_out",
    }
    if any(not isinstance(command, dict) or set(command) != command_keys or
           not isinstance(command.get("argv"), list) or
           not command["argv"] or
           any(not isinstance(item, str) or not item
               for item in command["argv"])
           for command in commands):
        raise CrossoverError("controlled build command shape is invalid")
    try:
        source = Path(
            record["build_metadata"]["entries"]["CMAKE_HOME_DIRECTORY"]
        ).resolve()
    except (KeyError, TypeError, ValueError):
        raise CrossoverError("controlled build omits its source directory")
    current_metadata = cmake_build_metadata(executable)
    if canonical_bytes(current_metadata) != canonical_bytes(
            record["build_metadata"]):
        raise CrossoverError(
            "controlled build CMake/object/link provenance changed"
        )
    validate_build_source_binding(
        current_metadata, source, record["source_identity"], "avx2",
        require_fresh=True
    )
    cmake = commands[0]["argv"][0]
    if (not Path(cmake).is_absolute() or Path(cmake).name != "cmake" or
            commands[1]["argv"][0] != cmake):
        raise CrossoverError("controlled build did not use one absolute cmake")
    expected_argv = (
        [
            cmake, "-S", str(source), "-B", str(build_root), "-G", "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "-DLEO2_BACKEND_VARIANT=avx2",
            "-DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git",
            "-DLEO2_BUILD_BENCHMARKS=ON",
            "-DLEO2_BUILD_TESTS=ON",
            "-DLEO2_BUILD_FUZZERS=OFF",
            "-DLEO2_ENABLE_CUDA=OFF",
            "-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF",
            "-DLEOPARD_ENABLE_GF8=ON",
            "-DLEOPARD_ENABLE_GF16=ON",
        ],
        [
            cmake, "--build", str(build_root), "--target",
            "bench_leopard2_direct_encode", "--parallel",
            str(min(len(settings["isolation"]["housekeeping_cpu_set"]), 128)),
        ],
    )
    for index, command in enumerate(commands):
        if (command.get("argv") != expected_argv[index] or
                command.get("cwd") != str(source) or
                command.get("environment") != GIT_ENVIRONMENT or
                type(command.get("returncode")) is not int or
                command.get("returncode") != 0 or
                command.get("timed_out") is not False):
            raise CrossoverError(
                "controlled build command {} is invalid".format(index)
            )
        for stream in ("stdout", "stderr"):
            log_path = artifact_path(
                result_dir, command[stream + "_log"],
                "controlled build {} log".format(stream)
            )
            try:
                value = log_path.read_bytes()
            except OSError as error:
                raise CrossoverError(
                    "cannot read controlled build log: {}".format(error)
                )
            if digest_bytes(value) != command[stream + "_sha256"]:
                raise CrossoverError("controlled build log hash changed")


def validate_manifest(manifest, path):
    expected_keys = {
        "configuration_id", "evidence_contract", "executables", "jobs",
        "machine", "schema", "settings", "source_fingerprint",
    }
    if (not isinstance(manifest, dict) or set(manifest) != expected_keys or
            manifest.get("schema") != SCHEMA):
        raise CrossoverError(
            "{} has an unknown, legacy, or incomplete schema".format(path)
        )
    settings = manifest.get("settings")
    jobs = manifest.get("jobs")
    executables = manifest.get("executables")
    source_state = manifest.get("source_fingerprint")
    if (not isinstance(settings, dict) or not isinstance(jobs, list) or
            not jobs or not isinstance(executables, dict) or not executables or
            not isinstance(source_state, dict)):
        raise CrossoverError("{} omits required v2 manifest fields".format(path))
    expected_job_keys = {
        "build_metadata", "cell", "configuration_id", "executable",
        "executable_artifact", "executable_sha256", "invocation_order",
        "job_id", "seed", "source_identity",
    }
    if any(not isinstance(job, dict) or set(job) != expected_job_keys
           for job in jobs):
        raise CrossoverError("{} contains an invalid v2 job".format(path))
    validate_manifest_settings(settings, jobs, path)
    validate_controlled_build(settings, jobs, path)
    contract = evidence_contract(settings.get("frozen_executable_required") is True)
    if canonical_bytes(manifest.get("evidence_contract")) != canonical_bytes(contract):
        raise CrossoverError("{} has a stale or relabeled evidence contract".format(path))
    validate_source_state(
        source_state, str(path) + " source_fingerprint",
        settings.get("mode") in AUTHORITATIVE_COMMANDS
    )
    if manifest.get("configuration_id") != digest_value(manifest_identity(manifest)):
        raise CrossoverError("{} configuration ID does not match its content".format(path))

    seen = set()
    expected_executables = {}
    for job in jobs:
        if not isinstance(job, dict) or set(job) != expected_job_keys:
            raise CrossoverError("{} contains an invalid v2 job".format(path))
        validate_cell_value(job.get("cell"), "manifest job")
        identity = job_identity(
            job["cell"], job["executable"], job["executable_artifact"],
            job["executable_sha256"], job["build_metadata"],
            job["source_identity"],
            manifest["machine"], settings
        )
        configuration_id = digest_value(identity)
        job_id = configuration_id[:24]
        if (job.get("configuration_id") != configuration_id or
                job.get("job_id") != job_id or job_id in seen or
                canonical_bytes(job.get("source_identity")) !=
                canonical_bytes(source_state) or
                job.get("seed") != stable_seed(job["cell"]) or
                job.get("invocation_order") != list(invocation_order(
                    settings.get("mode"), job_id,
                    settings.get("abba_rounds", 0)
                ))):
            raise CrossoverError("{} contains a stale or duplicate job".format(path))
        frozen_required = settings.get("frozen_executable_required") is True
        if frozen_required != (job.get("executable_artifact") is not None):
            raise CrossoverError("{} has inconsistent frozen job policy".format(path))
        if frozen_required:
            validate_frozen_executable(
                job["executable_artifact"], job["build_metadata"], source_state
            )
        backend = job["cell"].get("backend")
        if backend not in expected_executables:
            expected_executables[backend] = {
                "build_metadata": job["build_metadata"],
                "execution_artifact": job["executable_artifact"],
                "origin_path": job["build_metadata"]["executable"]["path"],
                "origin_sha256": job["build_metadata"]["executable"]["sha256"],
            }
        elif canonical_bytes(expected_executables[backend]) != canonical_bytes({
                "build_metadata": job["build_metadata"],
                "execution_artifact": job["executable_artifact"],
                "origin_path": job["build_metadata"]["executable"]["path"],
                "origin_sha256": job["build_metadata"]["executable"]["sha256"],
        }):
            raise CrossoverError("{} has inconsistent backend jobs".format(path))
        seen.add(job_id)
    if canonical_bytes(executables) != canonical_bytes(expected_executables):
        raise CrossoverError("{} executable map differs from its jobs".format(path))
    return manifest


def load_manifest(result_dir):
    path = result_dir / "manifest.json"
    try:
        manifest_bytes = path.read_bytes()
    except OSError as error:
        raise CrossoverError("cannot read {}: {}".format(path, error))
    manifest = decode_json_bytes(manifest_bytes, str(path))
    return validate_manifest(manifest, path)


def require_compatible_result_dir(result_dir, manifest):
    path = result_dir / "manifest.json"
    if not path.exists():
        job_dir = result_dir / "jobs"
        if job_dir.is_dir() and any(job_dir.glob("*.json")):
            raise CrossoverError(
                "result directory has jobs but no manifest: {}".format(result_dir)
            )
        return
    previous = load_manifest(result_dir)
    if previous.get("configuration_id") != manifest.get("configuration_id"):
        raise CrossoverError(
            "result directory belongs to configuration {}; new configuration is {}; "
            "select a new --result-dir rather than mixing jobs".format(
                previous.get("configuration_id"), manifest.get("configuration_id")
            )
        )


def load_job_results(result_dir, manifest):
    job_dir = result_dir / "jobs"
    if not job_dir.is_dir():
        raise CrossoverError("job directory does not exist: {}".format(job_dir))
    expected = {}
    for job in manifest["jobs"]:
        job_id = job.get("job_id")
        configuration_id = job.get("configuration_id")
        if not job_id or not configuration_id or job_id in expected:
            raise CrossoverError("manifest contains a duplicate or incomplete job")
        expected[job_id] = job
    actual_paths = {path.stem: path for path in sorted(job_dir.glob("*.json"))}
    missing = sorted(set(expected) - set(actual_paths))
    extra = sorted(set(actual_paths) - set(expected))
    if missing or extra:
        raise CrossoverError(
            "job set does not match manifest (missing: {}; stale/extra: {})".format(
                ",".join(missing) or "none", ",".join(extra) or "none"
            )
        )
    results = []
    for job_id in sorted(expected):
        path = actual_paths[job_id]
        try:
            item_bytes = path.read_bytes()
        except OSError as error:
            raise CrossoverError("cannot read {}: {}".format(path, error))
        item = decode_json_bytes(item_bytes, str(path))
        if not isinstance(item, dict) or item.get("schema") != JOB_SCHEMA:
            raise CrossoverError("{} has a legacy or unknown job schema".format(path))
        if item.get("job_id") != job_id:
            raise CrossoverError("{} contains the wrong job ID".format(path))
        if item.get("configuration_id") != expected[job_id]["configuration_id"]:
            raise CrossoverError("{} belongs to a stale configuration".format(path))
        validate_job_artifacts(
            result_dir, item, expected[job_id], manifest["settings"]
        )
        results.append(item)
    return results


def write_merged(
        result_dir, manifest, results, source_end, promotion_percent,
        execution_input_error=""):
    ordered = sorted(results, key=lambda item: item["job_id"])
    source_changed = source_end is not None and (
        canonical_bytes(source_end) !=
        canonical_bytes(manifest["source_fingerprint"])
    )
    failed = any(item.get("status") != "passed" for item in ordered)
    run_status = "failed" if failed or source_changed or execution_input_error else "passed"
    analysis = analyze_results(
        ordered, promotion_percent, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    merged = {
        "analysis": analysis,
        "execution_input_error": execution_input_error or None,
        "jobs": ordered,
        "manifest_configuration_id": manifest["configuration_id"],
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "source_fingerprint": manifest["source_fingerprint"],
        "source_fingerprint_after": source_end,
        "status": run_status,
    }
    atomic_write_json(result_dir / "matrix.json", merged)
    atomic_write_json(result_dir / "analysis.json", analysis)
    return merged


def invalidate_authoritative_result_dir(result_dir, reason):
    """Ensure a failed authoritative guard can never leave retained PASS."""
    result_dir = Path(result_dir).resolve()
    job_dir = result_dir / "jobs"
    retained = {}
    if job_dir.is_dir():
        for path in sorted(job_dir.glob("*.json")):
            value = decode_json_bytes(path.read_bytes(), str(path))
            if (not isinstance(value, dict) or
                    value.get("schema") != JOB_SCHEMA):
                continue
            if value.get("status") == "passed":
                value.pop("parity_identity", None)
                value.pop("summary", None)
                value["status"] = "failed"
                value["reason"] = (
                    "authoritative guard validation failed: {}".format(reason)
                )
                atomic_write_json(path, value)
            retained[value.get("job_id")] = value

    matrix_path = result_dir / "matrix.json"
    manifest_path = result_dir / "manifest.json"
    if not matrix_path.is_file():
        return
    matrix = decode_json_bytes(matrix_path.read_bytes(), str(matrix_path))
    try:
        manifest = decode_json_bytes(
            manifest_path.read_bytes(), str(manifest_path)
        )
        ordered = [retained[job["job_id"]] for job in manifest["jobs"]]
        promotion = required_finite_number(
            matrix["analysis"]["promotion_percent"],
            "retained failure promotion_percent"
        )
        write_merged(
            result_dir, manifest, ordered,
            matrix.get("source_fingerprint_after"), promotion,
            matrix.get("execution_input_error") or ""
        )
    except Exception:
        # Even if older/partial artifacts cannot be recomputed, never leave a
        # human- or machine-visible PASS after this process observed guard loss.
        if isinstance(matrix, dict):
            matrix["status"] = "failed"
            analysis = matrix.get("analysis")
            if isinstance(analysis, dict):
                analysis["run_status"] = "failed"
            atomic_write_json(matrix_path, matrix)
            if isinstance(analysis, dict):
                atomic_write_json(result_dir / "analysis.json", analysis)


def print_analysis(analysis):
    candidate = analysis["candidate"]
    print(
        "candidate cells: {}/{} passed; gain min={}; median={}; "
        "regressions={}; >= {:.1f}%={}; 95% lower min={}; "
        "confident >= threshold={}/{}; promotion gate={}".format(
            candidate["cell_count"] - candidate["failed_count"],
            candidate["cell_count"],
            ("n/a" if candidate["gain_min_percent"] is None else
             "{:.2f}%".format(candidate["gain_min_percent"])),
            ("n/a" if candidate["gain_median_percent"] is None else
             "{:.2f}%".format(candidate["gain_median_percent"])),
            candidate["regression_count"], analysis["promotion_percent"],
            candidate["promotion_count"],
            ("n/a" if candidate["confidence_worst_lower_gain_percent"] is None
             else "{:.2f}%".format(
                 candidate["confidence_worst_lower_gain_percent"])),
            candidate["confidence_promotion_count"],
            candidate["cell_count"],
            ("passed" if candidate[
                "all_cells_confidently_meet_promotion_threshold"]
             else "not passed"),
        )
    )
    neighbors = analysis["excluded_neighbors"]
    print(
        "excluded neighbors: {}/{} passed; gain min={}; median={}; regressions={}".format(
            neighbors["cell_count"] - neighbors["failed_count"],
            neighbors["cell_count"],
            ("n/a" if neighbors["gain_min_percent"] is None else
             "{:.2f}%".format(neighbors["gain_min_percent"])),
            ("n/a" if neighbors["gain_median_percent"] is None else
             "{:.2f}%".format(neighbors["gain_median_percent"])),
            neighbors["regression_count"],
        )
    )


def resolve_path(source, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (source / path).resolve()


def run_matrix(arguments):
    if arguments.command not in AUTHORITATIVE_COMMANDS:
        return run_matrix_body(arguments, None)
    if arguments.cpu is None or arguments.sibling is None:
        raise CrossoverError(
            "authoritative v2 requires explicit --cpu and --sibling"
        )
    source = Path(arguments.source).resolve()
    support = load_isolation_support(source)
    try:
        allowed, housekeeping = support.validate_topology(
            arguments.cpu, arguments.sibling
        )
    except Exception as error:
        raise CrossoverError(
            "authoritative CPU topology validation failed: {}".format(error)
        )
    original_affinity = set(os.sched_getaffinity(0))
    if original_affinity != set(allowed):
        raise CrossoverError(
            "launch affinity changed during topology validation"
        )
    result_dir = resolve_path(source, arguments.result_dir)
    campaign_state = {"evidence_started": False}
    try:
        with canonical_authoritative_lock() as lock_guard:
            pair_guard = AuthoritativePairGuard(
                support.PairLease(arguments.cpu, arguments.sibling)
            )
            try:
                with pair_guard as pair_identity:
                    validate_authoritative_guards(lock_guard, pair_guard)
                    lock_identity = lock_guard.identity
                    os.sched_setaffinity(0, housekeeping)
                    isolation = {
                        "allowed": sorted(allowed),
                        "canonical_lock": lock_identity,
                        "canonical_guard": lock_guard,
                        "campaign_state": campaign_state,
                        "cpu": arguments.cpu,
                        "housekeeping": sorted(housekeeping),
                        "pair_guard": pair_guard,
                        "pair_lease": pair_identity,
                        "sibling": arguments.sibling,
                        "support": support,
                    }
                    result = run_matrix_body(arguments, isolation)
                    validate_authoritative_guards(lock_guard, pair_guard)
                    return result
            finally:
                os.sched_setaffinity(0, original_affinity)
    except Exception as error:
        if (isinstance(error, AuthoritativeGuardError) and
                campaign_state["evidence_started"]):
            try:
                invalidate_authoritative_result_dir(result_dir, error)
            except Exception as invalidation_error:
                raise CrossoverError(
                    "authoritative isolation failed: {}; retained evidence "
                    "invalidation also failed: {}".format(
                        error, invalidation_error
                    )
                )
        if isinstance(error, CrossoverError):
            raise
        raise CrossoverError(
            "authoritative isolation failed: {}".format(error)
        )


def run_matrix_body(arguments, isolation):
    source = Path(arguments.source).resolve()
    executable_root = resolve_path(
        source, arguments.executable_root or arguments.build_root
    )
    result_dir = resolve_path(source, arguments.result_dir)
    backends = parse_backends(arguments.backends)
    authoritative = arguments.command in AUTHORITATIVE_COMMANDS

    def validate_campaign_guards():
        if isolation is not None:
            validate_authoritative_guards(
                isolation["canonical_guard"], isolation["pair_guard"]
            )

    validate_campaign_guards()
    if arguments.command == "historical-avx2" and backends != ["avx2"]:
        raise CrossoverError(
            "historical-avx2 requires exactly --backends avx2"
        )
    cpus = (
        list(isolation["allowed"]) if isolation is not None
        else allowed_cpus()
    )
    pin_cpu = None
    taskset = None
    if authoritative:
        pin_cpu = arguments.cpu
        if pin_cpu not in cpus:
            raise CrossoverError(
                "pinned CPU {} is outside allowed affinity {}".format(
                    pin_cpu, compact_cpu_list(cpus)
                )
            )
        taskset = shutil.which(
            arguments.taskset, path=BENCHMARK_ENVIRONMENT["PATH"]
        )
        if not taskset or Path(taskset).resolve() != Path("/usr/bin/taskset"):
            raise CrossoverError(
                "authoritative mode requires canonical /usr/bin/taskset"
            )
        if arguments.workers != 1:
            raise CrossoverError("pinned ABBA measurements require --workers 1")
    source_state = source_identity(source, require_clean=authoritative)
    machine = machine_identity(cpus)
    if pin_cpu is not None:
        machine["pinned_cpu_topology"] = cpu_topology(pin_cpu)
        machine["reserved_sibling_topology"] = cpu_topology(arguments.sibling)
    controlled_build = None
    if arguments.command == "historical-avx2":
        validate_campaign_guards()
        executable, controlled_build = controlled_avx2_build(
            source, result_dir, source_state,
            min(len(isolation["housekeeping"]), 128),
            validate_campaign_guards,
            {
                "canonical_lock": isolation["canonical_lock"],
                "pair_lease": isolation["pair_lease"],
            }
        )
        validate_campaign_guards()
        executables = {"avx2": executable}
    else:
        executables = {
            backend: find_executable(executable_root, backend)
            for backend in backends
        }
    build_metadata = {
        backend: cmake_build_metadata(executable)
        for backend, executable in executables.items()
    }
    for backend in backends:
        validate_build_source_binding(
            build_metadata[backend], source, source_state, backend,
            require_fresh=authoritative
        )
    settings = {
        "abba_rounds": arguments.abba_rounds if authoritative else 0,
        "benchmark": {
            "batch": arguments.batch,
            "iterations": arguments.iterations,
            "reuse": arguments.reuse,
            "warmups": arguments.warmups,
        },
        "mode": arguments.command,
        "frozen_executable_required": authoritative,
        "pin_cpu": pin_cpu,
        "placement_policy": (
            "external taskset single-CPU pin"
            if authoritative else
            "unpinned discovery jobs inherit allowed affinity; OS schedules workers"
        ),
        "taskset": str(Path(taskset).resolve()) if taskset else None,
        "taskset_sha256": (
            digest_bytes(Path(taskset).read_bytes()) if taskset else None
        ),
        "timeout_seconds_per_invocation": arguments.timeout,
        "workers": arguments.workers,
    }
    if authoritative:
        settings["isolation"] = {
            "allowed_cpu_set_at_launch": isolation["allowed"],
            "benchmark_cpu": isolation["cpu"],
            "canonical_lock": isolation["canonical_lock"],
            "child_environment": dict(BENCHMARK_ENVIRONMENT),
            "housekeeping_cpu_set": isolation["housekeeping"],
            "pair_lease": isolation["pair_lease"],
            "reserved_sibling": isolation["sibling"],
        }
    if arguments.command == "historical-avx2":
        cells = historical_avx2_grid()
        settings["campaign"] = {
            "cell_count": len(cells),
            "cells_sha256": digest_value(cells),
            "name": "corrected-high-profile-explicit-avx2",
        }
        settings["controlled_build"] = controlled_build
    else:
        r_values = parse_r_values(arguments.r)
        cells = sorted_grid(backends, r_values, arguments.full)

    executable_artifacts = {}
    execution_executables = dict(executables)
    if authoritative:
        for backend in backends:
            validate_campaign_guards()
            artifact = freeze_executable(
                result_dir, backend, executables[backend],
                build_metadata[backend], source_state
            )
            validate_campaign_guards()
            executable_artifacts[backend] = artifact
            execution_executables[backend] = Path(artifact["executable"])
        source_after_freeze = source_identity(source, require_clean=True)
        if canonical_bytes(source_after_freeze) != canonical_bytes(source_state):
            raise CrossoverError("source identity changed while freezing executables")

    jobs = make_jobs(
        cells, execution_executables, build_metadata, source_state, machine,
        settings, executable_artifacts
    )
    contract = evidence_contract(authoritative)
    executable_manifest = {
        backend: {
            "build_metadata": build_metadata[backend],
            "execution_artifact": executable_artifacts.get(backend),
            "origin_path": str(path),
            "origin_sha256": digest_bytes(path.read_bytes()),
        } for backend, path in sorted(executables.items())
    }
    manifest = {
        "configuration_id": None,
        "evidence_contract": contract,
        "executables": executable_manifest,
        "jobs": jobs,
        "machine": machine,
        "schema": SCHEMA,
        "settings": settings,
        "source_fingerprint": source_state,
    }
    manifest["configuration_id"] = digest_value(manifest_identity(manifest))
    validate_manifest(manifest, result_dir / "manifest.json")
    require_compatible_result_dir(result_dir, manifest)
    validate_campaign_guards()
    atomic_write_json(result_dir / "manifest.json", manifest)
    if isolation is not None:
        isolation["campaign_state"]["evidence_started"] = True
    validate_campaign_guards()
    print(
        "direct-encode crossover: mode={} cells={} backends={} cpus={}{}".format(
            arguments.command, len(jobs), ",".join(backends),
            machine["allowed_cpu_list"],
            " pinned={}".format(pin_cpu) if pin_cpu is not None else "",
        ), flush=True
    )
    context = {
        "isolation": isolation,
        "result_dir": result_dir,
        "resume": not arguments.no_resume,
        "settings": settings,
        "source": source,
        "timeout": arguments.timeout,
    }
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = {executor.submit(run_job, job, context): job for job in jobs}
        for future in concurrent.futures.as_completed(futures):
            job = futures[future]
            result = future.result()
            results.append(result)
            print(
                "{} {}: {}".format(
                    job["cell"]["backend"], job["job_id"], result["status"]
                ), flush=True
            )
    source_end = source_identity(source, require_clean=False)
    validate_campaign_guards()
    execution_input_error = ""
    try:
        checked = set()
        for job in jobs:
            key = (job["executable"], job["executable_sha256"])
            if key not in checked:
                checked.add(key)
                validate_execution_inputs(job)
    except CrossoverError as error:
        execution_input_error = str(error)
    validate_campaign_guards()
    merged = write_merged(
        result_dir, manifest, results, source_end, arguments.promotion_percent,
        execution_input_error
    )
    validate_campaign_guards()
    print_analysis(merged["analysis"])
    promotion_passed = merged["analysis"]["candidate"][
        "all_cells_confidently_meet_promotion_threshold"
    ]
    print(
        "matrix integrity: {}; promotion gate: {} ({})".format(
            merged["status"], "passed" if promotion_passed else "not passed",
            result_dir / "matrix.json"
        )
    )
    if merged["status"] != "passed":
        return 1
    if authoritative and not promotion_passed:
        return 2
    return 0


def analyze_command(arguments):
    result_dir = Path(arguments.result_dir).resolve()
    manifest = load_manifest(result_dir)
    results = load_job_results(result_dir, manifest)
    matrix_path = result_dir / "matrix.json"
    try:
        matrix = decode_json_bytes(matrix_path.read_bytes(), str(matrix_path))
    except (CrossoverError, OSError) as error:
        raise CrossoverError("cannot read completed matrix {}: {}".format(
            matrix_path, error
        ))
    matrix_keys = {
        "analysis", "execution_input_error", "jobs",
        "manifest_configuration_id", "schema", "source_changed_during_run",
        "source_fingerprint", "source_fingerprint_after", "status",
    }
    if (not isinstance(matrix, dict) or set(matrix) != matrix_keys or
            matrix.get("schema") != SCHEMA or
            matrix.get("manifest_configuration_id") != manifest["configuration_id"]):
        raise CrossoverError("matrix does not match the current manifest")
    if canonical_bytes(matrix.get("source_fingerprint")) != canonical_bytes(
            manifest.get("source_fingerprint")):
        raise CrossoverError("matrix source fingerprint does not match the manifest")
    source_after = validate_source_state(
        matrix.get("source_fingerprint_after"),
        str(matrix_path) + " source_fingerprint_after", False
    )
    matrix_jobs_value = matrix.get("jobs")
    if (not isinstance(matrix_jobs_value, list) or
            any(not isinstance(item, dict) for item in matrix_jobs_value)):
        raise CrossoverError("matrix jobs snapshot has an invalid shape")
    matrix_jobs = sorted(
        matrix_jobs_value, key=lambda item: item.get("job_id", "")
    )
    ordered_results = sorted(results, key=lambda item: item.get("job_id", ""))
    if canonical_bytes(matrix_jobs) != canonical_bytes(ordered_results):
        raise CrossoverError("matrix job snapshot differs from the validated job files")
    source_changed = canonical_bytes(source_after) != canonical_bytes(
        manifest["source_fingerprint"]
    )
    if matrix.get("source_changed_during_run") is not source_changed:
        raise CrossoverError("matrix source-change flag was not derived")
    execution_input_error = ""
    try:
        checked = set()
        for job in manifest["jobs"]:
            key = (job["executable"], job["executable_sha256"])
            if key not in checked:
                checked.add(key)
                validate_execution_inputs(job)
    except CrossoverError as error:
        execution_input_error = str(error)
    retained_input_error = matrix.get("execution_input_error")
    if retained_input_error != (execution_input_error or None):
        raise CrossoverError("matrix execution-input error was not derived")
    failed = any(item.get("status") != "passed" for item in results)
    run_status = (
        "failed"
        if failed or source_changed or execution_input_error
        else "passed"
    )
    if matrix.get("status") != run_status:
        raise CrossoverError("matrix status was not derived from retained evidence")
    retained_analysis = required_mapping(matrix.get("analysis"), "matrix.analysis")
    retained_promotion = required_finite_number(
        retained_analysis.get("promotion_percent"),
        "matrix.analysis.promotion_percent"
    )
    if retained_promotion < 0:
        raise CrossoverError("matrix retained a negative promotion threshold")
    recomputed_retained_analysis = analyze_results(
        results, retained_promotion, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    if canonical_bytes(retained_analysis) != canonical_bytes(
            recomputed_retained_analysis):
        raise CrossoverError("matrix analysis was not derived from retained jobs")
    analysis = analyze_results(
        results, arguments.promotion_percent, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    output = Path(arguments.output).resolve() if arguments.output else result_dir / "analysis.json"
    atomic_write_json(output, analysis)
    print_analysis(analysis)
    print("analysis: {}".format(output))
    if analysis["jobs_failed"] != 0 or run_status != "passed":
        return 1
    if (manifest["settings"]["mode"] in AUTHORITATIVE_COMMANDS and
            not analysis["candidate"][
                "all_cells_confidently_meet_promotion_threshold"]):
        return 2
    return 0


def self_test():
    def check(condition, message):
        if not condition:
            raise CrossoverError("self-test failed: {}".format(message))

    check(compact_cpu_list([3, 2, 1, 7, 9, 8]) == "1-3,7-9",
          "CPU-list compaction")
    check(digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2}),
          "canonical JSON digest")
    check(stable_seed({"cell": 1}) == stable_seed({"cell": 1}),
          "stable cell seed")
    check(parse_backends("avx512") == ["avx512"], "AVX-512 backend parsing")
    check(parser().parse_args(["screen"]).backends == "scalar,ssse3,avx2",
          "default backend list")
    repository = Path(__file__).resolve().parents[1]
    repository_state = source_identity(repository, require_clean=False)
    validate_source_state(repository_state, "self-test repository", False)
    with tempfile.TemporaryDirectory(
            prefix="leo2-crossover-path-shim-") as shim_directory:
        shim = Path(shim_directory) / "git"
        shim.write_text("#!/bin/sh\nexit 97\n", encoding="utf-8")
        shim.chmod(0o755)
        old_path = os.environ.get("PATH")
        os.environ["PATH"] = shim_directory
        try:
            shim_state = source_identity(repository, require_clean=False)
        finally:
            if old_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = old_path
        check(
            canonical_bytes(shim_state) == canonical_bytes(repository_state),
            "source identity ignored the canonical absolute Git executable"
        )
    repository_gitlink = repository_state["files"].get("sse2neon")
    check(
        isinstance(repository_gitlink, dict) and
        repository_gitlink.get("index_mode") == "160000" and
        repository_gitlink.get("submodule", {}).get("head") ==
            repository_gitlink.get("index_object"),
        "real tracked sse2neon gitlink closure"
    )
    grid = sorted_grid(["scalar", "ssse3", "avx2"], [16], False)
    candidates = [item for item in grid if item["region"] == "candidate"]
    check(bool(candidates), "compact grid has candidate cells")
    for item in candidates:
        check(item["profile"] == "low" and item["q"] == 1,
              "candidate profile and Q")
        check(item["shard_bytes"] >= 1024 and item["shard_bytes"] % 64 == 0,
              "candidate byte region")
        check(item["k"] >= threshold_k(item["backend"]),
              "candidate backend K threshold")
    regions = {item["region"] for item in grid}
    check({
        "excluded_scalar_k2", "excluded_q2", "excluded_high_profile",
        "excluded_ragged_tail",
    }.issubset(regions), "excluded-neighbor regions")
    check({item["field"] for item in candidates} == {"gf8", "gf16"},
          "candidate field coverage")
    check(any(item["region"] == "candidate_sparse_output" for item in grid),
          "sparse-output candidate")
    check(any(item["region"] == "excluded_q2_holey" for item in grid),
          "holey-Q2 neighbor")
    r1 = sorted_grid(["scalar"], [1], False)
    check(all(item["r"] == 1 and item["q"] == 1 for item in r1),
          "R=1 grid bounds Q")
    check(parse_r_values("1,16,1") == [1, 16], "R-list normalization")
    check(parse_r_values("all") == list(range(1, 17)), "all-R expansion")
    check(requested_indices("0,2-3", 4) == [0, 2, 3],
          "parity-mask expansion")
    check(invocation_order("pinned", "00000000", 2) == (
        "direct", "transform", "transform", "direct",
        "direct", "transform", "transform", "direct",
    ), "pinned ABBA order")

    historical = historical_avx2_grid()
    exact = [
        item for item in historical
        if item["region"] == "candidate_historical_exact"
    ]
    check(len(exact) == 1, "one exact historical cell")
    check(exact[0] == cell(
        "candidate_historical_exact", "avx2", 2, 16, "high",
        "gf8", 4096, 1, "0"
    ), "exact K=2,R=16,Q=1,4096-byte historical cell")
    check({
        "historical_neighbor_k", "historical_neighbor_r",
        "historical_neighbor_bytes_lower", "historical_neighbor_bytes_upper",
        "historical_neighbor_sparse_output", "historical_neighbor_q2",
    }.issubset({item["region"] for item in historical}),
          "historical one-coordinate neighbors")

    raw_job = {
        "build_metadata": {
            "entries": {
                "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_CXX_COMPILER": "/self-test/c++",
                "CMAKE_CXX_FLAGS": "",
                "CMAKE_CXX_FLAGS_DEBUG": "-g",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
                "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
                "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
                "ENABLE_OPENMP": "ON",
                "LEOPARD_ENABLE_GF8": "ON",
                "LEOPARD_ENABLE_GF16": "ON",
                "LEO2_BACKEND_VARIANT": "avx2",
                "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
                "LEO2_BUILD_BENCHMARKS": "ON",
                "LEO2_BUILD_TESTS": "ON",
                "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            },
        },
        "cell": exact[0],
        "seed": 12345,
        "source_identity": {
            "git": {
                "head": "1" * 40, "tree": "2" * 40,
                "worktree_clean": True,
            },
        },
    }
    raw_settings = {
        "benchmark": {
            "batch": 1, "iterations": 15, "reuse": 64, "warmups": 4,
        },
    }
    raw_fixture = {
        "schema": BENCHMARK_SCHEMA,
        "build": {
            "backend_variant": "avx2",
            "build_type": "Release",
            "build_configuration_sha256": "",
            "compiler": "self-test", "compiler_version": "1",
            "cplusplus": 201103, "test_hooks": True,
            "source_commit": "1" * 40,
            "source_tree": "2" * 40,
            "source_tracked_dirty": False,
        },
        "parameters": {
            "K": 2, "R": 16, "Q": 1,
            "forced_mode": "force_direct",
            "shard_bytes": 4096,
            "seed": 12345,
            "requested_parity_indices": [0],
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8",
            "batch": 1, "reuse": 64, "iterations": 15, "warmups": 4,
            "thread_count": 1,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": "avx2",
            "direct_capable": True, "thread_count": 1,
            "parent_count": 32, "padded_side": 16,
            "timed_path_is_direct": True,
        },
        "correctness": {
            "selected_transform_reference_parity_match": True,
            "direct_transform_parity_match": True,
            "unrequested_outputs_untouched": True,
            "parity_checksum_fnv1a64": "0x0123456789abcdef",
        },
        "memory": {},
        "operation_model": {
            "direct_row_terms": 2,
            "direct_output_initializations": 1,
            "direct_output_accumulations": 1,
            "fixed_coefficient_symbol_terms_before_unit_specialization":
                8192,
            "xor_accumulation_symbols": 4096,
            "modeled_source_bytes_read": 8192,
            "modeled_output_bytes_read": 4096,
            "modeled_output_bytes_written": 8192,
            "model_scope": (
                "direct streaming kernels before cache effects; unit "
                "coefficients specialize to copy/XOR"
            ),
            "direct_model_applies_to_timed_path": True,
            "transform_operation_counts": None,
            "hardware_counters": None,
        },
        "metrics": {
            # This deliberate distractor is the value selected by the invalid
            # recursive report-20 extractor.
            "codec_setup": {"median_us": 987654.0, "mad_us": 0.0},
            "encode_execution": {
                "median_us_per_batch_call": 7.25,
                "mad_us_per_batch_call": 0.5,
                "minimum_us_per_batch_call": 7.0,
                "maximum_us_per_batch_call": 8.0,
                "logical_input_GB_per_s": 1.0,
                "requested_parity_output_GB_per_s": 1.0,
                "logical_input_plus_output_GB_per_s": 2.0,
            },
        },
        "methodology": {},
    }
    raw_fixture["build"]["build_configuration_sha256"] = (
        build_configuration_digest(raw_job["build_metadata"]["entries"])
    )
    median_us, mad_us, parity_identity = validate_raw(
        raw_fixture, raw_job, "direct", raw_settings
    )
    check(median_us == 7.25 and mad_us == 0.5,
          "exact encode-execution metric path")
    check(parity_identity == {
        "algorithm": "fnv1a64",
        "digest": "0x0123456789abcdef",
        "hashed_bytes": 4096,
        "requested_parity_indices": [0],
    }, "non-vacuous parity identity")
    for label, encoded in (
            ("duplicate key", b'{"metric": 1, "metric": 2}'),
            ("non-standard NaN", b'{"metric": NaN}'),
            ("non-standard infinity", b'{"metric": Infinity}'),
            ("overflowing exponent", b'{"metric": 1e9999}'),
            ("huge integer", b'{"metric": 100000000000000000000}')):
        try:
            decode_json_bytes(encoded, "self-test {}".format(label))
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} JSON was accepted".format(label)
            )

    def reject_raw_mutation(name, mutator):
        changed = json.loads(json.dumps(raw_fixture))
        mutator(changed)
        try:
            validate_raw(changed, raw_job, "direct", raw_settings)
        except CrossoverError:
            return
        raise CrossoverError(
            "self-test failed: raw mutation was accepted: {}".format(name)
        )

    reject_raw_mutation(
        "wrong schema", lambda value: value.update({"schema": "unknown"})
    )
    reject_raw_mutation(
        "missing exact metric object",
        lambda value: value["metrics"].pop("encode_execution")
    )
    reject_raw_mutation(
        "codec-setup substitution",
        lambda value: value["metrics"]["encode_execution"].pop(
            "median_us_per_batch_call"
        )
    )
    reject_raw_mutation(
        "numeric string",
        lambda value: value["metrics"]["encode_execution"].update(
            {"median_us_per_batch_call": "7.25"}
        )
    )
    reject_raw_mutation(
        "boolean median",
        lambda value: value["metrics"]["encode_execution"].update(
            {"median_us_per_batch_call": True}
        )
    )
    for label, invalid in (
            ("NaN median", float("nan")),
            ("positive-infinity median", float("inf")),
            ("negative-infinity median", float("-inf")),
            ("overflowing integer median", 10 ** 10000)):
        reject_raw_mutation(
            label,
            lambda value, invalid=invalid:
                value["metrics"]["encode_execution"].update(
                    {"median_us_per_batch_call": invalid}
                )
        )
    for label, invalid in (
            ("negative MAD", -1.0),
            ("NaN MAD", float("nan")),
            ("infinite MAD", float("inf"))):
        reject_raw_mutation(
            label,
            lambda value, invalid=invalid:
                value["metrics"]["encode_execution"].update(
                    {"mad_us_per_batch_call": invalid}
                )
        )
    reject_raw_mutation(
        "empty parity digest string",
        lambda value: value["correctness"].update(
            {"parity_checksum_fnv1a64": ""}
        )
    )
    reject_raw_mutation(
        "empty parity digest object",
        lambda value: value["correctness"].update(
            {"parity_checksum_fnv1a64": {}}
        )
    )
    reject_raw_mutation(
        "empty correctness object",
        lambda value: value.update({"correctness": {}})
    )
    reject_raw_mutation(
        "extra top-level field",
        lambda value: value.update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra build field",
        lambda value: value["build"].update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra parameter field",
        lambda value: value["parameters"].update({"unversioned_extra": 1})
    )
    reject_raw_mutation(
        "extra resolved field",
        lambda value: value["resolved"].update({"unversioned_extra": 1})
    )
    reject_raw_mutation(
        "extra correctness field",
        lambda value: value["correctness"].update({"unversioned_extra": True})
    )
    reject_raw_mutation(
        "extra operation-model field",
        lambda value: value["operation_model"].update(
            {"unversioned_extra": 1}
        )
    )
    reject_raw_mutation(
        "extra encode-execution field",
        lambda value: value["metrics"]["encode_execution"].update(
            {"unversioned_extra": 1}
        )
    )
    reject_raw_mutation(
        "missing operation model",
        lambda value: value.pop("operation_model")
    )
    reject_raw_mutation(
        "flipped operation-model applicability",
        lambda value: value["operation_model"].update(
            {"direct_model_applies_to_timed_path": False}
        )
    )
    reject_raw_mutation(
        "wrong direct-row term count",
        lambda value: value["operation_model"].update(
            {"direct_row_terms": 3}
        )
    )
    reject_raw_mutation(
        "boolean K parameter",
        lambda value: value["parameters"].update({"K": True})
    )
    reject_raw_mutation(
        "wrong embedded backend variant",
        lambda value: value["build"].update({"backend_variant": "auto"})
    )
    reject_raw_mutation(
        "wrong embedded source commit",
        lambda value: value["build"].update({"source_commit": "3" * 40})
    )

    transform_fixture = json.loads(json.dumps(raw_fixture))
    transform_fixture["parameters"]["forced_mode"] = "force_transform"
    transform_fixture["resolved"]["timed_path_is_direct"] = False
    transform_fixture["correctness"]["direct_transform_parity_match"] = None
    transform_fixture["operation_model"][
        "direct_model_applies_to_timed_path"
    ] = False
    validate_raw(transform_fixture, raw_job, "transform", raw_settings)
    transform_fixture["correctness"]["direct_transform_parity_match"] = True
    try:
        validate_raw(transform_fixture, raw_job, "transform", raw_settings)
    except CrossoverError:
        pass
    else:
        raise CrossoverError(
            "self-test failed: transform invocation claimed direct parity"
        )

    measurements = [
        {"timed_mode": "direct", "median_us": 8.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "direct", "median_us": 8.0},
    ]
    summary = summarize_measurements(measurements)
    check(summary["gain_percent"] == 25.0, "gain summary")
    check(math.isclose(
        summary["rounds"][0]["gain_percent"], 25.0,
        rel_tol=1e-12, abs_tol=1e-12
    ), "ABBA round summary")
    inference_summary = summarize_measurements(measurements * 3)
    inference = inference_summary["paired_log_inference"]
    check(
        inference["degrees_of_freedom"] == 2 and
        abs(inference["gain_percent"] - 25.0) < 1e-12 and
        abs(inference["gain_percent_student_t_interval"][0] - 25.0) < 1e-12,
        "three-round paired-log Student-t inference"
    )
    asymmetric = summarize_measurements([
        {"timed_mode": "direct", "median_us": 4.0},
        {"timed_mode": "transform", "median_us": 9.0},
        {"timed_mode": "transform", "median_us": 25.0},
        {"timed_mode": "direct", "median_us": 16.0},
    ])
    check(math.isclose(
        asymmetric["rounds"][0]["direct_geometric_mean_us"], 8.0,
        rel_tol=1e-15, abs_tol=1e-15
    ) and math.isclose(
        asymmetric["rounds"][0]["transform_geometric_mean_us"], 15.0,
        rel_tol=1e-15, abs_tol=1e-15
    ), "ABBA round uses geometric representatives")
    for label, invalid in (
            ("string", "8.0"), ("boolean", True),
            ("NaN", float("nan")), ("infinity", float("inf"))):
        changed = json.loads(json.dumps(measurements))
        changed[0]["median_us"] = invalid
        try:
            summarize_measurements(changed)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: {} summary metric was accepted".format(label)
            )
    fake = {
        "cell": candidates[0], "job_id": "a", "status": "passed",
        "summary": summary,
    }
    analysis = analyze_results([fake], 5.0)
    check(analysis["candidate"]["gain_min_percent"] == 25.0,
          "candidate minimum gain")
    check(analysis["candidate"]["regression_count"] == 0,
          "candidate regression count")
    with tempfile.TemporaryDirectory(prefix="leo2-crossover-self-test-") as directory:
        root = Path(directory)
        path = root / "stable.json"
        atomic_write_json(path, {"z": [3, 2, 1], "a": "stable"})
        check(json.loads(path.read_text(encoding="utf-8"))["a"] == "stable",
              "atomic JSON output")
        record = run_command(
            [sys.executable, "-c", "print('ok')"], root,
            root / "stdout.log", root / "stderr.log", 5, os.environ.copy()
        )
        check(record["returncode"] == 0, "self-test subprocess exit")
        check((root / "stdout.log").read_bytes() == b"ok\n",
              "self-test subprocess output")
        result_root = root / "results"
        manifest_path = result_root / "manifest.json"
        for label, invalid_manifest in (
                ("v1 manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema":
                        "leopard2-direct-encode-crossover/v1",
                    "settings": {},
                }),
                ("relabeled v1 manifest", {
                    "configuration_id": "legacy",
                    "jobs": [], "schema": SCHEMA, "settings": {},
                })):
            atomic_write_json(manifest_path, invalid_manifest)
            try:
                load_manifest(result_root)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} was accepted".format(label)
                )
        submodule = root / "submodule"
        superproject = root / "superproject"

        def fixture_git(cwd, arguments):
            completed = subprocess.run(
                [str(CANONICAL_GIT)] + list(arguments), cwd=str(cwd),
                env=GIT_ENVIRONMENT,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=20, check=False
            )
            if completed.returncode != 0:
                raise CrossoverError(
                    "self-test git fixture failed: {}".format(
                        normalized_output(completed.stderr).decode(
                            "utf-8", errors="replace"
                        ).strip()
                    )
                )
            return normalized_output(completed.stdout)

        submodule.mkdir()
        fixture_git(submodule, ("init", "-q"))
        fixture_git(submodule, ("config", "user.name", "Self Test"))
        fixture_git(submodule, ("config", "user.email", "self@test.invalid"))
        (submodule / "tracked.txt").write_text("clean\n", encoding="utf-8")
        fixture_git(submodule, ("add", "tracked.txt"))
        fixture_git(submodule, ("commit", "-q", "-m", "initial"))
        (submodule / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )
        fixture_git(submodule, ("commit", "-q", "-am", "current"))
        superproject.mkdir()
        fixture_git(superproject, ("init", "-q"))
        fixture_git(superproject, ("config", "user.name", "Self Test"))
        fixture_git(superproject, ("config", "user.email", "self@test.invalid"))
        fixture_git(superproject, (
            "-c", "protocol.file.allow=always", "submodule", "add", "-q",
            str(submodule.resolve()), "module",
        ))
        fixture_git(superproject, ("commit", "-q", "-am", "gitlink"))
        (superproject / "module" / "tracked.txt").write_text(
            "worktree-only change\n", encoding="utf-8"
        )
        fixture_state = source_identity(superproject, require_clean=False)
        validate_source_state(fixture_state, "self-test gitlink fixture", False)
        fixture_submodule = fixture_state["files"]["module"]["submodule"]
        check(
            fixture_submodule["worktree_clean"] is False and
            fixture_submodule["status"] == [" M tracked.txt"],
            "gitlink retains exact worktree-only porcelain status"
        )
        module = superproject / "module"
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )
        clean_fixture = source_identity(
            superproject, require_clean=True
        )

        def reject_hidden_source(label):
            try:
                source_identity(superproject, require_clean=True)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "self-test failed: {} was accepted".format(label)
                )

        fixture_git(module, (
            "update-index", "--assume-unchanged", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "assume-unchanged mutation\n", encoding="utf-8"
        )
        reject_hidden_source("assume-unchanged dirty submodule")
        fixture_git(module, (
            "update-index", "--no-assume-unchanged", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )

        fixture_git(module, (
            "update-index", "--skip-worktree", "tracked.txt",
        ))
        (module / "tracked.txt").unlink()
        reject_hidden_source("skip-worktree missing submodule file")
        fixture_git(module, (
            "update-index", "--no-skip-worktree", "tracked.txt",
        ))
        (module / "tracked.txt").write_text(
            "clean current\n", encoding="utf-8"
        )

        fixture_git(module, ("config", "core.filemode", "false"))
        (module / "tracked.txt").chmod(0o755)
        reject_hidden_source("status-hidden executable-bit mutation")
        (module / "tracked.txt").chmod(0o644)

        fixture_git(module, ("replace", "HEAD", "HEAD^"))
        reject_hidden_source("Git replacement ref")
        fixture_git(module, ("replace", "-d", "HEAD"))

        hostile_config = root / "hostile.gitconfig"
        hostile_config.write_text(
            "[core]\n\tbare = true\n", encoding="utf-8"
        )
        inherited_names = ("GIT_CONFIG_GLOBAL", "HOME", "PATH")
        inherited_before = {
            name: os.environ.get(name) for name in inherited_names
        }
        os.environ["GIT_CONFIG_GLOBAL"] = str(hostile_config)
        os.environ["HOME"] = str(root / "hostile-home")
        os.environ["PATH"] = str(root)
        try:
            hostile_fixture = source_identity(
                superproject, require_clean=True
            )
        finally:
            for name, old_value in inherited_before.items():
                if old_value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = old_value
        check(
            canonical_bytes(hostile_fixture) ==
                canonical_bytes(clean_fixture),
            "hostile inherited Git config or PATH changed source identity"
        )

        fixture_git(module, ("config", "user.name", "Self Test"))
        fixture_git(module, ("config", "user.email", "self@test.invalid"))
        (module / "tracked.txt").write_text(
            "different submodule commit\n", encoding="utf-8"
        )
        fixture_git(module, ("commit", "-q", "-am", "different head"))
        try:
            source_identity(superproject, require_clean=False)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: mismatched gitlink HEAD was accepted"
            )

        fixture_git(superproject, (
            "submodule", "deinit", "-q", "-f", "--", "module",
        ))
        try:
            source_identity(superproject, require_clean=False)
        except CrossoverError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: uninitialized gitlink was accepted"
            )
    for invalid in ("1", "4"):
        try:
            exactly_three_abba_rounds(invalid)
        except argparse.ArgumentTypeError:
            pass
        else:
            raise CrossoverError(
                "self-test failed: authoritative ABBA rounds {} accepted".format(
                    invalid
                )
            )
    print(
        "leopard2 direct-encode crossover self-test passed "
        "(strict JSON/schema, provenance, parity, metric, summary, "
        "Git-closure, and ABBA mutations; no codec required)"
    )
    return 0


def finite_percentage(text):
    try:
        value = float(text)
    except (OverflowError, TypeError, ValueError):
        raise argparse.ArgumentTypeError("must be a finite nonnegative number")
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be a finite nonnegative number")
    return value


def exactly_three_abba_rounds(text):
    try:
        value = int(text)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            "authoritative v2 requires exactly 3 ABBA rounds"
        )
    if value != 3:
        raise argparse.ArgumentTypeError(
            "authoritative v2 requires exactly 3 ABBA rounds"
        )
    return value


def add_run_arguments(parser, default_result):
    default_source = str(Path(__file__).resolve().parents[1])
    parser.add_argument("--source", default=default_source)
    parser.add_argument("--build-root", default="build")
    parser.add_argument(
        "--executable-root", default=None,
        help="executable root; may contain a literal {backend} placeholder"
    )
    parser.add_argument("--result-dir", default=default_result)
    parser.add_argument("--backends", default="scalar,ssse3,avx2")
    parser.add_argument(
        "--r", default="16",
        help="R value, comma-separated R values, or 'all' (default 16)"
    )
    parser.add_argument("--batch", type=int, default=1)
    parser.add_argument("--reuse", type=int, default=32)
    parser.add_argument("--iterations", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--promotion-percent", type=finite_percentage, default=5.0)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--no-resume", action="store_true")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command")
    screen = subparsers.add_parser(
        "screen", help="run the compact grid without authoritative pinning"
    )
    add_run_arguments(
        screen, "results/leopard2/direct-encode-crossover/screen"
    )
    screen.set_defaults(workers=min(len(allowed_cpus()), 128))
    pinned = subparsers.add_parser(
        "pinned", help="run isolated, externally pinned ABBA measurements"
    )
    add_run_arguments(
        pinned, "results/leopard2/direct-encode-crossover/pinned"
    )
    pinned.set_defaults(iterations=15, warmups=4, reuse=64, workers=1)
    pinned.add_argument("--cpu", type=int, default=None)
    pinned.add_argument("--sibling", type=int, default=None)
    pinned.add_argument("--taskset", default="taskset")
    pinned.add_argument(
        "--abba-rounds", type=exactly_three_abba_rounds, default=3
    )
    historical = subparsers.add_parser(
        "historical-avx2",
        help=(
            "run the corrected high-profile cells with a frozen explicit-AVX2 "
            "binary and pinned ABBA order"
        ),
    )
    add_run_arguments(
        historical,
        "results/leopard2/direct-encode-crossover/historical-avx2",
    )
    historical.set_defaults(
        backends="avx2", iterations=15, warmups=4, reuse=64, workers=1
    )
    historical.add_argument("--cpu", type=int, default=None)
    historical.add_argument("--sibling", type=int, default=None)
    historical.add_argument("--taskset", default="taskset")
    historical.add_argument(
        "--abba-rounds", type=exactly_three_abba_rounds, default=3
    )
    analyze = subparsers.add_parser(
        "analyze", help="deterministically reanalyze completed job JSON"
    )
    analyze.add_argument(
        "--result-dir",
        default="results/leopard2/direct-encode-crossover/pinned",
    )
    analyze.add_argument(
        "--promotion-percent", type=finite_percentage, default=5.0
    )
    analyze.add_argument("--output", default=None)
    subparsers.add_parser("self-test", help="test the runner without a built codec")
    return result


def main():
    arguments = parser().parse_args()
    try:
        if arguments.command == "self-test":
            return self_test()
        if arguments.command == "analyze":
            return analyze_command(arguments)
        if arguments.command in ("screen", "pinned", "historical-avx2"):
            numeric = (
                arguments.batch, arguments.reuse, arguments.iterations,
                arguments.warmups, arguments.timeout, arguments.workers,
            )
            if any(value <= 0 for value in numeric):
                raise CrossoverError(
                    "batch, reuse, iterations, warmups, timeout, and workers "
                    "must be positive"
                )
            required_finite_number(
                arguments.promotion_percent, "promotion_percent"
            )
            if (arguments.command in AUTHORITATIVE_COMMANDS and
                    arguments.abba_rounds != 3):
                raise CrossoverError(
                    "authoritative paired-log inference requires --abba-rounds 3"
                )
            return run_matrix(arguments)
        parser().print_help()
        return 2
    except (CrossoverError, OSError, subprocess.SubprocessError) as error:
        print("direct-encode crossover error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
