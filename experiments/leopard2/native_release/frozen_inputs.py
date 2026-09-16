"""Small exact-inventory, read-only artifact guard for the native release lane."""
import os
from pathlib import Path
import re
import stat

from check_encode import require, sha256

TIMING_ARTIFACTS = frozenset({
    "main-steady", "main-synthetic", "main-abort", "current-steady", "current-synthetic",
    "current-abort", "timing-unit", "timing-unit-sanitized", "main.a", "current.a"})


class FrozenInputs:
    def __init__(self, directory, names, expected_manifest_sha256=None):
        self.directory = Path(directory)
        self.manifest = self.directory / "SHA256SUMS"
        self.manifest_sha256 = sha256(self.manifest)
        require(expected_manifest_sha256 is None or self.manifest_sha256 == expected_manifest_sha256,
                "manifest differs from preregistration")
        self.pins = {}
        for line in self.manifest.read_text().splitlines():
            words = line.split()
            require(len(words) == 2, "invalid manifest line")
            digest, name = words
            require(re.fullmatch("[0-9a-f]{64}", digest) and name not in self.pins and
                    re.fullmatch("[A-Za-z0-9][A-Za-z0-9._-]*", name), "invalid or duplicate artifact")
            self.pins[name] = digest
        require(set(self.pins) == set(names), "artifact inventory mismatch")
        self.identities = {}
        self.verify()

    def verify(self):
        for name in ["SHA256SUMS", *self.pins]:
            path = self.directory / name
            before = path.lstat()
            require(stat.S_ISREG(before.st_mode) and before.st_uid == os.getuid() and
                    not before.st_mode & 0o222, "unfrozen regular file required: " + name)
            digest = self.manifest_sha256 if name == "SHA256SUMS" else self.pins[name]
            require(sha256(path) == digest, "artifact bytes changed: " + name)
            after = path.lstat()
            identity = lambda value: (value.st_dev, value.st_ino, value.st_size,
                                      value.st_mode, value.st_ctime_ns, value.st_mtime_ns)
            require(identity(before) == identity(after), "artifact changed during hashing")
            require(name not in self.identities or self.identities[name] == identity(after),
                    "artifact metadata changed: " + name)
            self.identities[name] = identity(after)

    def executable(self, name):
        require(name in self.pins, "executed artifact is not pinned")
        self.verify()
        return str(self.directory / name)
