#!/usr/bin/python3
"""Explicit sealed-loader compiler dispatch; leopard-79h.38.5.4.8.2.2.2.3.

No complete-build integration, timing or broad runtime-closure claim. Caller
retains host/resources, sources, exact recipes and all non-executable inputs.
"""
from __future__ import annotations
from contextlib import ExitStack
import copy
import importlib.util
import os
from pathlib import Path
import re
import secrets
import struct

HERE = Path(__file__).resolve().parent
dependency = HERE / "v19_compiler_headers.py"
if dependency.resolve(strict=True) != dependency:
    raise RuntimeError("runtime dispatch dependency is not canonical")
spec = importlib.util.spec_from_file_location("v19_dispatch_headers", dependency)
header_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(header_module)
runtime = header_module.runtime
compiler, builder, provenance, require = runtime.compiler, runtime.builder, runtime.provenance, runtime.require
TEMPLATE = HERE / "v19_runtime_trampoline.S"
MAX_ARTIFACT = 65536


def fd_path(fd):
    require(type(fd) is int and 3 <= fd < 65536, "dispatcher descriptor is outside profile")
    return f"/proc/self/fd/{fd}"


def live_bindings(inventory):
    inventory.validate_current()
    phase = inventory.phase
    return {"loader_label": str(inventory.loader_path),
            "loader_descriptor": fd_path(inventory._files[inventory.loader_key].executable_descriptor),
            "library_prefix": fd_path(inventory.prefix.descriptor),
            "plugin_original": inventory._roots["gcc-link-plugin"],
            "plugin_descriptor": fd_path(inventory._files[inventory._roots["gcc-link-plugin"]].executable_descriptor),
            "helpers": {role: fd_path(tool.executable_descriptor) for role, tool in phase._tools.items() if role != "driver"}}


def source_bytes(template, bindings):
    require(type(template) is bytes and 0 < len(template) <= 16384 and template.endswith(b"\n"), "dispatcher template exceeds bound")
    require(type(bindings) is dict and set(bindings) == {"loader_label", "loader_descriptor", "library_prefix",
            "plugin_original", "plugin_descriptor", "helpers"}, "dispatcher bindings differ")
    helpers = bindings["helpers"]
    require(type(helpers) is dict and set(helpers) in ({"cc1", "as", "collect2", "ld"}, {"cc1plus", "as", "collect2", "ld"}),
            "dispatcher helper roles differ")
    for label, value in bindings.items():
        if label == "helpers": continue
        require(type(value) is str and len(value) < 4096 and re.fullmatch(r"/[A-Za-z0-9_./+-]+", value), "unsafe dispatcher binding string")
        builder.host.canonical_path(value, nonroot=True)
    for value in (bindings["loader_descriptor"], bindings["plugin_descriptor"], bindings["library_prefix"], *helpers.values()):
        require(type(value) is str and re.fullmatch(r"/proc/self/fd/[0-9]+", value) and fd_path(int(value.rsplit("/", 1)[1])) == value,
                "invalid dispatcher fd binding")
    lines = [".section .rodata"]
    for label in sorted(set(bindings) - {"helpers"}): lines.append(f'{label}: .asciz "{bindings[label]}"')
    roles = sorted(helpers)
    for index, role in enumerate(roles):
        lines += [f'role_{index}: .asciz "{role}"', f'target_{index}: .asciz "{helpers[role]}"']
    lines += [".balign 8", "helper_table:"] + [f".quad role_{index}, target_{index}" for index in range(4)]
    result = template + ("\n".join(lines) + "\n").encode("ascii")
    require(len(result) <= 32768, "dispatcher source exceeds bound")
    return result


def validate_static_elf(data):
    require(type(data) is bytes and 120 <= len(data) <= MAX_ARTIFACT, "dispatcher ELF exceeds bound")
    h = struct.unpack_from("<16sHHIQQQIHHHHHH", data)
    require(h[0][:7] == b"\x7fELF\x02\x01\x01" and h[1] == 2 and h[2] == 62 and h[3] == 1 and
            h[8] == 64 and h[9] == 56 and 0 < h[10] <= 16 and 64 <= h[5] <= len(data) - 56*h[10], "dispatcher ELF profile differs")
    segments = [struct.unpack_from("<IIQQQQQQ", data, h[5] + index*56) for index in range(h[10])]
    require(all(row[0] not in (2, 3) for row in segments), "dispatcher has a dynamic runtime dependency")
    loads = [row for row in segments if row[0] == 1]
    require(loads and any(row[1] & 1 and row[3] <= h[4] < row[3] + row[5] for row in loads), "dispatcher entry is not in executable file bytes")
    require(all(row[5] <= row[6] and row[2] <= len(data) and row[5] <= len(data)-row[2] and row[1] & 3 != 3 for row in loads),
            "dispatcher load segment is unsafe")
    stacks = [row for row in segments if row[0] == 0x6474e551]
    require(len(stacks) == 1 and stacks[0][1] == 6, "dispatcher requires exactly one non-executable stack declaration")


class RuntimeDispatch:
    def __init__(self, inventory, *, headers=None, _runner=None):
        require(type(inventory) is runtime.RuntimeInventory or _runner is not None, "dispatch requires a live runtime inventory")
        self.inventory, self.phase = inventory, inventory.phase
        require(headers is None or (type(headers) is header_module.CompilerHeaders and headers.inventory is inventory),
                "dispatch header owner differs from inventory")
        self.headers = headers
        self._runner = provenance._run if _runner is None else _runner
        self._stack, self._state = ExitStack(), "new"
        self._snapshots, self._commands = {}, []
        self._shim = None
        self.prefix = None

    def _hold(self, path):
        snapshot = self._stack.enter_context(provenance._RetainedFileSnapshot(path, "v19 dispatcher artifact", maximum_bytes=MAX_ARTIFACT))
        self._snapshots[str(path)] = snapshot
        return snapshot

    def _validate_inputs(self):
        self.inventory.validate_current()
        if self.headers is not None: self.headers.validate_current()
        require(live_bindings(self.inventory) == self.bindings, "dispatcher fd bindings changed")
        self._guard.verify()
        require(not os.get_inheritable(self._root_fd) and builder.streamed._directory_identity(os.fstat(self._root_fd)) ==
                self._root_identity == builder.streamed._directory_identity(self.root.lstat()), "dispatcher root changed")
        for snapshot in self._snapshots.values():
            builder.owners.verify_current_bytes(snapshot, snapshot.identity["sha256"])
        if self._shim is not None: self._shim.validate_current()
        if self.prefix is not None:
            self.prefix.verify()
            fd_path(self.prefix.descriptor)
            require(not os.get_inheritable(self.prefix.descriptor) and not os.get_inheritable(self._prefix_fd) and
                    provenance._stable_fields(os.fstat(self._prefix_fd)) ==
                    provenance._stable_fields(os.fstat(self.prefix.descriptor)), "dispatcher prefix descriptor changed")
        self._guard.verify()

    def _invoke(self, role, argv, *, logical_argv=None):
        self._validate_inputs()
        tool = self.phase._tools[role]
        loader = self.inventory._files[self.inventory.loader_key]
        effective = [self.bindings["loader_label"], "--inhibit-cache", "--glibc-hwcaps-mask", "", "--library-path",
                     self.bindings["library_prefix"], "--argv0", argv[0], fd_path(tool.executable_descriptor), *argv[1:]]
        record = {"role": role, "logical_argv": list(argv if logical_argv is None else logical_argv),
                  "effective_argv": list(argv), "loader_argv": effective,
                  "umask_policy": "caller" if role == "driver" else "0022", "status": "running"}
        self._commands.append(record)
        descriptors = [self.inventory.prefix.descriptor, *(tool.executable_descriptor for tool in self.phase._tools.values()),
                       *(tool.executable_descriptor for tool in self.inventory._files.values())]
        if self.prefix is not None: descriptors += [self.prefix.descriptor, self._shim.executable_descriptor]
        if self.headers is not None: descriptors += list(self.headers.descriptors())
        try:
            output = self._runner(effective, "v19 runtime dispatch " + role, maximum_bytes=1 << 20, timeout=600,
                executable_descriptor=loader.executable_descriptor, inherited_descriptors=tuple(descriptors), environment_overrides=builder.ENVIRONMENT)
            self._validate_inputs()
            record.update(status="exit-zero", stdout_sha256=builder.hashlib.sha256(output).hexdigest())
            return output
        except BaseException as error:
            record.update(status="failed", failure=type(error).__name__ + ": " + str(error))
            self._state = "failed"
            raise

    def __enter__(self):
        require(self._state == "new", "dispatcher cannot be reused")
        self._state = "entering"
        try:
            self.bindings = live_bindings(self.inventory)
            for role in self.inventory._roots:
                if role not in self.inventory._lists: self.inventory.list_dependencies(role)
            name = "v19-dispatch-" + secrets.token_hex(16)
            os.mkdir(name, mode=0o700, dir_fd=self.phase._parent_fd)
            self.root = self.phase.parent / name
            self._root_fd = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=self.phase._parent_fd)
            self._stack.callback(os.close, self._root_fd)
            self.phase.validate_current()
            value = os.fstat(self._root_fd)
            require(value.st_uid == os.getuid() and value.st_gid == os.getgid() and value.st_mode & 0o7777 == 0o700, "dispatcher root is not private")
            self._root_identity = builder.streamed._directory_identity(value)
            self._guard = self._stack.enter_context(provenance._InotifyMutationGuard("v19 dispatcher roots"))
            self._guard.add_directory_path(self.root)
            self._guard._add_watch(self.root, self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR | provenance.IN_ATTRIB, {b""})
            template = self._hold(TEMPLATE)
            require(template.resolved == TEMPLATE, "dispatcher template is not canonical")
            generated = source_bytes(template.content, self.bindings)
            source = self.root / "dispatch.S"
            fd = os.open(source.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600, dir_fd=self._root_fd)
            stream = None
            try:
                stream = os.fdopen(fd, "wb")
                with stream:
                    require(stream.write(generated) == len(generated), "dispatcher source write was incomplete")
            except BaseException:
                if stream is None: os.close(fd)
                else: stream.close()
                raise
            require(self._hold(source).content == generated, "generated dispatcher source differs")
            obj, binary = self.root / "dispatch.o", self.root / "dispatch"
            previous_umask = os.umask(0o022)
            try:
                self._invoke("as", [self.phase.pins["as"]["path"], "--64", "-o", str(obj), str(source)])
                self._hold(obj)
                self._invoke("ld", [self.phase.pins["ld"]["path"], "-static", "--build-id=none", "-z", "noexecstack", "-o", str(binary), str(obj)])
            finally:
                os.umask(previous_umask)
            final = self._hold(binary)
            validate_static_elf(final.content)
            self._shim = builder._StreamedTool(binary, maximum_bytes=MAX_ARTIFACT, _trusted_owner=(os.getuid(), os.getgid()))
            self._stack.callback(self._shim.close)
            require(self._shim.sha256 == final.identity["sha256"] and os.fstat(self._shim.fd).st_nlink == 1, "dispatcher output differs")
            descriptor = self._shim.executable_descriptor
            os.mkdir("helpers", mode=0o700, dir_fd=self._root_fd)
            prefix_fd = os.open("helpers", os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=self._root_fd)
            self._stack.callback(os.close, prefix_fd)
            self._prefix_fd = prefix_fd
            mappings = {role: fd_path(descriptor) for role in self.bindings["helpers"]}
            for role, target in mappings.items(): os.symlink(target, role, dir_fd=prefix_fd)
            os.fchmod(prefix_fd, 0o500)
            self._guard.verify()
            self._guard.add_directory_path(self.root / "helpers")
            self._guard._add_watch(self.root / "helpers", self._guard._DIRECTORY_MASK | provenance.IN_ONLYDIR | provenance.IN_ATTRIB, None)
            self.prefix = self._stack.enter_context(provenance._retain_exact_symlink_directory(self.root / "helpers", mappings, "v19 static helper prefix"))
            require(provenance._stable_fields(os.fstat(prefix_fd)) == provenance._stable_fields(os.fstat(self.prefix.descriptor)), "dispatcher prefix was replaced")
            self._state = "held"
            self.validate_current()
            return self
        except BaseException:
            self._state = "failed"
            self._stack.close()
            raise

    def validate_current(self):
        require(self._state == "held", "dispatcher is not held")
        try: self._validate_inputs()
        except BaseException:
            self._state = "failed"
            raise

    def run(self, argv):
        self.validate_current()
        try:
            require(type(argv) is list and len(argv) <= 511, "dispatcher argument count exceeds bound")
            selected = argv if self.headers is None else self.headers.arguments(argv)
            effective = compiler.driver_arguments(selected, self.phase.logical_driver, self.prefix.descriptor)
            require(len(effective) <= 512, "effective dispatcher argument count exceeds bound")
            return self._invoke("driver", effective, logical_argv=argv)
        except BaseException:
            self._state = "failed"
            raise

    def record(self):
        self.validate_current()
        return copy.deepcopy({"schema": "leopard2-v19-runtime-dispatch/v1", "root": str(self.root), "bindings": self.bindings,
            "helper_prefix": str(self.root / "helpers"), "commands": self._commands,
            "artifacts": {path: {"sha256": snap.identity["sha256"], "size": len(snap.content)} for path, snap in self._snapshots.items()},
            "sealed_dispatcher": self._shim.executable_record(), "compiler_data_owned": False,
            "headers": None if self.headers is None else self.headers.record(),
            "full_runtime_execution_owned": False, "fresh_build_recipe_integrated": False,
            "atomic_snapshot": False, "live_acquisition_armed": False, "benchmark_executed": False})

    def __exit__(self, kind, value, traceback):
        try:
            if self._state == "held": self.validate_current()
            elif value is None: raise builder.host.PreflightError("failed dispatcher cannot complete")
        finally:
            self._state = "closed"
            self._stack.__exit__(kind, value, traceback)
