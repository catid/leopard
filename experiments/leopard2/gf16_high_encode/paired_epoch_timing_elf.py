"""Variant-specific ELF oracle; does not relax the clock-free ELF verifier."""
import struct
from verify_paired_metadata import require

CLOCK = '_ZNSt6chrono3_V212steady_clock3nowEv'
VARIANTS = ('steady', 'abort', 'synthetic', 'fake-steady')


def bindings(defined, undefined, native, variant):
    require(variant in VARIANTS, 'clock variant')
    observed = {name for name in defined if name.startswith('__wrap_leo')}
    public = {'__wrap_leo_encode'} if native else {'__wrap_leo2_encode','__wrap_leo2_encode_batch'}
    witnessed = variant in ('synthetic','fake-steady')
    require(observed == (public if witnessed else set()), 'public wrapper inventory')
    require(('LeoPairedWitnessCalls' in defined) == witnessed, 'public call witness binding')
    require('LeoPairedWitnessMark' in defined and 'LeoPairedClockKind' in defined, 'clock/mark definitions')
    require((CLOCK in undefined) == (variant == 'steady'), 'real steady import binding')
    require(CLOCK not in defined, 'unexpected local steady clock implementation')
    require(('__wrap_'+CLOCK in defined) == (variant != 'steady'), 'clock wrapper binding')


def elf(path, native, variant):
    with path.open('rb') as stream:
        length = path.stat().st_size
        def read(offset, size):
            require(0 <= offset <= length and 0 <= size <= min(8*1024**2,length-offset), 'ELF range')
            stream.seek(offset); data = stream.read(size)
            require(len(data) == size, 'short ELF read')
            return data
        header = struct.unpack('<16sHHIQQQIHHHHHH', read(0,64))
        require(header[0][:7] == b'\x7fELF\x02\x01\x01' and header[1] in (2,3) and header[2] == 62,
                'ELF64 x86-64')
        phoff,shoff,phsize,phnum,shsize,shnum = header[5],header[6],header[9],header[10],header[11],header[12]
        require(phsize == 56 and shsize == 64 and 0 < phnum <= 64 and 0 < shnum <= 256, 'ELF tables')
        segments = []
        for i in range(phnum):
            p = struct.unpack('<IIQQQQQQ', read(phoff+i*phsize,phsize))
            if p[0] == 1:
                segments.append(dict(address=p[3],bytes=p[6],flags=p[1]))
        require(0 < len(segments) <= 16, 'ELF load count')
        sections = [struct.unpack('<IIQQQQIIQQ', read(shoff+i*shsize,shsize)) for i in range(shnum)]
        tables = [s for s in sections if s[1] == 2]
        require(len(tables) == 1, 'ELF static symbol table')
        table = tables[0]
        require(table[9] == 24 and table[5]%24 == 0 and table[6] < shnum, 'ELF symbols')
        strings = sections[table[6]]; require(strings[1] == 3, 'ELF string table')
        names,data = read(strings[4],strings[5]),read(table[4],table[5])
        functions,defined,undefined,addresses = {},set(),set(),{}
        desired = {'main':'main', 'leo_encode' if native else 'leo2_encode':'encode'}
        if not native: desired['leo2_encode_batch'] = 'batch'
        for i in range(0,len(data),24):
            name,info,_,index,value,size = struct.unpack_from('<IBBHQQ',data,i)
            require(name < len(names), 'ELF symbol name')
            stop = names.find(b'\0',name); require(stop >= name, 'ELF symbol terminator')
            name = names[name:stop].decode('ascii').split('@',1)[0]
            (defined if index else undefined).add(name)
            if index and (name in desired or name.startswith('__real_leo')):
                require(name not in addresses, 'duplicate public address symbol')
                addresses[name] = value
            label = desired.get(name)
            if 'paired_metadata' in name and name.endswith('6AnchorEv'): label = 'metadata_anchor'
            if label:
                require(label not in functions and info&15 == 2 and index != 0 and size > 0, 'function symbol')
                functions[label] = dict(value=value,size=size)
        require(set(functions) == set(desired.values())|{'metadata_anchor'}, 'function inventory')
        bindings(defined,undefined,native,variant)
        if variant in ('steady','abort'):
            for name in desired:
                if name != 'main':
                    require(addresses.get('__real_'+name) == addresses[name], 'plain metadata function alias')
        return dict(type=header[1],segments=segments,functions=functions)
