"""Independent structural counts for current production; not CPU-time attribution."""
from verify_gf16_callback_probe import validate_shape_counts

SHAPES = ((1000,200,32768,3,32768,1),(1000,199,65536,3,32768,2),
          (1000,200,65536,3,32768,2),(4096,512,4096,3,4096,1),
          (1000,199,32768,3,32768,1),(1000,200,32768,6,32768,1),
          (1000,200,32768,6,32768,1))


def validate(record, cell):
    if type(cell) is not int or not 0 <= cell <= 7:
        raise ValueError('cell')
    if cell == 7:
        if record != dict(schema='gf16-callback-counts/v1',timed=False,calls=0,passes=[],buckets=[]):
            raise ValueError('GF8 must bypass the GF16 observer')
        # Explicit type checks reject False masquerading as a numeric count.
        if type(record['calls']) is not int or type(record['timed']) is not bool:
            raise ValueError('GF8 record types')
        return dict(callbacks=0,operations={},affected=False,forward_pairs=0,
                    forward_blocks=0,accumulating_pairs=0,accumulating_blocks=0)
    operations = validate_shape_counts(record, SHAPES[cell])
    affected = SHAPES[cell][3] == 3
    forward = accumulating = forward_blocks = accumulating_blocks = 0
    if affected:
        for bucket in record['buckets']:
            op = bucket['op']
            pairs = 0
            if op == 'fft2':
                if bucket['zero_mask']:
                    raise ValueError('ordinary pair cannot receive zero sentinel')
                pairs = bucket['calls']
            elif op == 'fft4_range':
                if bucket['prefer_fused']:
                    raise ValueError('unexpected fused AVX2 range')
                zero = bucket['zero_mask']
                edges = 4 - int(bool(zero & 1)) - int(bool(zero & 2)) - 2 * int(bool(zero & 4))
                pairs = edges * bucket['distance'] * bucket['calls']
            forward += pairs
            forward_blocks += pairs * (bucket['bytes'] // 64)
            if op == 'ifft2_xor':
                if bucket['zero_mask']:
                    raise ValueError('ordinary accumulating pair cannot receive zero sentinel')
                accumulating += bucket['calls']
                accumulating_blocks += bucket['calls'] * (bucket['bytes'] // 64)
    return dict(callbacks=record['calls'],operations=operations,affected=affected,
                forward_pairs=forward,forward_blocks=forward_blocks,
                accumulating_pairs=accumulating,accumulating_blocks=accumulating_blocks)
