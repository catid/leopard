# 26 — LOW side-256 schedule-shape census

**Disposition: INCONCLUSIVE — no timing campaign armed**

The earlier LOW side-256 diagnostic measured a 1.185x decode ratio for one
nine-erasure mask.  Before spending another authoritative timing attempt, this
checkpoint sampled five deterministic mask shapes (clustered, spread,
boundary-heavy, interior, and gap-extreme) at `K=200`, `R=800`, GF16/AVX2,
64 KiB shards.  The benchmark was run once per case only to expose immutable
plan metadata; these samples are not performance measurements.

| Shape | Losses | Pruned operations | Full butterflies | Fused-4 groups |
| --- | ---: | ---: | ---: | ---: |
| clustered | 4 / 8 / 9 / 16 | 1568 / 1586 / 1596 / 1706 | 3072 | 207 / 205 / 204 / 205 |
| spread | 4 / 8 / 9 / 16 | 1604 / 1697 / 1735 / 1823 | 3072 | 211 / 266 / 265 / 265 |
| boundary-heavy | 4 / 8 / 9 / 16 | 1554 / 1582 / 1596 / 1798 | 3072 | 207 / 205 / 204 / 287 |
| interior | 4 / 8 / 9 / 16 | 1368 / 1666 / 1704 / 1738 | 3072 | 207 / 205 / 204 / 206 |
| gap-extreme | 4 / 8 / 9 / 16 | 1268 / 1556 / 1566 / 1824 | 3072 | 207 / 220 / 219 / 283 |

Every case selected `workspace_tiled` with 512 work slots.  Changing shard size
from 32 to 64 to 128 KiB changed only the aligned pass prefix and scratch
geometry (32, 64, and 128 KiB prefixes); it did not change the selected path or
plan shape.  The loss-9 operation range is 1566–1735 (10.8%), and the fused
group range is 204–265 (29.9%).  This is enough structural variation that the
single-mask decode gain cannot be promoted or pooled across patterns.

The side-256 candidate therefore remains diagnostic-only.  A future timing
campaign would need a fresh preregistration covering these finite shape classes
and neighboring loss/size/cache cells, with immutable binaries, the canonical
lock, zero-sibling checks, and the existing 5% gain / 2% regression gates.

Evidence was collected from the clean Release benchmark built from the current
master tree; raw one-sample JSON is retained in `/tmp/l52-census-cost.txt` for
the session and is intentionally not presented as timing evidence.
