# Low encoder coefficient-copy A/B evidence

`run_abba.py` is the authoritative evidence collector for the LOW_V1 encoder
change between control commit `4070e4e527935026fb87593567587558f0a08d51`
and candidate commit `6d3afee213b94d486cf5f8145ac18078883ebc20`.
It does not build either side.  It accepts only the clean, tests-disabled,
Release production builds described in the reproduction guide, verifies their
source/object/archive/link/runtime closure, and then runs the fixed matrix.

The evidence exit statuses are deliberately distinct:

- `0`: valid evidence and every encode promotion/no-regression gate passed;
- `2`: valid, reviewable evidence whose performance policy failed;
- `1`: malformed input, invalid evidence, unsafe isolation, child failure, or
  verified failed-run diagnostics (never a performance result).

The `self-test` command uses only `mock_benchmark.py`.  It exercises mock
A1/B1/B2/A2 collection, raw-sample validation, pass and policy-failure replay,
bounded JSON, mutation rejection, reservation replacement, pair-lease
interoperability, and no-replace publication.  It does not execute either real
benchmark binary and must not be cited as timing evidence.

See `docs/reproduction/leopard2_low_encode_copy.md` for the exact build,
reservation, collection, and replay commands.
