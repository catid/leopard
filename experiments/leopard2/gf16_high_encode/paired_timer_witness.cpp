// Reuse the qualified public-call witness; expose its count to the fake clock.
#include "paired_public_witness.cpp"
extern "C" unsigned LeoPairedWitnessCalls() { return witness.calls; }
