// Exercise the new frontend's actual Buffer implementation, never a codec.
#define main PairedUnusedDriverMain
#include "paired_r19932.cpp"
#undef main

int main(int argc, char** argv)
{
    try {
        Require(argc == 2, "guard test needs one case");
        Buffer empty(0), buffer(64);
        empty.Check(); buffer.Check();
        if (!std::strcmp(argv[1], "canary")) {
            // A normal boundary check must fail even without ASan.
            for (unsigned side = 0; side < 2; ++side) {
                buffer.Poison(false);
                uint8_t* byte = side ? buffer.data + buffer.bytes : buffer.raw;
                *byte ^= 1;
                bool refused = false;
                try { buffer.Check(); } catch (const std::runtime_error&) { refused = true; }
                Require(refused, "canary corruption accepted");
                buffer.Poison(false); *byte ^= 1; buffer.Poison(true);
                buffer.Check();
            }
            std::puts("both canaries rejected; zero-size and restored buffers pass");
            return 0;
        }
        Require(!std::strcmp(argv[1], "underflow") || !std::strcmp(argv[1], "overflow"), "guard test case");
        volatile uint8_t* data = buffer.data;
        const uint8_t value = !std::strcmp(argv[1], "underflow") ? data[-1] : data[buffer.bytes];
        std::printf("unexpected poisoned read: %u\n", static_cast<unsigned>(value));
        return 2;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
