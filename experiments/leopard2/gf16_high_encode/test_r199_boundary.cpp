// Clock-free exact-cell/boundary/both-field qualification; leopard-79h.38.5.4.19.
#define LEO_BOUNDARY_GUARD_NO_MAIN 1
#include "test_gfni_boundary.cpp"

int main(int argc, char** argv)
{
    try
    {
        const bool original = argc == 3 && !std::strcmp(argv[1], "--original");
        Require(argc == 2 || original, "usage: test_r199_boundary [--original] cell[0..7]");
        const char* selector = argv[original ? 2 : 1];
        Require(std::strlen(selector) == 1 && selector[0] >= '0' && selector[0] <= '7', "cell");
        const unsigned index = static_cast<unsigned>(selector[0] - '0');
        if (original)
        {
            Check(index);
            return 0;
        }
        const size_t bytes[] = {32768, 32768, 32770, 32766, 64, 66, 65, 66};
        const size_t offsets[] = {0, 1, 2, 1, 1, 2, 1, 1};
        CheckShape(index, index < 6 ? 1000 : 17, index < 6 ? 199 : 7,
                   bytes[index], offsets[index],
                   index == 6 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16,
                   LEO2_BACKEND_GFNI, LEO2_BACKEND_AUTO);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
