#include "../../../include/tbgal/using_Eigen.hpp" //TODO [DEBUG]
#include "../../../include/tbgal/assuming_Euclidean3.hpp" //TODO [DEBUG]

using namespace tbgal;
using namespace tbgal::Euclidean3;

int main(int argc, char *argv[]) {
    auto a = vector(0.5, 0.0, 0.5);    // a = 0.5 * (e1 + e3)

    auto M = 3.0 * (e1 ^ e2);          // M = 3.0 * e1^e2
    auto v = dual(M);                  // v = 3.0 * e3

    auto a_ = -gp(v, a, inverse(v));   // a_ = 0.5 * (e1 - e3)

    std::cout << "       a = " << a << std::endl;
    std::cout << "       M = " << M << std::endl;
    std::cout << "       v = " << v << std::endl;
    std::cout << "  inv(v) = " << inverse(v) << std::endl;
    std::cout << "gp(v, a) = " << gp(v, a) << std::endl;
    std::cout << "      a_ = " << a_ << std::endl; //TODO {DEBUG} There is a bug in the geometric product!

    return EXIT_SUCCESS;
}
