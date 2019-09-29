#include "../../../include/tbgal/using_Eigen.hpp" //TODO [DEBUG]
#include "../../../include/tbgal/assuming_Homogeneous3.hpp" //TODO [DEBUG]

#include <gatl/ga3h.hpp>
using namespace ga3h;

using namespace tbgal;

int main(int argc, char *argv[]) {
    double x, y, z;

    std::cout << "-- Input" << std::endl;
    std::cout << std::endl;
    std::cout << "x = "; std::cin >> x;
    std::cout << "y = "; std::cin >> y;
    std::cout << "z = "; std::cin >> z;
    std::cout << std::endl;

    auto p = Homogeneous3::point(x, y, z);
    auto d = Homogeneous3::direction(x, y, z);
    auto l = p ^ d;

    std::cout << "-- Result" << std::endl;
    std::cout << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "l = p ^ d = " << l << std::endl;
    std::cout << std::endl;

    {
        std::cout << "p = " << (p.scalar() ^ ga3h::vector(p.factors()(0, 0), p.factors()(1, 0), p.factors()(2, 0), p.factors()(3, 0))) << std::endl;
        std::cout << "d = " << (d.scalar() ^ ga3h::vector(d.factors()(0, 0), d.factors()(1, 0), d.factors()(2, 0), d.factors()(3, 0))) << std::endl;
        std::cout << "l = " << (l.scalar() ^ ga3h::vector(l.factors()(0, 0), l.factors()(1, 0), l.factors()(2, 0), l.factors()(3, 0)) ^ ga3h::vector(l.factors()(0, 1), l.factors()(1, 1), l.factors()(2, 1), l.factors()(3, 1))) << std::endl;

        std::cout << std::endl;
        std::cout << l.scalar() << std::endl; 
    }

    return EXIT_SUCCESS;
}
