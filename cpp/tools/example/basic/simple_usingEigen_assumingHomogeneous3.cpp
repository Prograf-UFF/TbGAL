#include "../../../include/tbgal/using_Eigen.hpp" //TODO [DEBUG]
#include "../../../include/tbgal/assuming_Homogeneous3.hpp" //TODO [DEBUG]

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

    return EXIT_SUCCESS;
}
