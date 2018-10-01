#include <multivector.cu>
#include <operations.cu>
#include <metric.cu>
#include <math.h>

using namespace std;
using namespace MultivectorOperations;

#define PI acos(-1)

int main() {

	Multivector::set_N(5);
	auto metric = EuclideanMetric();
	generate_T<EuclideanMetric>(metric);

	auto M = (e(1)^e(2)) + (e(1)^e(3));
	auto m = MultivectorOperations::FACT_BLADE<std::vector<Multivector>>(M);


	auto output = MultivectorOperations::GP(getElementFromContainer(m, 0), getElementFromContainer(m, 1));

	std::cout << "output" << output << std::endl;
	std::cout << "M" << M << std::endl;
	std::cout << "Is output equals to M ? " << (output == M) << std::endl;

	auto B = (e(1)^e(2))+(e(1)^e(3)+(e(2)^e(3)));
	B = B * (1.0 / MultivectorOperations::NORM(B));

	auto V = cos(PI/4) - ((+sin(PI/4.0)*B));
	auto v = MultivectorOperations::FACT_VERSOR<std::vector<Multivector>>(V);
	output = MultivectorOperations::GP(MultivectorOperations::GP(getElementFromContainer(v, 0), getElementFromContainer(v, 1)), getElementFromContainer(v, 2));
	std::cout << "output" << output << std::endl;
	std::cout << "M" << V << std::endl;
	std::cout << "Is output equals to V ? " << (output == V) << std::endl;



    return 0;
}
