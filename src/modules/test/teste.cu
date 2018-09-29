#include <multivector.cu>
#include <operations.cu>
#include <metric.cu>
// #include <multivector.cu>

int main() {

	Multivector::set_N(200);
	auto metric = EuclideanMetric();
	generate_T<EuclideanMetric>(metric);

	// auto t1 = T(e(1));
	// auto t2 = T(e(2));

	// std::cout << (~f(1)).to_string() << std::endl;
	// auto k = e(1);
	auto k = MultivectorOperations::ADD(e(1), MultivectorOperations::OP(e(1), e(100)));
	std::cout << k.to_string() << std::endl;
	std::cout << (~k).to_string() << std::endl;
	std::cout << k.to_string() << std::endl;

    return 0;
}
