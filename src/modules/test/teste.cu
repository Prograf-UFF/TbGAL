#include <multivector.cu>
#include <operations.cu>
#include <metric.cu>
#include <iostream>

int main() {
    Multivector::set_N(200);
    auto metric = EuclideanMetric();
    generate_T<EuclideanMetric>(metric);
    auto *m = MultivectorOperations::GP(e(1), e(2));
    std::cout << m->to_string() << std::endl;
    return 0;
}
