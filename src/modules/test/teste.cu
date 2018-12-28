#include <multivector.cu>
#include <operations.cu>
#include <metric.cu>
#include <math.h>
#include "benchmark_utils.cu"
#include <algorithm>
#include <iterator>
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;
using namespace MultivectorOperations;

#define PI acos(-1)

void run_loops(std::vector<Multivector> &R, std::vector<Multivector> &S, const int &RUNS, double &avg_time, double &std_dev) {
	std::vector<double> times(RUNS);
	for (int i = 0; i < RUNS; i++) {
		auto start = std::chrono::high_resolution_clock::now();
		GP(R, S);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> elapsed_seconds = finish - start;
		times[i] = elapsed_seconds.count();
		avg_time += times[i];
	}
	avg_time /= RUNS;
	for (int i = 0; i < RUNS; i++) {
		std_dev += (times[i] - avg_time)*(times[i] - avg_time);
	}
	std_dev /= (RUNS - 1);
	std_dev = sqrt(std_dev);
}

void run_iterations(const int &DIMS, const int &RUNS) {
	std::vector<Multivector> R;
	std::vector<Multivector> S;
	auto results = ofstream("results.csv", std::ios::app);
	for (IndexType r = 3; r <= DIMS; r++) {
		process_lines(read_lines("R", DIMS, r), R);
		for (IndexType s = 3; s <= DIMS; s++) {
			process_lines(read_lines("S", DIMS, s), S);

			double avg_time = 0;
			double std_dev = 0;

			run_loops(R, S, RUNS, avg_time, std_dev);
			S.clear();
			S.resize(0);
			S.shrink_to_fit();
			results << std::to_string(DIMS) << "," << std::to_string(r) << "," << std::to_string(s) << "," << std::to_string(avg_time) << "," << std::to_string(std_dev) << std::endl;
		}
		R.clear();
		R.resize(0);
		R.shrink_to_fit();
	}
	results.close();
}

int main() {

	IndexType MAX_DIMS = 100;
	int RUNS = 10;
	ofstream results ("results.csv");
	results << "N,R,S,avg_time,std_dev_time" << std::endl;
	results.close();
	for (IndexType DIMS = 3; DIMS <= MAX_DIMS; DIMS++) {
		Multivector::set_N(DIMS);
		generate_T<EuclideanMetric>(EuclideanMetric());
		run_iterations(DIMS, RUNS);
		delete_T();
	}


	// auto M = (e(1)^e(2)) + (e(1)^e(3));
	// auto m = MultivectorOperations::FACT_BLADE<std::vector<Multivector>>(M);

	// std::vector<string> lines = read_lines("R", 150, 6);


	// for (auto &i : k) {
	// 	std::cout << i << std::endl;
	// }
	// std::copy(v.begin(), v.end(), ostream_iterator<Multivector>(std::cout, "\n"));


	// auto output = MultivectorOperations::GP(getElementFromContainer(m, 0), getElementFromContainer(m, 1));
	//
	// std::cout << "output" << output << std::endl;
	// std::cout << "M" << M << std::endl;
	// std::cout << "Is output equals to M ? " << (output == M) << std::endl;
	//
	// auto B = (e(1)^e(2))+(e(1)^e(3)+(e(2)^e(3)));
	// B = B * (1.0 / MultivectorOperations::NORM(B));
	//
	// auto V = cos(PI/4) - ((+sin(PI/4.0)*B));
	// auto v = MultivectorOperations::FACT_VERSOR<std::vector<Multivector>>(V);
	// output = MultivectorOperations::GP(MultivectorOperations::GP(getElementFromContainer(v, 0), getElementFromContainer(v, 1)), getElementFromContainer(v, 2));
	// std::cout << "output" << output << std::endl;
	// std::cout << "M" << V << std::endl;
	// std::cout << "Is output equals to V ? " << (output == V) << std::endl;



    return 0;
}
