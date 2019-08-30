#include <multivector.cu>
#include <operations.cu>
#include <metric.cu>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace MultivectorOperations;

vector<string> read_lines(const string &op, const int &N, const int &grade) {
  string filename = "/home/eduardovera/Workspace/TbGAl/src/modules/test/test_data/" + op + "/" + to_string(N) + ".csv";
  int start = (grade * (grade + 1) / 2) + grade;
  std::vector<string> lines;
  ifstream file(filename);
  if (file.is_open()) {
    int size = 0;
    int idx = 1;
    string line;
    while (idx != start) {
      getline(file, line);
      idx++;
    }
    size = stoi(line) + 1;
    while (idx != start + size) {
      getline(file, line);
      lines.push_back(line);
      idx++;
    }
  }
  file.close();
  return lines;
}

// void process_lines(const std::vector<string> &factors, vector<Multivector> &ret) {
//   for (int i = 1; i < factors.size(); i++) {
//     3*e(i);
//   }
// }

void process_lines(const std::vector<string> &factors, vector<Multivector> &ret) {
  for (const string &f : factors) {
    std::vector<string> terms;
    boost::split(terms, f, boost::is_any_of(","));
    Multivector mv = e(1)^e(1);
    for (string &t : terms) {
      std::vector<string> basis_coeff;
      boost::split(basis_coeff, t, boost::is_any_of(":"));
      int basis = stoi(basis_coeff[0]);
      double coeff = stod(basis_coeff[1]);
      mv = mv + (coeff * e(basis));
    }
    ret.push_back(mv);
  }
}
