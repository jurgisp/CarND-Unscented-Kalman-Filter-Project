#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    auto n = estimations.size();
    VectorXd sum = VectorXd::Zero(4);
    for (size_t i = 0; i < n; i++) {
        VectorXd diff = estimations.at(i) - ground_truth.at(i);
        sum = sum + diff.cwiseProduct(diff);
    }

    sum = sum / n;
    sum = sum.cwiseSqrt();
    return sum;
}