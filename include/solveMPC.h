#pragma once

#include <chrono>
#include <sstream>

#include "OsqpEigen/OsqpEigen.h"

using namespace std::chrono;

std::stringstream ss;
bool verbose = false;

int count;

// // controller input and QPSolution vector
// Eigen::Matrix<double, 1, 1> U;
// Eigen::VectorXd QPSolution;

// set timinig variables
auto start = high_resolution_clock::now(), stop = high_resolution_clock::now();