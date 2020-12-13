#pragma once

#include <chrono>
#include <sstream>
#include <algorithm>

#include "OsqpEigen/OsqpEigen.h"

using namespace std::chrono;

std::stringstream ss;
bool verbose = false;

int count;

// set timinig variables
auto start = high_resolution_clock::now(), stop = high_resolution_clock::now();

/**
 * get command based on command line option
 */
char* getCmdOption(char ** begin, char ** end, const std::string & option);

/**
 * check if option exists
 */ 
bool cmdOptionExists(char** begin, char** end, const std::string& option);