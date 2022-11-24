#pragma once
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#endif
#include <iostream>
#include <fstream>
#include <string>

void readParameter();
struct Parameter {
    int iterations;
    double  dt, Kb, Kv, Ka,gamma,tolerance;
    float reduced_volume;
    std::string meshFile,outFile,resFile;
};
