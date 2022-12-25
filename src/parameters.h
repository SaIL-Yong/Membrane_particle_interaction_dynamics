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
    float reduced_volume,particle_radious,adhesion_strength,potential_range;
    std::string meshFile,outFile,resFile;
};
