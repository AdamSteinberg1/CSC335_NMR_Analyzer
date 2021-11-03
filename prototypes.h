//All the prototypes for functions used in main
#pragma once

configuration readConfig(std::string fileName);
std::vector<std::pair<double, double>> filter(std::vector<std::pair<double, double>> data, int filterType, int filterSize, int numPasses);
std::vector<std::pair<double, double>> readData(std::string fileName);
std::vector<std::pair<double, double>> baselineAdjustment(std::vector<std::pair<double, double>> data, double baseline, double& shift);
std::vector<peak> calculatePeaks(CubicSpline c, int integrationTechnique, double tolerance);
void outputResult(std::vector<peak> peaks, configuration config, double shift, double runtime);
