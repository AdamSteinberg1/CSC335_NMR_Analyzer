#pragma once
#include <string>
struct configuration
{
  std::string inputFile, outputFile;
  double baseline, tolerance;
  int filterType, filterSize, numPasses, integrationTechnique;
};

struct peak
{
  double begin, end, location, area;
  int numHydrogens;
};
