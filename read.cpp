//functions for reading in files
#include "structs.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>

//reads in the configuration file and returns all the options
configuration readConfig(std::string fileName)
{
  //data type we're going to return
  configuration result;

  auto max = std::numeric_limits<std::streamsize>::max();

  std::ifstream configFile(fileName);

  configFile >> result.inputFile;
  configFile.ignore(max, '\n'); //ignore the rest of the line
  configFile >> result.baseline;
  configFile.ignore(max, '\n');
  configFile >> result.tolerance;
  configFile.ignore(max, '\n');
  configFile >> result.filterType;
  configFile.ignore(max, '\n');
  configFile >> result.filterSize;
  configFile.ignore(max, '\n');
  configFile >> result.numPasses;
  configFile.ignore(max, '\n');
  configFile >> result.integrationTechnique;
  configFile.ignore(max, '\n');
  configFile >> result.outputFile;

  if(!configFile)
  {
    std::cerr << "Error reading in configuration file " << fileName << std::endl;
    exit(1);
  }

  return result;
}


//reads in the data and stores it in a vector of pairs
std::vector<std::pair<double, double>> readData(std::string fileName)
{
  std::ifstream file(fileName.c_str());
  if(!file)
  {
    std::cerr << "Error reading data file" << std::endl;
    exit(1);
  }

  std::vector<std::pair<double, double>> data;
  double x, y;
  while(file >> x >> y)
  {
    data.push_back({x,y});
  }
  return data;
}
