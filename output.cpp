#include "structs.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <chrono>

std::string printOptions(configuration config, double shift)
{
  std::stringstream out;
  out << "Program Options" << std::endl;
  out << "===============================" << std::endl;
  out << "Baseline Adjustment\t:\t" << config.baseline << std::endl;
  out << "Tolerance\t\t:\t" << config.tolerance << std::endl;
  switch(config.filterType)
  {
    case 0:
      out << "No Filtering" << std::endl;
      break;
    case 1:
      out << "Boxcar Filtering" << std::endl;
      out << "Boxcar Size (Cyclic)\t:\t" << config.filterSize << std::endl;
      out << "Boxcar Passes\t\t:\t" << config.numPasses << std::endl;
      break;
    case 2:
      out << "Savitzky-Golay Filtering" << std::endl;
      out << "SG Filter Size\t\t:\t" << config.filterSize << std::endl;
      out << "SG Filter Passes\t:\t" << config.numPasses << std::endl;
      break;
  }
  out << std::endl;
  out << "Integration Method" << std::endl;
  out << "===============================" << std::endl;
  const std::string methods[] = {"Newton-Cotes", "Romberg", "Adaptive Quadrature", "Quadrature"};
  out << methods[config.integrationTechnique] << std::endl << std::endl;
  out << "Plot File Data" << std::endl;
  out << "===============================" << std::endl;
  out << "File:\t" << config.inputFile << std::endl;
  out << "Plot shifted " << shift << " ppm for TMS calibration" << std::endl;
  out << std::endl << std::endl << std::endl << std::endl;
  return out.str();
}

std::string printPeaks(std::vector<peak> peaks)
{
  std::stringstream out;
  out.precision(10);
  int width = 17;
  out <<  std::left << std::setw(8) << "Peak" << std::setw(width) << "Begin" << std::setw(width) << "End" << std::setw(width) << "Location" << std::setw(width) << "Area" << std::setw(width) << "Hydrogens"<< std::endl;
  out << "======= ================ ================ ================ ================ ================"<< std::endl;
  int i = 1;
  width = 16;
  out << std::right;
  for(auto & p: peaks)
  {
    out << std::setw(7) << i++ << " ";
    out << std::setw(width) << p.begin << " ";
    out << std::setw(width) << p.end << " ";
    out << std::setw(width) << p.location << " ";
    out << std::setw(width) << p.area << " ";
    out << std::setw(width) << p.numHydrogens << std::endl;
  }
  out << std::endl;
  return out.str();
}

void outputResult(std::vector<peak> peaks, configuration config, double shift, double runtime)
{
  //a stringstream that will be our formatted result
  std::stringstream out;
  out << "                              -=> NMR ANALYSIS <=-" << std::endl << std::endl << std::endl;
  out << printOptions(config, shift);
  out << printPeaks(peaks);
  out << "Analysis took " << runtime << " seconds." << std::endl;


  //print out to stdout
  std::cout << out.str();

  //print out to an output text file
  std::ofstream outFile(config.outputFile.c_str());
  outFile << out.str();
  outFile << "test";
  outFile.close();
}
