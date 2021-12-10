#include "CubicSpline.h"
#include "structs.h"
#include "prototypes.h"

int main()
{
  auto startTime = std::chrono::high_resolution_clock::now(); //start timer
  auto config = readConfig("nmr.in"); //read in "nmr.in"
  auto data = readData(config.inputFile);   //read in the nmr data
  std::sort(data.rbegin(), data.rend());  //sort the data from most positive to most negative
  double shift = 0;
  data = baselineAdjustment(data, config.baseline, shift); //shift the data based on TMS and baseline
  data = filter(data, config.filterType, config.filterSize, config.numPasses);
  CubicSpline spline(data); //construct a cubic spline from the data
  auto peaks = calculatePeaks(spline, config.integrationTechnique, config.tolerance); //calculate the peak values

  auto endTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = endTime - startTime; //calculate elapsed time

  outputResult(peaks, config, shift, runtime.count());
  return 0;
}
