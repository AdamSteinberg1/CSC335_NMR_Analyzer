#include "CubicSpline.h"
#include "structs.h"
#include "prototypes.h"

int main()
{
  auto startTime = std::chrono::high_resolution_clock::now();
  auto config = readConfig("nmr.in");
  auto data = readData(config.inputFile);
  double shift = 0;
  data = baselineAdjustment(data, config.baseline, shift);
  data = filter(data, config.filterType, config.filterSize, config.numPasses);
  CubicSpline spline(data);
  auto peaks = calculatePeaks(spline, config.integrationTechnique, config.tolerance);

  auto endTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = endTime - startTime;

  outputResult(peaks, config, shift, runtime.count());
  return 0;
}
