#include <vector>
#include <utility>

//shifts all the data so that the TMS peak is at x=0
std::vector<std::pair<double, double>> baselineAdjustment(std::vector<std::pair<double, double>> data, double baseline, double& shift)
{

  //find the TMS peak
  //assumes the data is sorted from greatest to least by x-value
  for(auto & point : data)
  {
    if(point.second >= baseline)
    {
      shift = point.first;
      break;
    }
  }

  //shift all the data horiztonally so the TMS peak is at x=0
  //additionally shift all data down so that the baseline is at y=0
  //we do this so that the integrals will be the area between the spline and the baseline
  for(auto & point : data)
  {
    point.first -= shift;
    point.second -= baseline;
  }

  return data;
}
