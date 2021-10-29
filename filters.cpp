#include <vector>
#include <utility>
#include <iostream>

//applies a boxcar filter to data
std::vector<std::pair<double, double>> boxcarFilter(std::vector<std::pair<double, double>> data, int filterSize)
{
  int n = data.size();
  std::vector<std::pair<double, double>> result;
  result.reserve(n);
  for(int i = 0; i < n; i++)
  {
    int startIndex = i - (filterSize-1)/2;
    int endIndex = i + (filterSize-1)/2;
    double sum = 0;
    for(int j = startIndex; j <= endIndex; j++)
    {
      sum += data[(j+n)%n].second;
    }
    result.push_back({data[i].first, sum/filterSize});
  }
  return result;
}

//applies a boxcar filter to data for multiple passes
std::vector<std::pair<double, double>> boxcarFilter(std::vector<std::pair<double, double>> data, int filterSize, int numPasses)
{
  if(filterSize >= data.size())
  {
    std::cerr << "Error: number of filter points must be less than the number of data points" << std::endl;
    exit(1);
  }

  for(int i = 0; i < numPasses; i++)
    data = boxcarFilter(data, filterSize);

  return data;
}


//applies a Savitzky-Golay filter to the data
//filterSize must be 5, 11, or 17
std::vector<std::pair<double, double>> savitzkyGolayFilter(std::vector<std::pair<double, double>> data, int filterSize)
{
  int n = data.size();
  std::vector<std::pair<double, double>> result;
  result.reserve(n);

  std::vector<int> coefficients;   //the convoluting coefficients
  int norm; //the normalizing factor

  //these values are from Table 1 of the 1964 Savitzky-Golay paper
  switch (filterSize)
  {
    case 5:
      coefficients = {-3, 12, 17, 12, -3};
      norm = 35;
      break;
    case 11:
      coefficients = {-36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36};
      norm = 429;
      break;
    case 17:
      coefficients = {-21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21};
      norm = 323;
      break;
  }

  for(int i = 1+filterSize/2; i < n-1-filterSize/2; i++)
  {
    double sum = 0;
    for(int j = 0; j < filterSize; j++)
    {
      sum += coefficients[j] * data[i+j-filterSize/2].second;
    }
    result.push_back({data[i].first, sum/norm});
  }
  return result;
}

//applies a Savitzky-Golay filter to data for multiple passes
std::vector<std::pair<double, double>> savitzkyGolayFilter(std::vector<std::pair<double, double>> data, int filterSize, int numPasses)
{
  if(filterSize != 5 && filterSize != 11 && filterSize != 17)
  {
    std::cerr << "Error: Savitzky-Golay filter must have a size of 5, 11, or 17." << std::endl;
    exit(1);
  }

  for(int i = 0; i < numPasses; i++)
    data = savitzkyGolayFilter(data, filterSize);

  return data;
}

//filters the data according to the options specified
std::vector<std::pair<double, double>> filter(std::vector<std::pair<double, double>> data, int filterType, int filterSize, int numPasses)
{
  if(filterSize % 2 == 0)
  {
    std::cerr << "Error: filter size must be odd." << std::endl;
    exit(1);
  }

  switch (filterType)
  {
    case 0: //no filter
      return data;
    case 1: //boxcar
      return boxcarFilter(data, filterSize, numPasses);
    case 2: //Savitzky-Golay
      return savitzkyGolayFilter(data, filterSize, numPasses);
    default:
      std::cerr << "Error: filter type " << filterType << " is not valid." << std::endl;
      exit(1);
  }
}
