#include "CubicSpline.h"
#include "structs.h"
#include <fstream>
#include <sstream>
#include <string>
#include <limits>

std::string gnuPrint(Polynomial p)
{
  std::stringstream result;
  result.precision(std::numeric_limits<double>::max_digits10);
  int n = p.getDegree();
  for (int i =0; i <= n; i++)
  {
    if(p[i] == 0)
      continue;

    result << p[i];
    if (i != 0)
      result << "*x";

    if(i > 1)
      result << "**" << i;

    if(i < n)
      result << " + ";
  }
  return result.str();
}

//returns a string representing a cubic spline
std::string gnuPrint(CubicSpline c)
{
  std::stringstream result;
  for(int i = 0; i < c.cubics.size()-1; i++)
  {
    result << "x<" << c.xValues[i+1] << " ? " << gnuPrint(c.cubics[i]) << " : ";
  }
  result << gnuPrint(c.cubics.back());

  return result.str();
}


void graph(CubicSpline spline, std::vector<std::pair<double,double>> points)
{
  std::ofstream script("tmp.plt");

  script  << "set terminal pngcairo\n"
          << "set output 'graph.png'\n"
          //<< "set xrange [-6:1]\n"
          //<< "set yrange [-5000:25000]\n"
          //<< "set xrange [-3.6:-3.3]\n"
          //<< "set xrange [-2:-1.96]\n"
          //<< "set yrange [-800:1200]\n"
          //<< "set yrange [-500:500]\n"
          << "set samples 10000\n"
          << "p(x) = " << gnuPrint(spline) << "\n"
          << "plot p(x) title 'Spline', 0 title 'Baseline', '-' notitle\n";

  for(auto & point: points)
  {
    script << point.first << '\t' <<point.second << std::endl;
  }

  system("gnuplot tmp.plt");
  system("rm tmp.plt");
}
