#include "structs.h"
#include "CubicSpline.h"
#include "legendreConstants.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <gsl/gsl_poly.h>

#define MAX_ITERATIONS 10000

//finds all the x values on the interval (a,b] where p(x) = 0
//p is assumed to be a cubic polynomial
std::vector<double> findRoots(Polynomial p, double start, double end)
{
  //consider a cubic of the form x^3+ax^2+bx+c
  double a = p[2]/p[3];
  double b = p[1]/p[3];
  double c = p[0]/p[3];

  double x0, x1, x2;
  int numRoots =  gsl_poly_solve_cubic(a, b, c, &x0, &x1, &x2);

  std::vector<double> roots;
  switch (numRoots)
  {
    case 3:
      if(start<x2 && x2<=end)
        roots.push_back(x2);
    case 2:
      if(start<x1 && x1<=end)
        roots.push_back(x1);
    case 1:
      if(start<x0 && x0<=end)
        roots.push_back(x0);
      break;
  }

  //it is necessary for the roots to be sorted for when we iterate through them later
  std::sort(roots.begin(), roots.end());

  return roots;
}

//finds all the x-values at which the cubic spline intersects the x-axis
std::vector<double> findRoots(CubicSpline spline)
{
  std::vector<double> roots;
  //iterate through each cubic of the spline and accumulate their roots
  for(int i = 0; i < spline.getNumCubics(); i++)
  {
    Polynomial cubic = spline[i];
    std::pair<double,double> range = spline.getRange(i);
    std::vector<double> intermediateRoots = findRoots(cubic, range.first, range.second);
    roots.insert(roots.end(), intermediateRoots.begin(), intermediateRoots.end());
  }
  return roots;
}

//integrates f from a to b using Newton-Cotes
double newtonCotes(std::function<double(double)> f, double a, double b)
{
  //uses closed Newton-Cotes with n=4
  double h = (b-a)/4;
  double x0 = a;
  double x1 = a+h;
  double x2 = x1+h;
  double x3 = x2+h;
  double x4 = b;

  return (2*h/45)*(7*f(x0) + 32*f(x1) + 12*f(x2) + 32*f(x3) + 7*f(x4));
}

//performs Romberg integration over f from a to b
//computes until the error is less than tolerance or until MAX_ITERATIONS is exceeded, whichever comes first
double romberg(std::function<double(double)> f, double a, double b, double tolerance)
{
  double h = b-a;
  //we only keep two rows of the extrapolation table in memory at a time
  std::vector<double> currRow, lastRow;
  lastRow.push_back(0.5*h*(f(a)+f(b))); //R_1,1
  for(int i = 2; i <= MAX_ITERATIONS; i++)
  {
    currRow.clear();
    double sum = 0;
    for(int k = 1; k <= pow(2,i-2); k++)
      sum += f(a+(k-0.5)*h);
    currRow.push_back(0.5*(lastRow[0] + h*sum));
    for(int j = 1; j < i; j++)
      currRow.push_back(currRow[j-1] + (currRow[j-1]-lastRow[j-1])/(pow(4,j)-1));
    h *= 0.5;
    if(fabs(currRow.back() - lastRow.back()) < tolerance)
    {
      return currRow.back();
    }
    lastRow = currRow;
  }
  return currRow.back();
}

//integrates from a to b until error is less than tolerance
//performs adaptive quadrature with simpson's rule
//this version is slower, but simpler
//because it recalculates values that are used multiple times
double adaptiveQuadSlow(std::function<double(double)> f, double a, double b, double tolerance)
{
  //integrates f from x0 to x1 using simpson's rule
  std::function<double(double, double)> simpson = [&](double x0, double x1){
    double h = (x1-x0)/2;
    return (h/3)*(f(x0) + 4*f(x0+h) +f(x1));
  };

  double mid = (a+b)/2;
  double whole = simpson(a, b);
  double left = simpson(a, mid);
  double right = simpson(mid, b);
  double diff = left+right-whole;
  if(fabs(diff) < 10*tolerance) //10*tolerance is used because that's what is used in the Burden text
    return left + right;

  return adaptiveQuadSlow(f, a, mid, tolerance/2) + adaptiveQuadSlow(f, mid, b, tolerance/2);
}


//integrates from a to b until error is less than tolerance
//performs adaptive quadrature with simpson's rule
double adaptiveQuad(std::function<double(double)> f, double a, double b, double tolerance)
{
  //integrates f from x0 to x1 using simpson's rule
  //f_x0 is f evaluated at x0
  //f_x1 is f  evaluated at x1
  auto simpson = [&](double x0, double x1, double f_x0, double f_x1){
    double mid = (x0 + x1)/2;
    double f_mid = f(mid);
    double result =  (fabs(x1-x0)/6)*(f_x0 + 4*f_mid + f_x1);
    return std::make_tuple(mid, f_mid, result);
  };

  //performs adaptive quadrature recursively
  //function evaluations and intermediate integrals are passed through arguments to avoid recalcuating values
  std::function<double (double, double, double, double, double, double, double, double, int)> helper =
  [&](double x0, double x1, double f_x0, double f_x1, double tol, double whole, double mid, double f_mid, int currDepth){
    std::cout << "S(" << x0 << "," << x1 <<") tol = " << tol << " depth = " << currDepth << std::endl;
    if(currDepth > MAX_ITERATIONS)
    {
      std::cout << "Adaptive quadrature failed to reach tolerance within " << MAX_ITERATIONS << " iterations. Using best approximation obtained." << std::endl;
      return whole;
    }
    double left_mid, f_left_mid, left, right_mid, f_right_mid, right;
    std::tie(left_mid, f_left_mid, left) = simpson(x0, mid, f_x0, f_mid);
    std::tie(right_mid, f_right_mid, right) = simpson(mid, x1, f_mid, f_x1);
    double diff = left + right - whole;
    if (fabs(diff) < 10*tol)
      return left + right;
    return helper(x0, mid, f_x0, f_mid, tol/2, left, left_mid, f_left_mid, currDepth+1)
         + helper(mid, x1, f_mid, f_x1, tol/2, right, right_mid, f_right_mid, currDepth+1);
  };
  std::cout.precision(15);
  double f_a = f(a);
  double f_b = f(b);
  double mid, f_mid, whole;
  std::tie(mid, f_mid, whole) = simpson(a, b, f_a, f_b);
  return helper(a, b, f_a, f_b, tolerance, whole, mid, f_mid, 1);
}

//integrates f from a to b using Gaussian Quadrature with n=512
double gaussQuad(std::function<double(double)> f, double a, double b)
{
    //change of variable from x to t so we can integrate from -1 to 1
    std::function<double(double)> g = [&](double t){
      return f(((b-a)*t+b+a)/2)*(b-a)/2;
    };
    //the arrays coeff and roots are included from "legendreConstants.h"
    double sum = 0;
    for(int i = 0; i < 512; i++)
      sum += coeff[i]*g(roots[i]);
    return sum;
}

//integrates a cubic spline from a to b using the specified integration technique
double integrate(double a, double b, CubicSpline spline, int integrationTechnique, double tolerance)
{
  auto f = [&](double x) { return spline.evaluate(x); };  //lambda for evaluating the spline
  switch (integrationTechnique)
  {
    case 0: //Newton-Cotes
      return newtonCotes(f, a, b);
      break;
    case 1: //Romberg
      return romberg(f, a, b, tolerance);
      break;
    case 2: //Adaptive
      return adaptiveQuad(f, a, b, tolerance);
      break;
    case 3: //Quadrature
      return gaussQuad(f, a, b);
      break;
    default:
      std::cerr << "Error: integration technique " << integrationTechnique << " is not a valid option" << std::endl;
      exit(1);
      break;
  }
}

std::vector<peak> calculatePeaks(CubicSpline spline, int integrationTechnique, double tolerance)
{
  //find all the points that the cubic spline intersects the x-axis
  std::vector<double> roots = findRoots(spline);

  std::vector<peak> peaks; //what we will return
  peaks.reserve(roots.size()/2);
  //each pair of roots will enclose a peak
  for(int i=0; i < roots.size(); i+=2)
  {
    peak p;
    p.begin = roots[i];
    p.end = roots[i+1];
    p.location = (p.begin + p.end)/2;
    peaks.push_back(p);
  }

  double minArea = std::numeric_limits<double>::infinity(); //need a value that is bigger than all other values
  //calculate the area of each peak
  for(peak & p : peaks)
  {
    p.area = integrate(p.begin, p.end, spline, integrationTechnique, tolerance);
    minArea = std::min(p.area, minArea); //find the smallest area
  }

  //calculate the number of hydrogens each peak represents
  for(peak & p : peaks)
  {
    p.numHydrogens = int(std::round(p.area/minArea));
  }


  return peaks;
}
