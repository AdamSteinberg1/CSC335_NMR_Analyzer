//functions to calculate the peaks of the cubic spline
//finds their start and endpoints, their area, and their location
#include "structs.h" //peak struct is included here
#include "CubicSpline.h"
#include "legendreConstants.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <gsl/gsl_poly.h>

#define MAX_ITERATIONS 1000
#define MAX_RECURSION_DEPTH 10

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

//integrates f from a to b using composite Newton-Cotes
//performs n subdivisions. n must be even
double newtonCotes(std::function<double(double)> f, double a, double b, int n)
{
  //uses composite Newton-Cotes with Simpson's rule
  double h = (b-a)/n;
  double sum1 = 0;
  double sum2 = 0;
  for(int i = 1; i < n; i++)
  {
    double x = a + i*h;
    if(i%2==0)
      sum2 += f(x);
    else
      sum1 += f(x);
  }
  return h * (f(a) + 2*sum2 + 4*sum1 + f(b))/3;
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
      sum += f(a+(k-0.5)*h); //calculate value in first column of the extrapolation table
    currRow.push_back(0.5*(lastRow[0] + h*sum));
    for(int j = 1; j < i; j++)
      currRow.push_back(currRow[j-1] + (currRow[j-1]-lastRow[j-1])/(pow(4,j)-1)); //perform Richardson extrapolation
    h *= 0.5; //h halves for each row in the table
    if(fabs(currRow.back() - lastRow.back()) < tolerance) //estimate error and compare to tolerance
    {
      return currRow.back();
    }
    lastRow = currRow;
  }
  return currRow.back();
}

//recursively find the integral from a to b with simpson's method
//tolerance is halved at each recursive level
double adaptiveQuadhelper(std::function<double(double)> f, double a, double b, double tol, double whole, double f_a, double f_b, double f_mid, int recDepth) {
    double mid = (a + b)/2;
    double h = (b - a)/2;
    double left_mid  = (a + mid)/2;
    double right_mid  = (mid + b)/2;
    //check for signs it's not going to converge
    //if tolerance cannot be halved anymore
    //or if the average of a and m equals a => a and m are too close together
    if ((tol/2 == tol) || (a == left_mid))
      return whole;
    double f_left_mid = f(left_mid);
    double f_right_mid = f(right_mid);
    //simpson's method
    double left  = (h/6) * (f_a + 4*f_left_mid + f_mid);
    double right = (h/6) * (f_mid + 4*f_right_mid + f_b);
    double diff = whole - left - right;

    if (recDepth <= 0 || fabs(diff) <= 10*tol) //10*tolerance is used because that's what is used in the Burden text
        return left + right;
    return adaptiveQuadhelper(f, a, mid, tol/2, left,  f_a, f_mid, f_left_mid, recDepth-1) +
           adaptiveQuadhelper(f, mid, b, tol/2, right, f_mid, f_b, f_right_mid, recDepth-1);
}

//integrates from a to b until error is less than tolerance
//performs adaptive quadrature with simpson's rule
double adaptiveQuad(std::function<double(double)> f, double a, double b, double tolerance)
{
    if(a==b)
      return 0.0;
    double h = b - a;
    double f_a = f(a);
    double f_b = f(b);
    double f_m = f((a + b)/2);
    double simpsons = (h/6)*(f_a + 4*f_m + f_b);
    return adaptiveQuadhelper(f, a, b, tolerance, simpsons, f_a, f_b, f_m, MAX_RECURSION_DEPTH);
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
    case 0: //Adaptive
      return adaptiveQuad(f, a, b, tolerance);
      break;
    case 1: //Romberg
      return romberg(f, a, b, tolerance);
      break;
    case 2: //Composite Newton-Cotes with 20 subintervals
      return newtonCotes(f, a, b, 20);
      break;
    case 3: //Gaussian Quadrature
      return gaussQuad(f, a, b);
      break;
    default:
      std::cerr << "Error: integration technique " << integrationTechnique << " is not a valid option" << std::endl;
      exit(1);
      break;
  }
}

//calculate a vector of peak structs
//finds start and endpoints, area, and location
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
