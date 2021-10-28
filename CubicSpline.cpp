#include "CubicSpline.h"

//constructs a natural cubic spline from points
//each cubic we make is defined by four constants: a_i, b_i, c_i, d_i
CubicSpline::CubicSpline(std::vector<std::pair<double, double>> points)
{
  std::sort(points.begin(), points.end()); //points must be in order for this algorithm

  //how many cubics we're going to make
  int n = points.size()-1;

  xValues.reserve(n+1);
  cubics.reserve(n);

  for(int i = 0; i <= n; i++)
    xValues.push_back(points[i].first);

  //h contains the difference between consecutive x-values
  arma::vec h(n, arma::fill ::zeros);
  for(int i = 0; i < n; i++)
    h(i) = points[i+1].first - points[i].first;

  //alpha is a vector representing the constants on the right side of our system of linear equations
  arma::vec alpha(n+1, arma::fill ::zeros);
  for(int i = 1; i < n; i++)
    alpha(i) = 3/h[i]*(points[i+1].second - points[i].second) - 3/h[i-1]*(points[i].second - points[i-1].second);

  //A is a matrix representing the coefficients in our system of linear equations
  arma::mat A(n+1, n+1, arma::fill ::zeros);
  A(0,0) = 1.0;
  A(n,n) = 1.0;
  for(int i = 1; i < n; i++)
  {
    A(i,i-1) = h[i-1];
    A(i,i) = 2*(h[i-1]+h[i]);
    A(i,i+1) = h[i];
  }

  //solve A*c = alpha for c
  //c contains all our c_i constants
  arma::vec c = solve(A, alpha);

  for(int i = 0; i < n; i++)
  {
    //calculate all the constants that define the ith cubic
    double a_i = points[i].second;
    double b_i = (points[i+1].second - points[i].second)/h[i] - h[i]*(c(i+1)+2*c(i))/3;
    double c_i = c(i);
    double d_i = (c(i+1)-c(i))/(3*h[i]);
    //the x-values of the input points
    double x_i = points[i].first;

    //construct the ith cubic p
    Polynomial diff = {-x_i, 1}; // (x - x_i)
    Polynomial p = a_i + b_i*diff + c_i*diff.power(2) + d_i*diff.power(3);
    cubics.push_back(p);
  }
}

//gets how many cubic have been stitched together
int CubicSpline::getNumCubics() const
{
  return cubics.size();
}

//get the ith cubic polynomial
Polynomial CubicSpline::operator[](int i) const
{
  return cubics[i];
}

//get the range of x values that the ith cubic is valid over
std::pair<double, double> CubicSpline:: getRange(int i) const
{
  double inf = std::numeric_limits<double>::infinity();
  if(i == 0)
  {
    return {-inf, xValues[1]};
  }
  int lastIndex = cubics.size();
  if(i==lastIndex)
  {
    return {xValues[lastIndex-1], inf};
  }

  return {xValues[i], xValues[i+1]};
}

//returns the index of the cubic that x would be evaluated with
int CubicSpline::findIndex(double x) const
{
  //check edge cases where x is not between any two x-values
  if(x < xValues.front())
    return 0;
  if(x > xValues.back())
    return cubics.size()-1;

  //otherwise, peform binary search
  return findIndex(x, 0, xValues.size()-1);
}

//performs binary search for x's corresponding index
//searches between left and right indices
int CubicSpline::findIndex(double x, int left, int right) const
{
  if (right >= left)
  {
    int mid = (left + right)/ 2;

    if (xValues[mid] <= x && x <= xValues[mid+1])
        return mid;

    if (xValues[mid] > x)
        return findIndex(x, left, mid - 1);

    return findIndex(x, mid + 1, right);
  }

  //should never get here
  return -1;
}

//evaluate the cubic spline at x
double CubicSpline::evaluate(double x) const
{
  Polynomial p = cubics[findIndex(x)];
  return p.evaluate(x);
}
