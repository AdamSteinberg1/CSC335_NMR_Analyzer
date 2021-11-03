//class for a cubic spline
#include "Polynomial.h"
#include <limits>
#include <algorithm>
#include <armadillo>
#pragma once

class CubicSpline
{
  private:
    //the cubics that make up the spline
    std::vector<Polynomial> cubics;
    //the x-values at which the cubics are stitched together
    std::vector<double> xValues;

    //returns the index of the cubic that x would be evaluated with
    int findIndex(double x) const;
    //helper function to perform binary search
    int findIndex(double x, int left, int right) const;
  public:
    //constructs a natural cubic spline from points
    CubicSpline(std::vector<std::pair<double, double>> points);
    //gets how many cubics have been stitched together
    int getNumCubics() const;
    //get the ith cubic polynomial
    Polynomial operator[](int i) const;
    //get the range of x values that the ith cubic is valid over
    std::pair<double, double> getRange(int i) const;
    //evaluate the cubic spline at x
    double evaluate(double x) const;

    friend std::string gnuPrint(CubicSpline c);
};
