//class to represent a polynomial
#pragma once
#include <vector>
#include <cmath>
#include <ostream>

class Polynomial
{
  private:
    //the coefficients of each term are stored in order of ascending degree
    //for example ax^3+bx^2+cx+d is stored as [d,c,b,a]
    std::vector<double> coefficients;

    //returns the result of distributing ax^n to this polynomial
    Polynomial distribute(double a, int n) const;
  public:
    //default constructor creates the polynomial 0x^0
    Polynomial();

    //constructs a polynomial where args are the coefficients listed from lowest to highest degree
    Polynomial(std::initializer_list<double> args);

    //constructs a polynomial where args are the coefficients listed from lowest to highest degree
    Polynomial(std::vector<double> args);

    std::vector<double> getCoefficients() const;
    double operator[](int i) const;
    int getDegree() const;

    //returns the polynomial evaluated at a specific x value
    double evaluate(double x) const;
    //returns the derivative of the polynomial
    Polynomial derivative() const;

    //finds a root of the polynomial using Newton's Method with an initial approximation of p0
    double root(double p0) const;
    Polynomial power(int n) const;

    friend Polynomial operator*(Polynomial a, Polynomial b);
};
Polynomial operator+(Polynomial a, Polynomial b);
Polynomial operator*(double scalar, Polynomial p);
Polynomial operator+(double scalar, Polynomial p);
Polynomial operator/(Polynomial p, double scalar);
Polynomial operator-(Polynomial a, Polynomial b);
std::ostream& operator<<(std::ostream& os, const Polynomial& p);
