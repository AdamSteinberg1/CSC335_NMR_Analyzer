#include <armadillo>
#include <vector>
#include <complex>
#include <iostream>

arma::cx_mat makeZ(int n)
{
  arma::cx_mat Z(n, n); //Z is an nxn matrix
  for(int j = 0; j<n; j++)
  {
    for(int k = 0; k<n; k++)
    {
      Z(j,k) = std::pow(std::polar(1.0, -2*M_PI/n), j*k) / sqrt(n);
    }
  }
  return Z;
}

arma::cx_mat makeG(int n)
{
  //G is an nxn matrix initialized to all zeros
  arma::cx_mat G(n,n, arma::fill::zeros);
  //the Kronecker delta function in formula for the elements of G makes it so
  //that only the diagonal of G is nonzero
  for(int i = 0; i<n; i++)
  {
    G(i,i) = exp((-4*M_LN2*i*i) / pow(n, 1.5));
  }
  return G;
}

arma::cx_mat makeZ_inverse(const arma::cx_mat& Z)
{
  int n = Z.n_rows;
  arma::cx_mat result(n,n);
  for(int j = 0; j<n; j++)
  {
    for(int k = 0; k<n; k++)
    {
      result(j,k) = conj(Z(j,k));
    }
  }
  return result;
}

//solves the equation c = Z*y and returns y
arma::cx_vec recoverDFT(const arma::cx_mat& Z, const arma::cx_vec& c, int method)
{
  switch(method)
  {
    case 0: //inverse
    {
      arma::cx_mat Z_inv = makeZ_inverse(Z);
      return Z_inv*c;
    }
    case 1: //direct
      return arma::solve(Z, c);
    case 2: //iterative
      return arma::cx_vec(); //TODO implement iterative recovery
    default:
      std::cerr << "Error: DFT recovery method " << method << " is not valid." << std::endl;
      exit(1);
  }
}

//applies the Discrete Fourier Transform Filter to an arma::cx_vec
arma::cx_vec dftFilter(const arma::cx_vec& y, int method)
{
  int n = y.size(); //n is the dimension of y
  arma::cx_mat Z = makeZ(n);
  arma::cx_vec c = Z*y;
  arma::cx_mat G = makeG(n);
  c = G*c;
  return recoverDFT(Z, c, method);
}


//applies the Discrete Fourier Transform Filter to a std::vector
std::vector<std::pair<double, double>> dftFilter(std::vector<std::pair<double, double>> data, int method)
{
  int n = data.size();
  //put all the y-values of the elements in data into an armadillo vector
  arma::cx_vec y(n);
  std::cout << "unfiltered y-values:" << std::endl;
  for (int i = 0; i < n; i++)
  {
    y[i] = data[i].second;
    std::cout << y[i] << std::endl; //debug output
  }

  y = dftFilter(y, method);


  std::cout << "filtered y-values:" << std::endl;
  //put all the filtered y values back into the std::vector data
  for (int i = 0; i < n; i++)
  {
    std::cout << y[i] << std::endl; //debug output
    data[i].second = y[i].real(); //unsure if this should be done this way
  }

  return data;
}
