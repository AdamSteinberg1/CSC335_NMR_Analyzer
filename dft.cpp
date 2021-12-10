#include <armadillo>
#include <vector>
#include <complex>

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

//applies the Discrete Fourier Transform Filter to an arma::cx_vec
arma::cx_vec dftFilter(const arma::cx_vec& y)
{
  int n = y.size(); //n is the dimension of y
  arma::cx_mat Z = makeZ(n);
  arma::cx_vec c = Z*y;
  arma::cx_mat G = makeG(n);
  c = G*c;
  return conj(Z)*c;
}


//applies the Discrete Fourier Transform Filter to a std::vector
std::vector<std::pair<double, double>> dftFilter(std::vector<std::pair<double, double>> data)
{
  int n = data.size();
  //put all the y-values of the elements in data into an armadillo vector
  arma::cx_vec y(n);
  for (int i = 0; i < n; i++)
  {
    y[i] = data[i].second;
  }

  y = dftFilter(y);

  //put all the filtered y values back into the std::vector data
  for (int i = 0; i < n; i++)
  {
    data[i].second = y[i].real(); //discard imaginary part
  }

  return data;
}
