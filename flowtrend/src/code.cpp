#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// Generated from _main.Rmd: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <numeric>

using namespace arma;
using namespace Eigen;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision

//' Solve "barebones" sylvester equation that takes upper triangular matrices as coefficients.
//'
//' @param TA Upper-triangular matrix
//' @param TB Upper-triangular matrix
//' @param C matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd matrix_function_solve_triangular_sylvester_barebonesC2(const Eigen::MatrixXd & TA, 
								     const Eigen::MatrixXd & TB,
								     const Eigen::MatrixXd & C){
  // Eigen::eigen_assert(TA.rows() == TA.cols());
  // Eigen::eigen_assert(TA.Eigen::isUpperTriangular());
  // Eigen::eigen_assert(TB.rows() == TB.cols());
  // Eigen::eigen_assert(TB.Eigen::isUpperTriangular());
  // Eigen::eigen_assert(C.rows() == TA.rows());
  // Eigen::eigen_assert(C.cols() == TB.rows());

  // typedef typename MatrixType::Index Index; 
  // typedef typename MatrixType::Scalar Scalar;

  int m = TA.rows();
  int n = TB.rows();
  Eigen::MatrixXd X(m, n);

  for (int i = m - 1; i >= 0; --i) {
    for (int j = 0; j < n; ++j) {

      // Compute T_A X = \sum_{k=i+1}^m T_A_{ik} X_{kj}
      double TAX;
      if (i == m - 1) {
      	TAX = 0;
      } else {
	MatrixXd TAXmatrix = TA.row(i).tail(m-1-i) * X.col(j).tail(m-1-i);
      	TAX = TAXmatrix(0,0);
      }

      // Compute X T_B = \sum_{k=1}^{j-1} X_{ik} T_B_{kj}
      double XTB;
      if (j == 0) {
      	XTB = 0;
      } else {
      	MatrixXd XTBmatrix = X.row(i).head(j) * TB.col(j).head(j);
      	XTB = XTBmatrix(0,0);
      }

      X(i,j) = (C(i,j) - TAX - XTB) / (TA(i,i) + TB(j,j));
    }
  }
  return X;
}

// [[Rcpp::depends(RcppArmadillo)]]

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


// [[Rcpp::export]]
arma::vec dmvnorm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
    using arma::uword;
    uword const n = x.n_rows,
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())),
                constants = -(double)xdim/2.0 * log2pi,
              other_terms = rootisum + constants;

    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }

    if (logd)
      return out;
    return exp(out);
}

// All credit goes to https://gallery.rcpp.org/articles/dmvnorm_arma/

// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <cmath> // For standard C++ functions if needed, though not strictly required here


using namespace std;

// Dynamic programming algorithm for the 1d fused lasso problem
// (Ryan Tibshirani's implementation of Nick Johnson's algorithm)

void prox_dp_core_cpp(
    int n,
    const double *y_ptr, // Use const since input y is not modified
    double lam,
    double *beta_ptr)
{
    // Take care of a few trivial cases
    if (n == 0) return;
    if (n == 1 || lam == 0) {
        for (int i = 0; i < n; i++) beta_ptr[i] = y_ptr[i];
        return;
    }

    // Use std::vector for automatic memory management (RAII) instead of malloc/free
    // These store the derivative knots (x) and coefficients (a, b)
    std::vector<double> x(2 * n);
    std::vector<double> a(2 * n);
    std::vector<double> b(2 * n);
    int l, r;

    // These are the knots of the back-pointers
    std::vector<double> tm(n - 1);
    std::vector<double> tp(n - 1);

    // Aliases to make the original logic clearer (using ptrs to vector data)
    double *x_data = x.data();
    double *a_data = a.data();
    double *b_data = b.data();
    double *tm_data = tm.data();
    double *tp_data = tp.data();

    double afirst, alast, bfirst, blast;

    // We step through the first iteration manually
    tm_data[0] = -lam + y_ptr[0];
    tp_data[0] = lam + y_ptr[0];
    l = n - 1;
    r = n;
    x_data[l] = tm_data[0];
    x_data[r] = tp_data[0];
    a_data[l] = 1;
    b_data[l] = -y_ptr[0] + lam;
    a_data[r] = -1;
    b_data[r] = y_ptr[0] + lam;
    afirst = 1;
    bfirst = -lam - y_ptr[1];
    alast = -1;
    blast = -lam + y_ptr[1];

    // Now iterations 2 through n-1
    int lo, hi;
    double alo, blo, ahi, bhi;
    for (int k = 1; k < n - 1; k++) {
        // Compute lo: step up from l until the
        // derivative is greater than -lam
        alo = afirst;
        blo = bfirst;
        for (lo = l; lo <= r; lo++) {
            if (alo * x_data[lo] + blo > -lam) break;
            alo += a_data[lo];
            blo += b_data[lo];
        }

        // Compute the negative knot
        tm_data[k] = (-lam - blo) / alo;
        l = lo - 1;
        x_data[l] = tm_data[k];

        // Compute hi: step down from r until the
        // derivative is less than lam
        ahi = alast;
        bhi = blast;
        for (hi = r; hi >= l; hi--) {
            if (-ahi * x_data[hi] - bhi < lam) break;
            ahi += a_data[hi];
            bhi += b_data[hi];
        }

        // Compute the positive knot
        tp_data[k] = (lam + bhi) / (-ahi);
        r = hi + 1;
        x_data[r] = tp_data[k];

        // Update a and b
        a_data[l] = alo;
        b_data[l] = blo + lam;
        a_data[r] = ahi;
        b_data[r] = bhi + lam;
        afirst = 1;
        bfirst = -lam - y_ptr[k + 1];
        alast = -1;
        blast = -lam + y_ptr[k + 1];
    }

    // Compute the last coefficient: this is where
    // the function has zero derivative

    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
        if (alo * x_data[lo] + blo > 0) break;
        alo += a_data[lo];
        blo += b_data[lo];
    }
    beta_ptr[n - 1] = -blo / alo;

    // Compute the rest of the coefficients, by the
    // back-pointers
    for (int k = n - 2; k >= 0; k--) {
        if (beta_ptr[k + 1] > tp_data[k]) beta_ptr[k] = tp_data[k];
        else if (beta_ptr[k + 1] < tm_data[k]) beta_ptr[k] = tm_data[k];
        else beta_ptr[k] = beta_ptr[k + 1];
    }

    // Memory cleanup is now automatic because std::vector is used.
}


//' Implements the dynamic programming algorithm for the 1D Fused Lasso proximal operator.
//'
//' @param z A numeric vector of input values ($\mathbf{y}$).
//' @param lam A numeric scalar, the regularization parameter ($\lambda$).
//' @return A numeric vector, the result of the proximal operation ($\boldsymbol{\beta}$).
//' @export
// [[Rcpp::export]]
arma::vec prox_dp(arma::vec z, double lam) {

    // Set up the output vector
    int n = z.n_elem;
    arma::vec beta = arma::zeros<arma::vec>(n);

    // Call the core C++ logic
    // z.memptr() provides a pointer to the internal data (const for input)
    // beta.memptr() provides a pointer to the internal data (non-const for output)
    prox_dp_core_cpp(n, z.memptr(), lam, beta.memptr());

    return beta;
}
