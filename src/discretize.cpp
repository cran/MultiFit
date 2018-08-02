#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List discretizeCpp( arma::mat a,
                          arma::mat b,
                          arma::rowvec w,
                          arma::uvec mask,
                          arma::umat ij,
                          int d1,
                          int d2
                       ) {
  // # Examine only data points from parent window:
  arma::uvec msk = find(mask>0);
  a = a.rows(msk);
  b = b.rows(msk);

  int m = a.n_rows;
  int o = accu(mask);

  rowvec l1 = w.subvec( d1+d2, d1+d2+d1-1 );
  rowvec k1 = w.subvec( 0, d1-1 );
  rowvec p2k1(d1);
  for (int i=0; i<d1; i++) p2k1(i) = pow(2, -k1(i));
  rowvec xl = (l1-1.0) % p2k1;
  rowvec xh = l1 % p2k1;

  rowvec l2 = w.subvec( d1+d2+d1, d1+d2+d1+d2-1 );
  rowvec k2 = w.subvec( d1, d1+d2-1 );
  rowvec p2k2(d2);
  for (int j=0; j<d2; j++) p2k2(j) = pow(2, -k2(j));
  rowvec yl = (l2-1.0) % p2k2;
  rowvec yh = l2 % p2k2;

  // # For each pair of margins, discretize to a 2x2 contingecy table:
  rowvec xm = (xl+xh)/2.0;
  mat x_mid = repmat(xm, m, 1);
  rowvec ym = (yl+yh)/2.0;
  mat y_mid = repmat(ym, m, 1);
  umat x0 = a < x_mid ;
  umat y0 = b < y_mid ;
  umat tables(ij.n_rows, 8);
  tables.cols(0, 1) = ij;
  for (unsigned int c=0; c<ij.n_rows; c++) {
    int i = ij(c, 0) - 1;
    int j = ij(c, 1) - 1;
    tables(c,2) = accu(x0.col(i) % y0.col(j));
    tables(c,3) = accu(x0.col(i) % (1-y0.col(j)));
    tables(c,4) = accu((1-x0.col(i)) % y0.col(j));
    tables(c,5) = o - tables(c,2) - tables(c,3) - tables(c,4);
  }
  tables.col(6).fill(w(w.n_cols - 1));
  tables.col(7).fill(0);

  return Rcpp::List::create(
    Rcpp::Named( "tables" ) = tables,
    Rcpp::Named( "mask" ) = mask,
    Rcpp::Named( "x0.sub.mask" ) = x0,
    Rcpp::Named( "y0.sub.mask" ) = y0 );
}
