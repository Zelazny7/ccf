#include <Rcpp.h>
#include <stdlib.h>

/***
 * Fast centering function with inplace option
 *
 * @param m Numeric matrix to center
 * @param inplace Logical vector indicating whether to center the matrix inplace *
 *
 */

// [[Rcpp::export]]
SEXP center_cpp(Rcpp::NumericMatrix m, Rcpp::LogicalVector inplace) {

  size_t nrows = m.nrow(); size_t ncols = m.ncol();

  // calculate the means
  double* means = (double*) std::calloc(ncols, sizeof(double));

  for (size_t i = 0; i < ncols; i++) {
    for (size_t j = 0; j < nrows; j++) {

      means[i] += m[i * nrows + j];

    }
    means[i] /= nrows;
  }


  if (inplace[0] == TRUE) {

    for (size_t i = 0; i < ncols; i++) {
      for (size_t j = 0; j < nrows; j++) {

        m[i * nrows + j] -= means[i];

      }
    }

    free(means);
    return R_NilValue;

  } else {

    Rcpp::NumericMatrix out(Rcpp::clone(m));

    for (size_t i = 0; i < ncols; i++) {
      for (size_t j = 0; j < nrows; j++) {
        out(j, i) -= means[i];
      }
    }

    free(means);
    return out;
  }

}
