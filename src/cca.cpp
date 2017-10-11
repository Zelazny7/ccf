#include <RcppEigen.h>
#include "util.h"
#include <stdlib.h>
//#include <vector.h>

// [[Rcpp::depends(RcppEigen)]]

// let Eigen center the matrices
#define center(mat) mat.rowwise() - mat.colwise().mean()

// return the decom information in this struct
typedef struct  {
  //Eigen::MatrixXd Q;
  Eigen::MatrixXd R;
  Eigen::MatrixXf::Index* perm;
  int rk;
} qr_result;

// mostly used for debugging and verifying results in R
Rcpp::List qr_result_to_list(qr_result result) {

  // copy perm to numeric vector
  Rcpp::NumericVector perm(result.rk);
  for (int i = 0; i < perm.size(); i++) {
    perm[i] = result.perm[i] + 1;
  }

  // free the allocated mem for perm
  std::free(result.perm);

  return Rcpp::List::create(
    //Rcpp::Named("Q") = Rcpp::wrap(result.Q),
    Rcpp::Named("R") = Rcpp::wrap(result.R),
    Rcpp::Named("perm") = perm,
    Rcpp::Named("rank") = Rcpp::wrap(result.rk));

}

qr_result qr_decomp_cpp(Eigen::MatrixXd mat) {

  const double threshold = 0.0001;

  Eigen::ColPivHouseholderQR< Eigen::MatrixXd > qr;
  qr.setThreshold(threshold);
  qr.compute(mat);

  int rk = qr.rank();

  Eigen::MatrixXd perm = qr.colsPermutation();

  // find col-index positions of permutation matrix
  Eigen::MatrixXf::Index* indices = (Eigen::MatrixXf::Index*) std::calloc(perm.cols(), sizeof(indices));

  for(int i=0; i < perm.cols(); ++i) {
    perm.col(i).maxCoeff( &indices[i] );
  }

  //Eigen::ColPivHouseholderQR< Eigen::MatrixXd >::HouseholderSequenceType  seq = qr.householderQ();


  qr_result result = {
    //((Eigen::MatrixXd) qr.householderQ()).block(0, 0, mat.cols(), rk),
    qr.matrixQR().topLeftCorner(rk, rk).triangularView<Eigen::Upper>(),
    indices,
    rk
  };

  return result;

}



// [[Rcpp::export]]
Rcpp::List test_eigen(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y) {

  //MatrixXd centered = mat.rowwise() - mat.colwise().mean();

  qr_result qrX = qr_decomp_cpp(center(X));
  qr_result qrY = qr_decomp_cpp(center(Y));

  //Rcpp::List qr_result_to_list(qr_result result)

  return qr_result_to_list(qrX);


  // //  QR decomposition
  // // (https://cran.r-project.org/doc/contrib/Hiebeler-matlabR.pdf)
  // //  qrDecompX <- qr(x, tol = epsilon)
  // //  qX <- qr.Q(qrDecompX)
  // //  rX <- qr.R(qrDecompX)
  // //  pX <- qrDecompX$pivot
  // //  rankX <- qrDecompX$rank
  // //
  // //  qrDecomp <- qr(y, tol = epsilon)
  // //  qY <- qr.Q(qrDecomp)
  // //  rY <- qr.R(qrDecomp)
  // //  pY <- qrDecomp$pivot
  // //  rankY <- qrDecomp$rank
  // //
  //
  //
  // return Rcpp::List::create(
  //   Rcpp::Named("QR") = Rcpp::wrap(QR),
  //   Rcpp::Named("Q") = Rcpp::wrap(Q),
  //   Rcpp::Named("R") = Rcpp::wrap(R),
  //   Rcpp::Named("perm") = Rcpp::wrap(perm));
  //   //Rcpp::Named("rk") = Rcpp::as<NumericVector>(rk));
  //
  //return R_NilValue;

}
