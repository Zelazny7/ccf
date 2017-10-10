#include <RcppEigen.h>
#include "util.h"
#include <stdlib.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List test_eigen(Rcpp::NumericMatrix Xin, Rcpp::NumericMatrix Yin) {

  //Rcpp::NumericMatrix X_cen = center_cpp(Xin, Rcpp::LogicalVector(1));
  //Rcpp::NumericMatrix Y_cen = center_cpp(Yin, Rcpp::LogicalVector(1));

  // convert centered matrices to Eigen types
  Eigen::Map<Eigen::MatrixXd> X = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Xin);
  Eigen::Map<Eigen::MatrixXd> Y = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Yin);

  // QR Decomposition with column pivot for both centered matrices
  const double threshold = 0.0001;

  Eigen::ColPivHouseholderQR< Eigen::MatrixXd > qrX;
  //qrX.setThreshold(threshold);
  qrX.compute(X);

  // extract parts we need

  Eigen::ColPivHouseholderQR< Eigen::MatrixXd > qrY;
  qrY.setThreshold(threshold);
  qrY.compute(X);

  // extract and return parts and compare to R
  Eigen::Index rk = qrX.rank();


  Eigen::MatrixXd QR = qrX.matrixQR();
  Eigen::MatrixXd Q = ((Eigen::MatrixXd) qrX.householderQ()).block(0, 0, Xin.nrow(), rk);
  Eigen::MatrixXd R = qrX.matrixQR().topLeftCorner(rk, rk).triangularView<Eigen::Upper>();



  //qr.matrixQR().triangularView<Upper>()

  Eigen::MatrixXd perm = qrX.colsPermutation();

  // find the rank

  Rcpp::Rcout << "Rank of matrix" << rk << std::endl;

  //x = X.colPivHouseholderQr();


  // center the matrix
  //
  // Rcpp::NumericMatrix Y_cen = center_cpp(Yin, Rcpp::LogicalVector(1));


  //  QR decomposition
  // (https://cran.r-project.org/doc/contrib/Hiebeler-matlabR.pdf)
  //  qrDecompX <- qr(x, tol = epsilon)
  //  qX <- qr.Q(qrDecompX)
  //  rX <- qr.R(qrDecompX)
  //  pX <- qrDecompX$pivot
  //  rankX <- qrDecompX$rank
  //
  //  qrDecomp <- qr(y, tol = epsilon)
  //  qY <- qr.Q(qrDecomp)
  //  rY <- qr.R(qrDecomp)
  //  pY <- qrDecomp$pivot
  //  rankY <- qrDecomp$rank
  //


  return Rcpp::List::create(
    Rcpp::Named("QR") = Rcpp::wrap(QR),
    Rcpp::Named("Q") = Rcpp::wrap(Q),
    Rcpp::Named("R") = Rcpp::wrap(R),
    Rcpp::Named("perm") = Rcpp::wrap(perm));
    //Rcpp::Named("rk") = Rcpp::as<NumericVector>(rk));

  //return R_NilValue;

}
