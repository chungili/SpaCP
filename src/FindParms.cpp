#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
List ScoreCpp(const Eigen::SparseMatrix<double>& nrisk,
              const Eigen::MatrixXd& X,
              const Eigen::VectorXd& beta,
              const Eigen::VectorXd& ri,
              const Eigen::VectorXi& Mi,
              const Eigen::VectorXi& delta) {
  
  const double eps = 1e-12;
  
  int p = beta.size();
  int q = ri.size();
  const int n = X.rows();
  
  Eigen::VectorXd theta(p + q);
  theta.head(p) = beta;
  theta.tail(q) = ri;
  
  Eigen::VectorXd etai = X * theta;
  Eigen::VectorXd b = etai.array().exp();
  
  const int M = Mi.maxCoeff();
  if (M <= 0) stop("Mi must be positive (1..M).");
  
  Eigen::VectorXd sr0 = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd sr1 = Eigen::MatrixXd::Zero(M, n);
  
  // nrisk is column-compressed (dgCMatrix). Iterate by column j:
  for (int j = 0; j < n; ++j) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(nrisk, j); it; ++it) {
      const int i = it.row();           // row index of nonzero
      // if nrisk is lower triangular, typically i >= j; we don't need to enforce,
      // but you can uncomment next line as a safeguard:
      // if (i < j) continue;
      
      // value is 1 (or >0)
      const double bij = b[i];
      
      sr0[j] += bij;
      
      const int g = Mi[i] - 1;          // 1-based -> 0-based
      if (g < 0 || g >= M) stop("Mi out of range (expected 1..M).");
      sr1(g, j) += bij;
    }
  }
  
  // fbetar = sum( (b - log(sr0)) * delta )
  double fbetar =
    ((etai.array() - (sr0.array() + eps).log()) * delta.cast<double>().array()).sum();
  
  // ---- sr10, D ----
  // sr10 = sweep(sr1, 2, colSums(sr1), "/")
  Eigen::RowVectorXd colsum = sr1.colwise().sum();  // 1 x n
  colsum.array() += eps;
  
  Eigen::MatrixXd sr10 = sr1.array().rowwise() / colsum.array(); // M x n
  
  // sr10_w = sweep(sr10, 2, delta, "*")
  Eigen::RowVectorXd d = delta.cast<double>().transpose();       // 1 x n
  Eigen::MatrixXd sr10_w = sr10.array().rowwise() * d.array();   // M x n
  
  // D1 = sr10_w %*% t(sr10)   => (M x n)(n x M) = M x M
  Eigen::MatrixXd D1 = sr10_w * sr10.transpose();
  
  // D2 = diag(as.vector(sr10 %*% delta), nrow=M)
  Eigen::VectorXd v = sr10 * delta.cast<double>();               // M x 1
  Eigen::MatrixXd D2 = v.asDiagonal();                           // M x M
  
  Eigen::MatrixXd D = D1 - D2;
  // Eigen::MatrixXd sr = nrisk.cast<double>().array().colwise() * b.array();
  // 
  // Eigen::VectorXd sr0 = sr.colwise().sum().transpose();
  // 
  // 
  // if (M <= 0) stop("Mi must be positive (1..M).");
  // Eigen::MatrixXd sr1 = Eigen::MatrixXd::Zero(M, sr.cols());
  // 
  // for (int j = 0; j < sr.cols(); ++j) {
  //   for (int i = 0; i < sr.rows(); ++i) {
  //     const int g = Mi[i] - 1;          // 1-based -> 0-based
  //     if (g < 0 || g >= M) stop("Mi out of range (expected 1..M).");
  //     sr1(g, j) += sr(i, j);
  //   }
  // }
  // 
  // const double eps = 1e-12; // avoid log(0)
  // double fbetar =
  //   ( (b.array() - (sr0.array() + eps).log()) * delta.cast<double>().array() ).sum();
  // 
  // Eigen::RowVectorXd colsum = sr1.colwise().sum();
  // colsum.array() += eps;
  // Eigen::MatrixXd sr10 = sr1.array().rowwise() / colsum.array();
  // Eigen::RowVectorXd d = delta.cast<double>().transpose(); // 1 x n
  // Eigen::MatrixXd sr10_w = sr10.array().rowwise() * d.array();
  // Eigen::MatrixXd D1 = sr10_w * sr10.transpose();
  // Eigen::VectorXd v = sr10 * delta.cast<double>();
  // Eigen::MatrixXd D2 = v.asDiagonal();   // M x M
  // Eigen::MatrixXd D = D1 - D2;
  
  Eigen::VectorXd pfpr = Eigen::VectorXd::Zero(M);
  
  for (int j = 0; j < delta.size(); ++j) {
    const int g = Mi[j] - 1;     // 1-based -> 0-based
    if (g < 0 || g >= M) stop("Mi out of range (expected 1..M).");
    pfpr[g] += delta[j];
  }
  
  // 2) subtract sr10 %*% delta
  pfpr.noalias() -= sr10 * delta.cast<double>();
  
  return List::create(
    Named("b")         = b,
//    Named("sr")       = sr,
    Named("sr0")       = sr0,
    Named("sr1")       = sr1,
    Named("M")   = M,
    Named("delta")   = delta,
    Named("fbetar")   = fbetar,
    Named("pfpr")   = pfpr,
    Named("D")   = D
  );
}

// [[Rcpp::export]]
Eigen::VectorXd FindrCpp(const Eigen::SparseMatrix<double>& nrisk,
                      const Eigen::MatrixXd& X,
                      const Eigen::VectorXd& beta,
                      const Eigen::VectorXd& ri,
                      const Eigen::VectorXi& Mi,      // length N, 1..M
                      const Eigen::VectorXi& delta,   // length N
                      const double eta,
                      const double sigma2,
                        const Eigen::MatrixXd& d2,      // M x M
                        const double lambda = 1e-5) {
  
  const int M = Mi.maxCoeff();
  // ---- Call ScoreCpp to get pfpr ----
  List s = ScoreCpp(nrisk, X, beta, ri, Mi, delta);
  Eigen::VectorXd pfpr = as<Eigen::VectorXd>(s["pfpr"]);
  // ---- Q = sigma2 * exp(-d2/eta^2) + lambda*I ----
  const double eta2 = eta * eta;
  Eigen::MatrixXd Q = (-d2.array() / eta2).exp().matrix();
  Q *= sigma2;
  Q.diagonal().array() += lambda;
  
  // solve(Q, r)
  Eigen::VectorXd Qinv_r(M);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(Q);
  if (ldlt.info() == Eigen::Success) {
    Qinv_r = ldlt.solve(ri);
  } else {
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
    if (!lu.isInvertible()) stop("Q is not invertible (even after adding lambda).");
    Qinv_r = lu.solve(ri);
  }
  
  return pfpr - Qinv_r;
  //return List::create(
  //  Named("pfpr_Qinvr")         = pfpr - Qinv_r
  //);
}

// [[Rcpp::export]]
Eigen::MatrixXd MyjacCpp(const Eigen::SparseMatrix<double>& nrisk,
                         const Eigen::MatrixXd& X,
                         const Eigen::VectorXd& beta,
                         const Eigen::VectorXd& ri,
                         const Eigen::VectorXi& Mi,
                         const Eigen::VectorXi& delta,
                         const double eta,
                         const double sigma2,
                         const Eigen::MatrixXd& d2,
                         const double lambda = 1e-5) {
  
  if (eta <= 0.0) Rcpp::stop("eta must be > 0.");
  if (sigma2 < 0.0) Rcpp::stop("sigma2 must be >= 0.");
  
  const int n = X.rows();
  if (Mi.size() != n) Rcpp::stop("Mi length must equal X.rows().");
  if (delta.size() != n) Rcpp::stop("delta length must equal X.rows().");
  
  const int M = Mi.maxCoeff();
  if (M <= 0) Rcpp::stop("Mi must be positive (1..M).");
  if (d2.rows() != M || d2.cols() != M) Rcpp::stop("d2 must be M x M.");
  
  // 1) Call ScoreCpp and directly extract D
  Rcpp::List s = ScoreCpp(nrisk, X, beta, ri, Mi, delta);
  Eigen::MatrixXd D = Rcpp::as<Eigen::MatrixXd>(s["D"]);   // M x M
  
  // 2) Q = sigma2 * exp(-d2/eta^2) + lambda I
  const double eta2 = eta * eta;
  Eigen::MatrixXd Q = (-d2.array() / eta2).exp().matrix();
  Q *= sigma2;
  Q.diagonal().array() += lambda;
  
  // 3) Q^{-1} via LDLT
  Eigen::LDLT<Eigen::MatrixXd> ldlt(Q);
  if (ldlt.info() != Eigen::Success)
    Rcpp::stop("LDLT decomposition failed for Q.");
  
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(M, M);
  Eigen::MatrixXd Qinv = ldlt.solve(I);
  
  // 4) return D - solve(Q)
  return D - Qinv;
}

// [[Rcpp::export]]
double l21Cpp(const Eigen::VectorXd& parms,
              const Eigen::VectorXd& rhat,
              const Eigen::SparseMatrix<double>& nriskSp,
              const Eigen::MatrixXd& X,                // X.t from R
              const Eigen::VectorXi& Mi,             // Mi.t from R
              const Eigen::VectorXi& delta,          // delta.t from R
              const Eigen::MatrixXd& d2,               // dt$d2
              const double lambda = 1e-7,
              const double fail_value = 1e6) {
  
  const int npar = parms.size();
  if (npar < 3) return fail_value;
  
  // beta = parms[1:(npar-2)]
  const int p = npar - 2;
  Eigen::VectorXd beta = parms.head(p);
  
  // sigma2 = exp(tail(parms,2)[1]), eta = tail(parms,2)[2]
  const double sigma2 = std::exp(parms[npar - 2]);
  const double eta    = parms[npar - 1];
  if (!(sigma2 >= 0.0) || !(eta > 0.0)) return fail_value;
  
  const int M = rhat.size();
  if (M <= 0) return fail_value;
  
  // ---- dimension checks ----
  const int N = X.rows();
  // if (X.cols() != p + M) {
  //   // 你的模型是 X = [Xij | as.factor(Mi)]，theta = [beta | rhat]
  //   // 所以 X 需要有 p + M 欄
  //   return fail_value;
  // }
  // if (Mi.size() != N || delta.size() != N || Yij.size() != N) return fail_value;
  // 
  // if (nriskSp.rows() != N || nriskSp.cols() != N) return fail_value;
  // 
  // // Mi_t must be in 1..M
  // if (Mi.maxCoeff() != M) return fail_value;
  // if (Mi.minCoeff() < 1)  return fail_value;
  
  // d2 must be M x M
  // if (d2.rows() != M || d2.cols() != M) return fail_value;
  
  // ---- ScoreCpp ----
  Rcpp::List s = ScoreCpp(nriskSp, X, beta, rhat, Mi, delta);
  const double fbetar = Rcpp::as<double>(s["fbetar"]);
  Eigen::MatrixXd D = Rcpp::as<Eigen::MatrixXd>(s["D"]);   // M x M
  
  // if (D.rows() != M || D.cols() != M) return fail_value;
  
  // ---- Q = sigma2 * exp(-d2/eta^2) + lambda I ----
  const double eta2 = eta * eta;
  Eigen::MatrixXd Q = (-d2.array() / eta2).exp().matrix();
  Q *= sigma2;
  Q.diagonal().array() += lambda;
  
  // // ---- quad_term = r' * solve(Q, r) ----
  Eigen::LDLT<Eigen::MatrixXd> ldltQ(Q);
  if (ldltQ.info() != Eigen::Success) return fail_value;
  // 
  Eigen::VectorXd Qinv_r = ldltQ.solve(rhat);
  if (ldltQ.info() != Eigen::Success) return fail_value;
  // 
  const double quad_term = rhat.dot(Qinv_r);
  // 
  // // ---- det_term = det(I - Q %*% D) ----
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(M, M) - Q * D;
  
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
  if (!lu.isInvertible()) return fail_value;
  
  const Eigen::MatrixXd& LU = lu.matrixLU();
  
  double log_abs_det = 0.0;
  int sign_det = 1;
  
  for (int i = 0; i < M; ++i) {
    double di = LU(i, i);
    if (di == 0.0) return fail_value;
    if (di < 0.0) sign_det = -sign_det;
    log_abs_det += std::log(std::abs(di));
  }
  
  sign_det *= lu.permutationP().determinant();
  sign_det *= lu.permutationQ().determinant();
  
  if (sign_det <= 0) return fail_value;
  
  double logdet_term = log_abs_det;
  
  // ---- l2 = -(-0.5*logdet + fbetar - 0.5*quad) ----
  double value = -(-0.5 * logdet_term + fbetar - 0.5 * quad_term);
  
  if (!std::isfinite(value)) return fail_value;
  
  return value;
}