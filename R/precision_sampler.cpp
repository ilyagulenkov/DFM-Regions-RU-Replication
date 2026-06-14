// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SimplicialLLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>> CholSolver;

// Build G'*diag(invSig)*G for one AR block, appending triplets with global offset
static void addArPrecisionTriplets(
    const Eigen::VectorXd& psi, double sigma2, int TT,
    double sig2_init, int offset,
    std::vector<Triplet>& triplets
) {
  int Lu = psi.size();
  double invSigma2 = 1.0 / sigma2;
  double invSig2Init = 1.0 / sig2_init;

  for (int t = 0; t < TT; t++) {
    int bw = std::min(Lu, t);
    double w = (t < Lu) ? invSig2Init : invSigma2;
    int nk = bw + 1;
    for (int a = 0; a < nk; a++) {
      int col_a = t - a;
      double val_a = (a == 0) ? 1.0 : -psi(a - 1);
      for (int b = a; b < nk; b++) {
        int col_b = t - b;
        double val_b = (b == 0) ? 1.0 : -psi(b - 1);
        double entry = w * val_a * val_b;
        triplets.push_back(Triplet(offset + col_a, offset + col_b, entry));
        if (col_a != col_b) {
          triplets.push_back(Triplet(offset + col_b, offset + col_a, entry));
        }
      }
    }
  }
}

// Precision sampler using RcppEigen sparse Cholesky
// [[Rcpp::export]]
Eigen::VectorXd precisionSampleCpp(
    const Eigen::VectorXd& phi,
    const Eigen::MatrixXd& psi,
    const Eigen::VectorXd& sig2,
    double sig2f,
    double sig2_tau,
    bool hasTau,
    int N,
    int TT_ext,
    const Rcpp::List& obsCache,
    const Eigen::VectorXd& lamd,
    const Eigen::VectorXd& z,
    double sig2_init
) {
  int m = 1 + (int)hasTau + N;
  int mT = m * TT_ext;
  int Lu = psi.cols();

  // --- Build K (transition precision) ---
  std::vector<Triplet> triplets;
  triplets.reserve(m * TT_ext * (2 * Lu + 1) * 2);

  int offset = 0;
  addArPrecisionTriplets(phi, sig2f, TT_ext, sig2_init, offset, triplets);
  offset += TT_ext;
  if (hasTau) {
    Eigen::VectorXd tauPsi(1);
    tauPsi(0) = 1.0;
    addArPrecisionTriplets(tauPsi, sig2_tau, TT_ext, sig2_init, offset, triplets);
    offset += TT_ext;
  }
  for (int i = 0; i < N; i++) {
    Eigen::VectorXd psi_i = psi.row(i);
    addArPrecisionTriplets(psi_i, sig2(i), TT_ext, sig2_init, offset, triplets);
    offset += TT_ext;
  }

  // --- Build Omega + d (observation system) ---
  int offF = 0;
  int offTau = hasTau ? TT_ext : -1;
  int offUbase = TT_ext * (1 + (int)hasTau);
  Eigen::VectorXd d = Eigen::VectorXd::Zero(mT);

  for (int i = 0; i < N; i++) {
    Rcpp::List cc = obsCache[i];
    Rcpp::IntegerVector sp_i = cc["sp_i"];
    Rcpp::IntegerVector sp_j = cc["sp_j"];
    Rcpp::NumericVector sp_x = cc["sp_x"];
    Rcpp::NumericVector MtY_invR = cc["MtY_invR"];
    int nnz = sp_x.size();
    double li = lamd(i);
    double li2 = li * li;
    int offUi = offUbase + i * TT_ext;

    for (int k = 0; k < nnz; k++) {
      int ri = sp_i[k] - 1;
      int ci = sp_j[k] - 1;
      double xv = sp_x[k];
      triplets.push_back(Triplet(offF + ri, offF + ci, xv * li2));
      triplets.push_back(Triplet(offF + ri, offUi + ci, xv * li));
      triplets.push_back(Triplet(offUi + ri, offF + ci, xv * li));
      triplets.push_back(Triplet(offUi + ri, offUi + ci, xv));
    }
    for (int t = 0; t < TT_ext; t++) {
      d(offF + t) += li * MtY_invR[t];
      d(offUi + t) = MtY_invR[t];
    }
    if (hasTau && i == 0) {
      for (int k = 0; k < nnz; k++) {
        int ri = sp_i[k] - 1;
        int ci = sp_j[k] - 1;
        double xv = sp_x[k];
        triplets.push_back(Triplet(offTau + ri, offTau + ci, xv));
        triplets.push_back(Triplet(offF + ri, offTau + ci, xv * li));
        triplets.push_back(Triplet(offTau + ri, offF + ci, xv * li));
        triplets.push_back(Triplet(offTau + ri, offUi + ci, xv));
        triplets.push_back(Triplet(offUi + ri, offTau + ci, xv));
      }
      for (int t = 0; t < TT_ext; t++) {
        d(offTau + t) = MtY_invR[t];
      }
    }
  }

  // --- Assemble P, factorize, sample ---
  SpMat P(mT, mT);
  P.setFromTriplets(triplets.begin(), triplets.end());
  P.makeCompressed();

  CholSolver chol;
  chol.analyzePattern(P);
  chol.factorize(P);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky factorization failed");
  }

  Eigen::VectorXd mu = chol.solve(d);
  SpMat L = chol.matrixL();
  Eigen::VectorXd w = L.transpose().triangularView<Eigen::Upper>().solve(z);
  Eigen::VectorXd dev = chol.permutationP().transpose() * w;

  return mu + dev;
}
