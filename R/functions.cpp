#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat mvnrnd(arma::vec mu, arma::mat sigma, int cases) {
  // Get dimensions
  int p = mu.n_elem;
  int n_row = sigma.n_rows;
  int n_col = sigma.n_cols;
  
  // Input validation
  if(n_row != n_col) {
    stop("Sigma must be square");
  }
  if(n_row != p) {
    stop("The length of mu must equal the number of rows in sigma");
  }
  // Make sigma symmetric to handle numerical errors
  arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
  // Try Cholesky decomposition first (fastest method)
  arma::mat L;
  bool chol_success = arma::chol(L, sigma_sym, "lower");
  if(!chol_success) {
    // Cholesky failed, try with small jitter
    double jitter = 1e-10 * arma::trace(sigma_sym);
    if(jitter < 1e-12) jitter = 1e-12;  // Minimum jitter
    arma::mat sigma_jitter = sigma_sym + jitter * arma::eye<arma::mat>(p, p);
    chol_success = arma::chol(L, sigma_jitter, "lower");
    if(!chol_success) {
      // Still failed, use eigendecomposition (slow but robust)
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, sigma_sym);
      // Filter out negative/tiny eigenvalues
      double tol = arma::max(eigval) * p * arma::datum::eps;
      arma::uvec pos_idx = arma::find(eigval > tol);
      if(pos_idx.n_elem == 0) {
        // All eigenvalues are essentially zero - return mean
        arma::mat result(cases, p);
        result.each_row() = mu.t();
        return result;
      }
      arma::vec eigval_pos = eigval(pos_idx);
      arma::mat eigvec_pos = eigvec.cols(pos_idx);
      // Construct transformation matrix: L = eigvec * diag(sqrt(eigval))
      L = eigvec_pos * arma::diagmat(arma::sqrt(eigval_pos));
    }
  }
  // Generate standard normal random numbers
  arma::mat Z(cases, L.n_cols, arma::fill::randn);
  // Transform: X = Z * L' + mu
  arma::mat result = Z * L.t();
  result.each_row() += mu.t();
  return result;
}

// Durbin-Koopman (2002) simulation smoother
// [[Rcpp::export]]
Rcpp::List dk_simulation_smoother_cpp(const arma::mat& H,
                                             const arma::mat& F_mat,
                                             const arma::mat& Q,
                                             const arma::mat& dat,
                                             const arma::mat& idx,
                                             const arma::uvec& ENT_r,
                                             const arma::uvec& blockStarts,
                                             const arma::vec& a1) {
  int NS = H.n_cols;
  int Ny = H.n_rows;
  int tstar = dat.n_cols;
  arma::uvec km = ENT_r - 1;
  int nShocks = km.n_elem;
  int nBlocks = blockStarts.n_elem;

  arma::uvec blockEnds(nBlocks);
  for (int b = 0; b < nBlocks - 1; b++) blockEnds(b) = blockStarts(b + 1) - 1;
  blockEnds(nBlocks - 1) = NS - 1;
  std::vector<arma::mat> F_blocks(nBlocks);
  std::vector<arma::mat> Ft_blocks(nBlocks);
  for (int b = 0; b < nBlocks; b++) {
    F_blocks[b] = F_mat.submat(blockStarts(b), blockStarts(b),
                                blockEnds(b), blockEnds(b));
    Ft_blocks[b] = F_blocks[b].t();
  }

  arma::vec shock_sd(nShocks);
  for (int j = 0; j < nShocks; j++) shock_sd(j) = std::sqrt(Q(km(j), km(j)));

  // Precompute non-zero block indices for each row of H
  std::vector<std::vector<int>> varBlocks(Ny);
  for (int i = 0; i < Ny; i++) {
    for (int b = 0; b < nBlocks; b++) {
      for (arma::uword c = blockStarts(b); c <= blockEnds(b); c++) {
        if (H(i, c) != 0.0) { varBlocks[i].push_back(b); break; }
      }
    }
  }

  // Pass 1: KF on actual data, store P_t, K_t, F_inv_t
  arma::vec Sp = a1;
  arma::mat Pp = arma::eye(NS, NS);

  std::vector<arma::mat> Finv_store(tstar);
  std::vector<arma::mat> K_store(tstar);
  std::vector<arma::mat> P_store(tstar);
  std::vector<arma::uvec> obs_store(tstar);

  for (int t = 0; t < tstar; t++) {
    P_store[t] = Pp;
    arma::uvec obs = arma::find(idx.row(t) > 0.5);
    obs_store[t] = obs;
    arma::mat Hit = H.rows(obs);

    // Sparse HtP: each row of H touches only 2-3 blocks
    arma::mat HtP(NS, obs.n_elem, arma::fill::zeros);
    for (arma::uword j = 0; j < obs.n_elem; j++) {
      int vi = (int)obs(j);
      for (int bid : varBlocks[vi]) {
        arma::uword s = blockStarts(bid), e = blockEnds(bid);
        HtP.col(j) += Pp.cols(s, e) * H.submat(vi, s, vi, e).t();
      }
    }
    arma::mat Ft = Hit * HtP;
    arma::mat Ft_inv = arma::inv_sympd(Ft);
    // Block-sparse Kt: F_mat is block-diagonal
    arma::mat FHtP(NS, obs.n_elem, arma::fill::zeros);
    for (int b = 0; b < nBlocks; b++) {
      arma::uword s = blockStarts(b), e = blockEnds(b);
      FHtP.rows(s, e) = F_blocks[b] * HtP.rows(s, e);
    }
    arma::mat Kt = FHtP * Ft_inv;

    Finv_store[t] = Ft_inv;
    K_store[t] = Kt;

    arma::vec y_t(dat(obs, arma::uvec{(arma::uword)t}));
    arma::vec vt = y_t - Hit * Sp;
    arma::mat Kf = HtP * Ft_inv;
    arma::vec Stt = Sp + Kf * vt;
    // Restructured Ptt: compute Hit*Pp with sparse H, avoid O(NS^3)
    arma::mat HitPp(obs.n_elem, NS, arma::fill::zeros);
    for (arma::uword j = 0; j < obs.n_elem; j++) {
      int vi = (int)obs(j);
      for (int bid : varBlocks[vi]) {
        arma::uword s = blockStarts(bid), e = blockEnds(bid);
        HitPp.row(j) += H.submat(vi, s, vi, e) * Pp.rows(s, e);
      }
    }
    arma::mat Ptt = Pp - Kf * HitPp;
    Ptt = 0.5 * (Ptt + Ptt.t());

    if (t < tstar - 1) {
      Sp.zeros();
      for (int b = 0; b < nBlocks; b++) {
        int s = blockStarts(b), e = blockEnds(b);
        Sp.subvec(s, e) = F_blocks[b] * Stt.subvec(s, e);
      }
      Pp.zeros();
      for (int bi = 0; bi < nBlocks; bi++) {
        int si = blockStarts(bi), ei = blockEnds(bi);
        for (int bj = bi; bj < nBlocks; bj++) {
          int sj = blockStarts(bj), ej = blockEnds(bj);
          arma::mat result = F_blocks[bi] * Ptt.submat(si, sj, ei, ej) * F_blocks[bj].t();
          Pp.submat(si, sj, ei, ej) = result;
          if (bi != bj) Pp.submat(sj, si, ej, ei) = result.t();
        }
      }
      for (int j = 0; j < nShocks; j++) Pp(km(j), km(j)) += Q(km(j), km(j));
    }
  }

  // Forward simulation: draw alpha_plus
  arma::mat alpha_plus(tstar, NS, arma::fill::zeros);
  for (int t = 0; t < tstar; t++) {
    arma::vec state(NS, arma::fill::zeros);
    if (t > 0) {
      arma::vec prev = alpha_plus.row(t - 1).t();
      for (int b = 0; b < nBlocks; b++) {
        int s = blockStarts(b), e = blockEnds(b);
        state.subvec(s, e) = F_blocks[b] * prev.subvec(s, e);
      }
    }
    for (int j = 0; j < nShocks; j++) {
      state(km(j)) += shock_sd(j) * arma::randn();
    }
    alpha_plus.row(t) = state.t();
  }

  // Pass 2: cheap KF on y* = y - y_plus
  std::vector<arma::vec> v_star(tstar);
  std::vector<arma::vec> a_star(tstar);
  arma::vec Sp2(NS, arma::fill::zeros);

  for (int t = 0; t < tstar; t++) {
    a_star[t] = Sp2;
    arma::uvec& obs = obs_store[t];
    arma::mat Hit = H.rows(obs);
    arma::vec alpha_t = alpha_plus.row(t).t();
    arma::vec yplus_obs = Hit * alpha_t;
    arma::vec y_t(dat(obs, arma::uvec{(arma::uword)t}));
    arma::vec ys = y_t - yplus_obs;

    arma::vec vt = ys - Hit * Sp2;
    v_star[t] = vt;

    if (t < tstar - 1) {
      arma::vec Stt2 = Sp2 + P_store[t] * Hit.t() * Finv_store[t] * vt;
      Sp2.zeros();
      for (int b = 0; b < nBlocks; b++) {
        int s = blockStarts(b), e = blockEnds(b);
        Sp2.subvec(s, e) = F_blocks[b] * Stt2.subvec(s, e);
      }
    }
  }

  // Backward smoother
  arma::vec r(NS, arma::fill::zeros);
  arma::mat alpha_hat(tstar, NS, arma::fill::zeros);

  for (int t = tstar - 1; t >= 0; t--) {
    arma::mat Hit = H.rows(obs_store[t]);

    arma::vec w = Finv_store[t] * v_star[t] - K_store[t].t() * r;
    arma::vec Ftr(NS, arma::fill::zeros);
    for (int b = 0; b < nBlocks; b++) {
      int s = blockStarts(b), e = blockEnds(b);
      Ftr.subvec(s, e) = Ft_blocks[b] * r.subvec(s, e);
    }
    r = Hit.t() * w + Ftr;

    alpha_hat.row(t) = (a_star[t] + P_store[t] * r).t();
  }

  arma::mat Sdraw = alpha_plus + alpha_hat;
  arma::mat YHAT = H * Sdraw.t();
  arma::ivec km_r = arma::conv_to<arma::ivec>::from(km + 1);

  return Rcpp::List::create(
    Rcpp::Named("Z_MAT") = Sdraw,
    Rcpp::Named("YHAT") = YHAT,
    Rcpp::Named("km") = km_r
  );
}