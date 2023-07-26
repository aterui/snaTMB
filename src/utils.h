namespace snaTMB {

// non-spatial random effect template ------------------------------------------

template <class Type>
struct r_element {
  // Input from R
  int n_fix; // # fixed-effect terms with a random effect
  int n_group; // # group replicates within a random effect
  int n_param; // # parameters in each random effect

  matrix<Type> R; // correlation matrix for non-spatial random effects
};

template <class Type>
struct term_t : vector<r_element<Type> > {
  term_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i = 0; i < LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int n_fix = (int) REAL(getListElement(y, "n_fix", &isNumericScalar))[0];
      int n_group = (int) REAL(getListElement(y, "n_group", &isNumericScalar))[0];
      int n_param = (int) REAL(getListElement(y, "n_param", &isNumericScalar))[0];
      (*this)(i).n_fix = n_fix;
      (*this)(i).n_group = n_group;
      (*this)(i).n_param = n_param;
    }
  }
};

// non-spatial random effect distribution --------------------------------------

template<class Type>
Type dmarginal(array<Type> &V,
               vector<Type> log_theta,
               r_element<Type>& term) {
  Type ans = 0;
  int n = term.n_fix;

  if (n == int(1)) {

    Type sd = exp(log_theta(0));
    for(int i = 0; i < term.n_group; i++){
      ans -= dnorm(vector<Type>(V.col(i)), Type(0), sd, true).sum();
    }

  } else if (n > int(1)) {

    //sd_diag: diagonal SD elements
    //rho: off-diagonal rho elements
    vector<Type> sd_diag = exp(log_theta.head(n));
    vector<Type> rho = log_theta.tail(log_theta.size() - n);

    //dusmnorm: unstructured multivariate normal
    //sdusmnorm: scaled unstructured multivariate normal
    density::UNSTRUCTURED_CORR_t<Type> dusmnorm(rho);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > sdusmnorm = density::VECSCALE(dusmnorm, sd_diag);

    for(int i = 0; i < term.n_group; i++){
      // note: sdusmnorm returns negative log likelihood
      ans += sdusmnorm(V.col(i));
    }

    term.R = dusmnorm.cov();

  }

  return ans;
}

// negative log likelihood for non-spatial random effects ----------------------

template <class Type>
Type mnll(vector<Type> &v,
          vector<Type> log_theta,
          vector<r_element<Type> >& term) {
  Type ans = 0;

  // blsize: block size for each random effect
  // v_id: index for latent variable v
  // theta_id: index for theta vector
  int blsize;
  int v_id = 0;
  int theta_id = 0;

  for (int i = 0; i < term.size(); i++) {
    // block size for each random effect
    blsize = term(i).n_fix * term(i).n_group;

    // array for multi-variate normal
    vector<int> dim(2);
    dim << term(i).n_fix, term(i).n_group;
    array<Type> v_seg(&v(v_id), dim);

    // theta vector; segment length = # parameters for each random effect
    vector<Type> t_seg = log_theta.segment(theta_id, term(i).n_param);

    // likelihood
    ans += dmarginal(v_seg, t_seg, term(i));

    // update v and theta index
    v_id += blsize;
    theta_id += term(i).n_param;
  }

  return ans;
}


// negative log likelihood for spatial random effects --------------------------

template <class Type>
Type snll(vector<Type> &u,
          Type log_phi,
          Type log_lambda,
          matrix<Type> D,
          matrix<Type> W) {
  Type ans = 0;
  int n = u.size();

  // dimension check
  if(!(D.rows() == n && D.cols() == n ))
    error ("Dimension of distance matrix must equal blocksize.");

  Type phi0 = exp(log_phi);
  Type lambda0 = exp(log_lambda);

  // structured correlation matrix
  // W: weight matrix, D: distance matrix
  matrix<Type> Rd(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // ifelse: ensure diagonal elements equal one
      Rd(i, j) = (i == j ? Type(1) : W(i, j) * exp(-lambda0 * D(i, j) ) );
    }
  }

  // dsmnorm: structured multivariate normal
  // sdsmnorm: scaled structured multivariate normal
  density::MVNORM_t<Type> dsmnorm(Rd);
  density::SCALE_t<density::MVNORM_t<Type> > sdsmnorm = density::SCALE(dsmnorm, phi0);

  // note: sdsmnorm returns negative log likelihood
  ans += sdsmnorm(u);
  return ans;
}

}

