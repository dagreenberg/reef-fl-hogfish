functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
   int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_site; //number of sites in REEF
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_hab; //number of habitat classes - REEF
  int<lower=1,upper=N_hab> hab[N]; // vector of habitat class identities
  int<lower=0> N_mth; //number of months - REEF
  int<lower=1,upper=N_mth> mth[N]; // vector of month identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int<lower=0> N_my; //number of monthly clusters
  int<lower=1,upper=N_my> my[N]; // vector of monthly survey cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix for fixed effect covariates`
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
  int<lower=1,upper=N_site*TT> site_year_id[N]; // vector of year-site-id
 int<lower=1,upper=N_dv*TT> dv_year_id[N]; // vector of year-diver-id
}
parameters {
  ordered[K-1] cut; //cutpoints
  real x0; //initial globalpopn size
  row_vector[N_site] a0s; //initial site-level abundance
  row_vector[N_dv] a0d; //initial diver-level abundance

  //deviations from intercept
  vector[Z] beta; //effort coefficients - RVC
  vector[N_hab] a_hab; //deviation between habitats - RVC
  vector[N_dmy] a_dmy; //deviation between daily site survey clusters
  vector[N_my] a_my; //deviation between monthly site survey clusters
  vector[N_mth] a_mth; //deviation among months across years - RVC
  
  //variance on the deviance components
  real<lower = 0> sd_hab;
  real<lower = 0> sd_mth;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_my;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;

  //time-varying estimates by site/diver
  vector<lower=0>[N_site] sd_site_t; 
  vector<lower=0>[N_dv] sd_dv_t; 
  matrix[TT,N_site] z_site;  // process deviations by site
  matrix[TT,N_dv] z_dv;  // process deviations by diver

  //state-space parameters
  vector[TT-1] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
}

transformed parameters{
  matrix[TT,N_dv] a_dv_t;  //annual expected encounters by diver
  matrix[TT,N_site] a_site_t;  //process state expected encounters by site

  vector[TT] x;
  vector[N_yr] a_yr;
  
  a_site_t[1,]=a0s+z_site[1,].*to_row_vector(sd_site_t);
  a_dv_t[1,]=a0d+z_dv[1,].*to_row_vector(sd_dv_t);	
  x[1] = x0;
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t-1]*sd_q;
    a_site_t[t,]=z_site[t,].*to_row_vector(sd_site_t);
    a_dv_t[t,]=z_dv[t,].*to_row_vector(sd_dv_t);

  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]*sd_r; 
  }
 
 
}  

model{
  //priors
  cut ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates - rvc
  x0 ~ normal(0,5); //initial state - pop
  a0s ~ normal(0,5); //initial state - site
  a0d ~ normal(0,5); //initial state - diver
 
 
  //variance terms
  sd_hab ~ gamma(2,1);
  sd_mth ~ gamma(2,1);
  sd_q ~ gamma(2,2);
  sd_r ~ gamma(2,2);
  sd_site_t ~ normal(0,2);
  sd_dv_t ~ normal(0,2);
  sd_dmy ~ gamma(2,1);
  sd_my ~ gamma(2,1);
  
  to_vector(z_site) ~ std_normal();
  to_vector(z_dv) ~ std_normal();
  
  //varying intercepts
  a_hab ~ std_normal();
  a_dmy ~ std_normal();
  a_my ~ std_normal();
  a_mth ~ std_normal();
  
  obs_dev ~ std_normal();
  pro_dev ~ std_normal();
  
  
  y ~ ordered_logistic(a_yr[year_id]+to_vector(a_site_t)[site_year_id]+ to_vector(a_dv_t)[dv_year_id]  + a_hab[hab]*sd_hab+a_dmy[dmy]*sd_dmy+a_my[my]*sd_my+ a_mth[mth]*sd_mth+X*beta,cut); 
}
