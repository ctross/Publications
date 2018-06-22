functions{
  
  int[] IntColHead(int[,] GG, int II, int HH){
    int XX[II];
    for(i in 1:II){
      XX[i] = GG[i,HH];
    }
    return XX;
  }

  vector powZ(vector beta, real x){
    vector[rows(beta)] POW;
    for (i in 1:rows(beta)){
      POW[i] = beta[i]^x;
    }
    return POW;
  }

  vector gammaZ_rng(vector beta, real x) {
    vector[rows(beta)] SAMP;
    for (i in 1:rows(beta)){
      if( beta[i]>0){
        SAMP[i] = gamma_rng(beta[i],x);
        }else{
        SAMP[i] = 0;
      }
    }
    return SAMP;
  }

  vector neg_binomialZ_rng(vector beta, real x) {
    vector[rows(beta)] SAMP;
    for (i in 1:rows(beta)){
      if( beta[i]>0){
        SAMP[i] = neg_binomial_rng(beta[i],x);
        } else {
        SAMP[i] = 0;
      }
    }
    return SAMP;
  }

}

data{

  int P;
  int X;

  int N[P];
  int Proxy[P];

  real MaxExposure;

  int Wives[X,P];
  int RS[X,P];

  matrix[X,P] Rival_1;
  matrix[X,P] Rival_2;
  matrix[X,P] Exposure;

}

transformed data{

  matrix[X,P] Wives_r;

  for(j in 1:X){
    for(i in 1:P){
      Wives_r[j,i]  =  Wives[j,i];
    }
  }

}

parameters {

// Top-level Parameters
  vector[2] MU_Rival;
  vector<lower=0>[2] Sigma_Rival;
  cholesky_factor_corr[2] L_Rival;

  vector[3] MU_Wives;
  vector<lower=0>[3] Sigma_Wives;
  cholesky_factor_corr[3] L_Wives;

  vector[4] MU_RS;
  vector<lower=0>[4] Sigma_RS;
  cholesky_factor_corr[4] L_RS;

  vector<lower=0>[P] B_Rival_1;
  vector<lower=0>[P] B_Rival_2;
  vector<lower=0>[P] B_Wives;
  vector<lower=0>[P] B_RS;

  vector<lower=0, upper=1>[P] Trunc;

// Population Parameters
  vector[2] Beta_Rival_1_raw[P];
  vector[2] Beta_Rival_2_raw[P];
  vector[3] Beta_Wives_raw[P];
  vector[4] Beta_RS_raw[P];

// Shadow Price Parameters
  vector[P] ShadowPrice;

}

transformed parameters{

// Population Parameters
  vector[2] Beta_Rival_1[P];
  vector[2] Beta_Rival_2[P];

  vector[3] Beta_Wives[P];
  vector[4] Beta_RS[P];

// Regression Priors
  for(p in 1:P){
    Beta_Rival_1[p] = MU_Rival + diag_pre_multiply(Sigma_Rival, L_Rival)*Beta_Rival_1_raw[p];
    Beta_Rival_2[p] = MU_Rival + diag_pre_multiply(Sigma_Rival, L_Rival)*Beta_Rival_2_raw[p];
    Beta_Wives[p]   = MU_Wives + diag_pre_multiply(Sigma_Wives, L_Wives)*Beta_Wives_raw[p];
    Beta_RS[p]      = MU_RS    + diag_pre_multiply(Sigma_RS, L_RS)*Beta_RS_raw[p];
  }

}


model{

// Top-level Priors
  MU_Rival ~ normal(0,5);
  Sigma_Rival ~ cauchy(0,5);
  B_Rival_1 ~ cauchy(0,5);
  B_Rival_2 ~ cauchy(0,5);

  MU_Wives ~ normal(0,5);
  Sigma_Wives ~ cauchy(0,5);
  B_Wives ~ cauchy(0,5);

  MU_RS ~ normal(0,5);
  Sigma_RS ~ cauchy(0,5);
  B_RS ~ cauchy(0,5);

  L_Rival ~ lkj_corr_cholesky(2);
  L_Wives ~ lkj_corr_cholesky(2);
  L_RS ~ lkj_corr_cholesky(2);

  for(p in 1:P){
    Beta_Rival_1_raw[p] ~ normal(0, 1);
    Beta_Rival_2_raw[p] ~ normal(0, 1);
    Beta_Wives_raw[p]   ~ normal(0, 1);
    Beta_RS_raw[p]      ~ normal(0, 1);
  }

  ShadowPrice ~ normal(0,1);
  Trunc ~ beta(1,1);

  for( i in 1:P){
    {
      vector[N[i]] Rival;
      vector[N[i]] Mu_Rival_1;
      vector[N[i]] Mu_Rival_2;
      vector[N[i]] Mu_Wives;
      vector[N[i]] Mu_RS;

      vector[N[i]] s_Rival_1;
      vector[N[i]] s_Rival_2;

      vector[N[i]] WivesTrunc;

      int s_Wives [N[i]];
      int s_RS [N[i]];

      // Mean Vectors
      Mu_Rival_1 = exp(Beta_Rival_1[i,1] + Beta_Rival_1[i,2]*log(head(col(Exposure,i),N[i])) );

      if(Proxy[i]==1){
        Mu_Rival_2 = exp(Beta_Rival_2[i,1] + Beta_Rival_2[i,2]*log(head(col(Exposure,i),N[i])) );
      }

      Rival = head(col(Rival_1,i),N[i]) + exp(ShadowPrice[i])*head(col(Rival_2,i),N[i]);
      for(j in 1:N[i]){
        WivesTrunc[j] = if_else(Wives[j,i]==0, Trunc[i], Wives_r[j,i]) ;
      }

      Mu_Wives = exp(Beta_Wives[i,1] + Beta_Wives[i,2]*log(Rival) + Beta_Wives[i,3]*log(head(col(Exposure,i),N[i])) );
      Mu_RS = exp(Beta_RS[i,1]    + Beta_RS[i,2]*log(Rival)   + Beta_RS[i,3]*log(WivesTrunc) + Beta_RS[i,4]*log(head(col(Exposure,i),N[i])) );

      // Main Regressions
      s_Rival_1 = head(col(Rival_1,i),N[i]);
      s_Rival_2 = head(col(Rival_2,i),N[i]);

      s_Wives = IntColHead(Wives,N[i],i);
      s_RS = IntColHead(RS,N[i],i);

      s_Rival_1 ~ gamma(Mu_Rival_1*B_Rival_1[i],B_Rival_1[i]);
      
      if(Proxy[i]==1){
        s_Rival_2 ~ gamma(Mu_Rival_2*B_Rival_2[i],B_Rival_2[i]);
      }

      s_Wives ~ neg_binomial(Mu_Wives*B_Wives[i],B_Wives[i]);
      s_RS ~ neg_binomial(Mu_RS*B_RS[i],      B_RS[i]   );
    }
  }

}
