model_code <- vector("list",6)

model_code[[1]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + Beta_males[3]*log(NumberSpouses_males[i]) )*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i])  + Beta_females[3]*log(NumberSpouses_females[i]) )*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
 model_code[[2]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);


//################################################################################################################################
//################################################################################ Effective Spouses by Age - Males - Mate Quality
//################################################################################################################################
//  A little tricky. Men can have up to 5 wives. So we load five separate arrays.
//  If they are unmarried, effective spouse years, ESA, are not added.
//  If data is missing, then integrate over possible random effects using a weight vector.
//  Otherwise, insert the proper random effect.

 for (i in 1:N_males){
  for (a in 1:Age_males[i]){
   if(PartnersByAge1_males[i,a+11]==0){ESA[1] = 0;}
    else{ ESA[1] = 1;}

   if(PartnersByAge2_males[i,a+11]==0){ESA[2] = 0;}
    else{ ESA[2] = 1;}

   if(PartnersByAge3_males[i,a+11]==0){ESA[3] = 0;}
    else{ ESA[3] = 1;}

   if(PartnersByAge4_males[i,a+11]==0){ESA[4] = 0;}
    else{ ESA[4] = 1;}

   if(PartnersByAge5_males[i,a+11]==0){ESA[5] = 0;}
    else{ ESA[5] = 1;}

   EffectiveSpousesByAge_males[i,a] = ESA[1]+ESA[2]+ESA[3]+ESA[4]+ESA[5];
   }
   }

//################################################################################################################################
//############################################################################## Effective Spouses by Age - Females - Mate Quality
//################################################################################################################################
 for (i in 1:N_females){
  for (a in 1:Age_females[i]){
   if(PartnersByAge1_females[i,a+11]==0){EffectiveSpousesByAge_females[i,a] = 0;}
     else{EffectiveSpousesByAge_females[i,a] = 1;}
     }}

//################################################################################################################################
//################################################################################## Effective Partner Years - Males - Mate Timing
//################################################################################################################################
 for (i in 1:N_males){
   Ticker = 0;
 for (a in 1:Age_males[i]){
   Ticker = Ticker + 1*EffectiveSpousesByAge_males[i,a];
   }
   EffectiveSpouseYears_males[i] = Ticker;
   }


//################################################################################################################################
//################################################################################ Effective Partner Years - Females - Mate Timing
//################################################################################################################################
 for (i in 1:N_females){
   Ticker = 0;
 for (a in 1:Age_females[i]){
   Ticker = Ticker + 1*EffectiveSpousesByAge_females[i,a];
   }
   EffectiveSpouseYears_females[i] = Ticker;
   }

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + exp(Beta_males[4])*log(EffectiveSpouseYears_males[i]))*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i]) + exp(Beta_females[4])*log(EffectiveSpouseYears_females[i]))*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
model_code[[3]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);


//################################################################################################################################
//################################################################################ Effective Spouses by Age - Males - Mate Quality
//################################################################################################################################
//  A little tricky. Men can have up to 5 wives. So we load five separate arrays.
//  If they are unmarried, effective spouse years, ESA, are not added.
//  If data is missing, then integrate over possible random effects using a weight vector.
//  Otherwise, insert the proper random effect.

 for (i in 1:N_males){
  for (a in 1:Age_males[i]){
   if(PartnersByAge1_males[i,a+11]==0){ESA[1] = 0;}
    else{ ESA[1] = 1;}

   if(PartnersByAge2_males[i,a+11]==0){ESA[2] = 0;}
    else{ ESA[2] = 1;}

   if(PartnersByAge3_males[i,a+11]==0){ESA[3] = 0;}
    else{ ESA[3] = 1;}

   if(PartnersByAge4_males[i,a+11]==0){ESA[4] = 0;}
    else{ ESA[4] = 1;}

   if(PartnersByAge5_males[i,a+11]==0){ESA[5] = 0;}
    else{ ESA[5] = 1;}

   EffectiveSpousesByAge_males[i,a] = ESA[1]+ESA[2]+ESA[3]+ESA[4]+ESA[5];
   }
   }

//################################################################################################################################
//############################################################################## Effective Spouses by Age - Females - Mate Quality
//################################################################################################################################
 for (i in 1:N_females){
 for (a in 1:Age_females[i]){
   if(PartnersByAge1_females[i,a+11]==0){EffectiveSpousesByAge_females[i,a] = 0;}
     else{EffectiveSpousesByAge_females[i,a] = 1;}
     }
     }

//################################################################################################################################
//################################################################################## Effective Partner Years - Males - Mate Timing
//################################################################################################################################
 for (i in 1:N_males){
   Ticker = 0;
 for (a in 1:Age_males[i]){
   Ticker = Ticker + inv_logit(Delta_males[2]  + MarriageAgeEffect_males[a])*EffectiveSpousesByAge_males[i,a];
   }
   EffectiveSpouseYears_males[i] = Ticker;
   }


//################################################################################################################################
//################################################################################ Effective Partner Years - Females - Mate Timing
//################################################################################################################################
 for (i in 1:N_females){
   Ticker = 0;
 for (a in 1:Age_females[i]){
   Ticker = Ticker + inv_logit(Delta_females[2] + MarriageAgeEffect_females[a])*EffectiveSpousesByAge_females[i,a];
   }
   EffectiveSpouseYears_females[i] = Ticker;
   }

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + exp(Beta_males[4])*log(EffectiveSpouseYears_males[i]))*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i]) + exp(Beta_females[4])*log(EffectiveSpouseYears_females[i]))*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
model_code[[4]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);


//################################################################################################################################
//################################################################################ Effective Spouses by Age - Males - Mate Quality
//################################################################################################################################
//  A little tricky. Men can have up to 5 wives. So we load five separate arrays.
//  If they are unmarried, effective spouse years, ESA, are not added.
//  If data is missing, then integrate over possible random effects using a weight vector.
//  Otherwise, insert the proper random effect.

 for (i in 1:N_males){
  for (a in 1:Age_males[i]){
   if(PartnersByAge1_males[i,a+11]==0){ESA[1] = 0;}
    else{ if(PartnersByAge1_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[1] = sum(Scrap);}
    else{ ESA[1] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge1_males[i,a+11]]);}}

   if(PartnersByAge2_males[i,a+11]==0){ESA[2] = 0;}
    else{ if(PartnersByAge2_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[2] = sum(Scrap);}
    else{ ESA[2] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge2_males[i,a+11]]);}}

   if(PartnersByAge3_males[i,a+11]==0){ESA[3] = 0;}
    else{ if(PartnersByAge3_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[3] = sum(Scrap);}
    else{ ESA[3] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge3_males[i,a+11]]);}}

   if(PartnersByAge4_males[i,a+11]==0){ESA[4] = 0;}
    else{ if(PartnersByAge4_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[4] = sum(Scrap);}
    else{ ESA[4] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge4_males[i,a+11]]);}}

   if(PartnersByAge5_males[i,a+11]==0){ESA[5] = 0;}
    else{ if(PartnersByAge5_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[5] = sum(Scrap);}
    else{ ESA[5] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge5_males[i,a+11]]);}}

   EffectiveSpousesByAge_males[i,a] = ESA[1]+ESA[2]+ESA[3]+ESA[4]+ESA[5];
   }
   }

//################################################################################################################################
//############################################################################## Effective Spouses by Age - Females - Mate Quality
//################################################################################################################################
 for (i in 1:N_females){
  for (a in 1:Age_females[i]){
   if(PartnersByAge1_females[i,a+11]==0){EffectiveSpousesByAge_females[i,a] = 0;}
    else{ if(PartnersByAge1_females[i,a+11]==999){
         for(b in 1:MaxExposure){
         Scrap[b] = MissingDataWeights_females[b]*inv_logit(Phi_females[2] + SpouseAgeEffect_females[b]);}
         EffectiveSpousesByAge_females[i,a] = sum(Scrap);}
     else{EffectiveSpousesByAge_females[i,a] = inv_logit(Phi_females[2] + SpouseAgeEffect_females[PartnersByAge1_females[i,a+11]]);}}
     }}

//################################################################################################################################
//################################################################################## Effective Partner Years - Males - Mate Timing
//################################################################################################################################
 for (i in 1:N_males){
   Ticker = 0;
 for (a in 1:Age_males[i]){
   Ticker = Ticker + 1*EffectiveSpousesByAge_males[i,a];
   }
   EffectiveSpouseYears_males[i] = Ticker;
   }


//################################################################################################################################
//################################################################################ Effective Partner Years - Females - Mate Timing
//################################################################################################################################
 for (i in 1:N_females){
   Ticker = 0;
 for (a in 1:Age_females[i]){
   Ticker = Ticker + 1*EffectiveSpousesByAge_females[i,a];
   }
   EffectiveSpouseYears_females[i] = Ticker;
   }

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + exp(Beta_males[4])*log(EffectiveSpouseYears_males[i]))*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i]) + exp(Beta_females[4])*log(EffectiveSpouseYears_females[i]))*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
model_code[[5]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);


//################################################################################################################################
//################################################################################ Effective Spouses by Age - Males - Mate Quality
//################################################################################################################################
//  A little tricky. Men can have up to 5 wives. So we load five separate arrays.
//  If they are unmarried, effective spouse years, ESA, are not added.
//  If data is missing, then integrate over possible random effects using a weight vector.
//  Otherwise, insert the proper random effect.

 for (i in 1:N_males){
  for (a in 1:Age_males[i]){
   if(PartnersByAge1_males[i,a+11]==0){ESA[1] = 0;}
    else{ if(PartnersByAge1_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[1] = sum(Scrap);}
    else{ ESA[1] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge1_males[i,a+11]]);}}

   if(PartnersByAge2_males[i,a+11]==0){ESA[2] = 0;}
    else{ if(PartnersByAge2_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[2] = sum(Scrap);}
    else{ ESA[2] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge2_males[i,a+11]]);}}

   if(PartnersByAge3_males[i,a+11]==0){ESA[3] = 0;}
    else{ if(PartnersByAge3_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[3] = sum(Scrap);}
    else{ ESA[3] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge3_males[i,a+11]]);}}

   if(PartnersByAge4_males[i,a+11]==0){ESA[4] = 0;}
    else{ if(PartnersByAge4_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[4] = sum(Scrap);}
    else{ ESA[4] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge4_males[i,a+11]]);}}

   if(PartnersByAge5_males[i,a+11]==0){ESA[5] = 0;}
    else{ if(PartnersByAge5_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[5] = sum(Scrap);}
    else{ ESA[5] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge5_males[i,a+11]]);}}

   EffectiveSpousesByAge_males[i,a] = ESA[1]+ESA[2]+ESA[3]+ESA[4]+ESA[5];
   }
   }

//################################################################################################################################
//############################################################################## Effective Spouses by Age - Females - Mate Quality
//################################################################################################################################
 for (i in 1:N_females){
  for (a in 1:Age_females[i]){
   if(PartnersByAge1_females[i,a+11]==0){EffectiveSpousesByAge_females[i,a] = 0;}
    else{ if(PartnersByAge1_females[i,a+11]==999){
         for(b in 1:MaxExposure){
         Scrap[b] = MissingDataWeights_females[b]*inv_logit(Phi_females[2] + SpouseAgeEffect_females[b]);}
         EffectiveSpousesByAge_females[i,a] = sum(Scrap);}
     else{EffectiveSpousesByAge_females[i,a] = inv_logit(Phi_females[2] + SpouseAgeEffect_females[PartnersByAge1_females[i,a+11]]);}}
     }}

//################################################################################################################################
//################################################################################## Effective Partner Years - Males - Mate Timing
//################################################################################################################################
 for (i in 1:N_males){
   Ticker = 0;
 for (a in 1:Age_males[i]){
   Ticker = Ticker + inv_logit(Delta_males[2]  + MarriageAgeEffect_males[a])*EffectiveSpousesByAge_males[i,a];
   }
   EffectiveSpouseYears_males[i] = Ticker;
   }


//################################################################################################################################
//################################################################################ Effective Partner Years - Females - Mate Timing
//################################################################################################################################
 for (i in 1:N_females){
   Ticker = 0;
 for (a in 1:Age_females[i]){
   Ticker = Ticker + inv_logit(Delta_females[2] + MarriageAgeEffect_females[a])*EffectiveSpousesByAge_females[i,a];
   }
   EffectiveSpouseYears_females[i] = Ticker;
   }

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + exp(Beta_males[4])*log(EffectiveSpouseYears_males[i]))*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i]) + exp(Beta_females[4])*log(EffectiveSpouseYears_females[i]))*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
model_code[[6]] <- '
//##################################################################### Data Block
data {
//######################################################################## Indexes
  int N_males;                                    // Indexes
  int N_females;                                  //
  int N_years;                                    //
  real N_yearsreal;                               //
  int  MaxExposure;                               //
  real MaxExposurereal;                           //

  int RS_males[N_males];                          // Time invariant male data
  int Age_males[N_males];                         //
  int YOB_scaled_males[N_males];                  //
  int NumberSpouses_males[N_males];               //
  int Zeros_males[N_males];                       //
  vector[N_males] R_Age_males;                    //

  int PartnersByAge1_males [N_males,100];         // Time varying male data
  int PartnersByAge2_males [N_males,100];         //
  int PartnersByAge3_males [N_males,100];         //
  int PartnersByAge4_males [N_males,100];         //
  int PartnersByAge5_males [N_males,100];         //

  int RS_females[N_females];                      // Time invariant female data
  int Age_females[N_females];                     //
  int YOB_scaled_females[N_females];              //
  int NumberSpouses_females[N_females];           //
  int Zeros_females[N_females];                   //
  vector[N_females] R_Age_females;                //

  int PartnersByAge1_females [N_females,100];     // Time varying female data

  vector[MaxExposure] MissingDataWeights_males;   // Missing Data Weights
  vector[MaxExposure] MissingDataWeights_females; // Missing Data Weights

  vector[N_years] SecularEffectZeros;             // Additional Data
  vector[MaxExposure] MarriageAgeEffectZeros;     //
}

//############################################################### Parameters Block
parameters{
//################################################################ Male Parameters
  vector<lower=0, upper=1>[2] Theta_males;       // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_males;                 // Correlation Decays
  vector<lower=0>[2] Gamma_males;                // Correlation to Covariance Scalars

  vector[2] Delta_males;                         // Quality Model Parameters
  vector[2] Phi_males;                           // Spouse Quality Model Parameters
  vector[4] Beta_males;                          // RS Model Parameters - Married
  vector[2] Deta_males;                          // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_males;                // RS Model Dispersion Parameter
  vector[2] Ceta_males;                          // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_males;     // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_males;   // Random Effect on Age

//################################################################ Female Parameters
  vector<lower=0, upper=1>[2] Theta_females;     // Max Correlations (Spouse Age, Self Age)
  vector<lower=0>[2] Zeta_females;               // Correlation Decays
  vector<lower=0>[2] Gamma_females;              // Correlation to Covariance Scalars

  vector[2] Delta_females;                       // Quality Model Parameters
  vector[2] Phi_females;                         // Spouse Quality Model Parameters
  vector[4] Beta_females;                        // RS Model Parameters - Married
  vector[2] Deta_females;                        // RS Model Parameters - Unmarried
  vector<lower=0>[2] Sigma_females;              // RS Model Dispersion Parameter
  vector[2] Ceta_females;                        // Marriage Model Parameters

  vector[MaxExposure] SpouseAgeEffect_females;   // Random Effect on Spouse Age (Quality Measure)
  vector[MaxExposure] MarriageAgeEffect_females; // Random Effect on Age
}

//############################################################### Model Block
model{
//######################################################## Local Storage
  matrix[MaxExposure,MaxExposure] RhoMarriageAge;

  vector[N_males] EffectiveSpouseYears_males;
  matrix[N_males,MaxExposure] EffectiveSpousesByAge_males;

  vector[N_females] EffectiveSpouseYears_females;
  matrix[N_females,MaxExposure] EffectiveSpousesByAge_females;

  vector[N_males] A_males;
  vector[N_males] B_males;

  vector[N_females] A_females;
  vector[N_females] B_females;

  real G;
  real Ticker;
  vector[5] ESA;
  vector[MaxExposure] Scrap;

//######################################################## Male Priors
  Theta_males ~ beta(10,2);
  Zeta_males  ~ normal(0,5);
  Gamma_males ~ normal(0,5);

  Phi_males   ~ normal(0,5);
  Delta_males ~ normal(0,5);
  Beta_males  ~ normal(0,5);
  Deta_males  ~ normal(0,5);
  Ceta_males  ~ normal(0,5);
  Sigma_males ~ normal(0,10);

//######################################################## Female Priors
  Theta_females ~ beta(10,2);
  Zeta_females  ~ normal(0,5);
  Gamma_females ~ normal(0,5);

  Phi_females   ~ normal(0,5);
  Delta_females ~ normal(0,5);
  Beta_females  ~ normal(0,5);
  Deta_females  ~ normal(0,5);
  Ceta_females  ~ normal(0,5);
  Sigma_females ~ normal(0,10);

//################################################################################################################################
//####################################################################### Construct Gaussian Process for Spouse Age Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[1] * exp( -Zeta_males[1] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                       }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_males ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_males[1]*RhoMarriageAge);

//################################################################################################################################
//##################################################################### Construct Gaussian Process for Spouse Age Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[1] * exp( -Zeta_females[1] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  SpouseAgeEffect_females ~ multi_normal_cholesky( MarriageAgeEffectZeros , Gamma_females[1]*RhoMarriageAge);

//################################################################################################################################
//#################################################################### Construct Gaussian Process for Marriage Year Effect - Males
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_males[2] * exp( -Zeta_males[2] * G);     // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_males ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_males[2]*RhoMarriageAge);

//################################################################################################################################
//################################################################## Construct Gaussian Process for Marriage Year Effect - Females
//################################################################################################################################
 for (i in 1:(MaxExposure-1)){
 for (j in (i+1):MaxExposure){
   G = ((j-i) * (j-i))/(MaxExposurereal*MaxExposurereal);
                RhoMarriageAge[i,j] = Theta_females[2] * exp( -Zeta_females[2] * G); // Estimate Correlations
                RhoMarriageAge[j,i] = RhoMarriageAge[i,j];                           // Fill Other Triangle
                       }}

 for (i in 1:MaxExposure){
                RhoMarriageAge[i,i] = 1;                                             // Fill Diag
                   }

  RhoMarriageAge = cholesky_decompose(RhoMarriageAge);                               // Decompose Rho

  MarriageAgeEffect_females ~ multi_normal_cholesky(MarriageAgeEffectZeros, Gamma_females[2]*RhoMarriageAge);


//################################################################################################################################
//################################################################################ Effective Spouses by Age - Males - Mate Quality
//################################################################################################################################
//  A little tricky. Men can have up to 5 wives. So we load five separate arrays.
//  If they are unmarried, effective spouse years, ESA, are not added.
//  If data is missing, then integrate over possible random effects using a weight vector.
//  Otherwise, insert the proper random effect.

 for (i in 1:N_males){
  for (a in 1:Age_males[i]){
   if(PartnersByAge1_males[i,a+11]==0){ESA[1] = 0;}
    else{ if(PartnersByAge1_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[1] = sum(Scrap);}
    else{ ESA[1] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge1_males[i,a+11]]);}}

   if(PartnersByAge2_males[i,a+11]==0){ESA[2] = 0;}
    else{ if(PartnersByAge2_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[2] = sum(Scrap);}
    else{ ESA[2] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge2_males[i,a+11]]);}}

   if(PartnersByAge3_males[i,a+11]==0){ESA[3] = 0;}
    else{ if(PartnersByAge3_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[3] = sum(Scrap);}
    else{ ESA[3] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge3_males[i,a+11]]);}}

   if(PartnersByAge4_males[i,a+11]==0){ESA[4] = 0;}
    else{ if(PartnersByAge4_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[4] = sum(Scrap);}
    else{ ESA[4] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge4_males[i,a+11]]);}}

   if(PartnersByAge5_males[i,a+11]==0){ESA[5] = 0;}
    else{ if(PartnersByAge5_males[i,a+11]==999){
         for(b in 1:MaxExposure){
          Scrap[b] = MissingDataWeights_males[b]*inv_logit(Phi_males[2] + SpouseAgeEffect_males[b]);}
          ESA[5] = sum(Scrap);}
    else{ ESA[5] = inv_logit(Phi_males[2] + SpouseAgeEffect_males[PartnersByAge5_males[i,a+11]]);}}

   EffectiveSpousesByAge_males[i,a] = ESA[1]+ESA[2]+ESA[3]+ESA[4]+ESA[5];
   }
   }

//################################################################################################################################
//############################################################################## Effective Spouses by Age - Females - Mate Quality
//################################################################################################################################
 for (i in 1:N_females){
  for (a in 1:Age_females[i]){
   if(PartnersByAge1_females[i,a+11]==0){EffectiveSpousesByAge_females[i,a] = 0;}
    else{ if(PartnersByAge1_females[i,a+11]==999){
         for(b in 1:MaxExposure){
         Scrap[b] = MissingDataWeights_females[b]*inv_logit(Phi_females[2] + SpouseAgeEffect_females[b]);}
         EffectiveSpousesByAge_females[i,a] = sum(Scrap);}
     else{EffectiveSpousesByAge_females[i,a] = inv_logit(Phi_females[2] + SpouseAgeEffect_females[PartnersByAge1_females[i,a+11]]);}}
     }}

//################################################################################################################################
//################################################################################## Effective Partner Years - Males - Mate Timing
//################################################################################################################################
 for (i in 1:N_males){
   Ticker = 0;
 for (a in 1:Age_males[i]){
   Ticker = Ticker + inv_logit(Delta_males[2]  + MarriageAgeEffect_males[a])*EffectiveSpousesByAge_males[i,a];
   }
   EffectiveSpouseYears_males[i] = Ticker;
   }


//################################################################################################################################
//################################################################################ Effective Partner Years - Females - Mate Timing
//################################################################################################################################
 for (i in 1:N_females){
   Ticker = 0;
 for (a in 1:Age_females[i]){
   Ticker = Ticker + inv_logit(Delta_females[2] + MarriageAgeEffect_females[a])*EffectiveSpousesByAge_females[i,a];
   }
   EffectiveSpouseYears_females[i] = Ticker;
   }

//################################################################################################################################
//###################################################################################################### Main Regression Functions
//################################################################################################################################
//###### Male Link Functions
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i] = 1/Sigma_males[1];
   A_males[i] = exp(Beta_males[1] + Beta_males[2]*log(Age_males[i])  + Beta_males[3]*log(NumberSpouses_males[i]) + exp(Beta_males[4])*log(EffectiveSpouseYears_males[i]))*B_males[i];
   } else{
   B_males[i] = 1/Sigma_males[2];
   A_males[i] = exp(Deta_males[1] + Deta_males[2]*log(Age_males[i]))*B_males[i];
   }
   }

//###### Female Link Functions
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i] = 1/Sigma_females[1];
   A_females[i] = exp(Beta_females[1] + Beta_females[2]*log(Age_females[i])  + Beta_females[3]*log(NumberSpouses_females[i]) + exp(Beta_females[4])*log(EffectiveSpouseYears_females[i]))*B_females[i];
    } else{
   B_females[i] = 1/Sigma_females[2];
   A_females[i] = exp(Deta_females[1] + Deta_females[2]*log(Age_females[i]))*B_females[i];
   }
   }

//###### Model Outcomes
   RS_males ~ neg_binomial(A_males,B_males);
   RS_females ~ neg_binomial(A_females,B_females);

   Zeros_males ~ bernoulli_logit(1-(Ceta_males[1] + Ceta_males[2]*log(R_Age_males)));
   Zeros_females ~ bernoulli_logit(1-(Ceta_females[1] + Ceta_females[2]*log(R_Age_females)));

}
 '
 
 
 

 
 
 
 
 
 
 
 
 
 
 