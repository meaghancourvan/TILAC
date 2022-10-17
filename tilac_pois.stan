data {
  int NE; //num entries in data frame.
  int NF; //num features (e.g., transcripts)
  int TP[NE]; // type for each entry (experimental = 1, control = 0 for unlabeled)
  

  int FN[NE]; //feature number for each gene 
  
  int ML[NE]; // Are s4U labeled reads experimental? (yes = 1, no = 2) 
  int MT[NE]; // What type of mutation is this line in the data (1 = s4U, 2 = s6G)
  
  int num_mut[NE]; //  mutations in each entry
  int num_obs[NE]; // number of times
}

parameters {

  real log_lambda_unlabled_tc;  //mutations in unlabeled reads
  real log_lambda_unlabled_ga;  // mutations in unlabeled reads
  
  
  real<lower=0> TL_tc_mut_rate;  // induced TC mutations 
  real<lower=0> TL_ga_mut_rate;  // induced GA mutations 


  vector[NF] logit_frac_labeled_exp; // ie Heat Shock inferred fraction labeled 
  vector[NF] logit_frac_labeled_ctl; // ie control inferred fraction labeled 
}

transformed parameters {
  vector[NF] frac_labeled_exp = inv_logit(logit_frac_labeled_exp);  // unlogit to get frac sample labeled in exp
  vector[NF] frac_labeled_ctl = inv_logit(logit_frac_labeled_ctl);  // unlogit to get frac sample labeled in cntl
  vector[2] l_labeled_mat;  // matrix to hold lambdas for T and G labeled reads 
  vector[2] l_unlabled_mat; // matrix to hold lambdas for T and G unlabeled reads 
  real log_lambda_labeled_tc = log(exp(log_lambda_unlabled_tc) + TL_tc_mut_rate);  // calculated induced expected mutations TC
  real log_lambda_labeled_ga = log(exp(log_lambda_unlabled_ga) + TL_ga_mut_rate);  // calculated induced expected mutations GA
  matrix[NF,2] fn_mat;
  
  l_labeled_mat[1] = log_lambda_labeled_tc;   // put the tc rate into the matrix 
  l_labeled_mat[2] = log_lambda_labeled_ga;   // put the ga rate into the matrix 
  l_unlabled_mat[1] = log_lambda_unlabled_tc;  // put the tc rate into the matrix 
  l_unlabled_mat[2] = log_lambda_unlabled_ga;  // put the ga rate into the matrix 
  
  fn_mat[,1] = frac_labeled_exp;  // proportion of sequencing dataset from experimental sample 
  fn_mat[,2] = frac_labeled_ctl;  // proportion of sequencing dataset from control sample

}


model {
  
  log_lambda_unlabled_tc ~ normal(-2, 2);
  log_lambda_unlabled_ga ~ normal(-2, 2); 
  
  TL_tc_mut_rate ~ exponential(0.5);
  TL_ga_mut_rate ~ exponential(0.5);
  
  logit_frac_labeled_exp ~ normal(0, 1.5);
  logit_frac_labeled_ctl ~ normal(0, 1.5);
  
  for (i in 1:NE) {
    if (TP[i] == 1) 
        target += num_obs[i] * log_mix(fn_mat[FN[i],ML[i]], 
                                      poisson_log_lpmf(num_mut[i] | l_labeled_mat[MT[i]]),
                                       poisson_log_lpmf(num_mut[i] | l_unlabled_mat[MT[i]]));
    else
        target += num_obs[i] * poisson_log_lpmf(num_mut[i] | l_unlabled_mat[MT[i]]);
    } 
}

generated quantities {
  vector[NF] alpha; // 
  vector[NF] FC_alpha; //
  for (j in 1:NF){
     alpha[j] = fn_mat[j,1] / fn_mat[j,2];  // alpha or TILAC ratio = exp/cntl
     FC_alpha[j] = log2(alpha[j]);  // take the log2 to get the differential expression 
  }
}
