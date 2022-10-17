
# R script to take TILAC mutation data and calculate differential expression
# Used to get data in format for Stan model 


# Load packages ---------------------------------------------------------

  require(tidyverse)
  library(rstan)
  library(loo)
  require(extrafont)
  library(bayesplot)

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  date <- Sys.Date() %>% format('%y%m%d')
  

# Set up environment ------------------------------------------------------
  
  setwd('~/scratch60/2008_Exp119b/210126_stan_bounded_log2FC')
  
  # Which stan model are we using? 
  fl = 'tilac_pois.stan'
  

# Set up everything for the automatic save.  ------------------------------
  
  basename <- 'Exp119b_bounded' #don't add date, it will be added when files are made

  save_data <- paste(basename, '_data.rds', sep = '')                      # save data that goes into Stan
  save_fit <- paste(basename, '_fit.rds', sep = '')                        # save fit that comes out of Stan
  summary_save <- paste(basename, '_fitsummary.rds', sep = '')             # save the processed summary of the fit object 
  mutrates_save <- paste(basename, '_mutrates.rds', sep = '')              # save the expected rate of mutations, sub-element of summary 
  fraction_label_save <- paste(basename, '_fraction.rds', sep = '')        # save the fraction of signal from exp/cntl, sub-element of summary 
  alpha_save <- paste(basename, '_alpha.rds', sep = '')                    # save the TILAC ratio, sub-element of summary 
  FC_alpha_save <- paste(basename, '_FC_alpha.rds', sep = '')              # save the log2 of the TILAC ratio (log2FoldChange)  
  pairs_save_lambda <- paste(basename, '_pairs_lambda.pdf', sep = '')      # save pairs plots of expected rate of mutations
  pairs_save_alpha1 <- paste(basename, '_pairs_alpha1.pdf', sep = '')      # save pairs plots of frac_labele/unlabed and TILAC Ratio of randomly selected genes
  pairs_save_alpha2 <- paste(basename, '_pairs_alpha2.pdf', sep = '')      # save pairs plots of frac_labele/unlabed and Log2(TILAC Ratio) for randomly selected genes 
  rhats_save <- paste(basename, '_rhats.pdf', sep = '')                    # save histogram of rhats
  neff_save <- paste(basename, '_neff.pdf', sep = '')                      # save histogram of n effective 
  trace_save <- paste(basename, '_trace.pdf', sep = '')                    # save trace plots for expected mutation rates and randomly selected genes 





# Load the data -----------------------------------------------------------

  cB <- readRDS('~/scratch60/2008_Exp119b/cB.rds')
  

# Get reliable features ---------------------------------------------------
  # reliable features are those that do not have high background mutation rates 
  # and there is the option to set a read cut-off which 

  # for the paper I used all of the data

  reliableFeatures <- function(df = cB,
                               cv = c('MS200702_09', 'MS200702_10'),  # unfed samples 
                               meancut = 0.3, 
                               totcut = 200){
    y <- df %>%
      ungroup() %>%
      filter(sample %in% cv,
             !(grepl('__', XF))) %>%
      group_by(sample, XF) %>%
      mutate(totcounts = sum(n)) %>%
      summarize(mean_mut_TC = sum(TC*n)/sum(nT),
                mean_mut_GA = sum (GA*n)/sum(nG),
                totcounts = sum(n)) %>%
      filter(mean_mut_TC < meancut,
             mean_mut_GA < meancut,
             totcounts >= totcut) %>%
      ungroup() %>%
      select(XF) %>%
      unlist() %>% 
      unique()
    return(y)
  }


  keep <- reliableFeatures(meancut = 0.8, totcut = 0)


# What samples are being evaluated, and what type are they? 
  samp_list <- c( 'MS200702_01', 'MS200702_02', 'MS200702_05', 'MS200702_06', 'MS200702_09', 'MS200702_10')  
  type_list <- c(1, 1, 1, 1, 0, 0)  # experimental or unfed 
  label_list <- c(1, 2, 1, 2, 0, 0) # is it forward or reverse
  names(type_list) <- samp_list
  names(label_list) <- samp_list  



# Make a key that will relate gene name (XF) with a numerical identifier used in the model:
  ranked_features_df  <- cB %>%    
    ungroup() %>%
    filter(XF %in% keep) %>%  #keep only reliable features 
    group_by(XF) %>%          #set grouping for summary step below 
    summarize(n = sum(n)) %>%  # sum reads over all samples for genes in XF
    mutate(fnum = order(-n)) %>%  # make a column fnum, to give each XF a numeral identifier 
    arrange(fnum) %>%
    dplyr::select(XF, fnum) # select only XF and fnum, linking gene name to numeral identifier 

# Make a small dataframe that includes only the genes from previous step, that are linked 
# to numerical identifier 
  sdf <- cB %>%
    ungroup() %>%
    right_join(ranked_features_df, by = 'XF')

# calculate the # of reads with a given # of TC mutations
# rename nT and TC to be trials and muts
# add column identifying these rows as "TC" rows 
  tcdf <- sdf %>%
    mutate(trials = nT, muts = TC) %>%
    group_by(XF, fnum, sample, muts) %>%
    summarize(n = sum(n)) %>%
    mutate(mut_type = 1)

# calculate the # of reads with a given # of GA mutations
# rename nG and GA to be trials and muts
# add column identifying these rows as "GA" rows   
  gadf <- sdf %>%
    mutate(trials = nG, muts = GA) %>%
    group_by(XF, fnum, sample, muts) %>%
    summarize(n = sum(n)) %>%
    mutate(mut_type = 2)

# bind TC and GA dataframes into one long dataframe
# select out only samples to be analyzed 
  df <- bind_rows(tcdf, gadf) %>%
    ungroup() %>%
    filter(sample %in% slist)
  
  # slist = samp_list 
  # tlist = type_list
  # llist = label_list
  # keep = keep  

# Functions to map the sample name to the type (fed or unfed) and the label combination (fwd or rev)
  getType <- function(s) type_list[paste(s)]
  getLabel <- function(s) label_list[paste(s)]

  df$type <- paste(df$sample) %>% map_dbl(function(x) getType(x))
  df$type <- as.integer(df$type) 
  df$label <- paste(df$sample) %>% map_dbl(function(x) getLabel(x))
  df$label <- as.integer(df$label)
 
# rearrange df  
  df <- df %>%
    ungroup() %>%
    group_by(fnum) %>%
    arrange(type, .by_group = TRUE)
  
  
    
# Start enumerating parameters that Stan model will need   
  NE <- dim(df)[1]  # number of elements in the dataframe 
  NF <- length(keep)  # number of features (genes, XF, fnum, all the same)
  FN <- df$fnum
  TP <- df$type
  ML <- ifelse(df$mut_type == df$label, 1, 2)  # based on FWD/REV and TC or GA, should this line update exp or cntl fractions labeled
  MT <- df$mut_type
  num_mut <- df$muts # Number of mutations
  num_obs <- df$n # Number of times observed

  
  data_list <- list(
    NE = NE, # Number of reads
    NF = NF, # Number of features/transcripts
    TP = TP, # Experimental or unfed  
    FN = FN, # feature number (that correlates with transcript name)
    ML = ML, # "mutation label" ie do we think read is from experimental or control part of TILAC sample
    MT = MT, # mutation type 
    N_obs = N_obs, # Number of observed reads per feature
    N_ctl = N_ctl, # Number of control reads per feature
    FE = df$fnum, # Number of reads per feature
    num_mut = num_mut, # Number of mutations
    num_obs = num_obs # Number of times observed
    )
  
# save data that's going into stan model 
  data <- list(df, data_list)
  saveRDS(data, paste(date, save_data, sep='_'))
  

  
# Fit the data using stan -----------------------------------------------------------
  
# Model run parameters: 
  iter = 4000
  warmup = 1000
  chains = 4

  
  fit <- stan(file = fl,
              data = data_list,
              iter = iter,
              warmup = warmup,
              chains = chains)

   # save full data fit 
   saveRDS(fit, paste(date, save_fit, sep = '_'))





# Munge fit object into a smaller object 

   fit_summary <- summary(fit,
                          probs = c(0.1, 0.5, 0.9))$summary
   fit_summary <- round(fit_summary, 4)

   idf <- as.data.frame(fit_summary) #real usable dataframe

# This takes the fit summary and gets only the desired parameters
    
  #get alpha df 
   adf <- idf %>%
      rownames_to_column(var = 'term') %>%
      filter(str_detect(term, "alpha"),
             !str_detect(term, "logit_frac_"),
	     !str_detect(term, "FC")) %>%
      mutate(fnum = str_remove(term, 'alpha\\[')) %>%
      mutate(fnum = str_remove(fnum, '\\]'))
 
  #get FC_alpha df 
   adf <- idf %>%
      rownames_to_column(var = 'term') %>%
      filter(str_detect(term, "FC_alpha"),
             !str_detect(term, "logit_frac_")) %>%
      mutate(fnum = str_remove(term, 'FC_alpha\\[')) %>%
      mutate(fnum = str_remove(fnum, '\\]'))

   #get fraction df
    fdf <- idf %>%
      rownames_to_column(var = 'term') %>%
      filter(str_detect(term, "logit_frac_"), 
             !str_detect(term, "alpha")) %>% 
      mutate(fnum = str_remove(term, 'logit_frac_\\[')) %>%
      mutate(fnum = str_remove(fnum, '\\]'))
   
   #get mutation rates df 
     mdf <- idf %>%
      rownames_to_column(var = 'term') %>%
      filter(str_detect(term, "log"),
             !str_detect(term, "logit_frac_")) %>%
      mutate(fnum = str_remove(term, 'log\\[')) %>%
      mutate(fnum = str_remove(fnum, '\\]'))


# Grab key of FN with gene name ---

    key <- data[[1]] %>% # This gets the fnum and gene names
      dplyr::select(XF, fnum) %>% 
      unique()

    key$fnum <- as.character(key$fnum) #fnum needs to be as a character in order to do the join 



# Join small df with key to put model results with appropriate gene name ---

    mut_rates <- left_join(mdf, key, by = 'fnum')
    saveRDS(mut_rates, paste(date, mutrates_save, sep = ''))

    fractions <- left_join(fdf, key, by = 'fnum')
    saveRDS(fractions, paste(date, fraction_label_save, sep = ''))

    alpha <- left_join(adf, key, by = 'fnum')
    saveRDS(alpha, paste0(date, alpha_save, sep = ''))
    
    FC_alpha <- left_join(fadf, key, by = 'fnum')
    saveRDS(FC_alpha, paste0(date, FC_alpha_save, sep = ''))




# Plotting model validation -----------------------------------------------------

# Make pairs plot 
  pdf(file = paste(date, sep = "_", pairs_save_lambda))
        pplot <- pairs(fit, pars = c('log_lambda_o_tc', 'log_lambda_o_ga',
                                     'log_lambda_n_tc', 'log_lambda_n_ga' ))
        dev.off()


  pdf(file = paste(date, sep = "_", pairs_save_alpha1))
        pplot <- pairs(fit, pars = c('logit_frac_new_exp[1]', 'logit_frac_new_ctl[1]',
                                     'alpha[1]',
                                     'logit_frac_new_exp[2]', 'logit_frac_new_ctl[2]',
                                     'alpha[2]',
                                     'logit_frac_new_exp[3]', 'logit_frac_new_ctl[3]',
                                     'alpha[3]',
                                     'logit_frac_new_exp[4]', 'logit_frac_new_ctl[4]',
                                      'alpha[4]'))
        dev.off()

  pdf(file = paste(date, sep = "_", pairs_save_alpha2))
        pplot <- pairs(fit, pars = c('logit_frac_new_exp[100]', 'logit_frac_new_ctl[100]',
                                     'FC_alpha[100]',
                                     'logit_frac_new_exp[200]', 'logit_frac_new_ctl[200]',
                                     'FC_alpha[200]',
                                     'logit_frac_new_exp[300]', 'logit_frac_new_ctl[300]',
                                     'FC_alpha[300]',
                                     'logit_frac_new_exp[400]', 'logit_frac_new_ctl[400]',
                                      'FC_alpha[400]'))
        dev.off()


# Trace plots 
  pdf(file = paste(date, sep = "_", trace_save))
        traceplot(fit, pars = c('log_lambda_o_tc', 'log_lambda_o_ga',
                                     'log_lambda_n_tc', 'log_lambda_n_ga',
                                     'logit_frac_new_exp[1]', 'logit_frac_new_ctl[1]',
                                     'logit_frac_new_exp[2]', 'logit_frac_new_ctl[2]',
                                     'alpha[1]', 'alpha[2]', 'FC_alpha[1]', 'FC_alpha[2]'), 
				inc_warmup = FALSE, nrow = 2) #add others

        dev.off()

# N effective 
  n_eff <- neff_ratio(fit)
    pdf(paste(date, neff_save, sep = '_'))
    ggrat <- hist(n_eff)
  dev.off()


# Rhats 
  rhats <- rhat(fit)
    pdf(paste(date, rhats_save, sep = '_'))
    ggrhat <- hist(rhats)
  dev.off()

