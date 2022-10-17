

# Load libraries ----------------------------------------------------------

  library(tidyverse)
  library(rstan)
  library(loo)
  require(extrafont)
  library(bayesplot)
  library(dplyr)



# Set work environment ----------------------------------------------------

  date <- Sys.Date() %>% format('%y%m%d')
    # used when saving graphs and df's
  
  # Load data:
  
  setwd('~/Box Sync/TILAC/Current Version/')
    # set working directory
  save_dir <- '~/Box Sync/TILAC/Current Version/Figures/HeatShock/'
    # used when saving graphs and df's 
  
  hs_df <- readRDS('~/Box Sync/TILAC/Current Version/Data/HeatShock/cB.rds')
    # read in heat shock data, raw counts
  FC_alpha <- readRDS('~/Box Sync/TILAC/Current Version/Data/HeatShock/Stan/210128HS_boundedFC_alpha.rds')
    # read in log2 of the TILAC values indicating changes in RNA levels 
  mutrates <- readRDS('~/Box Sync/TILAC/Current Version/Data/HeatShock/Stan/210128HS_bounded_mutrates.rds')
    # read in the mutation rates
  
    

# Calculate average reads - for MA plot -----------------------------------

    samp_list_reads <- c('MS1908_01', 'MS1908_02', 'MS1908_03', 'MS1908_04') 
    
    # Get the average read count for plotting:
    log_reads  <- hs_df %>%    
      ungroup() %>%
      filter(sample %in% samp_list_reads) %>%
      group_by(XF) %>%
      filter(!grepl('28S', XF)) %>% 
      filter(!grepl('rRNA', XF)) %>% 
      summarize(log_reads = log10(sum(n)/length(samp_list_reads))) %>% 
      dplyr::rename(gene_name = 'XF')


# Select genes with sufficient sequencing depth  --------------------------
  # Filter out genes that have too low coverage

    tilac_df <- hs_df %>% 
      ungroup() %>% 
      filter(sample %in% c('MS1908_01', 'MS1908_02', 'MS1908_03', 'MS1908_04')) %>% 
      group_by(sample, XF) %>% 
      summarize(reads = sum(n))
  
    filter_df <- tilac_df %>% 
      mutate(logical = ifelse(reads > 200, 1, 0)) %>% 
      dplyr::select(XF, logical) %>% 
      spread(sample, logical) 
   
    filter_df[c("MS1908_01", "MS1908_02", "MS1908_03", "MS1908_04")][is.na(filter_df[c("MS1908_01", "MS1908_02", "MS1908_03", "MS1908_04")])] <- 0
    
    genes_filter <- filter_df %>% 
      mutate(rowsum = MS1908_01 + MS1908_02 + MS1908_03 + MS1908_04) %>% 
      filter(rowsum >= 2) 
      
      

      
#########################################    
######  Work done by Isaac Vock  ########
#########################################  
    

# Calculate a p-value -----------------------------------------------------

# in this case, we test a null-hypothesis that the fold change in expression 
# is < a conservative value of 0.5, following the test-statistic defined in 
# TREAT. We adjust that calculated p-value for multiple testing using the 
# Benjamini-Hochberg procedure and a false discovery rate of 0.05
    
  nconds <- 1
    #Number of experimental conditions (i.e. only hs)
  
  FDR_control <- 0.05
    #Set desired FDR control
  
  pop_sd <- 0.2 #Only relevant if mu_cutoff == 0
    #Population standard deviation, inferred by comparing 
    #replicates with Stan model
  
  mu_cutoff <- 0.5
    #Define absolute value of significance cutoff if you
    #don't have a population standard deviation estimate
    #from comparing replicates.
  
  
  nreps <- rep(1, times=nconds) 
    #Number of replicates for each experimental condition
    #Right now TILAC does not use replicates, it pools all 
    #replicates into one analysis 
  
  eff_gauss <- FC_alpha[, c("50%","mean", "sd")]
  eff_gauss <- as.data.frame(eff_gauss)
    #Pull out mean and standard deviation of parameter estimate    
    #Convert to data frame:
  
  ngs <- nrow(eff_gauss)/nconds
    #Calculate number of features
  
  pvals <- matrix(0, nrow=ngs*nconds, ncol=5)
  pvals[,2] <- rep(seq(from=1, to=ngs), each=nconds)
  pvals[,3] <- rep(seq(from=1, to=nconds), times=ngs)
  pvals[,5] <- eff_gauss$mean
    #Initialize pvalue matrix
      #The first column stores the pvalue
      #The second column stores the feature number ID
      #The third column stores the experimental condition ID
      #The fourth column stores the signifiance status
      #The fifth column stores the effect size mean
  
  
  
  #Calculate p-value assuming Gaussian z-stat
  for (i in 1:nrow(eff_gauss)){
    # use ROPE method if using a fold-change value as null
    if (mu_cutoff > 0){ 
      tstat <- eff_gauss$mean[i]/(eff_gauss$sd[i])
      if(abs(eff_gauss$mean[i]) > mu_cutoff){
        if(eff_gauss$mean[i] > 0){
          delta <- mu_cutoff/(eff_gauss$sd[i])
        }else{
          delta <- -mu_cutoff/(eff_gauss$sd[i])
        }
      }else{
        delta <- eff_gauss$mean[i]/(eff_gauss$sd[i])
      }
      pvals[i,1] <- pnorm(tstat + delta, lower.tail=ifelse(eff_gauss$mean[i] >0, FALSE, TRUE)) + pnorm(tstat - delta, lower.tail=ifelse(eff_gauss$mean[i] >0, FALSE, TRUE))
    }
    
    # or use Null-Hypothesis Statistical Testing-ish method
    else{ 
      if(nreps[pvals[i,3]] < 3){
        Zscore <- (eff_gauss$mean[i])/(eff_gauss$sd[i] + pop_sd)
        pvals[i,1] <- (1 - pnorm(abs(Zscore)))*2
      }else{
        Zscore <- (eff_gauss$mean[i])/(eff_gauss$sd[i] + pop_sd/sqrt(nreps[pvals[i,3]]))
        pvals[i,1] <- (1 - pnorm(abs(Zscore)))*2
      }
    }
  }
  


  
# Benjamini-Hochberg Significance Calling ---------------------------------
  
  #Order the p-values to prep for BH correction
  pvals_ordered <- pvals[order(pvals[,1]),] 
  pvals_ordered <- pvals_ordered[!is.na(pvals_ordered[,1]),]
  
  #Loop to call significance 
  for(k in 1:nconds){ 
    #for each condition 
    test_pvals <- pvals_ordered[pvals_ordered[,3]==k,]
      for(i in 1:ngs){
        #for each gene
        if(test_pvals[i,1] > ((i/ngs)*FDR_control)){
          crit_val <- i-1
          break
        }
      }
      
      if (crit_val == 1){
        Significant <- c()
        Not_sig <- seq(from=1, to = ngs)
      } else {
        Significant <- test_pvals[1:crit_val,2]
        Not_sig <- test_pvals[(crit_val+1):ngs,2]
      }
      
    #If significant, value in 4th column is 1
    pvals_ordered[(pvals_ordered[,3]==k) & (pvals_ordered[,2] %in% Significant),4] <- 1
    #If not significant, value in 4th column is 0
    pvals_ordered[(pvals_ordered[,3]==k) & (pvals_ordered[,2] %in% Not_sig),4] <- 0
    
  }
  
  #FDR adjust p-values
  pvals_adj <- pvals_ordered
  pvals_adj <- pvals_adj[order(pvals_adj[,3]),]
  for(k in 1:nconds){
    test_pvals <- pvals_ordered[pvals_ordered[,3]==k,1]
    pvals_adj[((k-1)*ngs + 1): (ngs*k),1] <- test_pvals*ngs/seq(from=1, to=ngs, by=1)
    for(i in 1:ngs){
      if(i > 1){
        if (pvals_adj[i + (k-1)*ngs,1] < pvals_adj[(i-1) + (k-1)*ngs,1]){
          pvals_adj[i + (k-1)*ngs,1] <- pvals_adj[(i-1) + (k-1)*ngs,1]
        }
      }
    }
  }



  

# Assemble all the information needed to plot into df ---------------------

  cond_choice <- 1
    # if you're working with multiple hypotheses 
  
  p_adj <- pvals_adj[pvals_adj[,3] == cond_choice,1] #FDR adjusted p-values (also called q values)
  gene_num <- pvals_adj[pvals_adj[,3] == cond_choice,2] #Feature/Gene ID number
  gene_id <- FC_alpha$XF[gene_num]
  condition_id <- pvals_adj[pvals_adj[,3] == cond_choice,3] #Experimental condition ID number
  significance <- pvals_adj[pvals_adj[,3] == cond_choice,4] #Significance numerical indicator
  L2FC <- pvals_adj[pvals_adj[,3] == cond_choice,5] #Log effect size

  toplot_df <- data.frame(p_adj, gene_id, condition_id, significance, L2FC)

  toplot_df <- toplot_df %>% 
    mutate(significance = ifelse(significance==1, ifelse(L2FC> 0, "Upregulated", "Downregulated") ,"Not Significant")) %>% 
    left_join(log_reads) %>% 
    filter(gene_id %in% genes_filter$XF) %>% 
    arrange(L2FC)



# Create annotations for the plot -----------------------------------------
  
  up <- toplot_df %>% filter(significance == "Upregulated")
  down <- toplot_df %>% filter(significance == "Downregulated")
  ns <- toplot_df %>% filter(significance == "Not Significant")

  hs_genes <- c("Hsp23", "Hsp68", "Hsp83", "Hsp67Bc", "Hsp70Bb") 

  hs_annotations <- toplot_df %>% filter(gene_id %in% hs_genes)
# label_ns <- volcano_df %>% filter(gene_id %in% c('Actn', "Set1", "AGO1"))
# label_down <- volcano_df %>% filter(gene_id %in% c('Ack', 'Arp5'))


# Make the MA plot --------------------------------------------------------

  pdf(file = paste(save_dir, date, '_HS_MA_contour.pdf', sep = ''), height = 2, width = 2)
    #open up pdf file in which to save the figure 
  
  ma_plot <- ggplot(ns, aes(x = log_reads, y = L2FC)) + 
    geom_point(size = 0.5,  alpha = 2/5, color = "gray85") +
    geom_point(data = up, aes(x = log_reads, y = L2FC), size = 0.75,  alpha = 3/5, color = "green4") +
    geom_point(data = down, aes(x = log_reads, y = L2FC), size = 0.75,  alpha = 3/5, color = "goldenrod3") +
    
    theme(panel.background = element_blank(),
          panel.grid.major.y = element_line(color = 'gray90'),
          panel.grid.major.x = element_blank(),
          axis.line.y = element_line(color = 'black'),
          axis.line.x = element_line(color = 'black'),
          axis.text = element_text(color = 'black', size = 6), 
          axis.title = element_text(size = 7), 
          plot.title = element_text(size=7, hjust = 0.5),
          text=element_text(family="Arial")) +
    
    geom_point(data = hs_annotations, aes(x = log_reads, y = L2FC), color = 'green4',
               size = 0.75, shape = 16) +
    geom_point(data = hs_annotations, aes(x = log_reads, y = L2FC), color = 'black',
               size = 0.75, shape = 1) +

    xlab('log10(read count)') + 
    ylab('TILAC Ratio (heat shocked vs untreated)') + 
    geom_abline(intercept = 0, slope = 0, color = "gray10", width = 30) + 
    
    geom_text(data = volcano_cntl,
              aes(label = gene_id),
              fontface = 'italic',
              color = 'black',
              size = 2,
              check_overlap = T,
              vjust = -0.5) +
    
    ylim(-5, 5) + 
    ggtitle("Heat Shock Response")
  ma_plot + geom_density_2d(data = volcano_df, aes(x = log_reads, y = L2FC), color = 'gray60', bins = 5)
    # add contour lines onto the graph to show where the majority of the data resides 
  dev.off()
    #close pdf file with figure 

