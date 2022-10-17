

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
  mutrates <- readRDS('~/Box Sync/TILAC/Current Version/Data/HeatShock/Stan/210128HS_bounded_mutrates.rds')
  
    

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


###########VOLCANO PLOT CODE################
#Goal: generate volcano plots based on previously calculated adjusted p-values 

####ONE USER INPUT FOR VOLCANO PLOTS####
#Which of the experimental conditions will data be plotted for?
cond_choice <- 1
####END OF USER INPUT###################

p_adj <- pvals_adj[pvals_adj[,3] == cond_choice,1] #FDR adjusted p-values (also called q values)
gene_num <- pvals_adj[pvals_adj[,3] == cond_choice,2] #Feature/Gene ID number
gene_id <- FC_alpha$XF[gene_num]
condition_id <- pvals_adj[pvals_adj[,3] == cond_choice,3] #Experimental condition ID number
significance <- pvals_adj[pvals_adj[,3] == cond_choice,4] #Significance numerical indicator
L2FC <- pvals_adj[pvals_adj[,3] == cond_choice,5] #Log effect size

volcano_df <- data.frame(p_adj, gene_id, condition_id, significance, L2FC)

volcano_df <- volcano_df %>% 
  mutate(significance = ifelse(significance==1, ifelse(L2FC> 0, "Upregulated", "Downregulated") ,"Not Significant")) %>% 
  left_join(log_reads) %>% 
  filter(gene_id %in% genes_filter$XF) %>% 
  arrange(L2FC)

setwd('/Users/meaghansullivan/Documents/Yale/Simon_Lab/TILAC/S2HeatShock/1908_HeatShock_2/1910_Pipeline2/Revisions')
write_csv(volcano_df, 'TILAC_fulldataset_HeatShock.csv')


# MA Plot -----------------------------------------------------------------

#min of log_reads == 1.761552

up <- volcano_df %>% filter(significance == "Upregulated")
down <- volcano_df %>% filter(significance == "Downregulated")
ns <- volcano_df %>% filter(significance == "Not Significant")

# write_csv(up, paste(save_dir, date, '_HS_up.csv', sep = ''))
# write_csv(down, paste(save_dir, date, '_HS_down.csv', sep = ''))
# write_csv(up, paste(save_dir, date, '_HS_ns.csv', sep = ''))

# Get things to plot on top: 
hs_genes <- c("Hsp23", "Hsp68", "Hsp83", "Hsp67Bc", "Hsp70Bb") # 

volcano_cntl <- volcano_df %>% filter(gene_id %in% hs_genes)
label_ns <- volcano_df %>% filter(gene_id %in% c('Actn', "Set1", "AGO1"))
label_down <- volcano_df %>% filter(gene_id %in% c('Ack', 'Arp5'))

pdf(file = paste(save_dir, date, '_HS_MA_contour.pdf', sep = ''), height = 2, width = 2)
ma_plot <- ggplot(ns, aes(x = log_reads, y = L2FC)) + 
  geom_point(size = 0.5,  alpha = 2/5, color = "gray85") +
  geom_point(data = up, aes(x = log_reads, y = L2FC), size = 0.75,  alpha = 3/5, color = "green4") +
  geom_point(data = down, aes(x = log_reads, y = L2FC), size = 0.75,  alpha = 3/5, color = "goldenrod3") +
  
  scale_color_manual(values=c("goldenrod3", "gray50", "green4"))  + #"burlywood3", "gray55", 
  # scale_color_gradient2(low = "darkgreen",
  #                       mid = "gray",
  #                       high = "turquoise4",
  #                       midpoint = 0) + 
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 6), 
        axis.title = element_text(size = 7), 
        plot.title = element_text(size=7, hjust = 0.5),
        text=element_text(family="Arial")) +
  
  geom_point(data = volcano_cntl, aes(x = log_reads, y = L2FC), color = 'green4',
             size = 0.75, shape = 16) +
  geom_point(data = volcano_cntl, aes(x = log_reads, y = L2FC), color = 'black',
             size = 0.75, shape = 1) +
  # geom_point(data = label_ns, aes(x = log_reads, y = L2FC, color = significance),
  #            size = 0) +
  # geom_point(data = label_down, aes(x = log_reads, y = L2FC, color = significance),
  #            size = 0) +
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
  
  # geom_text(data = label_ns,
  #           aes(label = gene_id),
  #           fontface = 'italic',
  #           color = 'black',
  
  
  #           size = 4.5,
  #           check_overlap = T,
  #           vjust = -0.5) +
  
  ylim(-5, 5) + 
  ggtitle("Heat Shock Response")
ma_plot + geom_density_2d(data = volcano_df, aes(x = log_reads, y = L2FC), color = 'gray60', bins = 5)
dev.off()




# TILAC Plot --------------------------------------------------------------

# TILAC PLOT --------------------------------------------------------------


hs_genes <- c("Hsp26", "Hsp22", "Hsp23", "Hsp68", "Hsp70Ba", "Hsp27", "Hsp83", "stv",  "CG6770", "DnaJ-1", 
              "Hsp67Bc", "Hsp70Bc", "Hsp70Aa", "Hsp70Ab", "Hsp70Ba", "Hsp70Bbb", "Hsp70Bb", "Hsromega", 
              "HSP70Ba", "Hsp67Ba") # 

   # Full dataset ---------------------------------------------------------
      water_full_df <- volcano_df %>% 
        #filter(significance == "Upregulated") %>% 
        arrange(desc(L2FC)) %>% 
        mutate(HS = (gene_id %in% hs_genes)) %>% 
        #filter(L2FC > 2) %>% 
        filter(!gene_id == 'CG12071') #intron of these gene overlaps exon of another, signal in intron



    water_full_df$gene_id <- factor(water_full_df$gene_id, levels = unique(water_full_df$gene_id))  #Takes GF and makes each gene a factor. Saves factors in that order.


# To plot "everything"  -----------------------------------------------------

  ggplot(water_full_df, aes(x = gene_id, y = L2FC, color = HS)) +
    geom_point(size =0.1, shape = 20) +
    #coord_flip() +
    # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) + 
    
    #geom_bar(data = hs_dye_HS, aes(x = hs_dye_HS$GF, y = hs_dye_HS$dye1, fill = "red"), stat = "identity") + 
    scale_color_manual(values=c("#999999", "red3")) +
    #geom_abline(intercept = 0, slope = 0) +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_blank(), #line(color = 'gray90'),
          panel.grid.major.y = element_line(color = 'gray90'), #blank(),
          axis.line.y = element_line(color = 'black'),
          axis.line.x = element_line(color = 'black'),
          axis.text.x = element_blank(), #text(size=9, angle = 90),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size=7)) +
    geom_hline(yintercept = 0) +  
    ylab('TILAC ratio (heat shock vs no treatment)') +
    xlab('Transcript') + 
    ggtitle('Upregulated During Stress') + 
  ggsave(paste0(save_dir, date, 'waterfall_HS_L2FC_full_dataset.pdf'), height = 1.25, width = 3.5)



  
  water_df <- volcano_df %>% 
    filter(significance == "Upregulated") %>% 
    arrange(desc(L2FC)) %>% 
    mutate(HS = (gene_id %in% hs_genes)) %>% 
    filter(L2FC > 2) %>% 
    filter(!gene_id == 'CG12071') #intron of these gene overlaps exon of another, signal in intron
    #CG17834 looks real
    #CG12071 is driven by aberant mutations in last exon of one sample 
    #CG44251 looks real 
    #CG44250 overlaps with 51 but looks real 
    #CG12224 looks real 


water_df$gene_id <- factor(water_df$gene_id, levels = unique(water_df$gene_id))  #Takes GF and makes each gene a factor. Saves factors in that order.


# To plot "everything"  -----------------------------------------------------

ggplot(water_df, aes(x = gene_id, y = L2FC, color = HS)) +
  geom_point(size = 3) +
  #coord_flip() +
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) + 
  
  #geom_bar(data = hs_dye_HS, aes(x = hs_dye_HS$GF, y = hs_dye_HS$dye1, fill = "red"), stat = "identity") + 
  scale_color_manual(values=c("#999999", "red3")) +
  #geom_abline(intercept = 0, slope = 0) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = 'gray90'),
        panel.grid.major.y = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size=9, angle = 90),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size=13)) +
  geom_hline(yintercept = 0) +  
  ylab('Log2(Heat Shock vs No Treatment)') +
  xlab('Transcript') + 
  ylim(-4, 4) + 
  ggtitle('Upregulated During Stress') + 
  ggsave(paste0(save_dir, date, '_HS_L2FC_up4fold.pdf'), height = 2.25, width = 5)







# -------------------------------------------------------------------------
# Duarte Annotations  -----------------------------------------------------
# -------------------------------------------------------------------------


# Load Duarte 2016 HS Data ------------------------------------------------

duarte_genes <- read.csv('/Users/meaghansullivan/Documents/Yale/Simon_Lab/TILAC/S2HeatShock/Literature/Duarte2016_updownreg_table.csv')
duarte_genes <- duarte_genes %>% as.tibble()
#duarte_sig_genes <- duarte_genes %>% filter(DESeq2_padj < 0.001)
hs_genes_activated <- duarte_genes %>% as.tibble() %>% filter(Class =='Activated') %>%
  select(Gene_name) 


# Get things to plot on top: 
#volcano_df_label_down <- volcano_df %>% filter(gene_id %in% c('Actn', "Set1", "AGO1"))
volcano_df_label_hs <- volcano_df %>% filter(gene_id %in% hs_genes_activated$Gene_name)


p <- ggplot(data=volcano_df, aes(x=L2FC, y=-log10(p_adj), col=significance)) + geom_point(size = 0.5) + theme_minimal()+ 
  scale_color_manual(values=c("burlywood3", "black", "olivedrab4")) +
  #xlim(-11, 11) + 
  
  geom_point(data = volcano_df_label_hs, aes(x=L2FC, y=-log10(p_adj)), size = 0.5, color ='darkblue') + 
  geom_text(data = volcano_df_label_hs,
            aes(label = gene_id),
            fontface = 'italic',
            color = 'black',
            size = 3,
            check_overlap = T,
            vjust = -0.5) #+ 
  
  # geom_point(data = volcano_df_label_down, aes(x=L2FC, y=-log10(p_adj)), size = 0.5, color ='brown') + 
  # geom_text(data = volcano_df_label_down,
  #           aes(label = gene_id),
  #           fontface = 'italic',
  #           color = 'black',
  #           size = 3,
  #           check_overlap = T,
  #           vjust = -0.5) #+ 

p2 <- p + geom_vline(xintercept=c(-mu_cutoff, mu_cutoff), col="gray60") +
  geom_hline(yintercept=-log10(FDR_control), col="gray60")

##By default, it is assigned to the categories in an alphabetical order):
p3 <- p2 

#Plot theme
theme_mcs <-  theme(panel.background = element_blank(),
                      panel.grid.major.y = element_line(color = 'gray95'),
                      panel.grid.major.x = element_line(color = "gray95"),
                      axis.line.y = element_line(color = 'black'),
                      axis.line.x = element_line(color = 'black'),
                      axis.text = element_text(color = 'black', size = 5), 
                      axis.title = element_text(size = 7), 
                      plot.title = element_text(size=7, hjust = 0.5)) 

#Generate L2FC_kd Volcano Plot
p3 + theme(legend.text=element_text(size=20)) + theme_classic() +
  xlab("TILAC Log2FC(Heat Shock vs Control)") + 
  ylab(expression(-log[10]("p_adj"))) + 
  ggtitle("TILAC Analysis of Heat Shock Response") +
  theme_mcs + 
  ggsave(paste(save_dir, date, '_TILAC_FlavopiridolTreatment_overlapDuarte_activated.pdf', sep = ''), height = 3, width = 5)


# Saving TILAC Data  ------------------------------------------------------

up <- volcano_df %>% filter(significance == "Upregulated")
down <- volcano_df %>% filter(significance == "Downregulated")
ns <- volcano_df %>% filter(significance == "Not significant")

write_csv(up, paste(save_dir, date, '_HS_up.csv', sep = ''))
write_csv(down, paste(save_dir, date, '_HS_down.csv', sep = ''))
write_csv(up, paste(save_dir, date, '_HS_ns.csv', sep = ''))


# Overlap between TILAC and Activated Duarte ----------------------------------------

volcano_df_water <- volcano_df %>% 
  mutate(Duarte = (volcano_df$gene_id %in% hs_genes$Gene_name)) %>% 
  arrange(desc(L2FC)) %>% 
  filter(L2FC > 2)
volcano_df_water$gene_id <- factor(volcano_df_water$gene_id, levels = unique(volcano_df_water$gene_id)) 


ggplot(volcano_df_water, aes(x = gene_id, y = L2FC, color = Duarte)) +
  geom_point() +
  #coord_flip() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) + 
  scale_color_manual(values=c("#999999", "red3")) +
  #geom_abline(intercept = 1, slope = 0) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(), #line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  ylab('Log2 Fold Change') +
  xlab('Gene') + 
  ggtitle('Duarte Activated Genes in my Dataset') +
  ggsave(paste(save_dir, date, "_DuarteActivated_TILAC_aboveL2FC_2.pdf", sep = ''), height = 3, width = 5)



# Overlap between TILAC and ALL UP Duarte ----------------------------------------

hs_genes_up <- duarte_genes %>% as.tibble() %>% filter(DESeq2_FoldChange > 2, DESeq2_padj < 0.001) %>%
  select(Gene_name) 

volcano_df_water <- volcano_df %>% 
  mutate(Duarte = (volcano_df$gene_id %in% hs_genes_up$Gene_name)) %>% 
  arrange(desc(L2FC)) %>% 
  filter(L2FC > 2)
volcano_df_water$gene_id <- factor(volcano_df_water$gene_id, levels = unique(volcano_df_water$gene_id)) 


ggplot(volcano_df_water, aes(x = gene_id, y = L2FC, color = Duarte)) +
  geom_point() +
  #coord_flip() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) + 
  scale_color_manual(values=c("#999999", "red3")) +
  #geom_abline(intercept = 1, slope = 0) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(), #line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  ylab('Log2 Fold Change') +
  xlab('Gene') + 
  ggtitle('Duarte Activated Genes in my Dataset') +
  ggsave(paste(save_dir, date, "_AllDuarteUp_TILAC_aboveL2FC_2.pdf", sep = ''), height = 3, width = 5)




# Make MA plot:  ----------------------------------------------------------

cB <- readRDS('/Users/meaghansullivan/Documents/Yale/Simon_Lab/TILAC/S2HeatShock/1908_HeatShock_2/1910_Pipeline2/Data/cB.rds')

samp_list <- c('MS1908_01', 'MS1908_02', 'MS1908_03', 'MS1908_04') 

# Get only the desired features:
ranked_features_df  <- cB %>%    
  ungroup() %>%
  filter(sample %in% samp_list) %>%
  ungroup() %>% 
  group_by(XF) %>%
  summarize(nreads = mean(n)/4) %>% 
  dplyr::rename(gene_id = 'XF')

TILAC_MA <- volcano_df %>% left_join(ranked_features_df)

duarte_MA <- TILAC_MA %>% filter(gene_id %in% hs_genes_activated$Gene_name)


ma_plot <- ggplot(TILAC_MA, aes(x = log10(nreads), y = L2FC, color = significance)) + 
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values=c("burlywood3", "gray60", "olivedrab4")) +
  
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 5), 
        axis.title = element_text(size = 7), 
        plot.title = element_text(size=7, hjust = 0.5)) +
  
  xlab('Log10(read count in Heat Shock)') + 
  ylab('DESeq2 Log2FC(Transcription Inhibition / Control)') + 
  geom_abline(intercept = 0, slope = 0, color = "gray10", size = 0.5) + 
  
  geom_point(data = duarte_MA, aes(x = log10(nreads), y = L2FC), color = 'darkblue', alpha = 0.5, size = 1) +
  # geom_text(data = volcano_df_label_down,
  #           aes(label = gene_name),
  #           fontface = 'italic',
  #           color = 'black',
  #           size = 3,
  #           check_overlap = T,
  #           vjust = -0.5) +
  # 
  # geom_point(data = volcano_df_label_hs, aes(x = nreads, y = log2FoldChange), color = 'brown') + 
  # geom_text(data = volcano_df_label_hs,
  #           aes(label = gene_name),
  #           fontface = 'italic',
  #           color = 'black',
  #           size = 3,
  #           check_overlap = T,
  #           vjust = -0.5) + 
  
  
  #ylim(-5, 10) +
  ggtitle("HS TILAC MA Plot") + 
  ggsave(paste(save_dir, date, "_HS_TILAC_MA_DuarteActivated.pdf", sep = ""), height = 3, width = 5)
ma_plot


# Compare sig to sig  -----------------------------------------------------


volcano_df_water <- volcano_df %>% 
  mutate(Duarte = (volcano_df$gene_id %in% hs_genes_activated$Gene_name)) %>% 
  filter(significance == "Upregulated") %>% 
  left_join(ranked_features_df) %>% 
  arrange(desc(L2FC)) 
volcano_df_water$gene_id <- factor(volcano_df_water$gene_id, levels = unique(volcano_df_water$gene_id)) 


ggplot(volcano_df_water, aes(x = gene_id, y = L2FC, color = Duarte)) +
  geom_point() +
  #coord_flip() +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) + 
  scale_color_manual(values=c("#999999", "red3")) +
  #geom_abline(intercept = 1, slope = 0) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(), #line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  ylab('Log2 Fold Change') +
  xlab('Gene') + 
  ggtitle('Duarte Activated Genes in my Dataset') +
  ggsave(paste(save_dir, date, "_SigDuarteUp_SigTILAC_aboveL2FC_2.pdf", sep = ''), height = 3, width = 5)






# Last thing is to compare to mutation rates...  --------------------------

mutrates <- cB %>% group_by(sample, XF) %>% 
  summarize (tc_tot = sum(TC*n), t_tot = sum(nT*n),ga_tot = sum(GA*n), g_tot = sum(nG*n), tot_reads = sum(n)) %>% 
  mutate(TC_mut = tc_tot/t_tot, GA_mut = ga_tot/g_tot) %>%
  group_by(XF) %>% 
  summarize(TC_mut = median(TC_mut), 
            GA_mut = median(GA_mut)) %>% 
  dplyr::rename(gene_id = 'XF') %>% 
  drop_na() #mutate(mut_ratio = replace_na(mut_ratio,0)) %>% 

volcano_df_muts <- left_join(volcano_df_water, mutrates)



# High GA mut rates -------------------------------------------------------

high_ga_mut <- volcano_df_muts %>% filter(GA_mut > 0.005)

ma_plot <- ggplot(volcano_df_muts, aes(x = GA_mut, y = L2FC, color = Duarte)) + 
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values=c("burlywood3", "gray60", "olivedrab4")) +
  
  geom_point(data = high_ga_mut, aes(x = GA_mut, y = L2FC), color = 'darkblue', alpha = 0.5, size = 1) +
  geom_text(data = high_ga_mut,
            aes(label = gene_id),
            fontface = 'italic',
            color = 'black',
            size = 3,
            check_overlap = T,
            vjust = -0.5) +
  
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 5), 
        axis.title = element_text(size = 7), 
        plot.title = element_text(size=7, hjust = 0.5)) +
  
  xlab('GA Mutation Rate') + 
  ylab('TILAC Log2FC(Transcription Inhibition / Control)') + 
  geom_abline(intercept = 0, slope = 0, color = "gray10", size = 0.5) + 
  ggtitle("Correlation of GA Mutation Rate and L2FC") + 
  ggsave(paste(save_dir, date, "_L2FC_GAmutrate.pdf", sep = ''), height = 3, width = 5)
ma_plot


# High TC mut rate --------------------------------------------------------

high_tc_mut <- volcano_df_muts %>% filter(TC_mut > 0.005)


ma_plot <- ggplot(volcano_df_muts, aes(x = TC_mut, y = L2FC, color = Duarte)) + 
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values=c("burlywood3", "gray60", "olivedrab4")) +
  
  # geom_point(data = high_tc_mut, aes(x = TC_mut, y = L2FC), color = 'darkblue', alpha = 0.5, size = 1) +
  # geom_text(data = high_tc_mut,
  #           aes(label = gene_id),
  #           fontface = 'italic',
  #           color = 'black',
  #           size = 3,
  #           check_overlap = T,
  #           vjust = -0.5) +
  
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = 'gray90'),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 5), 
        axis.title = element_text(size = 7), 
        plot.title = element_text(size=7, hjust = 0.5)) +
  
  xlab('TC Mutation Rate') + 
  ylab('TILAC Log2FC(Transcription Inhibition / Control)') + 
  geom_abline(intercept = 0, slope = 0, color = "gray10", size = 0.5) + 
  ggtitle("Correlation of TC Mutation Rate and L2FC") + 
  ggsave(paste(save_dir, date, "_L2FC_TCmutrate.pdf", sep = ''), height = 3, width = 5)
ma_plot

