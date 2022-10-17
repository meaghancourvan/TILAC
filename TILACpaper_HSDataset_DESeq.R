
# Load libraries  ---------------------------------------------------------
  
  require(DESeq2)
  require(tidyverse)


# Set work environment  ---------------------------------------------------

  date <- Sys.Date() %>% format('%y%m%d')
    # used when saving graphs and df's

  # Load data:

  hs_df <- readRDS('~/Box Sync/TILAC/Current Version/Data/HeatShock/cB.rds')
    # read in heat shock data  
  setwd('/Users/meaghansullivan/Box Sync/TILAC/Current Version/Data/HeatShock/')
    # set working directory
  save_dir <- './Output/'
    # used when saving graphs and df's 


# Calculate average reads - for MA plot -----------------------------------

  #Select the TILAC samples that will be analyzed in this script
  samp_list_reads <- c('MS1908_01', 'MS1908_02', 'MS1908_03', 'MS1908_04')
    #MS1908_01 --> HS-s4U | cntl-s6G, rep1
    #MS1908_02 --> HS-s4U | cntl-s6G, rep2
    #MS1908_03 --> HS-s6G | cntl-s4U, rep1
    #MS1908_04 --> HS-s6G | cntl-s4U, rep2
  
  
  log_readcount  <- hs_df %>%    
    ungroup() %>%
    filter(sample %in% samp_list_reads) %>%
    group_by(XF) %>%
    filter(!grepl('28S', XF)) %>% 
    filter(!grepl('rRNA', XF)) %>% 
    summarize(log_readcount = log10(sum(n)/4)) %>% 
    dplyr::rename(gene_name = 'XF')


# Select genes with sufficient sequencing depth  --------------------------

  select_genes <- function(df = hs_df,
                           samples = c('MS1908_05', 'MS1908_06'),
                           totcut = 200){
    y <- df %>%
      ungroup() %>%
      filter(sample %in% samples,
             !(grepl('__', XF))) %>%
      group_by(sample, XF) %>%
      mutate(totcounts = n) %>%
      summarize(totcounts = sum(n)) %>%  # I could choose to average here.
      filter(totcounts >= totcut) %>%
      ungroup() %>%
      dplyr::select(XF) %>%
      unlist() %>% 
      unique()
    return(y)
  }


  genes_to_analyze <- select_genes(hs_df, samples = c('MS1908_05', 'MS1908_06',
                                                      'MS1908_07', 'MS1908_08',
                                                      'MS1908_09', 'MS1908_10'), 
                                   totcut = 200)


# Select and label samples according to treatment for DESeq2 --------------

  samp_list <- c('MS1908_05', 'MS1908_06', 'MS1908_07', 'MS1908_08', 'MS1908_09', 'MS1908_10') 
  type_list <- c(1, 0, 1, 0, 1, 0) # Heat shock or not heat shock. 


# Create df containing just genes and samples chosen above ----------------

  hs_df_selected_genes  <- hs_df %>%    
    ungroup() %>%
    filter(sample %in% samp_list) %>%
    filter(XF %in% genes_to_analyze) %>%
    group_by(sample, XF) %>%
    filter(!grepl('28S', XF)) %>% 
    filter(!grepl('rRNA', XF)) %>% 
    summarize(nreads = sum(n))


# Data now needs to be organized for DESeq2 formatting --------------------

  RNAseq_data <- spread(hs_df_selected_genes, sample, nreads) 
    #make each sample a column 
  RNAseq_data[is.na(RNAseq_data)] <- 0
    #replace na with 0 (not seen)
  RNAseq_data <- as.data.frame(RNAseq_data)
    #needs to be a data frame, not a tibble
  colnames(RNAseq_data) <- c('gene', 'HS1', 'Ctl1', 
                             'HS2', 'Ctl2', 'HS3', 'Ctl3')
    #set column names according to use by DESeq2
  rownames(RNAseq_data) <- RNAseq_data$gene
    #set row names according to use by DESeq2
  RNAseq_data <- RNAseq_data %>% dplyr::select('HS1', 'Ctl1', 'HS2', 'Ctl2', 'HS3', 'Ctl3')
    #drop the 'gene' column, for correct DESeq2 formmmatting 


# Make the object that contains metadata describing the experiment --------
  
  condition <- c('HS', 'Ctl', 'HS', 'Ctl', 'HS', 'Ctl')
  condition <- factor(condition, levels = c("Ctl", "HS"))
  colData <- data_frame(condition = condition)
  colData <- as.data.frame(colData)
  rownames(colData) <- colnames(RNAseq_data)




# Running DESeq2 ----------------------------------------------------------

  dds <- DESeqDataSetFromMatrix(countData = RNAseq_data,
                                    colData = colData,
                                    design = ~ condition)
    #create DESeq2-specific object of read counts 
  
  dds <- DESeq(dds)
    # run the differential expression pipeline 


# Checking/interpretting DESeq2 results -----------------------------------

  resultsNames(dds) 
    # lists the coefficients
  
  res_alpha <- results(dds, name="condition_HS_vs_Ctl", alpha=0.05)
    # extract results
    # alpha does not affect results 
    # alpha set a p-value of 0.05 used in summary function below
    # it will only report as significant results with padj < 0.05
    # adjusted with Benjamini-Hochberg procedure 
  
  summary(res_alpha)
    # print a summary of the results so that user can get an idea of 
    # how regulation is working, at a high level, before graphing
    # results in more detail



# Prepared DESeq2 results for plotting ------------------------------------

  res_alpha$gene_name <- row.names(res_alpha)
    #turn rownames back into a column, labeled gene_name
  
  diffexp <- as_tibble(res_alpha) %>% 
    filter(!grepl('__', gene_name), 
           !grepl('28S', gene_name), 
           !grepl('18S', gene_name), 
           !grepl('5.8S', gene_name), 
           !grepl('5S', gene_name)) 
    # Make object a tibble and filter out rRNA that is contamination

  diffexp_toplot <- diffexp %>% 
    mutate(significance = ifelse(padj < 0.05, ifelse(log2FoldChange > 0, "Upregulated", "Downregulated") ,"Not Significant")) %>% 
    left_join(log_readcount)
    # add information about significance based on padj
    # join with log_reads calculated above. This means that the data points 
    # will be on the same x-axis here as in the comparable TILAC graph
  
  write_csv(diffexp_toplot, paste0(date, '_HS_DESeq_padj05.csv', sep = ''))
    # save the processed data under today's date 
    # (to avoid any possibility of over-writing old files)

  

# Create annotations for plotting -----------------------------------------

  # First, separate out up, down, ns to plot in different colors 
    up <- diffexp_toplot %>% filter(significance == "Upregulated")
    down <- diffexp_toplot %>% filter(significance == "Downregulated")
    ns <- diffexp_toplot %>% filter(significance == "Not Significant")

  # Select out specific genes I'd like to label 
    hs_labels <- diffexp_toplot %>% filter(gene_name %in% c("Hsp67Bc", "Hsp70Bb", "Hsp23", "Hsp68", "Hsp83"))  
    ns_labels <- diffexp_toplot %>% filter(gene_name %in% c('Actn', "Set1", "AGO1"))
    dn_labels <- diffexp_toplot %>% filter(gene_name %in% c('Ack', 'Arp5'))


  
# Create MA plot ----------------------------------------------------------

  pdf(file = paste('./Output/', date, '_HS_DESeqMA_contour.pdf', sep = ''), height = 2, width = 2)
    # open a pdf object in which to save the graph 
  
  ma_plot <- ggplot(ns, aes(x = log_readcount, y = log2FoldChange)) + 
    geom_point(size = 0.5,  alpha = 2/5, color = "gray85") +
    # first plot no sign change genes in gray 
    geom_point(data = up, aes(x = log_readcount, y = log2FoldChange), size = 0.75,  alpha = 2/5, color = "green4", shape = 16) +
    # second plot positive change genes in green 
    geom_point(data = down, aes(x = log_readcount, y = log2FoldChange), size = 0.75,  alpha = 2/5, color = "goldenrod3", shape = 16) +
    # third, plot the negative changes in gold 
    
    theme(panel.background = element_blank(),
          panel.grid.major.y = element_line(color = 'gray90'),
          panel.grid.major.x = element_blank(),
          axis.line.y = element_line(color = 'black'),
          axis.line.x = element_line(color = 'black'),
          axis.text = element_text(color = 'black', size = 6), 
          axis.title = element_text(size = 7), 
          plot.title = element_text(size=7, hjust = 0.5), 
          text=element_text(family="Arial")) +
    
    geom_point(data = hs_labels, aes(x = log_readcount, y = log2FoldChange), color = 'green4',
               size = 0.75, shape = 16, alpha = 3/5) +
    geom_point(data = hs_labels, aes(x = log_readcount, y = log2FoldChange), color = 'black',
               size = 0.75, shape = 1, alpha = 3/5) +
    
    xlab('log10(read count)') + 
    ylab('log2(heat shock vs untreated)') + 
    geom_abline(intercept = 0, slope = 0, color = "gray10", width = 30) + 
    
    geom_text(data = hs_labels,
              aes(label = gene_name),
              fontface = 'italic',
              color = 'black',
              size = 2,
              check_overlap = T,
              vjust = -0.5) +
    
    ylim(-12, 12) + 
    ggtitle("Heat Shock Response") 
  
  ma_plot + geom_density_2d(data = diffexp_toplot, aes(x = log_readcount, y = log2FoldChange), color = 'gray50', bins = 4)
    # this adds contour lines indicating where the mass of data lies 
  
  dev.off()
    # Close the pdf document.

