# [TILAC](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac693/6677324?login=false)

To make accurate estiamtes of differences in RNA levels between experimental and control samples, measurement methods need to distinguish between biological variation and technical noise. We developed a method called TILAC, which achieves this using the metabolic labels s4U and s6G. In a TILAC experiment, RNA from one sample is labeled with s4U, and the other is labeled with s6G. Cells from each sample are mixed prior to any downstream biochemical steps, like fractionation or polysome purification. During sequencing library preparation, RNA is treated with TimeLapse chemistry to recode the hydrogen bonding pattern of s4U -> C and s6G -> A. After sequencing, the amount of RNA that was present in either sample can be inferred, using a statistical model, from the presence of T -> C or G -> A mutations in reads. 

- Courvan, M.C.S., Vock, I.W., Niederer, R.O., Kiefer, L., Gilbert, W.V., Simon, M.D. 2022. Internally normalized RNA sequencing
comparisons using nucleoside recoding chemistry. Nucleic Acids Res, doi: 10.1093/nar/gkac693.


![alt text](https://github.com/meaghancourvan/TILAC/blob/main/TILAC_Fig1.png?raw=true)

The statistical model used to analyze TILAC data is an extension of the model published by Schofield et al. 2018. 

Raw sequencing data was processed using a [pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by former members of the Simon Lab. 

- **hs_dataset.rds** - contains data from a heat shock experiment peformed in *Drosophila* S2 cells
- **TILACpaper_HSDataset_DESeq.R** - analyses the data using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- **TILAC_data_to_model.R**  - is a wrapper script that takes the heat shock data and feeds it into a Bayesian statistical model written in [Stan](https://mc-stan.org/) (see below) 
- **tilac_pois.stan** - is Poissan mixture model that calculates the differences in RNA levels between two conditions given mutation rates of the two independent RNA metabolic labels. It is written in the probabilistic programming language [Stan](https://mc-stan.org/)
- **TILACpaper_HSDataset_TILAC_graph.R** - is a script that takes the results from the statistical model and calculates a p-value, used to determine what genes are significantly different between conditions. Genes and their TILAC values are graphed on an MA plot and colored by significance.

This work was done in collaboration with **Isaac Vock** and Matthew Simon. It follows from two previous publications
- Schofield,J.A., Duffy,E.E., Kiefer,L., Sullivan,M.C. and Simon,M.D. (2018) TimeLapse-seq: adding a temporal dimension to RNA sequencing through nucleoside recoding. Nat. Methods, 15, 221–225. [[link]](https://www.nature.com/articles/nmeth.4582)
- Kiefer,L., Schofield,J.A. and Simon,M.D. (2018) Expanding the nucleoside recoding toolkit: revealing RNA population dynamics with 6-thioguanosine. J. Am. Chem. Soc., 140, 14567–14570. [[link]](https://pubs.acs.org/doi/abs/10.1021/jacs.8b08554)
