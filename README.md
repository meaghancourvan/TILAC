# TILAC

This repository contains data and scripts from my graduate school work developing [TILAC](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac693/6677324?login=false). TILAC is a method for internally controlling RNA-sequencing data, thereby allowing for more accurate comparisons between two conditions. 

- **hs_dataset.rds** - contains data from a heat shock experiment peformed in *Drosophila* S2 cells
- **TILACpaper_HSDataset_DESeq.R** - analyses the data using DESeq2
- **TILAC_data_to_model.R**  - is a wrapper script that takes the heat shock data and feeds it into a Bayesian statistical model written in Stan 
- **tilac_pois.stan** - is a mixed Poisson model that calculates the differences in RNA levels between two conditions given mutation rates of the two independent RNA metabolic labels
- **TILACpaper_HSDataset_TILAC_graph.R** - is a script that takes the results from the statistical model, calculates a p-value, and graphs the data on an MA plot

This work was done in collaboration with **Isaac Vock** and **Matthew Simon**. It follows from two previous publications

- Schofield,J.A., Duffy,E.E., Kiefer,L., Sullivan,M.C. and Simon,M.D. (2018) TimeLapse-seq: adding a temporal dimension to RNA sequencing through nucleoside recoding. Nat. Methods, 15, 221–225. [[link]](https://www.nature.com/articles/nmeth.4582)
- Kiefer,L., Schofield,J.A. and Simon,M.D. (2018) Expanding the nucleoside recoding toolkit: revealing RNA population dynamics with 6-thioguanosine. J. Am. Chem. Soc., 140, 14567–14570. [[link]](https://pubs.acs.org/doi/abs/10.1021/jacs.8b08554)
