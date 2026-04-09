## Embedding General Sparse Hypergraphs with Non-uniformity and Multiplicity

Code in this repository reproduces the results in  
[A general latent embedding approach for modeling non-uniform high-dimensional sparse hypergraphs with multiplicity](https://arxiv.org/abs/2410.12108) using R.

**Folders and files**:

- **`cocitation_network/`**  
  Contains the code and data for the co-citation hypergraph analysis. See the folder for details.  

- **`HGraphEmbed.zip`**  
  An R package implementing hypergraph embedding with Rcpp.  

- **`Theta_result.R`**  
  Code to reproduce the simulation results in Figures 2 and 3 of the paper.

- **`evaluate_Theta.R`**  
  Code to aggregate and evaluate the simulation results for Figures 2 and 3 of the paper.

- **`coverage_result.R`**  
  Code to reproduce the simulation inference results reported in Figures 4 and 5 of the paper.

- **`evaluate_coverage.R`**  
  Code to plot the simulation inference results reported in Figures 4 and 5 of the paper. 

- **`functions_limit.R`**  
  Utility functions used in the analysis.

**Reproducibility guide**:

- To reproduce Figures 2 and 3 in the paper, first run `Theta_result.R` with seeds 1, ..., 100 and save the results. Then aggregate the results using `evaluate_Theta.R`.

- To reproduce Figures 4 and 5 in the paper, first run `coverage_result.R` with seeds 1, ..., 1000 and save the results. Then aggregate the results using `evaluate_coverage.R`.

- To reproduce Figures 6-8 in the co-citation hypergraph analysis in the paper, run `./cocitation_network/learn_the_embeddings.R` to generate the embeddings and use `./cocitation_network/learn_the_crs.R` to produce the plots.

- To reproduce Figures 1-6 in the Supplmentary material, run `supp.R`.

**Author embedding plot from the co-citation hypergraph**:

<p align="center">
<img src="author_embed.png" alt="drawing" width="600"/>
</p>


