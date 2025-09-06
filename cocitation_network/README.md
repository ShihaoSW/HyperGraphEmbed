## Data Resource

The MADStat dataset is available at [ZhengTracyKe/MADStat](https://github.com/ZhengTracyKe/MADStat).  
If you use this data, please reference the repository and the paper:  
**"Co-citation and Co-authorship Networks of Statisticians"**  
([Journal of Business & Economic Statistics](https://www.tandfonline.com/doi/full/10.1080/07350015.2021.1978469)).  

In this repository, we make use of only a subset of the data from MADStat.  
For full details of the data, please see the [MADStat repository](https://github.com/ZhengTracyKe/MADStat).  

---

## Repository Files

- **`AuPapMat.csv`**  
  Author–paper incidence matrix summarizing the BibTeX data.  

- **`AuthorPaperInfo.RData`**  
  R data object containing metadata on authors and papers.  

- **`BibtexInfo.RData`**  
  R data object with bibliographic information.  

- **`PapPapMat.csv`**  
  Paper–paper incidence matrix summarizing the citation data.  

- **`author_name.txt`**  
  Text file listing author names.  

- **`functions_limit.R`**  
  R script with utility functions used in the analyses.  

- **`hg1991-2000_3000_final_embeddings.RData`**  
  R data object containing the final hypergraph embeddings of 3,000 authors trained on papers published during 1991–2000.  

- **`learn_the_crs.R`**  
  R script for constructing confidence intervals and generating plots from the embeddings.  

- **`learn_the_embeddings.R`**  
  R script for learning embeddings from the hypergraph data.  
