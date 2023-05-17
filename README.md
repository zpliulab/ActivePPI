#  ActivePPI: Quantifying Protein-Protein Interaction Network Activity with Markov Random Fields 

**We have developed a novel framework, ActivePPI, which effectively measures the activity of PPIs network, the activity of pathways, and the optimal PPIs network architecture from mass spectrometry abundance data under different cellular conditions.**

**Please cite our paper if it helps you (under reviewing).**

![workfolw](https://github.com/zpliulab/ActivePPI/figure/Fig1.png)

## Environment requirements
- Download anaconda/miniconda if needed
- `R==4.2.2`
- `igraph`
- `ks`
- `parallel`
- `doParallel`
- `doSNOW`
- `ggplot2`
- `KEGGREST`

## Files:
- **main_BRCA.R: Perform pathway activity analysis on BRCA dataset (PDC000120)**
- **main_SARS-CoV-2.R: Perform pathway activity analysis on SARS-CoV-2 dataset**
- **main_simu.R: Perform PPIs network activity analysis on a toy dataset**
- **code**: subfunction used in the ActivePPI process*
- **BRCA_TCGA**: Related documents of BRCA
- **SARS-CoV-2**: Related documents of BRCA
- **simu**: Related documents of simulation examples

## Run the code
Run the following code in the R environment

- main_BRCA.R
- main_SARS-CoV-2.R
- main_simu.R

## Use ActivePPI on a new dataset

####Function 1 (estimating pathway activity):
```result <- activePPI_pathway(Protein.Pathway, protein, Energe.single, corr = prob.pro, parallel = FALSE, verbose = FALSE, rtxt = FALSE, rtxtfile='SARS_result.csv')```

	Input:
	  Protein.Pathway : a list containing gene set networks. It must contain "from" and "to" columns, eg.
	       		form     to
	         protein_A protein_B
	         protein_A protein_C
	         protein_B protein_D
	  protein:   a (n*m) protein mass spectrometry information matrix
	  Energe.single:   Kernel density of each protein (n*m)
	  corr :   distance matrix between proteins (n*n)
	  parallel (default = FALSE):   whether to run in parallel (unstable,it is recommended to run one by one)
	  verbose (default = FALSE):   whether to load the progress bar
	  rtxt (default=FALSE):    whether to output documentation
	  rtxtfile: output file name
	Output:
	 result: a list, including pvalue (result$pvalue),
	        joint probability density of the pathway (result$dis),
	        energy of each clique in each network (result$re),
	        optimal network architecture (result$gra)	
####Function 2 (estimating network activity):
```result <- activePPI_network(gra, protein, Energe. single, prob. pro)```

	Input:
	  gra: a network of igraph structures
	  The rest of the settings are consistent with the function one
	
	Output:
	  result: a list, including pvalue (result$pvalue),
	         joint probability density of each network (result$prob),
	         energy of each clique in each network (result$prob.list),
	         optimal network architecture (result$ gra)

Please write to [zpliu@sdu.edu.cn](mailto:zpliu@sdu.edu.cn) if you have any questions.
