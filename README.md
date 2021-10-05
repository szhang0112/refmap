# refmap (REgional Fine-MAPping)
RafMap is a Bayesian network to prioritize disease-associated genes by integrating GWAS summary statistics with funtional annotations (e.g. ATAC-seq, ChIP-seq, etc.). The core code of RefMap was written in MATLAB. This package is in active development.
# Inputs
1. GWAS summary statistics: CHR, ID, POS, BETA, SE, and P.
2. LD matrix. Use in-sample estimation if you have access to the individual-level data. Otherwise, you can perform out-sample estimation based on external cohort with matched population, e.g., the 1000 Genomes.
3. Genomic annotations: epigenetic profiling (e.g. ATAC-seq and ChIP-seq) or evolutionary features (e.g. positive selection scores).
4. RNA-seq and HiC (both optional).
# Reference
S. Zhang, J. Cooper-Knock, A.K. Weimer, M. Shi, T. Moll, C. Harvey, H.G. Nezhad, J. Franklin, C.D.S. Souza, C. Wang, J. Li, C. Eitan, E. Hornstein, K.P. Kenna, Project MinE Sequencing Consortium, J. Veldink, L. Ferraiuolo, P.J. Shaw, and M.P. Snyder. Genome-wide identification of the genetic basis of amyotrophic lateral sclerosis. bioRxiv, doi: 10.1101/2020.11.14.382606, 2020.
