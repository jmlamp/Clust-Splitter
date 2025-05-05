# Clust-Splitter $-$ Nonsmooth Optimization Based Incremental Clustering Software Using LMBM
Clust-Splitter is a nonsmooth optimization (NSO) based clustering algorithm for solving the minimum sum-of-squares clustering (MSSC) problem in large data sets. The method comprises two main components: an incremental algorithm and the Limited memory bundle method (LMBM), which is applied at each iteration to solve both the main and auxiliary clustering problems. Due to its incremental nature, Clust-Splitter not only solves the $k$-partition problem but also all intermediate $l$-partition problems for $l=1,\ldots,k−1$.
## Files included
* clustsplitter.f95
  - Program file (Fortran95 code).
* Makefile
  - For compiling the Clust-Splitter method. Requires a Fortran compiler (gfortran) to be installed.
## Program usage
To use the code:
1)	Select the data, name of output file, numbers of records and features, the maximum number of clusters "max_cl", and other parameters (for example n_outlier, ncenter1, and ncenter2) at the end of clustsplitter.f95 file.
2)	Run Makefile by typing "make".
3)	Finally, just type "ajo.exe".
   
The algorithm outputs a txt-file containing clustering function values, the Davies-Bouldin and Dunn validity indices, the number of function and subgradient evaluations, and the elapsed CPU times for up to max_cl clusters. The file also includes the data reading time.
## References
* Clust-Splitter:
  - Lampainen, J., Joki, K., Karmitsa, N., & Mäkelä, M. M. (2025). Clust-Splitter − an efficient nonsmooth optimization-based algorithm for clustering large datasets. arXiv:2505.????1 [math.OC].
* LMBM:
  - Haarala, N., Miettinen, K., & Mäkelä, M. M. (2007). Globally convergent limited memory bundle method for large-scale nonsmooth optimization. Mathematical Programming, 109, 181-205.
  - Haarala, M., Miettinen, K., & Mäkelä, M. M. (2004). New limited memory bundle method for large-scale nonsmooth optimization. Optimization Methods and Software, 19(6), 673-692.
## Acknowledgements
The work was financially supported by the Research Council of Finland (projects no. #345804 and #345805 led by Prof. Tapio Pahikkala and Prof. Antti Airola, respectively), and Jenny and Antti Wihuri Foundation.
