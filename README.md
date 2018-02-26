# Pca

The point of this package is to be able to do permutation/bootstrap analysis of pca results to remove inconsistent components and loadings respectively. See componentSig.m to use these scripts and check out the paper this is based on (contained in the help section of componentSig). The trick to this approach is to figure out which components bootstrapped/permutted components match up with. For this we use a couple of different approaches, but the default is a specific procrustes rotation that differs slightly from matlab's version (again, see McIntosh paper and read me in componentSig). PCA/SVD itself is performed with scripts native to matlab.

To make this script friendly for larger datasets, it runs its permutations/bootstraps in batches and saves them out in batches. This helps in the case that something crashes or if you want to resume the analysis later. To make sure you can leave off right where you started, parameters are saved and can be loaded in to run the analysis from the point you left off.

The script can also project data into your pca space, and do some other cool stuff. You may encounter some bugs if you use the script outside the way I have in the past (not on purpose!). Let me know if you run into bugs. Some functionality needs to be finished: for example, parallelization. This is noted in the documentation where applicable. 
