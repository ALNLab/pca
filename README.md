# Pca

# The point of this package is to be able to do permutation/bootstrap analysis of pca results to remove inconsistent components. See componentSig.m to start and check out the paper this is based on. The trick to this approach is to figure out which components bootstrapped/permutted components match up with. For this we use a couple of different approaches, but the default is a specific procrustes rotation that differs slightly from matlab's version (again, see McIntosh paper and read me in componentSig).

# To make this script big data friendly it runs its permutations/bootstraps in batches and saves them out in batches. Significance is assessed for you. It can also project data into your pca space, and do some other cool stuff. You may encounter some bugs if you use the script outside the way I have in the past. Let me know if you run into bugs. Some functionality needs to be finished: for example, parallelization. 
