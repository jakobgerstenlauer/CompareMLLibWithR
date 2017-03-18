# CompareMLLibWithR
In this project, we compared two implementations of the Lasso (L1) regression algorithm:
1) the R function cv.glmnet() in the glmnet package, 
2) the function LassoWithSGD() in MLlib. 

Our objective is to compare the two algorithms in terms of:
1) usability,
2) performance,
3) accuracy.

The following files are part of this project:

CC.R              The main R script that both creates the input files and analyses each file with the cv.glmnet() function.
properties.R      A properties file that controls the latin hypercube sampling.
workingDir.R      This file is not submitted because it contains the hardcoded path variables of the local user. 
ScalaShell.scala  A scala script file that can be run in the Spark shell. 
CompareResultsLasso_CompareSparkvsR.csv 
                  A csv file with the combined results of the R and Spark simulation runs. 
analyzeResults.R  A script to analyze the above scv file.
