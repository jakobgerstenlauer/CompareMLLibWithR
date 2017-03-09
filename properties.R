#*******************************************************************************
# Here all parameters are set for the R script:
#*******************************************************************************

#Define number of replications for 10-fold Cross-validation!
#E.g 10 means that we use 10-times 10-fold cross-validation
numCVReplicates<-1

#Define number of replications for each  parameter combination
#sampled in the Latin Hyper Cube.
maxReplicatesLHC<-5

#TODO Check if the sample size is correct!
#number of samples from the LHC 
SampleSize<-16;#TODO 200
NumVariables<-4;   

#Now define the ranges for all four parameters of the LHC:
#V1: signal-to-noise ratio
low_V1= 0.01;
high_V1= 0.99;

#V2: number of observations N
low_V2  = 16;
high_V2 = 256;

#V3: ratio observations / number of variables D
low_V3  = 0.5;
high_V3 = 20.0;

#V4: polynomial degree of the inputs.
low_V4  = 1;
high_V4 = 2;
#*******************************************************************************
