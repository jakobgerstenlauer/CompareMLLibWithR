###############################################################################################
# Cloud Computing - MIRI Master
# Lecturers: 
# LEANDRO NAVARRO MOLDES, JORDI TORRES VIÑALS
# Term project
# February-March 2017  
# "Compare the Accuracy and Computational Efficiency of Centralized and
# Decentralized Algorithms for the Support Vector Machine 
# and the Lasso in Regression Analysis."
#
# First, an instance generater is set-up.
#
# Second, we define a latin-hypercube (LH) scheme for efficient sampling of parameter space.
#
# Third, for each parameter combination of the LH sampling scheme, one instance (data set) is generated.
#
# Fourth, for each instance we tune the parameters of the SVM and the kernel.
#   Based on the optimum parameters (both kernel and SVM parameters, criteria: cv error) 
#   we record generalization error, sparsity, and total computation time (including tuning) for the selected model.
#
# Sixth, for each instance we tune the parameter Lambda of the Lasso. 
#   Based on the optimum parameter of Lambda, we record generalization error, sparsity, and total computation time (including tuning) for the selected model.
# Seventh, we analyse the results of step 5 and 6 comparing SVM and Lasso across the following dimensions:
#   dim 1: signal-to-noise ratio (between 0: no signal only noise and 1: only signal no noise),
#   dim 2: number of observations N,
#   dim 3: number of parameters D,
#   dim 4: polynomial degree of the inputs.
#
# Date: 02.03.2017
# Jakob Gerstenlauer
# jakob.gerstenlauer@gjakob.de
###############################################################################################

#remove old objects for safety resons
rm(list=ls(all=TRUE))

#utility function
glue<-function(...){paste(...,sep="")}

#define path of standard directories
source("workingDir.R")
source("properties.R")

#################################################################################
#
# Create the synthetical data sets 
#
#################################################################################

#' @title Test validity of numeric values.
#'
#' @description
#' Checks if a numeric vector with given name exists, is of type numeric, and has the given minimum length.
#' This function should be used to test pre- and postconditions within functions.
#' 
#' @param x The numeric vector to be tested.
#' @param min_length The minimum number of elements of the vector.
check.numeric.values<-function(x,min_length){
  stopifnot(exists("x"))
  stopifnot(is.numeric(x))
  stopifnot(length(x)>=min_length)
}

#' @title Creates a new problem instance as a data frame.
#'
#' @description
#' This functions creates a new data set (instance) given four different parameters: 
#' signal-to-noise ratio, number of observations, number of features,  and polynomial degree.
#'
#' @param signal-to-noise_ratio Ratio between signal and noise in the data (between 0: no signal only noise and 1: only signal no noise)
#' @param N Number of observations (cases) in the new data set
#' @param D Number of features (inputs) in the new data set
#' @param polynomialDegree Highest power of the inputs affecting the output
#' @param isDebug Should debug statements be printed? Default is FALSE.
#' 
#' @examples d1<-instance.generator(signal_to_noise_ratio=0.5, N=100, D=10, polynomialDegree=4, isDebug=FALSE)
instance.generator<-function(signal_to_noise_ratio, N, D, polynomialDegree, isDebug=FALSE){
  #test preconditions
  check.numeric.values(signal_to_noise_ratio,1);
  check.numeric.values(N,1);
  check.numeric.values(D,1);
  check.numeric.values(polynomialDegree,1);
  stopifnot(signal_to_noise_ratio>=0)
  stopifnot(signal_to_noise_ratio<=1)
  stopifnot(D>=1)
  stopifnot(N>=1)  
  #Here the variability in the data is based on a hierarchical Gaussian model:
  #The mean and the variance of the inputs and the coefficients follow a normal distribution.
  
  #means of the inputs
  means<-rnorm(n=D, mean=0, sd=1)
  if(isDebug) print(means)
  
  #standard deviations of the inputs
  sds<-abs(rnorm(n=D, mean=1, sd=1))
  if(isDebug) print(sds)
  d<-data.frame("output"=rep(0,N))
  
  #the coefficients of the polynomials of the inputs on the output are normally distributed
  #with mean 0 and standard deviation 1:
  for(nrVar in 1:D){
    for(power in 1:polynomialDegree){
      eval(parse(text=
                   glue("coefs_",nrVar,"_",power,"<-rnorm(n=1, 
                        mean=0,
                        sd=1)")
      ));
    }
    }
  for(nrVar in 1:D){
    eval(parse(text=
                 glue("d$input_",nrVar,"<-rnorm(n=N, 
                      mean=means[",nrVar,"],
                      sd=abs(sds[",nrVar,"]))")
    ));
  }
  
  #additive linear component for each variable
  #and for each power function of the input (from 1 to polynomial degree)
  for(nrVar in 1:D){
    for(power in 1:polynomialDegree){
      eval(parse(text=glue("d$output <- d$output + coefs_",nrVar,"_",power," * d$input_",nrVar,"**",power
      )));
    }
  }
  #add Gaussian noise 
  d$output <- rnorm(n=N, mean=d$output, sd= 1-signal_to_noise_ratio)
  return(d)
  }


#Step 1: define the LH scheme 
require("lhs")


#set-up the Latin Hypercube sampling scheme
LHS<-improvedLHS(n=SampleSize, k=NumVariables, dup=1)

setwd(dataDir)
file.names<-scala.inputfile.names<-""

#V1: signal-to-noise ratio
signal.to.noise.ratio.grid<-rep(0,SampleSize*maxReplicatesLHC)

#V2: number of observations N
num.observations.grid<-rep(0,SampleSize*maxReplicatesLHC)

#V3: number of variables D
num.vars.grid<-rep(0,SampleSize*maxReplicatesLHC)

#V4: polynomial degree of the inputs.
polynomial.degree.grid<-rep(0,SampleSize*maxReplicatesLHC)

i<-1
for (simulation in seq(1,dim(LHS)[1]))
{
  for (arguments in seq(1,NumVariables))
  {   
    #Here we use the quantile function for the uniform distribution to "translate" from the standard uniform distribution to the respective trait range
    eval(parse(text=paste(
      'A',arguments,'<-round(qunif(LHS[simulation,',arguments,'], min=low_V',arguments,', max=high_V',arguments,'),digits=3)'
      ,sep="")));
  } 
  
  for(replicate in 1:maxReplicatesLHC){
    
    #calculate the number of variables based on the ratio observations / variables which is stored in A3
    D<-round(A2/A3)
    
    #Create the new data set with the specific variables
    d<-instance.generator(signal_to_noise_ratio=A1, N=round(A2), D, polynomialDegree=round(A4));
    
    signal.to.noise.ratio.grid[i]<-A1;
    num.observations.grid[i]<-round(A2);
    num.vars.grid[i]<-D;
    polynomial.degree.grid[i]<-round(A4);
    
    i<-i+1
    
    
    dump.file.name<-glue("data_signal_to_noise_", A1,
                         "_N_", round(A2),
                         "_D_", D,
                         "_poly_", round(A4),
                         "_replicate_",replicate,
                         ".RData");
    save(list="d", file=dump.file.name);
    file.names<-c(file.names, dump.file.name);
    dump.file.name<-glue("data_signal_to_noise_", A1,
                         "_N_", round(A2),
                         "_D_", D,
                         "_poly_", round(A4),
                         "_replicate_",replicate,
                         ".csv");
    
    scala.inputfile.names<-c(scala.inputfile.names, dump.file.name);
    
    #Scala specific output:
    #first separator is ";" all following separator are blanks
    seps <- c(";", rep(" ", (dim(d)[2]) - 1), "\n")
    
    apply(d, 1, function(x, s=seps) {
      cat(x, sep=s, file = dump.file.name, append = TRUE)
    })

    #normal csv output
    #write.table(d, file=dump.file.name, row.names = FALSE, sep=";")
  }
}     

file.names<-file.names[-1]
scala.inputfile.names<-scala.inputfile.names[-1]
cat(scala.inputfile.names, file = "fileNames.txt", fill = TRUE)

#################################################################################
#
# Lasso regression 
#
#################################################################################
require("glmnet")
set.seed(1365)
sample.size <- length(file.names) 

#Id of the parameter combination (has maxReplicatesLHC replicates!)
id_parameter_combination<-vector(mode="numeric",length=sample.size) 

#vectors for optimal parameters
lambda_setting<-vector(mode="numeric",length=sample.size)  

#vectors for sparsity 
sparsity<-vector(mode="numeric",length=sample.size)

#vector for R2
r2<-vector(mode="numeric",length=sample.size) 

#vector for mean squared error
mse<-vector(mode="numeric",length=sample.size) 

#vector for computation time
compu.time<-vector(mode="numeric",length=sample.size)

#index for files
i<-1

#index for parameter combination
j<-1

for(fileName in file.names){
  print(paste("read input file:", fileName))
  load(fileName)
  
  #For each parameter combination there are maxReplicatesLHC replicates.
  #Index j identifies the replicate.
  id_parameter_combination[i]<-j
  
  if((i %% maxReplicatesLHC)==0){
    j<-j+1
  }

  #the output / response variable
  y<-d$output
  #a matrix of inputs / predictor variables
  x<-as.matrix(subset(d,select=-output))
  #a range of possible values for lambda
  grid<-10^ seq (10,-2 ,length =100)
  
  #start recording computation time
  ptm <- proc.time();
  
  #fit the Lasso (alpha=1) for the best lambda as found by cross-validation
  #assuming a normal distribution of the outputs (family="gaussian")
  cv.out<-cv.glmnet(
    x,
    y,
    family="gaussian",
    alpha = 1,
    lambda = grid)

  lambda_setting[i] <- cv.out$lambda.min
  predictions    <- predict(cv.out, newx = x, s = cv.out$lambda.min)
  
  #coefficient of determination
  r2[i]<-cor(predictions,y)
  
  #mean squared error
  mse[i]<- mean((predictions-y)**2)
  
  #store the coefficients as an object (it is a sparse matrix)
  c.lasso <- coef(cv.out)
  #####################################
  #Internals of the object "c.lasso":
  #####################################
  #
  # It is an object of S4 class dgCMatrix,
  # that means it is a sparse matrix where only non-zero entries are stored internally.
  #
  # It contains 6 slots or internal data fields:
  #@i This slot stores the indices of the non-zero data entries.
  #@p This slot stores the dimensions of the hidden data field which contains all non-zero data entries.
  #@Dim This slot contains the "official" dimensions of the matrix, as seen from outside.
  #Here, only the columns with indices given in slot i are occupied.
  #@Dimnames A list with character vectors for the names of columns and rows of the matrix.
  #@x This slot contains the actual (non-zero) coefficients.
  #The matrix has dimension
  # See output of the command "str(c.lasso)"
  # Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  #..@ i       : int [1:27] 0 86 129 178 670 790 1024 1500 1817 2124 ...
  #..@ p       : int [1:2] 0 27
  #..@ Dim     : int [1:2] 10000 1
  #..@ Dimnames:List of 2
  #.. ..$ : chr [1:10000] "(Intercept)" "V1" "V2" "V3" ...
  #.. ..$ : chr "s0"
  #..@ x       : num [1:27] 45.92511 -3.21903 0.00372 0.19754 -0.99976 ...
  #..@ factors : list()
 
  #Here I determine the names of all variables whose coefficient is not zero.
  #c.lasso@Dimnames[[1]][c.lasso@i + 1] : I have to add 1 to the values in slot i,
  #because the indexing starts with 0 (C and Java style)instead of 1 (R style).
  
  names.positive.coefficients.lasso<-c.lasso@Dimnames[[1]][c.lasso@i + 1]
  sparsity[i]<-length(names.positive.coefficients.lasso)
  time.used <- proc.time() - ptm
  print(time.used[3])
  compu.time[i]<-time.used[3]
  i<-i+1
}

#Set up a data set with all the results of the simulation
d.results<-data.frame(
  compu.time,
  signal.to.noise.ratio=signal.to.noise.ratio.grid,
  num.vars=num.vars.grid,
  num.observations=num.observations.grid,
  polynomial.degree=polynomial.degree.grid,
  id_parameter_combination,
  #optimal parameters
  lambda_setting,
  r2,
  mse,
  sparsity)

setwd(dataDir)
#get current date and replace hyphens by underline
Date<-gsub(pattern="-", replacement="_",Sys.Date())
#paste new filename
fileName<-paste("Results_Simulation_Lasso_",Date,".csv",sep="")
write.table(d.results, file=fileName, append=FALSE, row.names = FALSE, sep = ";")
