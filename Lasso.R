
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
