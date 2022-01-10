###############################################
## A function to change letters to 0 or 1  ####
###############################################
G_to_binary <- function(G) {
  M <- matrix(0, nrow(G), ncol(G))
  for (j in 1:ncol(G)) {
    most.common <- names(which.max(table(G[,j])))
    M[G[,j]==most.common,j] = 1
  }
  return(M)
}

###############################################
## A function to compute relatedness       ####
###############################################
Relatedness <- function(Z, PLOT=FALSE) {
  m <- ncol(Z)
  G <- tcrossprod(scale(Z))/m
  if (PLOT) {
    par(mfrow=c(2,2))
    plot(diag(G), main="Inbreeding if values >1", ylab="Inbreeding coefficient")
    hist(as.numeric(G-diag(diag(G))), main="Relatedness Between Individuals", xlab="")
    image(G, main="Matrix of relationships")
  }
  cat("Average relatedness between all individuals ", round(mean(G-diag(diag(G))),3), "\n")
  cat("Maximum relatedness ", round(max(G-diag(diag(G))),3), "\n")
  return(G)
}

#################################
## GET THE PHENOTYPE - Mean  ####
#################################
get.phenotype <- function(pt, trait="sawGras") {
  trait.index <- which(names(pt) %in% trait)
  N <- nrow(pt)
  Phenotype <- numeric(N)
  for(j in 1:N)
  {
    m <- nchar(as.character(pt[[trait.index]][j]))
    vec <- as.numeric(unlist(strsplit(substr(pt[[trait.index]][j], 2, m-1), ',')))
    Phenotype[j] <- mean(vec) #variance
  }
  return(Phenotype)
}

#####################################
## GET THE PHENOTYPE - Variance  ####
#####################################
get.phenotype_new <- function(pt, trait="sawGras") {
  trait.index <- which(names(pt) %in% trait)
  N <- nrow(pt)
  Phenotype_new <- numeric(N)
  for(j in 1:N)
  {
    m <- nchar(as.character(pt[[trait.index]][j]))
    vec <- as.numeric(unlist(strsplit(substr(pt[[trait.index]][j], 2, m-1), ',')))
    Phenotype_new[j] <- var(vec) #variance
  }
  return(Phenotype_new)
}

############################################################
## COMPUTE THE HERITABILITY  ########
#####################################################
library(hglm)
compute.heritability <- function(Phenotype,G) {
  tmp <- svd(G)
  D <- tmp$d
  D[abs(D)<1e-8]=0
  L = tmp$u%*%diag(sqrt(tmp$d))
  N = length(Phenotype)
  X <- matrix(1, N, 1)
  mod1 <- hglm(y=Phenotype, X=matrix(1, N, 1), Z=L)
  h2 <- mod1$varRanef/(mod1$varFix + mod1$varRanef)
  return(h2) #Heritability
}


#########################################################################################
## PREDICT INDIVIDUAL PHENOTYPES BASED ON THE GENETIC DATA ONLY  ########
#########################################################################################
cross.predict <- function(Phenotype, G, proportion.to.predict=0.25) {
  tmp <- svd(G)
  D <- tmp$d
  D[abs(D)<1e-8]=0
  L = tmp$u%*%diag(sqrt(tmp$d))
  N = length(Phenotype)
  hidden <- sort(sample(1:N, round(N*proportion.to.predict)))
  X <- matrix(1, N-length(hidden),1)
  mod1.subset <- hglm(y=Phenotype[-hidden], X=X, Z=L[-hidden,])
  h2 <- mod1.subset$varRanef/(mod1.subset$varFix + mod1.subset$varRanef)
  BLUP <- L[hidden,]%*%mod1.subset$ranef
  ratio = as.numeric(sd(BLUP)/sd(Phenotype[-hidden]))
  predicted <- as.numeric(mod1.subset$fixef + BLUP/ratio)
  list(predicted=predicted, unobserved.phenotypes=Phenotype[hidden])
}



###################################################################
####### CREATING DATAFRAME FOR ALL PHENOTYPES - USING MEAN ########
###################################################################

Phenotype_Mean <- function(){
  
  ############### scoreIndividual ##################
  
  scoreIndividual <- get.phenotype(pt, trait="scoreIndividual")
  
  ############### scoreGroup ##################
  
  scoreGroup <- get.phenotype(pt, trait="scoreGroup")
  
  ############### score0 ##################
  
  score0 <- get.phenotype(pt, trait="score0")
  
  ############### score1 ##################
  
  score1 <- get.phenotype(pt, trait="score1")
  
  ############### score2 ##################
  
  score2 <- get.phenotype(pt, trait="score2")
  
  ############### score3 ##################
  
  score3 <- get.phenotype(pt, trait="score3")
  
  ############### beepSent ##################
  # Gives Zero
  
  beepSent <- get.phenotype(pt, trait="beepSent")
  
  ############### beepReceived ##################
  # Gives Zero
  
  beepReceived <- get.phenotype(pt, trait="beepReceived")
  
  ############### turnLeft ##################
  # Gives Zero
  
  turnsLeft <- get.phenotype(pt, trait="turnsLeft")
  
  ############### turnRight ##################
  # Gives Zero
  
  turnsRight <- get.phenotype(pt, trait="turnsRight")
  
  ############### didNothing ##################
  #Gives Zero
  
  didNothing <- get.phenotype(pt, trait="didNothing")
  
  ############### moveForward ##################
  # Gives Zero
  
  #movedForward <- get.phenotype(pt, trait="movedForward")
  
  ############### pickedGras ##################
  
  pickedGras <- get.phenotype(pt, trait="pickedGras")
  
  ############### droppedGras ##################
  # Gives Zero
  
  
  ############### handedGrasOver ##################
  # Gives Zero
  
  
  ############### sawNothing ##################
  
  sawNothing <- get.phenotype(pt, trait="sawNothing")
  
  ############### sawGras ##################
  
  sawGras <- get.phenotype(pt, trait="sawGras")
  
  ############### sawAgent ##################
  
  sawAgent <- get.phenotype(pt, trait="sawAgent")
  
  ############### sawWall ##################
  
  sawWall <- get.phenotype(pt, trait="sawWall")
  
  
  ############### Creating DF ##############
  
  phenotype_pairs_mean <- data.frame(sawWall, sawAgent, sawGras, sawNothing, 
                                pickedGras, didNothing, turnsRight, score3, 
                                score2, score1, score0, scoreGroup, scoreIndividual)
  
  return(phenotype_pairs_mean)
}






#######################################################################
####### CREATING DATAFRAME FOR ALL PHENOTYPES - USING VARIANCE ########
#######################################################################


Phenotype_Var <- function(){
  
  ############### scoreIndividual ##################
  
  scoreIndividual <- get.phenotype(pt, trait="scoreIndividual")
  
  ############### scoreGroup ##################
  
  scoreGroup <- get.phenotype(pt, trait="scoreGroup")
  
  ############### score0 ##################
  
  score0 <- get.phenotype(pt, trait="score0")
  
  ############### score1 ##################
  
  score1 <- get.phenotype(pt, trait="score1")
  
  ############### score2 ##################
  
  score2 <- get.phenotype(pt, trait="score2")
  
  ############### score3 ##################
  
  score3 <- get.phenotype(pt, trait="score3")
  
  ############### beepSent ##################
  # Gives Zero
  
  beepSent <- get.phenotype(pt, trait="beepSent")
  
  ############### beepReceived ##################
  # Gives Zero
  
  beepReceived <- get.phenotype(pt, trait="beepReceived")
  
  ############### turnLeft ##################
  # Gives Zero
  
  turnsLeft <- get.phenotype(pt, trait="turnsLeft")
  
  ############### turnRight ##################
  # Gives Zero
  
  turnsRight <- get.phenotype(pt, trait="turnsRight")
  
  ############### didNothing ##################
  #Gives Zero
  
  didNothing <- get.phenotype(pt, trait="didNothing")
  
  ############### moveForward ##################
  # Gives Zero
  
  #movedForward <- get.phenotype(pt, trait="movedForward")
  
  ############### pickedGras ##################
  
  pickedGras <- get.phenotype(pt, trait="pickedGras")
  
  ############### droppedGras ##################
  # Gives Zero
  
  
  ############### handedGrasOver ##################
  # Gives Zero
  
  
  ############### sawNothing ##################
  
  sawNothing <- get.phenotype(pt, trait="sawNothing")
  
  ############### sawGras ##################
  
  sawGras <- get.phenotype(pt, trait="sawGras")
  
  ############### sawAgent ##################
  
  sawAgent <- get.phenotype(pt, trait="sawAgent")
  
  ############### sawWall ##################
  
  sawWall <- get.phenotype(pt, trait="sawWall")
  
  
  ############### Creating DF ##############
  
  phenotype_pairs_var <- data.frame(sawWall, sawAgent, sawGras, sawNothing, 
                                pickedGras, didNothing, turnsRight, score3, 
                                score2, score1, score0, scoreGroup, scoreIndividual)
  
  return(phenotype_pairs_var)
}