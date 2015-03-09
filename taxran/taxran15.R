### The `taxran` function takes a dataframe where taxonomic level are supplied in columns.
### In addition, the dataframe contains continuous traits supplied in columns.
### The function returns the proportion of total variation explained for a trait at each taxonomic level. 
### The function also returns null expectation of explained variation for `n` randomizations of data. 

taxran15 <- function(data,trait){

  # Load `lme4` package
  require(lme4)  
  
  # Number of randomizations
  n = 10 

  # Create 7 x n output matrix of zeros
  output <- matrix(0, 7, n) 

  # Add taxa rownames
  rownames(output) <- c(colnames(data)[1:6], "Residual(Species)") 
  
  # Use linear mixed-effect model (`lmer`) with nested taxonomic heirarchy
  B <- lmer(data[,trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)

  # Calculate variance explained at each level
  vars <- as.numeric(lme4::VarCorr(B))

  # Exract residual variance (a.k.a., species-level variance)
  residvar <- var(resid(B))
  
  # Calculate percent variance explained at each taxonomic level
  pctvar <- c(vars, residvar) / sum(c(vars, residvar))
  
  # Add names to percent variance output
  names(pctvar) <- c(colnames(data)[1:6], "Residual(Species)")

  # Now, simulate the null expectation of no taxonomic structure 
  # Observed trait values are randomly assigned to a speices
 
  for(i in 1:n){
   
    # Permute trait values
    data[,trait] <- sample(data[, trait])
    
    # Fit the randomized data to the `lmer` model to calculate percent variance explained
    sim.B <- lmer(data[, trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)
    
    # Calculate variance explained at each level
    sim.vars <- as.numeric(VarCorr(sim.B))
    
    # Calculate residual variance (equivalent to species-level variance)
    sim.residvar <- var(resid(sim.B))
    
    # Calculate percent variance explained at each taxonomic level
    sim.pctvar <- c(sim.vars, sim.residvar) / sum(c(sim.vars, sim.residvar))
    
    # Create output
    output[,i] <- sim.pctvar    
  }
  
  # Defines list of observed and simulated percent variance explained
  return(list(act.pct = pctvar, sim.pct = output))
}