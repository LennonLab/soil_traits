### This `taxran` function takes a dataframe where taxonomic level are supplied in columns.
### Taxonomic levels range from Domain to Species. 
### In addition, the dataframe contains continuous traits for supplied in columns.
### The function returns the proportion of total variation explained for a trait at each taxonomic level. 
### The function also returns `n` simulations. 
### Simulations are used to test the hypothesis of no taxonomic structure.
### This is ccomplished by randomly assigning trait values to taxa and refitting the model.

taxran15 <- function(data,trait){

require(lme4)  
  
  #set up number of times to simulate the null hypothesis
  n = 10
  output <- matrix(0, 7, n)
  rownames(output) <- c(colnames(data)[1:6], "Residual(Species)")
  
  # Use lmer with stucture nested by taxonomic heirarchy
  B <- lmer(data[,trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)

  # Collect variance explained at each level
  vars <- as.numeric(lme4::VarCorr(B))

  # Fit residual variance
  residvar <- var(resid(B))
  
  #combine variances and divide by total to get variance explained at each level
  pctvar <- c(vars, residvar) / sum(c(vars, residvar))
  names(pctvar) <- c(colnames(data)[1:6], "Residual(Species)")

  #simulate the null hypothesis of no taxonomic structure 
  #for the trait by assuming any trait could be assigned 
  #to any speices independent of taxonomy
 
  for(i in 1:n){
   
    #permute trait values
    data[,trait] <- sample(data[, trait])
    #fit the simualted data to same model as actual data
    sim.B <- lmer(data[, trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)
    
  sim.vars <- as.numeric(VarCorr(sim.B))
  sim.residvar <- var(resid(sim.B))
  sim.pctvar <- c(sim.vars, sim.residvar) / sum(c(sim.vars, sim.residvar))
  output[,i] <- sim.pctvar    
  }
  return(list(act.pct = pctvar, sim.pct = output))
}