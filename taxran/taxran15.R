### The `taxran` function takes a dataframe where taxonomic level are supplied in columns.
### In addition, the dataframe contains continuous traits supplied in columns.
### The function returns the proportion of total variation explained for a trait at each taxonomic level. 
### The function also returns null expectation of explained variation for `n` randomizations of data. 

#taxran15 <- function(input = "", trait = ""){
taxran15 <- function(input, trait){
  # Load packages
  options(warn = -1)
  require(nlme)
  
  # Number of randomizations
  n = 10 # default = 1000 

  # Create 7 x n output matrix of zeros
  output <- matrix(NA, 7, n) 

  # Add taxa rownames
  rownames(output) <- c(colnames(input)[1:6], "Residual(Species)") 
  
  # Use non-linear mixed-effect model (`nlme`) with nested taxonomic heirarchy
  # B <- lmer(data[,trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)
  # B <- lme(Optimum ~ 1, data = data, random = ~1|Domain/Division/Class/Order/Family/Genus)
  
  # Identify trait of interest
  
  # trait <- "Optimum"
  #trait <- input[,trait]
  # dat1 <- data.frame(data, resp, sim.resp = rep(0, nrow(data)))
  # resp is a function of intercept; which is a random variable, which we are attempting to explain via tax
  #B <- lme(resp ~ 1, random = ~1|Domain/Division/Class/Order/Family/Genus,data = dat1)
  
  B <- lme(trait ~ 1, random = ~1|Domain/Division/Class/Order/Family/Genus, data = input)
  
  vc <- VarCorr(B)[,1]
  vars <- as.numeric(vc)[!is.na(as.numeric(vc))]
  names(vars) <- c("Domain","Division","Class","Order","Family","Genus","Residual(Species)")
  
  pctvar <- round(vars/sum(vars),3)
  
  # Beginning of permutation test
  trait.sim <- rep(NA, length(trait))
  
  for(i in 1:n){
  
    trait.sim <- sample(trait)
    sim.B <- lme(trait.sim ~ 1, random = ~1|Domain/Division/Class/Order/Family/Genus, data = input)
    sim.vc <- VarCorr(sim.B)[,1]
    sim.vars <- as.numeric(sim.vc)[!is.na(as.numeric(sim.vc))]
    names(sim.vars) <- c("Type","Division","Class","Order","Family","Genus","Residual(Species)")
    sim.pctvar <- round(sim.vars/sum(sim.vars),3)
    output[,i] <- sim.pctvar    
  }
  return(list(act.pct = pctvar, sim.pct = output))
}
    
  
  # Calculate variance explained at each level
  #vars <- as.numeric(lme4::VarCorr(B))
  
  #varcorr <- nlme::VarCorr(B)
  #vars <- as.numeric(varcorr[seq(2,nrow(varcorr), by = 2), 1]) 
    
  # Exract residual variance (a.k.a., species-level variance)
  #residvar <- var(resid(B))
  
  # Calculate percent variance explained at each taxonomic level
  #pctvar <- c(vars, residvar) / sum(c(vars, residvar))
  
  # Add names to percent variance output
  #names(pctvar) <- c(colnames(data)[1:6], "Residual(Species)")

  # Now, simulate the null expectation of no taxonomic structure 
  # Observed trait values are randomly assigned to a speices
 
#   #for(i in 1:n){
#    
#     # Permute trait values
#     data[,trait] <- sample(data[, trait])
#     
#     # Fit the randomized data to the `lmer` model to calculate percent variance explained
#     # sim.B <- lmer(data[, trait] ~ (1|Domain/Division/Class/Order/Family/Genus), data = data)
#     sim.B <- lme(Optimum ~ 1, data = data, random = ~1|Domain/Division/Class/Order/Family/Genus)
#     
#     # Calculate variance explained at each level
#     # sim.vars <- as.numeric(VarCorr(sim.B))
#     sim.varcorr <- nlme::VarCorr(sim.B)
#     sim.vars <- as.numeric(sim.varcorr[seq(2,nrow(sim.varcorr), by = 2), 1]) 
#     
#     # Calculate residual variance (equivalent to species-level variance)
#     sim.residvar <- var(resid(sim.B))
#     
#     # Calculate percent variance explained at each taxonomic level
#     sim.pctvar <- c(sim.vars, sim.residvar) / sum(c(sim.vars, sim.residvar))
#     
#     # Create output
#     output[,i] <- sim.pctvar    
#   }
#   
#   # Defines list of observed and simulated percent variance explained
#   return(list(act.pct = pctvar, sim.pct = output))
# }