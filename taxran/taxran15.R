### The `taxran` function takes a dataframe where taxonomic level are supplied in columns.
### In addition, the dataframe contains continuous traits supplied in columns.
### The function returns the proportion of total variation explained for a trait at each taxonomic level.
### The function also returns null expectation of explained variation for `n` randomizations of data.

taxran15 <- function(input, trait){
  # Load packages
  require(nlme)
  options(warn = -1)
  
  # Number of randomizations; default = 1000
  n = 10 

  # Create 7 x n output matrix of zeros (7 = levels of taxonomic hierarchy)
  output <- matrix(NA, 7, n)

  # Add taxa rownames
  rownames(output) <- c(colnames(input)[1:6], "Residual(Species)")

  # Identify trait values
  trait <- input[,trait]
  input <- cbind(input, trait)

  # Use linear mixed-effect model (`lme`) with nested taxonomic heirarchy
  # fitting intercepts as a random variable, which we are attempting to explain via taxonomy
  B <- lme(trait ~ 1, random = ~1|Domain/Division/Class/Order/Family/Genus, data = input)

  # Obtain observed variances
  vc <- VarCorr(B)[,1]
  vars <- as.numeric(vc)[!is.na(as.numeric(vc))]
  names(vars) <- c("Domain","Division","Class","Order","Family","Genus","Residual(Species)")
  pctvar <- round(vars/sum(vars),3)

  # BEGIN PERMUTATION TESTS
  
  # Create empty matrix for output
  trait.sim <- rep(NA, length(trait))
  
  input <- cbind(input, trait.sim)

  for(i in 1:n){
    # Resample the trait data
    input$trait.sim <- sample(trait)

    # Run `lme` model with randomized data
    sim.B <- lme(trait.sim ~ 1, random = ~1|Domain/Division/Class/Order/Family/Genus, data = input)
    
    # Calculate variance from randomized data
    sim.vc <- VarCorr(sim.B)[,1]
    sim.vars <- as.numeric(sim.vc)[!is.na(as.numeric(sim.vc))]
    names(sim.vars) <- c("Type","Division","Class","Order","Family","Genus","Residual(Species)")
    sim.pctvar <- round(sim.vars/sum(sim.vars),3)
    output[,i] <- sim.pctvar
  }
  return(list(act.pct = pctvar, sim.pct = output))
}