taxran15<-function(dat,trait){
##This function takes in a dataframe with a columns for each taxonomic
##level from Type to Species and the names of a trait that is also a column name in the
##dataframe. It returns the pct of total variation in the trait
##attributable to each taxonomic level. It aslo returns n simulations of the 
##hypothesis of no taxonomic structure by randomly assigning trait values to
##species and refitting the model.
require(lme4)  
  #set up number of times to simulated the null hypothesis
  n=999
 output<-matrix(0,7,n)
 rownames(output)<-c(colnames(dat)[1:6],"Residual(Species)")
 #Use lmer with stucture nested by taxonomic heirarchy
  B<-lmer(dat[,trait]~(1|Type/Division/Class/Order/Family/Genus),data=dat)

  #collect variance explained at each level
  vars<-as.numeric(VarCorr(B))

 #get residual variance
  residvar<-var(resid(B))
  
 #combine those variances and divide by total to get
  #percent of variance explained at each level
  pctvar<-c(vars,residvar)/sum(c(vars,residvar))
names(pctvar)<-c(colnames(dat)[1:6],"Residual(Species)")
#simulate the null hypothesis of no taxonomic structure 
#for the trait by assuming any trait could be assigned 
#to any speices independent of taxonomy
  for(i in 1:n){
   #permute trait values
    dat[,trait]<-sample(dat[,trait])
   #fit the simualted data to same model as actual data
    sim.B<-lmer(dat[,trait]~(1|Type/Division/Class/Order/Family/Genus),data=dat)
    
    
  sim.vars<-as.numeric(VarCorr(sim.B))
  sim.residvar<-var(resid(sim.B))
  sim.pctvar<-c(sim.vars,sim.residvar)/sum(c(sim.vars,sim.residvar))
  output[,i]<-sim.pctvar    
  }
  return(list(act.pct=pctvar,sim.pct=output))
}