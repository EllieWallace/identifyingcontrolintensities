####################################

# FULL SCRIPT FOR OUR MODEL 
# Ellie Wallace and Timothy E. Walsworth
# Department of Watershed Sciences
# The Ecology Center
# Utah State University
# as used for publication in the Journal of Applied Ecology
# final combined copy generated on: 08-22-2024
# article title: "Identifying control intensities to overcome aquatic invader resilience across a range of ecosystem sizes"
# Wallace and Walsworth, 2024

####################################
####################################
# GENERALIZED MODEL CODE
####################################
####################################

unlink(".RData")
rm(list=ls())

#===========================================================
# Load in data from the estimation model runs and set
#  output working directory
#===========================================================
setwd("C:/Users/A02353829/Documents/CarpCode")
load("VaryQ_boundSR2024_2024-03-06.RData") #load para estimates
setwd("C:/Users/A02353829/Documents/CarpCode/PB_update_2024") #set working directory

nysim<-50                                                   # How many years do you want to project the populations forward?
ny<-nysim+nyr                                                # How many total years in storage vectors (simulated plus observed)
niter<-1000                                                 #How many iterations to run
matage<-matureage                                            #age of sexual maturity 
set.seed(10)                                                 # Pick a seed to keep results consistent
sigrec<-1
pb<- c(0,1,2,3,4,5,6,7,8,9, 10,20,30,40, 50,60,70,80,90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,2000,3000,4000,5000,6000,7000,8000,9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000,17000,18000,19000,20000,30000,40000,50000) #poison bait quantities
pbx<-length(pb) #number of poison bait quantities.
ld<-c(0.005,0.01) #lethal dosage potencies
ldx<-length(ld) #number of lethal dosage potencies

amod<-c(.3,.5,1,2,10) #population productivity modifications
bmod<-c(1,10,100,1000,10000,100000) #used for carrying capacity modification for ecosystem sizes

ab.scen<-expand.grid(amod, bmod) #combine amod and bmod
xab.scen<-nrow(ab.scen) #number of scenarios for ecosystem size and productivities

NatPred<-array(data=NA,dim=c(nage,ny,niter,pbx,ldx, xab.scen))   # Storage array for number at age
poison.by.age<-NatPred #storage array for number of poisoned individuals at age
poison.biomass<-NatPred #storage array for biomass of poisoned individuals at age

###for calculating stable-state abundance for each ecosystem size/population productivity
mt<-apply(mcmcoutput[,grep("m",colnames(mcmcoutput))],
          MARGIN=2,median)
m<-as.numeric(mt)
arec<-median(mcmcoutput[,colnames(mcmcoutput)=="RecA"])
brec<-median(mcmcoutput[,colnames(mcmcoutput)=="RecB"])

a1<-arec
b1<-brec

b2<- c(ab.scen[,2])*b1 
a2<-c(ab.scen[,1])*a1

nage<-8
nyrr<-1000

Natage3<-array(1000,dim=c(nage,nyrr,length(b2)))

for(k in 1:length(b2)){
  for(i in 2:nyrr){
    for(j in nage:2){
      if(j==nage)  Natage3[j,i,k]<-Natage3[j-1,i-1,k]*(1-m[j-1])+Natage3[j,i-1,k]*(1-m[j])
      if(j!=nage)  Natage3[j,i,k]<-Natage3[j-1,i-1,k]*(1-m[j-1])
    }
    spns<-sum(Natage3[c(4:nage),i,k]*wgts[4:nage])
    Natage3[1,i,k]<-a1*spns*exp(-b2[k]*spns)
  }
}

Natage2<-array(1000,dim=c(nage,nyrr,length(b2))) #for carrying capacity with combination of ES and PP scens

xa<-array(NA, dim=c(xab.scen)) ###x as calculated in eq 4b.
for(k in 1:xab.scen){
  for(i in 2:nyrr){
    for(j in nage:2){
      if(j==nage)  Natage2[j,i,k]<-Natage2[j-1,i-1,k]*(1-m[j-1])+Natage2[j,i-1,k]*(1-m[j])
      if(j!=nage)  Natage2[j,i,k]<-Natage2[j-1,i-1,k]*(1-m[j-1])
    }
    spns<-sum(Natage2[c(4:nage),i,k]*wgts[4:nage])
    if(ab.scen[k,1]!=1) x<-(1-((log(arec/(ab.scen[k, 1]*arec)))/((ab.scen[k, 2]*brec)*(sum(Natage3[(matage:nage),1000,k]*wgts[(matage:nage)])))))
    if(ab.scen[k,1]==1) x<-1    #### factor in scen where alpha isn't modified
    Natage2[1,i,k]<-ab.scen[k,1]*arec*spns*exp(-brec*ab.scen[k,2]*x*spns)
  }
  xa[k]<-x 
}

xa
remove(m)
#====================================
#Simulation model:
#====================================


for(ww in 1:niter)
{
  mcuse<-sample(seq(1,nrow(mcmcoutput)),1,FALSE) #draws random row number from mcmcoutput without repeats
  print(ww) #shows progress when model is running
  arec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecA"]) #population productivity
  brec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecB"]) #density dependent factor for recruitment
  rec.anoms<-rnorm(nysim,-sigrec^2/2,sigrec) #allows recruitment anomolies 
  mort<-mcmcoutput[mcuse,grep("m",colnames(mcmcoutput))] #pulls mortality from mcmc estimates
  m<-as.numeric(mort[1:nage]) #format
  
  for(sld in 1:ldx){
    sigrec<-1
    #===========================================================
    # Run simulations forward
    #===========================================================
    for(pbi in 1:pbx){   #for each poison bait quantity
      for(ab in 1:xab.scen){ #for each PP/ES scenario
        NatPred[,nyr,ww,pbi,sld, ab]<-Natage2[,1000, ab]  #applies stable-state # of individuals as the starting values
        for(i in nyr:ny){   #for each simulated year
          if(i>nyr){ 
            
            for(j in 2:nage)                                         # For all ages greater than 0
            {
              if(j<nage) NatPred[j,i,ww,pbi,sld, ab]<-(NatPred[j-1,i-1,ww,pbi,sld, ab]*(1-m[j-1]))  # Calculate the number of non-plus-age individuals
              if(j==nage) NatPred[j,i,ww,pbi,sld, ab]<-(NatPred[j-1,i-1,ww,pbi,sld, ab]*(1-m[j-1]))+
                  (NatPred[j,i-1,ww,pbi,sld, ab]*(1-m[j])) # Calculate the number of plus-age individuals
            }
            NatPred[1,i,ww,pbi,sld, ab]<-exp(log(arec*ab.scen[ab, 1])+(((-brec*(xa[ab]*ab.scen[ab, 2]))*sum(NatPred[c((matage):nage),i,ww,pbi,sld, ab]*wgts[(matage):nage])+rec.anoms[i-nyr])))*sum(NatPred[c((matage):nage),i,ww,pbi,sld, ab]*wgts[(matage):nage]) #### added alpha and beta coeff
          }
          
          if(sum(NatPred[,i,ww,pbi,sld, ab])>0){ #for sim of poison bait application
            acomp<-(NatPred[,i,ww,pbi,sld, ab]/sum(NatPred[,i,ww,pbi,sld, ab])) #age composition
            pool<-rmultinom(1,100000000,acomp) #selects number of individuals based on acomp
            poisoned<-round(pool*(pb[pbi]/sum(pool*wgts*ld[sld]))) #calcs the number of individuals poisoned for each age class
            poison.by.age[,i,ww,pbi,sld, ab]<-poisoned #saves number of individuals poisoned
            poison.biomass[,i,ww,pbi,sld, ab]<-(poisoned*wgts) #saves poisoned biomass by age
            NatPred[,i,ww,pbi,sld, ab]<-(NatPred[,i,ww,pbi,sld, ab]-poisoned) #removes poisoned individuals from population
            acomp<-(NatPred[,i,ww,pbi,sld, ab]/sum(NatPred[,i,ww,pbi,sld, ab])) #recalcs age composition
            NatPred[NatPred[,i,ww,pbi,sld, ab]<0,i,ww,pbi,sld, ab]<-0 #for if poison bait is greater than number of individuals remaining, this is important for scenarios where eradication occurs
          }
        }
      }
    } 
  }
}



print(NatPred)

####################################
####################################
# CASE STUDY MODEL CODE
####################################
####################################


unlink(".RData")
rm(list=ls())

#===========================================================
# Load in data from the estimation model runs and set
#  output working directory
#===========================================================
setwd("C:/Users/A02353829/Documents/CarpCode")
load("VaryQ_boundSR2024_2024-03-06.RData") #loads parameter estimates

sim.lake.level<-function(ny){ #function to create a simulated lake level in our model
  deg.to.rad<-pi/360
  lev<-cos(seq(1,720)*deg.to.rad*50)+rnorm(720,0,.25)
  start.yr<-sample(seq(1,720-ny),1)
  stop.yr<-start.yr-1+ny
  return(lev[start.yr:stop.yr])
}


#===========================================================
# Set the harvest scenarios to be examined
#  - create a vector of annual harvest biomass to be
#     simulated
#===========================================================
effort.scenario<-c(0,1)*mean(effcomm) #######FOR PBM
nysim<-50 #change back to 100, possibly 50                                                   # How many years do you want to project the populations forward?
ny<-nysim+nyr                                                # How many total years in storage vectors (simulated plus observed)
n.scen.harv<-length(effort.scenario)                         # How many effortlevels to test
niter<-1                                                 # number of iterations
matage<-matureage                                            # age of sexual maturation
set.seed(10)                                                 # Pick a seed to keep results consistent
qt<-exp(median(mcmcoutput[,colnames(mcmcoutput)=="lnqc1"]))    # Read median catchability value from MCMC output
s.est<-apply(mcmcoutput[,grep("sel",colnames(mcmcoutput))],
             MARGIN=2,median)                                    # Read median selectivity by age values from MCMC output

s<-as.numeric(s.est)                                                # Selectivity scenario from observed data

sigrec<-1

Natqs<-apply(Natmcmc,MARGIN=2,median)                     #calcs starting abundance numbers
Natmed<-matrix(as.numeric(Natqs),nrow=8,ncol=nyr,byrow=T) #format starting abundance numbers


pb<- c(0,1,2,3,4,5,6,7,8,9, 10,20,30,40, 50,60,70,80,90, 100,125,150,175, 200,225,250,275, 300,325,350,375, 400,425,450,475, 500,550, 600,650, 700,750, 800,850, 900,950, 1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500, 10000, 11000, 12000, 13000, 14000, 15000, 16000,17000,18000,19000,20000,25000,30000,35000,40000,45000,50000) #poison bait quantities explored (kg)
pbx<-length(pb) #number of poison bait quantities explored
ld<-c(0.005,0.01) #lethal dosage potencies
ldx<-length(ld) #number of lethal dosage potencies


NatPred<-array(data=NA,dim=c(nage,ny,n.scen.harv,niter,pbx,ldx))   # Storage array for Numbers at age
poison.by.age<-NatPred                                        # Storage array for number of individuals poisoned at age
poison.biomass<-NatPred                                       # Storage array for biomass poisoned at age
catch.by.age<-NatPred                                        # Storage array for catch at age
cp.by.age<-NatPred                                           # Storage array for catch at age proportion
effort<-array(data=NA,dim=c(nysim+1,n.scen.harv,niter,1))    # Storage matrix for commercial effort by year
zlakedat<-(lakearea-mean(lakearea))/sd(lakearea)             #z-scores historic lake level data for application in our model


#===========================================================
# Create a function to calculate the commercial effort
#   required to harvest the set annual harvest level given
#   numbers at age and selectivity by age
#===========================================================
calc.eff<-function(eff,Naty,q,s,wgts,set.harvest){
  c.age<-rep(0,nage)                                         # Storage for catch at age
  h.age<-rep(0,nage)                                         # Storage for harvest at age
  for(z in 1:nage)                                           # For each age
  {
    c.age[z]<-Naty[z]*(1-exp(-1*q*s[z]*eff))                 # Calculate catch at age given effort
    h.age[z]<-c.age[z]*wgts[z]                               # Calculate harvest biomass at age
  }
  return(sum(h.age)-set.harvest)                             # Return the difference between harvest at given effort and the harvest goal
}

for(ww in 1:niter){ #for each iteration
  mcuse<-sample(seq(1,nrow(mcmcoutput)),1,FALSE) #draws random row number from mcmcoutput for each iteration without repeats
  print(ww) #shows progress when model is running
  arec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecA"]) #population productivity for recruitment
  brec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecB"]) #density dependent parameter for recruitment
  lakeb<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="LakeB"]) #lake effect for recruitment
  rec.anoms<-rnorm(nysim,-sigrec^2/2,sigrec) #recruitment anomolies 
  s.est<-mcmcoutput[mcuse,grep("sel",colnames(mcmcoutput))] #selectivity for all age classes for mechanical removal                         
  mort<-mcmcoutput[mcuse,grep("m",colnames(mcmcoutput))] #mortality for all age classes
  m<-as.numeric(mort[1:nage]) #transform to usable format
  s<-as.numeric(s.est) #transform to usable format
  
  zlakeareasim<-sim.lake.level(nysim) #creates lake-level fluctuation for each iteration
  
  
  ###set ld scenario for loop 
  for(sld in 1:ldx){ #for each lethal dosage potency
    
    sigrec<-1
    
    #===========================================================
    # Run simulations forward
    #===========================================================
    
    for(v in 1:length(effort.scenario)){                                      # For each harvest scenario
      for(pbi in 1:pbx){ #for each poison bait scenario
        
        catch.by.age[,1:(nyr-1),v,ww,pbi,sld]<-catch.out[,1:(nyr-1)]                      # Set initial catch at age from model estimates
        NatPred[,1:nyr,v,ww,pbi,sld]<-Natmed #set starting # of individuals to historical estimates   
        
        for(i in nyr:ny){   
          if(i>nyr){
            
            for(j in 2:nage){                                         # For all ages greater than 0
              if(j<nage) NatPred[j,i,v,ww,pbi,sld]<-(NatPred[j-1,i-1,v,ww,pbi,sld]*(1-cp.by.age[j-1,i-1,v,ww,pbi,sld]))*(1-m[j-1])  # Calculate the number of non-plus-age individuals
              if(j==nage) NatPred[j,i,v,ww,pbi,sld]<-((NatPred[j-1,i-1,v,ww,pbi,sld]*(1-cp.by.age[j-1,i-1,v,ww,pbi,sld]))*(1-m[j-1]))+
                  ((NatPred[j,i-1,v,ww,pbi,sld]*(1-cp.by.age[j,i-1,v,ww,pbi,sld]))*(1-m[j])) # Calculate the number of plus-age individuals
            }
            # recruitment function, eqn 7 (see supplement): 
            NatPred[1,i,v,ww,pbi,sld]<-exp(log(arec)+(((-brec)*sum(NatPred[c((matage):nage),i,v,ww,pbi,sld]*wgts[(matage):nage])+zlakeareasim[i-nyr]*lakeb+rec.anoms[i-nyr])))*sum(NatPred[c((matage):nage),i,v,ww,pbi,sld]*wgts[(matage):nage]) 
          }
          
          if(sum(NatPred[,i,v,ww,pbi,sld])>0){ #simulation for poison bait application
            acomp<-(NatPred[,i,v,ww,pbi,sld]/sum(NatPred[,i,v,ww,pbi,sld])) #calculate age composition
            pool<-rmultinom(1,100000000,acomp) #selects number of individuals based on acomp
            poisoned<-round(pool*(pb[pbi]/sum(pool*wgts*ld[sld]))) #calcs the number of individuals poisoned for each age class
            poison.by.age[,i,v,ww,pbi,sld]<-poisoned #saves number of individuals poisoned
            poison.biomass[,i,v,ww,pbi,sld]<-(poisoned*wgts) #saves poisoned biomass by age
            NatPred[,i,v,ww,pbi,sld]<-(NatPred[,i,v,ww,pbi,sld]-poisoned) #removes poisoned individuals from population
            acomp<-(NatPred[,i,v,ww,pbi,sld]/sum(NatPred[,i,v,ww,pbi,sld])) #recalcs age composition
            NatPred[NatPred[,i,v,ww,pbi,sld]<0,i,v,ww,pbi,sld]<-0 #for if poison bait is greater than number of individuals remaining, this is important for effort scenarios where eradication might occur
          }
          
          
          if(i>nyr) q<-exp((mcmcoutput[mcuse,colnames(mcmcoutput)=="lnqc1"]+zlakeareasim[i-nyr]*mcmcoutput[mcuse,colnames(mcmcoutput)=="qslope"]))    # Calc catchability for simulation
          if(i<=nyr) q<-exp((mcmcoutput[mcuse,colnames(mcmcoutput)=="lnqc1"]+zlakedat[i]*mcmcoutput[mcuse,colnames(mcmcoutput)=="qslope"])) #calc catchability for prior years
          
          if(is.na(q)) stop("Q problem")
          
          for(k in 1:nage){                                           # For each age
            catch.by.age[k,i,v,ww,pbi,sld]<-NatPred[k,i,v,ww,pbi,sld]*(1-exp(-1*(q)*s[k]*effort.scenario[v])) # Calculate catch at age
            cp.by.age[k,i,v,ww,pbi,sld]<-catch.by.age[k,i,v,ww,pbi,sld]/NatPred[k,i,v,ww,pbi,sld] #put catch into proportion for application to NatPred calcs
            cp.by.age[cp.by.age[k,i,v,ww,pbi,sld]>1,i,v,ww,pbi,sld]<-1 #ensures high efforts can be applied without NaN creation
            if(is.nan(cp.by.age[k,i,v,ww,pbi,sld])) cp.by.age[k,i,v,ww,pbi,sld] <- 0 #allows eradication of an age-class/population without NaN creation
          }
        }
      }
    } 
  }
}


print(NatPred)