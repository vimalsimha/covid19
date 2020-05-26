#Functions for Computing Spread of Covid-19

newcases <- function(n,r,P,ntot)
{
  newcases = n*r*(1-ntot/P)
  newcases
}

gencases <- function(n0,oldcases,r,P,ndays,erec,emort,tau)
  # n0 - initial cases
  # oldcases - pre-existing cases
  # r - expansion rate
  # P - population
  # ndays - number of days
  # erec - recovery parameter
  # emort - mortality parameter
  # tau - time delay parameter for recovery or death
  
  # In principle, oldcases, r, erec, emort and tau are free parameters while n0, P and ndays must be supplied.
{
  
  N=vector()
  dN=vector()
  RN=vector() #Recovered Cases
  MN=vector() #Deaths
  NACT=vector() #Active Cases
  
  t=1 #time in days
  n=n0 #number of cases (total)
  nactive=n0 # number of active cases
  deadcases=0 # Dead Cases
  reccases=0 # Recovered Cases
  elimcases=0 # Eliminated Cases
  tmin=14 # minimum number of days before elimination
  
  while (t <= ndays)
  {
    dN[t] = newcases(nactive,r,P,n) #New Cases
    
    if(t>tau) #After tmin, find old cases
    {
      oldcases=N[t-tau] #Old Cases = Total Cases at t-tau - eliminated cases
    }
    
    #Update old cases and recovered/dead cases
    
    oldcases = oldcases-elimcases #remove eliminated cases. Do it outside loop because it can happen before t=tau too!
    #if(oldcases < 0) oldcases=0
    
    deadcases = deadcases + oldcases*emort
    reccases = reccases + oldcases*erec
    elimcases = deadcases+reccases
  
    #Update Variables
    N[t]=n+dN[t]
    n = N[t] # Cumulative number of infected people
    nactive=N[t]-elimcases
    #if(nactive<0) nactive=0
        
    NACT[t]=nactive
    RN[t]=reccases
    MN[t]=deadcases
    
    t=t+1
  }
  
  d <- 1:length(N)
  COVID <- data.frame(d,N,dN,NACT,MN,RN,RN+MN+NACT)
  cnames <- c("Days","Cases","New","Active","Deaths","Recovered","Total")
  names(COVID) <- cnames
  
  COVID  
}

# Likelihood Functions

chisqs <- function(data,model) 
  #return chi-squares
{
  errs = data^0.5 #Assume Poisson errors because number counts
  errs = 3*errs
  chisq = sum((data-model)^2/errs^2)
  chisq
}

covidchisqs <- function(data,model)
{
  
#  dN <- vector()
  
#  for (i in 1:length(data$Cases)-1)
#  {
#    dN[i]=data$Cases[i+1]-data$Cases[i]
#  }
 
  dN = as.numeric(data$dailyconfirmed)
  dDeaths = as.numeric(data$totaldeceased)
  dRecs = as.numeric(data$totalrecovered)
   
  chisq=vector()
  chisq[1] = chisqs(dN,model$New[1:length(dN)])
  #chisq[1]=  chisqs(data$Cases[2:items],model$Cases[1:items-1])
  chisq[2] = chisqs(dDeaths,model$Deaths[1:length(dDeaths)])
  chisq[3] = chisqs(dRecs,model$Recovered[1:length(dRecs)])
  
  
  
  covidchisqs=sum(chisq)
  covidchisqs
}

#Other models

simpexp <- function(n0,e)
{
  newcases = n0*e
  newcases
}

imexp <- function(nactive,e,P,ntot)
{
  newcases = nactive*e*(1-ntot/P)
  newcases
}

fullmod <- function(n0,e,P,el,oldcases)
{
  newcases = n0*e*(1-n0/P) - el*oldcases
  newcases
}

fillfullmod <- function(n0,e,P,nelems,el)
  # n0 - initial cases
  # e - expansion rate
  # P - population
  # n - number of vector elements
  #el - elimination rate
{
  
  N=vector()
  dN=vector()
  RN=vector() #Recovered Cases
  NACT=vector() #Active Cases
  t=1 #time in days
  n=n0 #number of cases (total)
  nactive=n0 # number of active cases
  oldcases=0 # number of old cases
  elimcases=0 # number of eliminated cases
  tmin=14 # minimum number of days before elimination
  
  while (t < nelems)
  {
    dN[t] = imexp(nactive,e,P,n) #New Cases
    
    if(t>tmin) #After tmin, eliminate cases from active pool
      {
      oldcases=N[t-tmin]-elimcases # All cases existing before time (t-tmin) minus all cases already inactive
      elimcases=elimcases + (oldcases*el) # total number of cases eliminated
      }
    
    #Update Variables
    N[t]=n+dN[t]
    n = N[t] # Cumulative number of infected people
    nactive=N[t]-elimcases
    
    NACT[t]=nactive
    RN[t]=elimcases
    
    t=t+1
  }

  d <- 1:length(N)
  COVID <- data.frame(d,N,dN,RN,NACT)
  #cnames <- c("Days","N","Nnew")
  
  COVID  
}

fillimexp <- function(n0,e,P,nelems)
  # n0 - initial cases
  # e - expansion rate
  # P - population
  # n - number of vector elements
{
  
  N=vector()
  dN=vector()
  t=1 #time in days
  n=n0
  
  while (t < nelems)
  {
    dN[t] = imexp(n,e,P,n)
    
    #Update Variables
    N[t] = n+dN[t]
    n=N[t] #Update number of active cases and number of total cases
    t=t+1
  }
  
  d <- 1:length(N)
  COVID <- data.frame(d,N,dN)
  #cnames <- c("Days","N","Nnew")
  
  COVID  
}

#To compute growth rate over last x days
grater <- function(caselist,ndays)
{
grate = vector()
  for (i in 1:length(caselist))
  {
    if (i<=ndays) grate[i]=0
    if(i>ndays)
    {
      newcases <- caselist[i]-caselist[i-ndays]
      r <- newcases/caselist[i-ndays]
      grate[i] <- (1+r)^(1/ndays)
    }
    
  }
  grate
}

