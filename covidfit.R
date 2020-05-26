#Covid Max Likelihood Code
#Read in data
#loop through parameters - r(0.05:0.4::0.01), e(0:85:5), tau (range 0:t) - given n_0 and 
#generate model predictions
#compute likelihood for each set of parameters
#create data frame with parameter values, likelihood,dN derived parameter values

#Source scripts containing functions for reading data and computing likelihoods
source("covidmodels.R")
source("readdatacovid.R")

#Read in data
#data <- readcoviddata("India")
data <- readtrackerdata(100)

#Set parameter ranges

#Fixed parameters
n0 = 100 # n0 - initial cases
P = 1.37e9 # P - population
ndays = 300 # ndays - number of days

#Variable parameters with range and step size
# oldcases - pre-existing cases
# r - expansion rate
# erec - recovery parameter
# emort - mortality parameter
# tau - time delay parameter for recovery or death

#oldcases <- seq(0,0.9*n0,0.1*n0)
oldcases <- seq(0,40,5)
r <- seq(0.02,0.25,0.01)
erec <- seq(0.4,0.8,0.1)
emort <- seq(0.01,0.05,0.01)
tau <- seq(5,20,5)
#erec <- seq(0.5,0.8,0.1)
#emort <- seq(0.01,0.05,0.01)
#tau <- seq(0,25,5)
#erec <- seq(0.1,0.9,0.1)
#emort <- seq(0.01,0.2,0.02)
#tau <- seq(0,21,3)

likes <- data.frame()


for(i in 1:length(oldcases))
{
  for(j in 1:length(r))
  {
    for(k in 1:length(erec))
    {
      for(l in 1:length(emort))
      {
        for(m in 1:length(tau))
        {
          mod <- gencases(n0,oldcases[i],r[j],P,ndays,erec[k],emort[l],tau[m])
          chisq <- covidchisqs(data,mod)
          
          
          df <- rbind(c(oldcases[i],r[j],erec[k],emort[l],tau[m],chisq))
          likes <- rbind(likes,df)
          
          
                    
        }
      }
    }
    print(c(i,j))
  }
}

cnames <- c("OldCases","R","E(rec)","E(mort)","tau","chisq") 
names(likes) <- cnames

