source("covidmodels.R")

#Parameters

n<-100 #Initial number of active cases
e<-(2^(1/2)-1) #Number of encounters
P<-1e10 #total population

N=vector()
dN=vector()
t=1 #time in days
print(c(n,P,t))
      
while (t < 1000)
  {
    dN[t] = imexp(n,e,P)
    
    #Update Variables
    N[t] = n+dN[t]
    n=N[t] #Update number of active cases
    t=t+1
  }

d <- 1:length(N)
COVID <- data.frame(d,N,dN)
cnames <- c("Days","N","Nnew")

#Make Plot

g <- ggplot(COVID, aes(x=d, y=N)) + geom_line()
g <- g + scale_x_continuous(name="Number of Days", limits=c(0, length(N))) + scale_y_continuous(name="Number of Active Cases",trans='log10',labels=comma)
plot(g)


