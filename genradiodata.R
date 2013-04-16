# Generate data
radiodata <- function(N,tau,dt) {
  no_data<-floor(2*tau/dt)
  radiodataarray<-array(N,c(no_data))
  for (i in 1:no_data) {
    if (i==1) {
      prevN<-N
    } else {
      prevN<-radiodataarray[i-1]
    }
    deltaN<-rpois(1,prevN*dt/tau)
    #print(deltaN)
    if (deltaN<=prevN) {
      radiodataarray[i]<-prevN-deltaN
    } else {
      radiodataarray[i]<-0
    }
  }
  return(radiodataarray)
}

# Linear Regression, return lifetime
radiolnregres <- function(N,tau,dt,radiodataarray) {
  tarray<-seq(dt,2*tau,dt)
  lnradiodataarray=log(radiodataarray)
  plot(tarray,lnradiodataarray,xlab="t",ylab="ln(N)")
  linregresult<-lm(lnradiodataarray~tarray)
  abline(linregresult)
  lifetime<-(-1/linregresult$coefficients[2])
  return(lifetime)
}

# Pearson's chi-squared test, return p-value
radiochitest <- function(N,tau,dt,radiodataarray) {
  tarray<-seq(dt,2*tau,dt)
  thN<-N*exp(-tarray/tau)
  chisqresult<-chisq.test(radiodataarray,p=thN/sum(thN))
  #critchisq<-qchisq(1-alpha,df=length(tarray)-1)
  return(chisqresult)
}

# Average lifetime, return p-value of t-test
radioavgtau <- function(N,tau,dt,radiodataarray) {
  no_data<-floor(2*tau/dt)
  tauarray<-array(tau,c(no_data-1))
  for (i in 1:(no_data-1)) {
    deltaN<-radiodataarray[i+1]-radiodataarray[i]
    if (deltaN==0) {
      if (i==1) {
        tauarray[i]<-tau
      } else {
        tauarray[i]<-tauarray[i-1] 
      }
    } else {
      tauarray[i]<-N*dt/deltaN
    }
  }
  meantau<-mean(tauarray)
  stdtau<-sqrt(var(tauarray))
  t<-(meantau-tau)/(stdtau/sqrt(no_data))
  return(pt(t,df=no_data-1))
}
