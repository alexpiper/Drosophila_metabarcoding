#########################################################################################################
#########################################################################################################
### THIS SCRIPT PROVIDES FUNCTIONS FOR SIMULATING AND ANALYSING E-DNA DATA (2-LEVEL SAMPLING STRUCTURE), 
### IN COMBINATION WITH CALIBRATION DATA (equipment and/or PCR level), AND DATA FROM AN UNAMBIGUOUS 
### DETECTION METHOD
### 
### Guillera-Arroita et al
### v 2016/12/04
#########################################################################################################
#########################################################################################################

#########################################################################################################
### GENERATE SIMULATED DATA
### Input parameters
###       S: number of sites
###       V: number of water samples per site 
###          (can be one number if the same in all sites; a vector if site specific)
###       L: number of qPCR tests in each water sample 
###          (can be one number if the same in all water samples; a matrix if sample specific)
###       K: number of repeats per site of the unambigous method 
###          (can be one number if the same in all sites; a vector if site specific)
###       mypars: (psi,theta11,theta10,p11,p10,r11); as defined in manuscript
###       nCal1: number of trials of calibration experiment (equipment level)
###              (2 numbers, number of water samples, and number of PCRs per water sample)
###       nCal2: number of trials of calibration experiment (PCR level)
###              (1 number, number of PCR blanks)
###
###  e.g. mydata<-simdataFPMS(S=20,V=3,L=2,K=2,mypars=c(0.6,0.7,0.1,0.3,0.05,0.3),nCal1=c(50,4),nCal2=50)
#########################################################################################################
simdataFPMS<-function(S,V,L,K,mypars,nCal1,nCal2){
  
  ## simulate field data: eDNA
  if (length(V)==1) V<-rep(V[1],S)    #water samples in each site (vector) - if one value given, use that for all sites
  if (length(L)==1) L<-matrix(L,nrow=S,ncol=max(V))    #PCR samples in each water sample (matrix)
  Z2 <- matrix(NA,nrow=S,ncol=max(V)) #to hold eDNA status in each water sample
  D  <- matrix(NA,nrow=S,ncol=max(V)) #to hold number of positive PCR samples in each water sample                                         
  Z1 <- 1*(runif(S)<mypars[1])   #simulate occupancy status of sites (eDNA present or not)
  for (ii in 1:S){
    thetaU<-mypars[2]*Z1[ii]+mypars[3]*(1-Z1[ii])
    for (vv in 1:V[ii]){
      Z2[ii,vv]<-1*(runif(1)<thetaU)   #simulate presence of eDNA in each water sample
      pU<-mypars[4]*Z2[ii,vv]+mypars[5]*(1-Z2[ii,vv])
      D[ii,vv]<-rbinom(n=1,L,pU)       #simulate number of PCR detections per water sample
    }
  }
  
  ## simulate field data: unambigous method (one replication level)
  if (length(K)==1) K<-rep(K[1],S)    #surveys unambiguous method in each site (vector)
  Du  <- rep(NA,nrow=S)  #number of unambigous detections in each site  
  for (ii in 1:S){
     rU<-mypars[6]*Z1[ii]
     Du[ii]<-rbinom(n=1,K[ii],rU)       #simulate number of detections per site
  }  
  
  ## simulate calibration data: equipment level (one or more PCRs per water sample blank)
  CalL1<-rep(nCal1[2],nCal1[1])
  CalD1<-rep(NA,nCal1[1])
  Cal1_Z2<-1*(runif(nCal1[1])<mypars[3])  #simulate contamination of eDNA in each water sample
  for (vv in 1:length(Cal1_Z2)){
    pU<-mypars[4]*Cal1_Z2[vv]+mypars[5]*(1-Cal1_Z2[vv])
    CalD1[vv]<-rbinom(n=1,CalL1[vv],pU)       #simulate number of PCR detections per water sample
  }
     
  ## simulate calibration data: PCR-level 
  CalD2<-rbinom(n=1,nCal2,mypars[5])
  
  #write dataset out
  mydata<-list(V=V,L=L,K=K,D=D,Du=Du,Z1=Z1,Z2=Z2,Cal1_Z2=Cal1_Z2,caldata1=list(CalL1=CalL1,CalD1=CalD1),caldata2=c(nCal2,CalD2))  
}

#########################################################################################################
### NEGATIVE LOGLIKELIHOOD FUNCTION
### Input parameters
###       params: values of the parameters (logit scale) for which the loglik is computed
###       mydata: list with the different bits of data (format as in simulation function output)
###       fixpar: can be used to fix some parameter values (NA for parameters to be estimated, 
###               otherwise fixed to the value given in probability scale) ** this is useful to test
###               special cases of the model, and required to construct profile likelihood functions **
#########################################################################################################
loglikf_FPMS <- function(params, mydata, fixpar) { 
  
  V<-mydata$V; L<-mydata$L; D<-mydata$D; #eDNA survey data 
  K<-mydata$K; Du<-mydata$Du; caldata1<-mydata$caldata1; caldata2<-mydata$caldata2; #additional data 
  
  npar<-1
  if (is.na(fixpar[1])) {psi <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{psi<-fixpar[1]}
  if (is.na(fixpar[2])) {theta11 <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{theta11<-fixpar[2]} 
  if (is.na(fixpar[3])) {theta10 <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{theta10<-fixpar[3]} 
  if (is.na(fixpar[4])) {p11 <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{p11<-fixpar[4]} 
  if (is.na(fixpar[5])) {p10 <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{p10<-fixpar[5]} 
  if (is.na(fixpar[6])) {r11 <- 1/(1+exp(-params[npar])); npar<-npar+1;}else{r11<-fixpar[6]} 
  
  # field data component (ambiguous and unambiguous method)
  ND <- L-D   # number of non detections (ambiguous method)
  tmpZ2_1 <- p11^D * (1-p11)^ND
  tmpZ2_0 <- p10^D * (1-p10)^ND
  tmp <- theta11*tmpZ2_1 + (1-theta11)*tmpZ2_0
  tmpZ1_1 <- apply(tmp, 1, prod)
  tmp <- theta10*tmpZ2_1 + (1-theta10)*tmpZ2_0
  tmpZ1_0 <- apply(tmp, 1, prod)
  NDu <- K-Du #number of non detections (unambiguous method)
  tmpZ1_1u <- r11^Du * (1-r11)^NDu
  tmpZ1_0u <- 0^Du * 1^NDu
  loglik<-sum(log(psi * tmpZ1_1 * tmpZ1_1u + (1-psi) * tmpZ1_0 * tmpZ1_0u))
  
  # calibration data component (equipment level)
  caldata1.U<-caldata1$CalL1-caldata1$CalD1
  tmp <- theta10 * p11^caldata1$CalD1 * (1-p11)^caldata1.U +
         (1-theta10) * p10^caldata1$CalD1 * (1-p10)^caldata1.U
  loglik<-loglik + sum(log(tmp))
  
  #calibration data component (PCR level)
  if (caldata2[1]>0) loglik<-loglik + caldata2[2]*log(p10)+(caldata2[1]-caldata2[2])*log(1-p10)

  #change sign for minimization
  loglik <- -loglik

}


#################################################################################################
### MODEL FITTING (find MLEs via numerical optimization, run several times with different starts)
### Input parameters
###       mydata: list with the different bits of data (format as in simulation function output)
###       fixpar: can be used to fix some parameter values (NA for the param to be estimated, 
###               otherwise fixed to the value given)
###       nstarts: number of times the optimization is run (with different starting values)
###       meanstarts: mean value used to generate starting values for parameters (vector, logit scale)
###       sdstart: spread in generation of starting values (normally distributed)
###       method: optimization algorithm used
#################################################################################################
fit_model<-function(mydata,fixpar,nstarts=10,meanstarts,sdstart=1,method="Nelder-Mead",dohess=F){
  
  parshat_m <- matrix(fixpar,nrow=nstarts,ncol=6,byrow=T) #to hold estimates from each optimization
  theLLvals <- rep(NA,nstarts) #to hold the loglikelihood value from each optimization
  optimres  <- list() #to store complete model output from optim in each optimization
  
  for (ii in 1:nstarts){
    mystarts<-rnorm(n=length(meanstarts),mean=meanstarts,sd=sdstart)  # init to given values + noise
    tmp <- optim(mystarts,loglikf_FPMS,mydata=mydata,fixpar=fixpar,method=method,hessian=dohess);
    theLLvals[ii]<-tmp$value;parshat_m[ii,is.na(fixpar)]<-round(plogis(tmp$par),3)
    optimres[[ii]]<-tmp
  }
  
  diff<-(theLLvals-min(theLLvals))
  allres<-cbind(parshat_m,theLLvals,diff=round(diff,2))
  return(list(myMLEs=optimres[[which.min(theLLvals)]],allres=allres))
}

#################################################################################################
### CONSTRUCT PROFILE LOGLIKELIHOOD FOR ALL MODEL PARAMETERS
### Input parameters
###       mydata: list with the different bits of data (format as in simulation function output)
###       fixpar: can be used to fix some parameter values (NA for the param to be estimated, 
###               otherwise fixed to the value given)
###       nstarts: number of times the optimization is run (with different starting values)
###       meanstarts: mean value used to generate starting values for parameters (vector, logit scale)
###       sdstart: spread in generation of starting values (normally distributed)
###       method: optimization algorithm used
###       thestep: step used to construct the profile loglikelihood function
#################################################################################################
get_profCIall<-function(mydata,fixpar,nstarts=10,sdstart=2,method="Nelder-Mead",thestep=0.02){
  
  parnames <- c("psi","theta11","theta10","p11","p10","r11")
  npar <- sum(is.na(fixpar))   # total number of parameters in the model (those not fixed)
  meanstarts <- rep(0,npar-1)  # optim has one dimension less than the total number of parameters
  theseq <- seq(0.01,0.99,by=thestep) # parameter values to sweep through
  profliks <- matrix(NA,ncol=length(theseq),nrow=6) # to hold the profile loglik function values
  for (ii in 1:6){
    if (is.na(fixpar[ii])){   #this is a parameter for which to calculate the prof loglik
        for (jj in 1:length(theseq)){
          print(paste0(parnames[ii],"=",theseq[jj]))
          fixpar2<-fixpar; fixpar2[ii]<-theseq[jj]  #fix the parameter to the corresponding value
          m1<-fit_model(mydata,fixpar2,nstarts=nstarts,meanstarts,sdstart=sdstart)  #max ll wrt remaining params
          profliks[ii,jj]<-m1$myMLEs$value  #keep value of ll function at the maximum 
      }
    }
  }
  return(list(profliks=profliks,theseq=theseq))
}

#################################################################################################
### PLOT PROFILE LOGLIKELIHOOD FOR ALL MODEL PARAMETERS
### Input parameters
###       profliks: the results of getting profile likelihoods with the previous function
###       theseq: points used to sweep the values for each parameter
###       fixpar: whether any of the parameters was fixed (a special case of the full model)
###       ylimZoom: whether to force the range to be 10 units (T), or whether display full Y range of function
###       parnames: labels with parameter names to be used in the plots
###
#################################################################################################
plot_profCIall<-function(profliks,theseq,mypars,fixpar,ylimZoom=T,parnames=c("psi","theta11","theta10","p11","p10","r11")){
  
  for (ii in 1:6){
    if (is.na(fixpar[ii])){ 
      proflik<-profliks[ii,]
      themin<-min(proflik)
      proflik<-proflik-themin
      if (ylimZoom){
        myYlim<-c(-10,0)  # fix axes to this range
      }else{
        myYlim<-c(-max(max(proflik),3.84),0)  # adjust for each case, but at least show 3.84 units
      }
      plot(theseq,-proflik,type="l",ylab="",xlab="",ylim=myYlim,xaxt="n")
      axis(1,at=c(0,0.5,1))
      mtext(text=parnames[ii],side=1,line=2.5,cex=0.7)
      points(theseq,-proflik,pch=20)
      IDin<-which((proflik)<3.84/2)   
      points(theseq[IDin],-proflik[IDin],pch=20,col="green")
      abline(v=mypars[ii],col="gray")
    }
  }
}


#################################################################################################
### FUNCTION TO COMPUTE THE PROBABILITY OF SPECIES PRESENCE AT A SITE, CONDITIONAL TO THE SURVEY DATA
### Input parameters
###       mypar: estimated parameter values (probability scale)
###       mydata: list with containing the survey data (format as in simulation function output)
### Example of use
###      load("litewi.RData"); mypar<-c(0.57,0.44,0.07,0.88,0.03,0.19); psicond_FPMS(mypar,mydata_litewi)
###
#################################################################################################
psicond_FPMS <- function(mypar, mydata) { 
  
  psi<-mypar[1];theta11<-mypar[2];theta10<-mypar[3];p11<-mypar[4];p10<-mypar[5];r11<-mypar[6]
  L<-mydata$L;D<-mydata$D;K<-mydata$K;Du<-mydata$Du;
  #K<-Du<-rep(0,S) #set to zero if only want to consider eDNA survey data
  
  #probability calculations (as in function computing the likelihood)
  ND <- L-D   # number of non detections (ambiguous method)
  tmpZ2_1 <- p11^D * (1-p11)^ND
  tmpZ2_0 <- p10^D * (1-p10)^ND
  tmp <- theta11*tmpZ2_1 + (1-theta11)*tmpZ2_0
  tmpZ1_1 <- apply(tmp, 1, prod)
  tmp <- theta10*tmpZ2_1 + (1-theta10)*tmpZ2_0
  tmpZ1_0 <- apply(tmp, 1, prod)
  NDu <- K-Du #number of non detections (unambiguous method)
  tmpZ1_1u <- r11^Du * (1-r11)^NDu
  tmpZ1_0u <- 0^Du * 1^NDu
  
  #numerator and denominator for the calculation of the conditional probability of presence
  mynum<-psi * tmpZ1_1 * tmpZ1_1u   # Pr(Z1=1 & survey data)
  myden<-psi * tmpZ1_1 * tmpZ1_1u + (1-psi) * tmpZ1_0 * tmpZ1_0u    # Pr(survey data)
  
  #conditional probability of presence
  psicond<-mynum/myden  # Pr(Z1=1|survey data)
  print(round(psicond,3))  #print the probability for each of the sampled sites
}

