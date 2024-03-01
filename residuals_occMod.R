
library(mgcv)

residuals.occMod<-function(object,is.detect.constant=FALSE)
{
  det.data<-object$model$data()[["y"]];     # Binary matrix of detections, sites in rows, visits in columns.
  nSites<-dim(object$model$data()[["y"]])[1];           # Number of sites.
  psi<- object$mean$psi        # Occupancy probabilities, assuming a vector.
  
  ## Get number of visits (ni) and probability of non-detection (prob0).
  
  ni<-apply(is.na(det.data)==FALSE,1,sum);    # Number of site visits, a vector across sites.
  xi<-apply(det.data,1,sum,na.rm=TRUE);       # Number of detections.
  xiOcc<-pmin(xi,1);                          # Make it binary to study occupancy "only".
  
  pi<-matrix(object$mean$p,nrow=nSites);  # Detection probabilities.
  
  ##### need to include some code to compute is.detect.constant from pi #####
  
  ## Below assumes equal probabilities of detection for each sampling time!!!
  
  if(is.detect.constant==TRUE)
  {
    pi<-object$mean$p[1:nSites];   # Detection probabilities, assuming equal across sites.
    prob0<-pbinom(0,ni,pi);            # Probability of no detections when present, a site vector.
    
    ## Get cdf's for detection residuals - as a function of sum of detection events
    
    xi[xi==0]<-NA;                               # Ignore sites with no detections here.
    pdet<-(pbinom(xi,ni,pi)-prob0)/(1-prob0);   # CDF for number of detections xi, positive binomial.
    pdetMinus<-(pbinom(xi-1,ni,pi)-prob0)/(1-prob0); # Previous value of the cdf of xi.
  }
  
  if(is.detect.constant==FALSE)
  {
    prob0<-apply(1-pi,1,prod);    # Probability of no detections when present, a site vector.
    
    ## Define a function to get the pdf under unequal detections.
    
    hetpdf<-function(xiSite,niSite,piSite)
    {
      ind<-combn(niSite,xiSite);
      piMat<-matrix(piSite[ind],nrow=xiSite);
      
      return(sum(apply(piMat/(1-piMat),2,prod))*prod(1-piSite));
    }
    
    hetcdf<-function(xiSite,niSite,piSite)
    {
      if(xiSite==0){cdf<-0;}
      else
      {
        detiSite<-rep(NA,xiSite);
        
        for(iX in 1:xiSite)
        {
          detiSite[iX]<-hetpdf(iX,niSite,piSite);
        }
        
        cdf<-sum(detiSite);       
      }
      
      return(cdf);
    }
    
    ## Get cdf's for detection residuals - as a function of sum of detection events.
    
    isDetected<-xi>0;
    xi[isDetected==FALSE]<-NA;  # Ignores sites with no detections.
    
    pdet=pdetMinus=rep(NA,nSites);
    
    for(iSite in which(isDetected))
    {
      xiSite<-xi[iSite];
      niSite<-ni[iSite];
      piSite<-pi[iSite,];
      pdetMinus[iSite]<-hetcdf(xiSite-1,niSite,piSite);
      pdet[iSite]<-pdetMinus[iSite]+hetpdf(xiSite,niSite,piSite);
    }
    
    pdet<-pdet/(1-prob0);   # 'CDF' for number of detections xi in heterogeneous case.
    pdetMinus<-pdetMinus/(1-prob0); # Previous value of the cdf of xi.
  }
  
  ## Get cdf's for occupancy residuals - as a function of binary detected/not.
  
  probOcc<-psi*(1-prob0);                     # Probability of occupancy.
  pOcc<-1-probOcc+xiOcc*probOcc;              # CDF for occupancy, Bernoulli variable with param probOcc.
  pOccMinus<-xiOcc*(1-probOcc);               # Previous value of the cdf of occupancy.
  
  ## Jitter and get occupancy residuals.
  
  uOcc<-runif(nSites);                             # Standard uniform value to "jitter" the cdf.
  residOcc<-qnorm(pOcc*uOcc+pOccMinus*(1-uOcc));   # Dunn-Smyth residual, standard normal if cdf correct.
  
  ## Jitter and get detection residuals.
  
  u<-runif(nSites);                             # Standard uniform value to "jitter" the cdf.
  residDet<-qnorm(pdet*u+pdetMinus*(1-u));      # Dunn-Smyth residual, standard normal if cdf correct.
  
  residuals<-list(occ=residOcc,det=residDet);  
  
  return(residuals) # Return output (i.e., a list with occupancy residuals (occ) and detection residuals (det)).
}


DS.resid.plot<-function(x,y,ylim=c(-1,1)*max(abs(y)),alpha=0.05,k=5)
{
  plot(x,y,pch=16,cex=1.2,col="blue",cex.axis=1.4,ylim=ylim,cex.main=0.9,ylab="",xlab="");
  
  lsmod<-gam(y~s(x,k=k));
  lsmod.p<-predict(lsmod,se.fit=TRUE);
  z.crit<-qnorm(1-alpha/2)
  upr<-lsmod.p$fit+(z.crit*lsmod.p$se.fit);
  lwr<-lsmod.p$fit-(z.crit*lsmod.p$se.fit);
  polygon(c(rev(sort(x)),sort(x)),c(rev(upr[order(x)]),(lwr[order(x)])),col='grey80',border=NA);
  points(x,y,pch=16,cex=1.2,col="blue");
  lines(sort(x),fitted(lsmod)[order(x)],lwd=2);
  abline(h=0,col="black",lwd=1);
}
