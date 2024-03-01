
# load data ---------------------------------------------------------------

load("GAPS_detection_array.RData")
load("catchment_fact.RData")
load("covs_data.RData")


# packages ----------------------------------------------------------------
library(rjags)
library(jagsUI)
library(AHMbook)


# standardise -------------------------------------------------------------

hd = as.numeric(covs_data$human_dist) #human disturbance
erosion = (as.numeric(covs_data$erosion)) #bank erosion
b_veg = (as.numeric(covs_data$b_veg)) #bank vegetation
v_veg =  (as.numeric(covs_data$v_veg)) #verge vegetation
b_bur =  (as.numeric(covs_data$b_bur)) #banks for burrowing 
s_complx = (as.numeric (covs_data$stream_complx)) #stream complexity
c_complx =  (as.numeric(covs_data$chan_complx)) #Chanel complexity
lu3 = standardize(covs_data$LU_3_result1) # Agriculture
lu4 = standardize(covs_data$LU_4_result1) # Shrubs and grasslands
lu5 = standardize(covs_data$LU_5_result1) # Forest
ro = standardize(covs_data$runoff_catch1) # runoff 



# model -------------------------------------------------------------------


str(whole_state_output <- list(y = GAPS_detection_array,site=1:dim(GAPS_detection_array)[1],
                          n.site = dim(GAPS_detection_array)[1],
                          n.samples = dim(GAPS_detection_array)[2], 
                          n.pcr = dim(GAPS_detection_array)[3],
                          catchment_fact=catchment_fact,
                          n.catch=length(unique(catchment_fact)),
                          hd=hd,
                          erosion=erosion,
                          b_veg=b_veg,
                          b_bur=b_bur,
                          lu3=lu3,
                          lu4=lu4,
                          lu5=lu5,
                          runoff=ro
                          
))


sink("whole_state_output.txt")
cat("
    model {
    
    # Priors and model for params
    int.psi ~ dunif(0,1)         # Intercept of occupancy probability
    int.theta ~ dunif(0,1) # Intercepts availability probability
    int.p ~ dunif(0,1)     # Intercepts detection probability (1-PCR error)
    
    sigma_site~dunif(0,1)
    sigma_catch~dunif(0,1)
    
    # convert it to a precision (1 / variance)
    
    tau_site<-pow(sigma_site, -2)
    tau_catch<-pow(sigma_site, -2)
    
    # define the random effect; one value for each site, all share the same precision
    
    for (i in 1:n.site){
    gamma_site[i]~dnorm(0,tau_site)
    }
    for (i in 1:n.catch){
    gamma_catch[i]~dnorm(0,tau_catch)
    }
    
    
    beta_hd~ dnorm(0, 0.1)  
    
    beta_erosion[1] <- 0
    for (k in 2:5) {
    beta_erosion[k] ~ dnorm(0, 0.1)
    }
    
    beta_b_veg[1] <- 0
    for (k in 2:5) {
    beta_b_veg[k] ~ dnorm(0, 0.1)
    }
    
   
    
    beta_b_bur[1] <- 0
    for (k in 2:5) {
    beta_b_bur[k] ~ dnorm(0, 0.1)
    }
    
 
    
 
    
    beta_lu3 ~ dnorm(0, 0.1)
    
    beta_lu4 ~ dnorm(0, 0.1)
    
    beta_lu5 ~ dnorm(0, 0.1)
    
    beta_ro ~ dnorm(0, 0.1)
    
    
    #quad
    
    
  
    beta_lu4_quad~ dnorm(0, 0.1)
    beta_lu5_quad~ dnorm(0, 0.1)
    beta_runoff_quad~ dnorm(0, 0.1)
    
    
    
    # 'Likelihood' (or basic model structure)
    for (i in 1:n.site){
    # Occurrence in pond i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi)+
  
    beta_lu3 * lu3[i] +  
    beta_lu4 * lu4[i] +   
    beta_lu4 * lu4[i] +   
    beta_hd*hd[i]+
    beta_erosion[erosion[i]]+ 
    beta_b_veg[b_veg[i]]+
    beta_b_bur[b_bur[i]]+
    beta_ro*runoff[i]+
    
    beta_lu4_quad*pow(lu4[i],2) +
    beta_lu5_quad*pow(lu5[i],2) +
    beta_runoff_quad*pow(runoff[i],2) +
    
    
    gamma_site[site[i]]+ gamma_catch[catchment_fact[i]]
    
    
    for (j in 1:n.samples){
    # Occurrence in sample j
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- logit(int.theta) 
    
    for (k in 1:n.pcr){
    # PCR detection error process in sample k
    
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- logit(int.p) 
    }
    }
    }
    }# end model
    
    ",fill=TRUE)
sink()

zst <- apply(GAPS_detection_array, 1, max)        # inits for presence (z)
zst<-ifelse(is.na(zst),1,zst)
ast <- apply(GAPS_detection_array, c(1,2), max)   # inits for availability (a)
ast<-ifelse(is.na(ast),1,ast)
inits <- function() list(z = zst, a = ast,int.psi = 0.5)

params <- c("int.psi", "int.theta", "int.p",  
            "p", "theta","psi",
            "beta_hd","beta_erosion","beta_b_veg","beta_v_veg","beta_b_bur",
            "beta_s_complx","beta_c_complx",
            "beta_ro","beta_lu1","beta_lu3","beta_lu4","beta_lu5",
            "beta_runoff_quad","beta_lu3_quad","beta_lu4_quad","beta_lu5_quad")

ni <- 350000  ;   nt <- 1000  ;   nb <- 50000  ;   nc <- 3 
 

# Call WinBUGS and summarize posterior
out_whole_state_output<- jags(whole_state_output, inits, params, "whole_state_output.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)



# Figure 2 --------------------------------------------------------

 #human disturbance


all1<-as.matrix(out_whole_state_output$sims.list$beta_hd)
pm1 <- apply(all1, 2, mean)    # Get posterior means and 95% CRIs
cri1 <- apply(all1, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs

covs<-as.data.frame(cbind(pm1,cri1[1],cri1[2]))
covs$position<-"1"
covs$name<-"Human disturbance"
covs$overlap<-1
colnames(covs)<-c("mean", "cri_l", "cri_u","position", "Name","overlap")




# land use 3 
all2<-as.matrix(out_whole_state_output$sims.list$beta_lu3)
pm2 <- apply(all2, 2, mean)   # Get posterior means and 95% CRIs
cri2 <- apply(all2, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[2,]<-c(pm2,cri2[1],cri2[2],"2", "Agriculture",0)




# land use 4
all3<-as.matrix(out_whole_state_output$sims.list$beta_lu4)
pm3 <- apply(all3, 2, mean)   # Get posterior means and 95% CRIs
cri3 <- apply(all3, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[3,]<-c(pm3,cri3[1],cri3[2],"3", "Shrubs and grasslands",0)

# land use 4 quad
all4<-as.matrix(out_whole_state_output$sims.list$beta_lu4_quad)
pm4 <- apply(all4, 2, mean)   # Get posterior means and 95% CRIs
cri4 <- apply(all4, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[4,]<-c(pm4,cri4[1],cri4[2],"4", "Shrubs and grasslands quad",0)

# land use 5
all5<-as.matrix(out_whole_state_output$sims.list$beta_lu5)
pm5 <- apply(all5, 2, mean)   # Get posterior means and 95% CRIs
cri5 <- apply(all5, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[5,]<-c(pm5,cri5[1],cri5[2],"5", "Forest",1)

# land use 5 quad
all6<-as.matrix(out_whole_state_output$sims.list$beta_lu5_quad)
pm6 <- apply(all6, 2, mean)   # Get posterior means and 95% CRIs
cri6 <- apply(all6, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[6,]<-c(pm6,cri6[1],cri6[2],"6", "Forest quad",1)

# runoff 
all7<-as.matrix(out_whole_state_output$sims.list$beta_ro)
pm7 <- apply(all7, 2, mean)   # Get posterior means and 95% CRIs
cri7 <- apply(all7, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[7,]<-c(pm7,cri7[1],cri7[2],"7", "Runoff",0)


# runoff quad
all8<-as.matrix(out_whole_state_output$sims.list$beta_runoff_quad)
pm8 <- apply(all8, 2, mean)   # Get posterior means and 95% CRIs
cri8 <- apply(all8, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[8,]<-c(pm8,cri8[1],cri8[2],"8", "Runoff quad",1)


# erosion
all19<-as.matrix(out_whole_state_output$sims.list$beta_erosion[,1])
covs[9,]<-c(0,0,0,"9", "Erosion (very poor)",1)

all10<-as.matrix(out_whole_state_output$sims.list$beta_erosion[,2])
pm10 <- apply(all10, 2, mean)   # Get posterior means and 95% CRIs
cri10 <- apply(all10, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[10,]<-c(pm10,cri10[1],cri10[2],"10", "Erosion (poor)",1)

all11<-as.matrix(out_whole_state_output$sims.list$beta_erosion[,3])
pm11 <- apply(all11, 2, mean)   # Get posterior means and 95% CRIs
cri11 <- apply(all11, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[11,]<-c(pm11,cri11[1],cri11[2],"11", "Erosion (fair)",1)

all12<-as.matrix(out_whole_state_output$sims.list$beta_erosion[,4])
pm12 <- apply(all12, 2, mean)   # Get posterior means and 95% CRIs
cri12 <- apply(all12, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[12,]<-c(pm12,cri12[1],cri12[2],"12", "Erosion (good)",1)

all13<-as.matrix(out_whole_state_output$sims.list$beta_erosion[,5])
pm13 <- apply(all13, 2, mean)   # Get posterior means and 95% CRIs
cri13 <- apply(all13, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[13,]<-c(pm13,cri13[1],cri13[2],"13", "Erosion (excellent)",1)


# bank veg 
all14<-as.matrix(out_whole_state_output$sims.list$beta_b_veg[,1])
covs[14,]<-c(0,0,0,"14", "Bank vegetation (very poor)",1)

all15<-as.matrix(out_whole_state_output$sims.list$beta_b_veg[,2])
pm15 <- apply(all15, 2, mean)   # Get posterior means and 95% CRIs
cri15 <- apply(all15, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[15,]<-c(pm15,cri15[1],cri15[2],"15", "Bank vegetation (poor)",1)

all16<-as.matrix(out_whole_state_output$sims.list$beta_b_veg[,3])
pm16 <- apply(all16, 2, mean)   # Get posterior means and 95% CRIs
cri16 <- apply(all16, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[16,]<-c(pm16,cri16[1],cri16[2],"16", "Bank vegetation (fair)",1)

all17<-as.matrix(out_whole_state_output$sims.list$beta_b_veg[,4])
pm17 <- apply(all17, 2, mean)   # Get posterior means and 95% CRIs
cri17 <- apply(all17, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[17,]<-c(pm17,cri17[1],cri17[2],"17", "Bank vegetation (good)",1)

all18<-as.matrix(out_whole_state_output$sims.list$beta_b_veg[,5])
pm18 <- apply(all18, 2, mean)   # Get posterior means and 95% CRIs
cri18 <- apply(all18, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[18,]<-c(pm18,cri18[1],cri18[2],"18", "Bank vegetation (excellent)",1)


# burrowing banks 
all19<-as.matrix(out_whole_state_output$sims.list$beta_b_bur[,1])
covs[19,]<-c(0,0,0,"19", "Burrowing banks (very poor)",1)

all20<-as.matrix(out_whole_state_output$sims.list$beta_b_bur[,2])
pm20 <- apply(all20, 2, mean)   # Get posterior means and 95% CRIs
cri20 <- apply(all20, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[20,]<-c(pm20,cri20[1],cri20[2],"20", "Burrowing banks (poor)",1)

all21<-as.matrix(out_whole_state_output$sims.list$beta_b_bur[,3])
pm21 <- apply(all21, 2, mean)   # Get posterior means and 95% CRIs
cri21 <- apply(all21, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[21,]<-c(pm21,cri21[1],cri21[2],"21", "Burrowing banks (fair)",0)

all22<-as.matrix(out_whole_state_output$sims.list$beta_b_bur[,4])
pm22 <- apply(all22, 2, mean)   # Get posterior means and 95% CRIs
cri22 <- apply(all22, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[22,]<-c(pm22,cri22[1],cri22[2],"22", "Burrowing banks (good)",0)

all23<-as.matrix(out_whole_state_output$sims.list$beta_b_bur[,5])
pm23 <- apply(all23, 2, mean)   # Get posterior means and 95% CRIs
cri23 <- apply(all23, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
covs[23,]<-c(pm23,cri23[1],cri23[2],"23", "Burrowing banks (excellent)",1)



# plot 
covs$position<-as.numeric(covs$position)
covs<-covs[order(covs$position),]


covs1<-lapply(covs, as.numeric)
covs1$Name<-covs$Name
covs1<-as.data.frame(covs1)
covs1$position<-as.factor(covs1$position)

labels<-covs$Name
covs1$overlap<-as.factor(covs1$overlap)


plot1<-ggplot(covs1,aes(x=mean, y=position,xmin=cri_l, xmax=cri_u))+
  geom_rect(xmin=-13, xmax=13, ymin=0, ymax=8.5, fill="#f0f0f0", alpha=0.05)+
  geom_rect(xmin=-13, xmax=13, ymin=8.5, ymax=13.5, fill="#bdbdbd", alpha=0.05)+
  geom_rect(xmin=-13, xmax=13, ymin=13.5, ymax=18.5, fill="#969696", alpha=0.05)+
  geom_rect(xmin=-13, xmax=13, ymin=18.5, ymax=23.5, fill="#737373", alpha=0.05)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  geom_point(aes(colour=overlap),size=2)+
  geom_errorbarh(height=.1, aes(colour=overlap),size=1)+
  scale_colour_manual(values = c("#1b9e77", "#d95f02"))+
  xlim(-7.5, 7.5)+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14,angle = 15)) +
  
  geom_vline(xintercept = 0, 
             color = "gray65", size=1)+
  xlab("Parameter estimate")+
  ylab(element_blank())+
  scale_y_discrete(breaks=c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                            "14","15","16","17","18","19","20","21","22","23"),
                   labels=labels)


# Figure 3 --------------------------------------------

ro = covs_data$runoff_catch1
o.ro<- seq(min(ro), max(ro),,500)
ro.pred <- standardize(o.ro)

lu4<- covs_data$LU_4_result1
o.lu4<- seq(min(lu4), max(lu4),,500)
lu4.pred <-standardize(o.lu4)

ag<-covs_data$LU_3_result1
o.ag<- seq(min(ag), max(ag),,500)
ag.pred <- standardize(o.ag)


str( tmp <- out_whole_state_output$sims.list )              # grab MCMC samples
nsamp <- length(tmp[[1]])    # number of mcmc samples
predC <- array(NA, dim = c(500, nsamp, 3))
for(i in 1:nsamp){
  
  predC[,i,1] <- plogis(tmp$int.psi[i] + tmp$beta_ro[i] * ro.pred+
                          tmp$beta_runoff_quad[i]*ro.pred^2)
  
  predC[,i,2] <- plogis(tmp$int.psi[i] + tmp$beta_lu4[i] * lu4.pred+
                          tmp$beta_lu4_quad[i]*lu4.pred^2 )
  
  predC[,i,3] <- plogis(tmp$int.psi[i] + tmp$beta_lu3[i] * ag.pred)
  
}


pmC <- apply(predC, c(1,3), mean)
criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

plot(o.ro, pmC[,1], col = "#fc8d59", lwd = 3, type = 'l', lty = 1, frame = F, 
     ylim = c(0, 1),
     xlab = "Runoff", ylab = "Mean occupancy",
     cex.lab=1.5)
rug(ro)
matlines(o.ro, t(criC[,,1]), col = "grey", lty = 1)
mtext("A", side=3,adj = 0)


plot(o.lu4, pmC[,2], col = "#fc8d59", lwd = 3, type = 'l', lty = 1, frame = F, 
     ylim = c(0, 1),
     xlab = "Shrubs and grasslands", ylab = "Mean occupancy",
     cex.lab=1.5)
rug(lu4)
matlines(o.lu4, t(criC[,,2]), col = "grey", lty = 1)
mtext("B", side=3,adj = 0)

plot(o.ag, pmC[,3], col = "#fc8d59", lwd = 3, type = 'l', lty = 1, frame = F, 
     ylim = c(0, 1),
     xlab = "Agriculture", ylab = "Mean occupancy",
     cex.lab=1.5)
rug(ag)
matlines(o.ag, t(criC[,,3]), col = "grey", lty = 1)
mtext("C", side=3,adj = 0)

boxplot(plogis(tmp$int.psi + tmp$beta_b_bur ), xlab="Burrowing banks",
        ylim=c(0,1), col="#fc8d59",
        names=c("V.poor","Poor","Fair","Good", "Excellent"),
        cex.lab=1.5)
mtext("D", side=3,adj = 0)




# example residual plot ---------------------------------------------------
source("residuals_occMod.R")
res1<-residuals.occMod(out_whole_state_output) 

m1_res_occ<-res1$occ # Stores residuals for occupancy.
m1_res_det<-res1$det[-which(is.na(res1$det))] # Stores residuals for detection (and remove NAvalues)

DS.resid.plot( standardize(covs_data$runoff_catch1)[-which(is.na(res1$det))],m1_res_det)
title(main="covariate occupancy Dunn-Smyth residuals",cex.main=2)
mtext(text="Runoff",side=1,las=1,line=2.7,cex=1.2)


