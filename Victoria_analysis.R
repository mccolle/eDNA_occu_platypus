
# load data ---------------------------------------------------------------
load("melb_vic_platy_data.RData")
covs_data<-melb_vic_platy_data

#data array
load("sampled_array_vic.RData")

# standardize -------------------------------------------------------------
lu1 = melb_vic_platy_data$LU_1_result1
lu1<-ifelse(lu1>0,1,0)

runoff=standardize(melb_vic_platy_data$runoff_catch1)
z_flow=standardize(melb_vic_platy_data$AnnlZrF)

# model -------------------------------------------------------------------

site<-seq(289)
str( vic_output <- list(y = sampled_array_vic,site=site,n.site = 289, n.samples = 2, n.pcr = 3, 
                               LU1=lu1,
                               runoff=runoff,
                               z_flow=z_flow))

sink("vic_output.txt")
cat("
    model {
    
    # Priors and model for params
    int.psi ~ dunif(0,1)         # Intercept of occupancy probability
        
    int.theta ~ dunif(0,1) # Intercepts availability probability
   
    int.p ~ dunif(0,1)     # Intercepts detection probability (1-PCR error)
    
    
    
    sigma_site~dunif(0,1)
    
    # convert it to a precision (1 / variance)
    
    tau_site<-pow(sigma_site, -2)
    
    # define the random effect; one value for each site, all share the same precision
    
    for (i in 1:n.site){
    gamma_site[i]~dnorm(0,tau_site)
    }
    
    
    beta_runoff~dnorm(0, 0.1)
    
    beta_flow ~ dnorm(0, 0.1)
  
    beta_LU1 ~ dnorm(0, 0.1)
    

    # 'Likelihood' (or basic model structure)
    for (i in 1:n.site){
    # Occurrence in pond i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi)+
     beta_LU1*LU1[i]+

    beta_runoff * runoff[i] + 
    beta_flow*z_flow[i]+
    gamma_site[site[i]]
    
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

zst <- apply(sampled_array_vic, 1, max)        # inits for presence (z)
zst<-ifelse(is.na(zst),1,zst)
ast <- apply(sampled_array_vic, c(1,2), max)   # inits for availability (a)
ast<-ifelse(is.na(ast),1,ast)
inits <- function() list(z = zst, a = ast,int.psi = 0.5)

# Parameters monitored
params <- c("int.psi", "int.theta", "int.p",  "sigma_site", 
            "p","theta","psi",
            "beta_LU1","beta_runoff",
            "beta_flow")

ni <- 350000  ;   nt <- 1000  ;   nb <- 50000  ;   nc <- 3 

out_vic_output<- jags(vic_output, inits, params, "vic_output.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
