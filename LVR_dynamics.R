## Adapted from: Nguyen, V. K., & Hernandez-Vargas, E. A. (2018). 
## Parameter Estimation in Mathematical Models of Viral Infections Using R. 
## Methods in molecular biology (Clifton, N.J.), 1836, 531-549. 
## https://doi.org/10.1007/978-1-4939-8678-1_25

library("deSolve")
library("DEoptim")

## filter by subject
# set working directory to directory with data
setwd("~/Desktop/simulation/param_estimation")

## read csv file as myData
df1 <-  read.csv("D21_random_effect_QVOA_subj.csv")
head(df1)
str(df1)


df1 <- df1[,c(1,2,5)]
head(df2)
str(df2)

df3 <- df2[!is.na(df2[5]),]
head(df3)
str(df3)


## create tab
final_params_all = c("beta"=0, "alpha"=0, "mu"=0)
final_aic_all <- c(0,0)

## Loop through LVR measurements for each individual
for (i in unique(df3$subject)) {
  x <-data.frame(df3[df3$subject==i,][1:dim(df3[df3$subject==i,])[1],])
  if (dim(x)[1]>=2) {
    a=1
    while (a<dim(x)[1]) {
      reald <- (x[a:(a+1),])
      
      file_name <- paste(i,a, sep="_")
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/All/Log_all/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      plot(reald$time, log(reald$L), type="b", main=paste("Subject", i, sep=" "), xlab="Time", ylab="log L size")
      dev.off()
      
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/All/L_all/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      plot(reald$time, reald$L, type="b", main=paste("Subject", i, sep=" "), xlab="Time", ylab="L size")
      dev.off()
      
      # writing the model and the cost function
      myModel <- function(t,state, parameters) {
        with(as.list(c(state, parameters)), {
          dL = beta*L - alpha*L - mu*L
          list(dL)
        })
      }
      
      myStates <- c(L=reald[1,3])
      myParams <- c('beta','alpha','mu')
      myL <- 'L' # component in the model observed in data
      modelTime <- seq(from=0, to=10, by=1)
      
      # cost function
      ## the root mean square errors (RMSE)
      ### measures magnitude of difference between model output and experimental data
      
      myCostFn <- function(x) {
        parms <- x[1:length(myParams)]
        names(parms) <- myParams
        yhat <- ode(myStates, modelTime, myModel, 10^parms)
        yMatch <- yhat[as.character(yhat[,1]) %in% as.character(reald$time), ]
        nm <- rle(reald$time)$lengths
        x <- reald[,myL] - rep(yMatch[,myL], times = nm)
        rmse <- sqrt(mean(x^2))
        return(rmse)
      }
      
      
      # Defining Parameter boundaries, optimizer conditions 
      ## parameters beta: clonal expansion rate, alpha: activation rate, mu: death rate
      lower = log(c(0.1,0.1,0.01))
      upper = log(c(2,2,0.5))
      
      ## increasing number of trials (itermax and septol) and decreasing thr relative tolerance (reltol: measurement of the error relative to the size of each solution component)
      ### this forces the optimizer to work more exhaustively
      myOptions <- DEoptim.control(itermax = 10000, steptol = 100, reltol = 1e-8)
      
      # fit the model
      fit <- do.call("DEoptim", list(myCostFn, lower, upper, myOptions))
      
      # Visualization
      (bestPar <- fit$optim$bestmem)
      names(bestPar) <- myParams
      bestPar
      
      nam_par <- paste("bestPar",i, a, sep="_")
      assign(nam_par, bestPar)
      
      out <- ode(myStates, modelTime, myModel, 10^bestPar)
      
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/sim/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      plot(out[,"time"], out[,"L"], type="l", main=paste("Subject", i, sep=" "), xlab="Time", ylab="L size")
      points(reald$time, reald$L)
      dev.off()
      
      # Model comparison with Akaike Information Criteria (AIC)
      ## AIC gives a penalty to the number of parametes to avoid overfitting
      ### The smaller the AIC, the better the model
      
      myAIC <- function(fit, np=NULL, rms=NULL, n=NULL) {
        if (is.null(n)) stop("How many observations were used? n=#")
        if (is.null(np)) np <- length(fit$optim$bestmem)
        if (is.null(rms)) rms <- fit$optim$bestval
        return(2*np + n*log(rms))
      }
      ## n-sample size used for fitting
      aic <- myAIC(fit, n = dim(reald[1]))
      nam_aic <- paste("aic",i,a, sep="_")
      assign(nam_aic, aic)
      
      # Likelihood profile of model parameters
      myProfile <- function(lower, upper, bestPar) {
        pro.ll <- NULL
        for (v in 1:length(bestPar)) {
          # creating parameter sequence
          tmpl <- seq(lower[v], bestPar[[v]], length.out = 100)
          tmpl <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]]
          tmpr <- seq(bestPar[[v]], upper[v], length.out = 100)
          tmpr <- tmpr[cumsum(1:13)]
          pars <- sort(unique(c(lower[v], tmpl, bestPar[[v]], tmpr, upper[v])))
          ppl <- NULL
          # Run optimization for each and record the parameters and RMSE
          for (p in pars) {
            DEargs <- list(myCostFn, replace(lower,v, p), replace(upper, v, p), myOptions)
            fit <- do.call("DEoptim", DEargs)
            ppl <- c(ppl, fit$optim$bestval)
          }
          pro.ll[[v]] <- cbind(pars, ppl)
        }
        return(pro.ll)
      }
      
      outProfiles <- myProfile(lower, upper, bestPar)
      
      ## plot profile of the first parameters
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/All/RMSE_all/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      par(mfrow = c(2,2))
      sapply(1:3, function(x) plot(outProfiles[[x]], xlab=myParams[x], ylab = 'RMSE',type="b", main=paste("Subject", i, sep=" ")))
      dev.off()
      
      ## Bootstrapping parameters
      myBoot <- function(numboot = 10, numpar = 3) {
        results <- matrix(NA, numboot, numpar)
        original <- reald
        sampling <- function(x) sample(original$L[original$time==x],
                                       length(original$L[original$time==x]),
                                       replace = 1)
        for (i in 1:numboot) {
          message("Bootstraping sample", i)
          tmp <- sapply(unique(original$time), sampling)
          reald <- cbind(original$time, as.vector(tmp))
          DEarguments <- list(myCostFn, lower, upper, myOptions)
          fit <- do.call("DEoptim", DEarguments)
          results[i,] <- fit$optim$bestmem
        }
        results <- as.data.frame(results)
        colnames(results) <- myParams
        return(results)
      }
      
      bootResults <- myBoot()
      
      ## save bootsrap histograms
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/All/boot_all/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      par(mfrow = c(2,2))
      sapply(1:3, function(x) hist(bootResults[,x], main = paste("Subject", i, sep=" ")))
      dev.off()
      
      ## save bootstrap line graphs
      file_path <- paste("~/Desktop/simulation/param_estimation/plots/All/boot_line_all/",file_name,".jpeg", sep="")
      jpeg(filename = file_path)
      par(mfrow = c(2,2))
      sapply(1:3, function(x) plot(bootResults[,x],  type="l", main = paste("Subject", i, sep=" ")))
      dev.off()
      
      
      ## calculate 95% confidence intervals from bootstrap sample with percentile method
      apply(bootResults, 2, quantile, probs = c(0.25,0.975))
      
      par(mfrow = c(1,1))
      
      final_params_all <- rbind(final_params_all, get(nam_par)) 
      final_aic_all <- rbind(final_aic_all ,get(nam_aic)) 
      
      a = a + 1
      
    } 
  }
}
