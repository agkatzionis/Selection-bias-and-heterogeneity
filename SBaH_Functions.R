
##########   TYPE 2 SELECTION BIAS IN MR   ##########

## This file contains R functions to implement the 
## simulations for the project "Selection Bias Due
## to Heterogeneity in IV Analyses". The simulations
## themselves are run elsewhere.

## Set this up.
library(ivreg)
library(OneSampleMR)
library(estimatr)
library(survival)
library(ivtools)
library(sandwich)

## Note: "dag" takes the values:
## 1: V modifies the X-Y effect.
## 2: U modifies the X-Y effect.
## 3: V modifies the G-X effect.
## 4: U modifies the G-X effect.

## Simulation generator.
Sim_gen <- function (iter = 1000, n = 1e5, y.distr = "normal", dag = 1,
                     em.strength, beta.sel, theta, seed = NULL, ipw = FALSE) {
  
  ## Set seed, if needed.
  if(!is.null(seed)) set.seed(seed)
  
  ## Store results here.
  results <- list("diagnostics" = matrix(0, iter, 4), "oracle" = 0, "cca" = 0, "int.oracle" = 0, "int.cca" = 0)
  colnames(results$diagnostics) <- c("y.mean", "s.prob", "f.stat", "r2")
  if (y.distr == "normal" | y.distr == "cox") {
    results$oracle <- matrix(0, iter, 5)
    colnames(results$oracle) <- c("est", "se", "ci.lower", "ci.upper", "pval")
  } else {
    results$oracle <- matrix(0, iter, 4)
    colnames(results$oracle) <- c("est", "se", "ci.lower", "ci.upper")
  }
  results$cca <- results$oracle
  results$int.oracle <- results$oracle
  results$int.cca <- results$oracle
  
  ## Store IPW results here.
  if (ipw == TRUE) {
    results$ipw <- results$oracle
    results$int.ipw <- results$oracle
  }
  
  ## Start looping.
  for (I in 1:iter) {
    
    ## ----- Generate data. ----- ##
    
    ## Instrument.
    G <- rnorm(n, 0, 1)
    
    ## Confounders and effect modifiers.
    U <- rnorm(n, 0, 1)
    if (dag == 1 | dag == 3 | dag == 4) V <- rnorm(n, 0, 1)
    
    ## Exposure.
    if (dag == 3) {
      X <- 0.4 * G + U + V + em.strength * G * V + rnorm(n, 0, 1)
    } else if (dag == 4) {
      X <- 0.4 * G + U + em.strength * G * U + rnorm(n, 0, 1)
    } else {
      X <- 0.4 * G + U + rnorm(n, 0, 1)
    }
    fit.gx <- summary(lm(X ~ G))
    results$diagnostics[I, 3] <- fit.gx$f[1]
    results$diagnostics[I, 4] <- fit.gx$r.sq
    
    ## Outcome.
    if (y.distr == "normal") {
      if (dag == 1) {
        Y <- theta * X + U + V + em.strength * X * V + rnorm(n, 0, 1)
      } else if (dag == 2) {
        Y <- theta * X + U + em.strength * X * U + rnorm(n, 0, 1)
      } else {
        Y <- theta * X + U + rnorm(n, 0, 1)
      }
      results$diagnostics[I, 1] <- mean(Y)
    } else if (y.distr == "binary") {
      if (dag == 1) {
        Y.prob <- plogis(theta * X + U + V + em.strength * X * V)
      } else if (dag == 2) {
        Y.prob <- plogis(theta * X + U + em.strength * X * U)
      } else {
        Y.prob <- plogis(theta * X + U)
      }
      Y <- rbinom(n, 1, Y.prob)
      results$diagnostics[I, 1] <- mean(Y)
    } else if (y.distr == "cox") {
      if (dag == 1) {
        lp <- theta * X + U + V + em.strength * X * V
      } else if (dag == 2) {
        lp <- theta * X + U + em.strength * X * U
      } else {
        lp <- theta * X + U
      }
      TTE <- -log(runif(n)) / exp(lp)   ## Survival time for Exp(1) baseline.
      E <- rep(1, n)   ## Event indicator (all events, no censoring).
      Y <- Surv(TTE, E)
      results$diagnostics[I, 1] <- mean(Y[, "time"])
    }
    
    ## Selection coefficient.
    if (dag == 2 | dag == 4) {
      S.prob <- plogis(-0.5 + beta.sel * U)
    } else {
      S.prob <- plogis(-0.5 + beta.sel * V)
    }
    S <- rbinom(n, 1, S.prob)
    results$diagnostics[I, 2] <- mean(S)
    
    ## ----- Fit models. ----- ##
    
    ## Oracle fit.
    if (y.distr == "normal") {
      oracle <- summary(ivreg(Y ~ X | G))
      results$oracle[I, c(1, 2, 5)] <- oracle$coefficients[2, c(1, 2, 4)]
      results$oracle[I, c(3, 4)] <- results$oracle[I, 1] + qnorm(c(0.025, 0.975)) * results$oracle[I, 2]
      results$int.oracle[I, c(1, 2, 5)] <- oracle$coefficients[1, c(1, 2, 4)]
      results$int.oracle[I, c(3, 4)] <- results$int.oracle[I, 1] + qnorm(c(0.025, 0.975)) * results$int.oracle[I, 2]
    } else if (y.distr == "binary") {
      oracle <- tsri(Y ~ X | G, link = "logit")
      results$oracle[I, c(1, 3, 4)] <- oracle$estci[4, ]
      results$oracle[I, 2] <- (oracle$estci[4, 1] - oracle$estci[4, 2]) / qnorm(0.975)
      results$int.oracle[I, c(1, 3, 4)] <- oracle$estci[3, ]
      results$int.oracle[I, 2] <- (oracle$estci[3, 1] - oracle$estci[3, 2]) / qnorm(0.975)
    } else if (y.distr == "cox") {
      cox_df <- data.frame("G" = G, "X" = X, "TTE" = TTE, "E" = E, "Y" = Y)
      fitX.LZ <- glm(X ~ G, data = cox_df)
      fitT.LX <- coxph(Surv(TTE, E) ~ X, data = cox_df)
      oracle <- ivcoxph(estmethod = "ts", fitX.LZ = fitX.LZ, fitT.LX = fitT.LX, data = cox_df, ctrl = TRUE)
      results$oracle[I, c(1, 2, 5)] <- summary(oracle)$coefficients[1, c(1, 2, 4)]
      results$oracle[I, c(3, 4)] <- results$oracle[I, 1] + qnorm(c(0.025, 0.975)) * results$oracle[I, 2]
      results$int.oracle[I, 1:5] <- NA
    }
    
    ## Complete-case analysis.
    if (y.distr == "normal") {
      cca <- summary(ivreg(Y ~ X | G, subset = which(S == 1)))
      results$cca[I, c(1, 2, 5)] <- cca$coefficients[2, c(1, 2, 4)]
      results$cca[I, c(3, 4)] <- results$cca[I, 1] + qnorm(c(0.025, 0.975)) * results$cca[I, 2]
      results$int.cca[I, c(1, 2, 5)] <- cca$coefficients[1, c(1, 2, 4)]
      results$int.cca[I, c(3, 4)] <- results$int.cca[I, 1] + qnorm(c(0.025, 0.975)) * results$int.cca[I, 2]
    } else if (y.distr == "binary") {
      cca <- tsri(Y ~ X | G, link = "logit", subset = which(S == 1))
      results$cca[I, c(1, 3, 4)] <- cca$estci[4, ]
      results$cca[I, 2] <- (cca$estci[4, 1] - cca$estci[4, 2]) / qnorm(0.975)
      results$int.cca[I, c(1, 3, 4)] <- cca$estci[3, ]
      results$int.cca[I, 2] <- (cca$estci[3, 1] - cca$estci[3, 2]) / qnorm(0.975)
    } else if (y.distr == "cox") {
      cox_sel <- data.frame("G" = G[S == 1], "X" = X[S == 1], "TTE" = TTE[S == 1], "E" = E[S == 1], "Y" = Y[S == 1])
      fitX.LZ <- glm(X ~ G, data = cox_sel)
      fitT.LX <- coxph(Surv(TTE, E) ~ X, data = cox_sel)
      cca <- ivcoxph(estmethod = "ts", fitX.LZ = fitX.LZ, fitT.LX = fitT.LX, data = cox_sel, ctrl = TRUE)
      results$cca[I, c(1, 2, 5)] <- summary(cca)$coefficients[1, c(1, 2, 4)]
      results$cca[I, c(3, 4)] <- results$cca[I, 1] + qnorm(c(0.025, 0.975)) * results$cca[I, 2]
      results$int.cca[I, 1:5] <- NA
    }
    
    ## ----- Perform IPW. ----- ##
    
    ## This is done only for a normal or binary outcome.
    if (ipw == TRUE) {
      
      ## Build the weighting model.
      if (dag == 2 | dag == 4) {
        ipw_fit <- glm(S ~ U, family = binomial(link = "logit"))
      } else {
        ipw_fit <- glm(S ~ V, family = binomial(link = "logit"))
      }
      ipw_weights <- 1 / ipw_fit$fitted.values
      
      if (y.distr == "normal") {
        
        ## Run the re-weighted IV analysis for a normal outcome.
        ipw_results <- summary(iv_robust(Y ~ X | G, subset = which(S == 1), weights = ipw_weights))
        results$ipw[I, 1:5] <- ipw_results$coefficients[2, c(1, 2, 5, 6, 4)]
        results$int.ipw[I, 1:5] <- ipw_results$coefficients[1, c(1, 2, 5, 6, 4)]
        
      } else if (y.distr == "binary") {
        
        ## Run the re-weighted IV analysis for a binary outcome.
        ## OneSampleMR does not support robust SEs, so do it manually.
        ipw_df <- data.frame("G" = G[S == 1], "X" = X[S == 1], "Y" = Y[S == 1], "ipw_weights" = ipw_weights[S == 1])
        ipw_first_stage <- glm(X ~ G, weights = ipw_weights, data = ipw_df)
        ipw_second_stage <- glm(Y ~ X + residuals(ipw_first_stage), family = binomial(link = "logit"), weights = ipw_weights, data = ipw_df)
        
        ## Calculate robust standard errors.
        se_robust <- sqrt(diag(vcovHC(ipw_second_stage, type = "HC1")))
        
        ## Store results.
        results$ipw[I, 1] <- summary(ipw_second_stage)$coef[2, 1]
        results$ipw[I, 2] <- se_robust[2]
        results$ipw[I, 3:4] <- results$ipw[I, 1] + qnorm(c(0.025, 0.975)) * results$ipw[I, 2]
        results$int.ipw[I, 1] <- summary(ipw_second_stage)$coef[1, 1]
        results$int.ipw[I, 2] <- se_robust[1]
        results$int.ipw[I, 3:4] <- results$int.ipw[I, 1] + qnorm(c(0.025, 0.975)) * results$int.ipw[I, 2]
        
      }
    } 
  }
  
  ## Goodbye.
  return(results)
}

## Function to compute empirical coverage.
cover <- function(sim, true.values = NULL, ipw = FALSE) {
  
  ## Initialize.
  ll <- length(sim)
  if (ipw) {
    res <- matrix(NA, ll, 3)
    colnames(res) <- c("Oracle.cover", "Cca.cover", "Ipw.cover")
  } else {
    res <- matrix(NA, ll, 2)
    colnames(res) <- c("Oracle.cover", "Cca.cover")
  }
  
  ## Loop through simulation scenarios.
  for (i in 1:ll) {
    
    ## Coverage for oracle analysis.
    mat <- as.data.frame(sim[[i]]$oracle)
    if (is.null(true.values)) current.value <-  mean(mat$est) else current.value <- true.values[i]
    res[i, 1] <- mean( (mat$ci.lower - current.value) * (mat$ci.upper - current.value) < 0) 
    
    ## Coverage for CCA.
    mat <- as.data.frame(sim[[i]]$cca)
    res[i, 2] <- mean( (mat$ci.lower - current.value) * (mat$ci.upper - current.value) < 0) 

    ## Coverage for IPW, if needed.
    if (ipw == TRUE) {
      mat <- as.data.frame(sim[[i]]$ipw)
      res[i, 3] <- mean( (mat$ci.lower - current.value) * (mat$ci.upper - current.value) < 0)
    }
  }
  
  ## Goodbye.
  return(res)
}

## Create a table for the results of one scenario, reporting
## bias in point estimates and empirical coverage.
analyze_results <- function (scenario, se = FALSE, true.values = NULL, ipw = FALSE) {
  
  ## Load the data.
  filename <- paste("SBaH_sim", scenario, ".RData", sep = "")
  load(filename)
  
  ## Get point estimates.
  if (ipw) {
    est <- t(sapply(all.sims, function(x) cbind(mean(x$oracle[, 1]), mean(x$cca[, 1]), mean(x$ipw[, 1]))))
    colnames(est) <- c("Oracle", "Cca", "Ipw")
  } else {
    est <- t(sapply(all.sims, function(x) cbind(mean(x$oracle[, 1]), mean(x$cca[, 1]))))
    colnames(est) <- c("Oracle", "Cca")
  }
  
  ## Get standard deviations, if necessary.
  if (se) {
    if (ipw) {
      sds <- t(sapply(all.sims, function(x) cbind(mean(x$oracle[, 2]), mean(x$cca[, 2]), mean(x$ipw[, 2]))))
      colnames(sds) <- c("Oracle", "Cca", "Ipw")
    } else {
      sds <- t(sapply(all.sims, function(x) cbind(mean(x$oracle[, 2]), mean(x$cca[, 2]))))
      colnames(sds) <- c("Oracle", "Cca")
    }
  }
  
  ## If true.values not specified, use oracle estimates.
  if (is.null(true.values)) true.values <- est[, 1]
  
  ## Get empirical coverage.
  coverage <- cover(all.sims, true.values = true.values, ipw = ipw)
  if (ipw) {
    colnames(coverage) <- c("Oracle", "Cca", "Ipw")
  } else {
    colnames(coverage) <- c("Oracle", "Cca")
  }
  
  ## Get other diagnostics.
  diagnostics <- t(sapply(all.sims, function(x) colMeans(x$diagnostics)))
  probs <- t(sapply(all.sims, function(x) c(min(x$diagnostics[, 2]), mean(x$diagnostics[, 2]), max(x$diagnostics[, 2]))))
  
  ## Create the table.
  if (se) {
    Table <- matrix(0, 6, 15)
    colnames(Table) <- c("e-05.bias", "e-05.se", "e-05.cov", "e-025.bias", "e-025.se", "e-025.cov", "e0.bias",
                         "e0.se", "e0.cov", "e025.bias", "e025.se", "e025.cov", "e05.bias", "e05.se", "e05.cov")
  } else {
    Table <- matrix(0, 6, 10)
    colnames(Table) <- c("e-05.bias", "e-05.cov", "e-025.bias", "e-025.cov", "e0.bias",
                         "e0.cov", "e025.bias", "e025.cov", "e05.bias", "e05.cov")
  }
  rownames(Table) <- c("th0.b02", "th0.b05", "th0.b1", "th05.b02", "th05.b05", "th05.b1")
  
  ## Fill in the values.
  if (se) {
    if (ipw) {
      for (i in 1:6) {
        Table[i, c(1, 4, 7, 10, 13)] <- est[(5*(i-1) + 1):(5*i), 3] - true.values[(5*(i-1) + 1):(5*i)]
        Table[i, c(2, 5, 8, 11, 14)] <- sds[(5*(i-1) + 1):(5*i), 3]
        Table[i, c(3, 6, 9, 12, 15)] <- coverage[(5*(i-1) + 1):(5*i), 3]
      }  
    } else {
      for (i in 1:6) {
        Table[i, c(1, 4, 7, 10, 13)] <- est[(5*(i-1) + 1):(5*i), 2] - true.values[(5*(i-1) + 1):(5*i)]
        Table[i, c(2, 5, 8, 11, 14)] <- sds[(5*(i-1) + 1):(5*i), 2]
        Table[i, c(3, 6, 9, 12, 15)] <- coverage[(5*(i-1) + 1):(5*i), 2]
      }  
    }
  } else {
    if (ipw) {
      for (i in 1:6) {
        Table[i, c(1, 3, 5, 7, 9)] <- est[(5*(i-1) + 1):(5*i), 3] - true.values[(5*(i-1) + 1):(5*i)]
        Table[i, c(2, 4, 6, 8, 10)] <- coverage[(5*(i-1) + 1):(5*i), 3]
      }  
    } else {
      for (i in 1:6) {
        Table[i, c(1, 3, 5, 7, 9)] <- est[(5*(i-1) + 1):(5*i), 2] - true.values[(5*(i-1) + 1):(5*i)]
        Table[i, c(2, 4, 6, 8, 10)] <- coverage[(5*(i-1) + 1):(5*i), 2]
      }  
    }
  }
  
  ## End.
  res.list <- list("est" = est, "coverage" = coverage, "Table" = Table, 
                 "true.values" = true.values, "diagnostics" = diagnostics, "probs" = probs)
  if (se) res.list$se <- sds
  return(res.list)
}

## Function to create boxplots of results.
## This was used in an older version of the code.
plot_res <- function (res, Scenarios, dag, theta, ylim = NULL, ipw = FALSE, title = "") {
  
  ## Extract data to be plotted.
  iter <- nrow(res[[1]]$cca)
  which.scen <- which(Scenarios$DAG == dag & Scenarios$theta == theta)
  ll <- length(which.scen)
  estimates <- matrix(NA, iter * ll, 3)
  colnames(estimates) <- c("beta.sel", "em.strength", "est")
  estimates[, 1] <- rep(Scenarios$beta.sel[which.scen], each = iter)
  estimates[, 2] <- rep(Scenarios$em.strength[which.scen], each = iter)
  if (ipw) {
    for (i in 1:ll) estimates[((i-1) * iter + 1):(i * iter), 3] <- res[[which.scen[i]]]$ipw[, 1]
  } else {
    for (i in 1:ll) estimates[((i-1) * iter + 1):(i * iter), 3] <- res[[which.scen[i]]]$cca[, 1]
  }
  estimates <- as.data.frame(estimates)
  
  ## Plot the results.
  plot1 <- ggplot(estimates, aes(x = factor(em.strength), y = est, fill = factor(beta.sel))) +
    geom_boxplot(outlier.size = 1) +
    labs(x = expression(beta[het]), y = "Causal Effect", fill = expression(delta[VS])) +
    geom_hline(yintercept = theta, color = "grey", linetype = "dashed", linewidth = 1) +
    ggtitle(title) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1)))
  if (!is.null(ylim)) plot1 <- plot1 +
    scale_y_continuous(limits = ylim)
  plot1
    
}

## Save results.
save(Sim_gen, cover, analyze_results, plot_res, file = "SBaH_Functions.RData")

##################################################
