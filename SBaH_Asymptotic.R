
##########   TYPE 2 SELECTION BIAS IN MR   ##########

## Take one scenario, investigate for a large n to show
## that there is indeed bias under G-X heterogeneity.

## Set this up.
options(digits = 4)
library(ivreg)
library(OneSampleMR)
library(estimatr)
library(survival)
library(ivtools)
library(sandwich)
library(flextable)
library(officer)

## Load functions to be used in the simulation.
load("SBaH_Functions.RData")

## Specify parameter values for each run.
Parameters.asy <- data.frame("y.distr" = rep(c("normal", "binary"), each = 20), 
                             "theta" = rep(rep(c(0, 0.5), each = 10), times = 2),
                             "beta.sel" = rep(rep(c(1, 2), each = 5), times = 4),
                             "em.strength" = rep(c(-1, -0.5, 0, 0.5, 1), times = 8))

## No need for a "Scenarios" table, as we only consider one scenario.

## Store all results here.
all.asy.sims <- vector("list", 40)

## Run it.
set.seed(3383)
system.time({
  for (J in 1:40) {
    sim <- Sim_gen(iter = 1, n = 1e7, y.distr = Parameters.asy$y.distr[J], dag = 3, 
                   em.strength = Parameters.asy$em.strength[J], 
                   beta.sel = Parameters.asy$beta.sel[J], 
                   theta = Parameters.asy$theta[J], ipw = FALSE)
    all.asy.sims[[J]] <- sim
    print(paste("----- Scenario", J, "done. -----", sep = " "))
  }
})

## Save results.
save(all.asy.sims, file = "SBaH_Asymptotic.RData")
#load("SBaH_Asymptotic.RData")

##################################################

##########   RESULTS TABLE   ##########

## Create a table for the paper.
Asy.table <- matrix(NA, 20, 15)
Asy.table[, 1:3] <- as.matrix(Parameters.asy[1:20, 2:4])
colnames(Asy.table) <- c("theta", "beta.sel", "em.strength", "No_mean", "No_l", "No_u", "Nc_mean", "Nc_l", "Nc_u", 
                         "Bo_mean", "Bo_l", "Bo_u", "Bc_mean", "Bc_l", "Bc_u")

## Fill in the values.
for (i in 1:20) {
  Asy.table[i, 4:6] <- c(all.asy.sims[[i]]$oracle[c(1, 3, 4)])
  Asy.table[i, 7:9] <- c(all.asy.sims[[i]]$cca[c(1, 3, 4)])
  Asy.table[i, 10:12] <- c(all.asy.sims[[20 + i]]$oracle[c(1, 3, 4)])
  Asy.table[i, 13:15] <- c(all.asy.sims[[20 + i]]$cca[c(1, 3, 4)])
}

## Use flextable to export to MS Word (with 3 decimals).
at1 <- flextable(as.data.frame(apply(Asy.table, 2, function(x) sprintf("%.3f", x))))

## Do the export.
doc <- read_docx()
doc <- body_add_par(doc, "Table S4", style = "heading 1")
doc <- body_add_flextable(doc, at1)
print(doc, target = "Asymptotic_results_tables.docx")

##################################################