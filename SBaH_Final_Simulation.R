
##########   TYPE 2 SELECTION BIAS IN MR   ##########

## Simulations on the cluster.

## Set this up.
options(digits = 4)
library(ivreg)
library(OneSampleMR)
library(estimatr)
library(survival)
library(ivtools)
library(sandwich)

## Load functions to be used in the simulation.
load("SBaH_Functions.RData")

## Generate array ID and file name.
ID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
filename <- paste("SBaH_sim", ID, ".RData", sep = "")

## Specify parameter values for each run.
Parameters <- data.frame("theta" = rep(c(0, 0.5), each = 15),
                        "beta.sel" = rep(rep(c(0.2, 0.5, 1), each = 5), times = 2),
                        "em.strength" = rep(c(-0.5, -0.25, 0, 0.25, 0.5), times = 6))

## Specify DAGs and models for different runs.
Scenarios <- data.frame("DAG" = rep(c(1, 2, 3, 4), times = 3),
                        "y.distr" = rep(c("normal", "binary", "cox"), each = 4),
                        "ipw" = rep(c(TRUE, TRUE, FALSE), each = 4))

## Store all results here.
all.sims <- vector("list", 30)

## Run it.
set.seed(3183 + ID)
system.time({
for (J in 1:30) {
  sim <- Sim_gen(iter = 1000, n = 1e4, y.distr = Scenarios$y.distr[ID], dag = Scenarios$DAG[ID], em.strength = Parameters$em.strength[J], beta.sel = Parameters$beta.sel[J], theta = Parameters$theta[J], ipw = Scenarios$ipw[ID])
  all.sims[[J]] <- sim
  print(paste("----- Scenario", J, "done. -----", sep = " "))
}
})

## Save results.
save(all.sims, file = filename)

