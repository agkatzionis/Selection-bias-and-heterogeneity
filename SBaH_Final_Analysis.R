
##########   TYPE 2 SELECTION BIAS IN MR   ##########

## Analyze the results of simulations run on the cluster.

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

## Set working directory.
setwd("...")

## Load functions to be used in the simulation.
load("SBaH_Functions.RData")

########################################

##########   RESULTS PER SCENARIO   ##########

## Assess the results for one of the 12 scenarios.
## Use function "analyze_results" for that.

## Reminder of what the different scenarios were.
Parameters <- data.frame("theta" = rep(c(0, 0.5), each = 15),
                         "beta.sel" = rep(rep(c(0.2, 0.5, 1), each = 5), times = 2),
                         "em.strength" = rep(c(-0.5, -0.25, 0, 0.25, 0.5), times = 6))
Scenarios <- data.frame("DAG" = rep(c(1, 2, 3, 4), times = 3),
                        "y.distr" = rep(c("normal", "binary", "cox"), each = 4),
                        "ipw" = rep(c(TRUE, TRUE, FALSE), each = 4))

## Scenario 1:
theta <- Parameters$theta
res <- analyze_results(scenario = 1, true.values = theta, ipw = FALSE)
res$est
res$coverage
res$diagnostics
res$Table

## Scenario 3 (G-X effect modification):
res <- analyze_results(scenario = 3, true.values = theta, ipw = FALSE)
res$Table

## Scenario 1 with IPW:
res <- analyze_results(scenario = 1, true.values = theta, ipw = TRUE)
res$Table

## Scenario 5 (binary Y):
res <- analyze_results(scenario = 5, ipw = FALSE)
res$Table

## Scenario 11 (exponential Y, G-X effect modification):
res <- analyze_results(scenario = 11, ipw = FALSE)
res$Table

########################################

##########   CREATE TABLES   ##########

## Create the tables for the manuscript.

## Main Table:
Table1 <- rbind(analyze_results(scenario = 1, se = TRUE, true.values = theta, ipw = FALSE)$Table, 
                analyze_results(scenario = 5, se = TRUE, ipw = FALSE)$Table, 
                analyze_results(scenario = 3, se = TRUE, true.values = theta, ipw = FALSE)$Table, 
                analyze_results(scenario = 7, se = TRUE, ipw = FALSE)$Table)

## Supplementary Table 1: effect modification by U.
STable1 <- rbind(analyze_results(scenario = 2, se = TRUE, true.values = theta, ipw = FALSE)$Table, 
                 analyze_results(scenario = 6, se = TRUE, ipw = FALSE)$Table, 
                 analyze_results(scenario = 4, se = TRUE, true.values = theta, ipw = FALSE)$Table, 
                 analyze_results(scenario = 8, se = TRUE, ipw = FALSE)$Table)

## Supplementary Table 2: survival outcome.
STable2 <- rbind(analyze_results(scenario = 9, se = TRUE, ipw = FALSE)$Table, 
                 analyze_results(scenario = 11, se = TRUE, ipw = FALSE)$Table)

## Supplementary Table 3: Inverse Probability Weighting.
STable3 <- rbind(analyze_results(scenario = 1, se = TRUE, true.values = theta, ipw = TRUE)$Table, 
                 analyze_results(scenario = 3, se = TRUE, true.values = theta, ipw = TRUE)$Table)

## Use flextable to export to MS Word (with 3 decimals).
t1 <- flextable(as.data.frame(apply(Table1, 2, function(x) sprintf("%.3f", x))))
st1 <- flextable(as.data.frame(apply(STable1, 2, function(x) sprintf("%.3f", x))))
st2 <- flextable(as.data.frame(apply(STable2, 2, function(x) sprintf("%.3f", x))))
st3 <- flextable(as.data.frame(apply(STable3, 2, function(x) sprintf("%.3f", x))))

## Do the export.
doc <- read_docx()
doc <- body_add_par(doc, "Table 1", style = "heading 1")
doc <- body_add_flextable(doc, t1)
doc <- body_add_par(doc, "Supplementary Table 1", style = "heading 1")
doc <- body_add_flextable(doc, st1)
doc <- body_add_par(doc, "Supplementary Table 2", style = "heading 1")
doc <- body_add_flextable(doc, st2)
doc <- body_add_par(doc, "Supplementary Table 3", style = "heading 1")
doc <- body_add_flextable(doc, st3)
print(doc, target = "Results_tables.docx")

########################################
