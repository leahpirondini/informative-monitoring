# ##################################################
# Script: G-computation adapted to incorporate monitoring
# Date: 01/06/2026
# Author: Leah Pirondini
# Notes: 
# ##################################################

library(tidyverse)
library(data.table)
library(boot)
library(gfoRmula)

dat.long.gcomp <- dat.long

# -------------------
# setup
# -------------------

dat.long.gcomp <- as.data.table(dat.long.gcomp)
id <- "id"
time_name <- "time"
time_points <- n.obs
n <- 5000

# -------------------
# covariates
# -------------------

covnames <- c("A","N","L")
covtypes <- c("binary","binary","normal")
basecovs <- NA
covparams <- list(covmodels = c(
  
  A ~ L + lag1_A + N,
  
  N ~ lag1_L + lag1_A + lag1_N,
  
  L ~ lag1_L + N + lag1_N + time +
    I(N*lag1_L) +  
    I(N*lag1_A) + 
    I(N*lag2_A) +
    I(N*lag3_A) +
    I(N*lag4_A) +
    I(N*lag1_N) +
    I(N*lag1_N*time) +
    I(N*lag1_N*lag2_A) + 
    I(N*lag1_N*lag1_L)
  
))

# -------------------
# outcome
# -------------------

outcome_name <- "event"
outcome_type <- "survival"

ymodel <- event ~ factor(time) + 
  L + lag1_L + lag2_L +
  A + lag1_A + lag2_A +
  N + lag1_N + lag2_N +
  L*N + lag1_L*lag1_N + lag2_L*lag2_N


histories <- c(lagged)
histvars <- list(c("A","N","L"))



# ------------------------------------
# ------------------------------------
# functions for treatment strategies
# ------------------------------------
# ------------------------------------

# ---------------------------
# EARLY TREAT
# ---------------------------

# set A0, A1, A2 = 1
# then allow treatment to follow natural course
# the default in gfoRmula is to let A take its natural (simulated) course
# so function only needs to set values for t<=2

int_fn_early <- function (newdf, pool, intvar, intvals, time_name, t){
  if(t<=2){
    newdf[t <= 2, (intvar) := 1]
  }  
}

# ---------------------------
# WAIT TO TREAT
# ---------------------------

# set A0, A1 = 0
# then require:
# - treatment initiation at either k=2 or k=3 AND
# - two consecutive A_k = 1 during A2, A3, A4
# this means:
# - A2 is a grace-period and can take natural value
# - if A2=0, require A3=1 and A4=1
# - if A2=1, require A3=1 and A4 is a grace-period and can take natural value

int_fn_wait <- function (newdf, pool, intvar, intvals, time_name, t){
  if(t %in% c(0,1)){
    newdf[, (intvar) := 0]
  } else if(t==3){
    newdf[, (intvar) := 1]
  } else if(t==4){
    # set A4 to 1 if A2 == 0, else let A4 take natural course
    # 'pool' contains all prior history - use that for lookback
    # pool is a list of length t + 1, with pool[[i + 1]] being the data at time i
    assign(paste0("pool_t", t), pool, envir = .GlobalEnv)
    A2_vals <- pool_t4$A[pool_t4$time==2]
    newvals <- ifelse(A2_vals == 0, 1, newdf[[intvar]])
    set(newdf, j = intvar, value = newvals)
    
  }
}


# variable(s) on which we wish to intervene
# list of length = no. regimes
intvars <- list("A","A","A","A")

# name of the intervention(s) (for reference)
int_descript <- c("Treat early","Wait to treat","Always treat","Never treat")

# interventions
interventions <- list(list(c(int_fn_early)),
                      list(c(int_fn_wait)),
                      list(c(static, rep(1, time_points))),
                      list(c(static, rep(0, time_points)))
)

# -------------------
# run gformula
# -------------------

gform_output <- gfoRmula::gformula(
  obs_data = dat.long.gcomp,
  id = id,
  time_points = time_points,
  time_name = time_name,
  basecovs = NA,
  covnames = covnames,
  covtypes = covtypes,
  covparams = covparams,
  outcome_name = outcome_name,
  outcome_type = outcome_type,
  ymodel = ymodel,
  intvars = intvars,
  interventions = interventions,
  int_descript = int_descript,
  histories = histories,
  histvars = histvars,
  nsimul = n,
  nsamples = 250,
  ci_method = "normal",
  boot_diag = TRUE,
  seed = 100, #+ seed.vector[a],
  sim_data_b = FALSE)

gform_results <- as.data.frame(gform_output$result[,c(1,2,4,5)])


# -----------
# early
# -----------

gcomp.surv.early <- gform_results %>%
  filter(Interv.==1) %>%
  dplyr::select(`g-form risk`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.surv.early <- as.data.frame(t(1-gcomp.surv.early))
colnames(gcomp.surv.early) <- paste0("surv.time",1:n.obs,".early")

gcomp.se.early <- gform_results %>%
  filter(Interv.==1) %>%
  dplyr::select(`Risk SE`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.se.early <- as.data.frame(t(gcomp.se.early))
colnames(gcomp.se.early) <- paste0("se.time",1:n.obs,".early")

# -----------
# wait
# -----------

gcomp.surv.wait <- gform_results %>%
  filter(Interv.==2) %>%
  dplyr::select(`g-form risk`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.surv.wait <- as.data.frame(t(1-gcomp.surv.wait))
colnames(gcomp.surv.wait) <- paste0("surv.time",1:n.obs,".wait")

gcomp.se.wait <- gform_results %>%
  filter(Interv.==1) %>%
  dplyr::select(`Risk SE`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.se.wait <- as.data.frame(t(gcomp.se.wait))
colnames(gcomp.se.wait) <- paste0("se.time",1:n.obs,".wait")

# -----------
# always
# -----------

gcomp.surv.always <- gform_results %>%
  filter(Interv.==3) %>%
  dplyr::select(`g-form risk`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.surv.always <- as.data.frame(t(1-gcomp.surv.always))
colnames(gcomp.surv.always) <- paste0("surv.time",1:n.obs,".always")

gcomp.se.always <- gform_results %>%
  filter(Interv.==1) %>%
  dplyr::select(`Risk SE`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.se.always <- as.data.frame(t(gcomp.se.always))
colnames(gcomp.se.always) <- paste0("se.time",1:n.obs,".always")

# -----------
# never
# -----------

gcomp.surv.never <- gform_results %>%
  filter(Interv.==4) %>%
  dplyr::select(`g-form risk`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.surv.never <- as.data.frame(t(1-gcomp.surv.never))
colnames(gcomp.surv.never) <- paste0("surv.time",1:n.obs,".never")

gcomp.se.never <- gform_results %>%
  filter(Interv.==1) %>%
  dplyr::select(`Risk SE`) %>%
  as.vector() %>%
  unlist() %>%
  unname()

gcomp.se.never <- as.data.frame(t(gcomp.se.never))
colnames(gcomp.se.never) <- paste0("se.time",1:n.obs,".never")

# -----------
# combine
# -----------

summary.gcomp <- cbind(gcomp.surv.early,
                       gcomp.surv.wait,
                       gcomp.surv.always,
                       gcomp.surv.never,
                       gcomp.se.early,
                       gcomp.se.wait,
                       gcomp.se.always,
                       gcomp.se.never)
