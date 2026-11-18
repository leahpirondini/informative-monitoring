# ##################################################
# Script: TMLE adapted to incorporate monitoring
# Date: 01/06/2026
# Author: Leah Pirondini
# Notes: 
# ##################################################

library(tidyverse)
library(lmtp)
library(SuperLearner)
library(gam)
library(nnet)
library(glmnet)

# ----------------------------
# ----------------------------
# setup for lmtp package
# ----------------------------
# ----------------------------

# ----------------------------
# data setup
# ----------------------------

dat_tmle <- dat.long %>%
  select(id,time,L,A,N,event) %>%
  rename(Y=event)

dat_tmle <- dat_tmle %>%
  pivot_wider(
    id_cols = id,
    names_from = time,
    values_from = c(L, A, N, Y),
    names_sep = ""
  )


dat_tmle <- dat_tmle %>%
  mutate(across(all_of(paste0("Y", 0:(n.obs-1))), ~ifelse(is.na(.),1,.))
  )

vars.time.dep <- c("N","L","A","Y")
timepoints <- 1:n.obs
neworder_timevars <- as.vector(sapply(timepoints,function(i) paste0(vars.time.dep,i-1)))
neworder <- c("id",neworder_timevars)
dat_tmle <- dat_tmle[,neworder]
dat_tmle <- dat_tmle %>%
  select(-N0)


# ----------------------------
# SuperLearner libraries
# ----------------------------

SL.lib.outcome <- c("SL.glm","SL.glm.interaction")
SL.lib.treat <- "SL.glm"


# ----------------------------
# covariates
# ----------------------------

Anodes <- paste0("A",0:{n.obs-1})
Ynodes <- paste0("Y",0:{n.obs-1})

for (i in 2:n.obs){
  assign(paste0("Lnodes.",i-1),
         c(paste0("N",i-1),
           paste0("L",i-1)
         ))
}

Lnodes.baseline <- "L0"
Lnodes.nonbaseline <- mget(paste0("Lnodes.", 1:{n.obs-1}))
Lnodes <- c(list(Lnodes.baseline), Lnodes.nonbaseline)
names(Lnodes)[1] <- "Lnodes.0"


# ----------------------------
# treatment strategy functions
# ----------------------------

treat_early <- function(data, trt){
  if (trt %in% c("A0","A1","A2")) 1
  else if (trt %in% c("A3","A4")) data[[trt]]
}

treat_wait <- function(data, trt){
  if (trt %in% c("A0","A1")) 0
  else if(trt == "A2") data[[trt]]
  else if(trt == "A3") 1
  else if(trt == "A4") ifelse(data$A2 == 0, 1, data[[trt]])
}

treat_always <- function(data, trt){
  1
}

treat_never <- function(data, trt){
  0
}


# ----------------------------
# ----------------------------
# run tmle
# ----------------------------
# ----------------------------

#progressr::handlers(global = TRUE)

run_lmtp_strategies <- function(data,
                                shifts,
                                trt,
                                outcome,
                                time_vary,
                                learners_outcome,
                                learners_trt,
                                folds,
                                estimator = "lmtp_tmle",
                                return_fits = FALSE,
                                ...) {
  
  stopifnot(is.list(shifts), !is.null(names(shifts)), all(names(shifts) != ""))
  
  # fit each strategy
  fits <- lapply(shifts, function(s) {
    lmtp_survival(
      data = data,
      trt = trt,
      outcome = outcome,
      baseline = NULL,
      time_vary = time_vary,
      cens = NULL,
      compete = NULL,
      shift = s,
      shifted = NULL,
      estimator = estimator,
      k = Inf,
      mtp = TRUE,
      id = "id",
      learners_outcome = learners_outcome,
      learners_trt = learners_trt,
      folds = folds,
      weights = NULL,
      control = lmtp_control(),
      ...
    )
  })
  
  # extract estimates and SEs for each time point, one block per strategy
  summ <- lapply(names(fits), function(nm) {
    tt  <- lapply(fits[[nm]], tidy)
    est <- vapply(tt, function(x) as.numeric(x$estimate),  numeric(1))
    se  <- vapply(tt, function(x) as.numeric(x$std.error), numeric(1))
    out <- as.data.frame(t(c(est, se)))
    colnames(out) <- c(paste0("surv.time", seq_along(est), ".", nm),
                       paste0("se.time",   seq_along(se),  ".", nm))
    out
  })
  
  results <- do.call(cbind, summ)
  
  if (return_fits) list(results = results, fits = fits) else results
}


summary.tmle <- run_lmtp_strategies(
   data = dat_tmle,
   shifts = list(early  = treat_early,
                 wait   = treat_wait,
                 always = treat_always,
                 never  = treat_never),
   trt = Anodes, outcome = Ynodes, time_vary = Lnodes,
   learners_outcome = SL.lib.outcome,
   learners_trt = SL.lib.treat,
   folds = 10
)
