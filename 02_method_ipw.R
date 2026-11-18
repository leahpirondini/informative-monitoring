# ##################################################
# Script: IPW adapted to incorporate monitoring
# Date: 01/06/2026
# Author: Leah Pirondini
# Notes: 
# ##################################################

library(tidyverse)
library(boot)

# ------------------
# data
# ------------------

dat.long.ipw <- dat.long

# ------------------
# get unique ids
# ------------------

id.data <- data.frame(id = unique(dat.long.ipw$id))

# ------------------
# calculate weights and estimate counterfactual risks
# within a function that can be bootstrapped
# ------------------

ipw.formula <- function(id.data, indices){
  
  boot.ids <- id.data$id[indices]
  
  d <- bind_rows(
    lapply(seq_along(boot.ids), function(i){
      
      tmp <- dat.long %>%
        filter(id == boot.ids[i])
      
      # renumber duplicated ids
      tmp$id <- i
      
      tmp
    })
  )
  
  dat.long.ipw <- d
  
  # ------------------
  # ------------------
  # run IPW
  # ------------------
  # ------------------
  
  # ------------------
  # stabilised weights
  # ------------------
  
  # denominator
  wt.mod.A0 <- glm(A~L, family="binomial", data=dat.long.ipw %>% filter(time==0))
  wt.mod.A1plus <- glm(A~L+Alag1+N, family="binomial", data=dat.long.ipw %>% filter(time>0))
  
  pred.wt.A0 <- predict(wt.mod.A0, type = "response", newdata = dat.long.ipw)
  pred.wt.A1plus <- predict(wt.mod.A1plus, type = "response", newdata = dat.long.ipw)
  pred.wt <- data.frame(id = dat.long.ipw$id,
                        time = dat.long.ipw$time,
                        pred.wt.A0,
                        pred.wt.A1plus)
  pred.wt <- pred.wt %>% mutate(pred.wt = ifelse(time==0,pred.wt.A0,pred.wt.A1plus))
  pred.wt <- pred.wt$pred.wt
  
  # numerator
  wt.mod.num <- glm(A~Alag1, family="binomial", data=dat.long.ipw)
  pred.wt.num <- predict(wt.mod.num, type = "response")
  
  # add weights to data
  dat.long.ipw <- dat.long.ipw %>%
    mutate(pred.wt = pred.wt,
           pred.wt.num = pred.wt.num)
  
  # -------------------------------------------------------------------------------------
  # construct the weights based on treatment strategy
  # -------------------------------------------------------------------------------------
  
  dat.long.ipw <- dat.long.ipw %>%
    mutate(
      
      # weights for the "early treat" strategy
      wt.early = case_when(time >= 3 ~ 1,
                           time %in% c(0,1,2) & A==1 ~ pred.wt,
                           time %in% c(0,1,2) & A==0 ~ 0),
      wt.early.cum = ave(wt.early, id, FUN=cumprod),
      wt.early.num = case_when(time >= 3 ~ 1,
                               time %in% c(0,1,2) & A==1 ~ pred.wt.num,
                               time %in% c(0,1,2) & A==0 ~ 0),
      wt.early.cum.num = ave(wt.early.num, id, FUN=cumprod),
      ipw.early = ifelse(wt.early.cum==0,0,1/wt.early.cum),
      ipw.early.stab = ifelse(wt.early.cum==0,0,wt.early.cum.num/wt.early.cum),
      
      # weights for the "wait to treat" strategy
      wt.wait = case_when(time %in% c(0,1) & A==0 ~ 1 - pred.wt,   # t=0 and t=1, A must be 0 => wt = normal IPW
                          time %in% c(0,1) & A==1 ~ 0,             # t=0 and t=1, A must be 0 => if not, censor
                          time==2 ~ 1,                             # t=2, grace-period => wt = 1
                          time==3 & A==1 ~ pred.wt,                # t=3, A must be 1 => wt = normal IPW
                          time==3 & A==0 ~ 0,                      # t=3, A must be 1 => if not, censor
                          time==4 & Alag2==0 & A==1 ~ pred.wt,     # t=4, if A2=0, then A4 must be 1 => wt = normal IPW
                          time==4 & Alag2==0 & A==0 ~ 0,           # t=4, if A2=0, then A4 must be 1 => if not, censor
                          time==4 & Alag2==1 ~ 1),                 # t=4, if A2=1, then A4 is a grace-period => wt=1
      wt.wait.cum = ave(wt.wait, id, FUN=cumprod),
      wt.wait.num = case_when(time %in% c(0,1) & A==0 ~ 1 - pred.wt.num,
                              time %in% c(0,1) & A==1 ~ 0,
                              time==2 ~ 1,        
                              time==3 & A==1 ~ pred.wt.num,
                              time==3 & A==0 ~ 0,
                              time==4 & Alag2==0 & A==1 ~ pred.wt.num, 
                              time==4 & Alag2==0 & A==0 ~ 0,
                              time==4 & Alag2==1 ~ 1),
      wt.wait.cum.num = ave(wt.wait.num, id, FUN=cumprod),
      ipw.wait = ifelse(wt.wait.cum==0,0,1/wt.wait.cum),
      ipw.wait.stab = ifelse(wt.wait.cum==0,0,wt.wait.cum.num/wt.wait.cum),
      
      # weights for the "always treat" strategy
      wt.always = ifelse(A==1, pred.wt, 0),
      wt.always.cum = ave(wt.always, id, FUN=cumprod),
      wt.always.num = ifelse(A==1, pred.wt.num, 0),
      wt.always.cum.num = ave(wt.always.num, id, FUN=cumprod),
      ipw.always = ifelse(wt.always.cum==0,0,1/wt.always.cum),
      ipw.always.stab = ifelse(wt.always.cum==0,0,wt.always.cum.num/wt.always.cum),
      
      # weights for the "never treat" strategy
      wt.never = ifelse(A==0, 1 - pred.wt, 0),
      wt.never.cum = ave(wt.never, id, FUN=cumprod),
      wt.never.num = ifelse(A==0, 1 - pred.wt.num, 0),
      wt.never.cum.num = ave(wt.never.num, id, FUN=cumprod),
      ipw.never = ifelse(wt.never.cum==0,0,1/wt.never.cum),
      ipw.never.stab = ifelse(wt.never.cum==0,0,wt.never.cum.num/wt.never.cum)
      
    )
  
  # ----------------------------
  # "early treat" strategy
  # ----------------------------
  
  ipw_results_early <- dat.long.ipw %>%
    group_by(time) %>%
    
    # hazard
    summarise(
      ipw_died_early =
        mean(event * ipw.early.stab, na.rm = TRUE) /
        mean(ipw.early.stab, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time) %>%
    mutate(
      ipw_died_early = ifelse(is.na(ipw_died_early), 0, ipw_died_early),
      
      # cumulative event
      ipw_cum_event_early =
        Reduce(
          function(prev, cur) prev + (1 - prev) * cur,
          ipw_died_early,
          accumulate = TRUE
        ),
      
      # survival
      ipw_surv_early = 1 - ipw_cum_event_early
    )
  
  # ----------------------------
  # "wait to treat" strategy
  # ----------------------------
  
  ipw_results_wait <- dat.long.ipw %>%
    group_by(time) %>%
    
    # hazard
    summarise(
      ipw_died_wait =
        mean(event * ipw.wait.stab, na.rm = TRUE) /
        mean(ipw.wait.stab, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time) %>%
    mutate(
      ipw_died_wait = ifelse(is.na(ipw_died_wait), 0, ipw_died_wait),
      
      # cumulative event
      ipw_cum_event_wait =
        Reduce(
          function(prev, cur) prev + (1 - prev) * cur,
          ipw_died_wait,
          accumulate = TRUE
        ),
      
      # survival
      ipw_surv_wait = 1 - ipw_cum_event_wait
    )
  
  # ----------------------------
  # "always treat" strategy
  # ----------------------------
  
  ipw_results_always <- dat.long.ipw %>%
    group_by(time) %>%
    
    # hazard
    summarise(
      ipw_died_always =
        mean(event * ipw.always.stab, na.rm = TRUE) /
        mean(ipw.always.stab, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time) %>%
    mutate(
      ipw_died_always = ifelse(is.na(ipw_died_always), 0, ipw_died_always),
      
      # cumulative event
      ipw_cum_event_always =
        Reduce(
          function(prev, cur) prev + (1 - prev) * cur,
          ipw_died_always,
          accumulate = TRUE
        ),
      
      # survival
      ipw_surv_always = 1 - ipw_cum_event_always
    )
  
  # ----------------------------
  # "never treat" strategy
  # ----------------------------
  
  ipw_results_never <- dat.long.ipw %>%
    group_by(time) %>%
    
    # hazard
    summarise(
      ipw_died_never =
        mean(event * ipw.never.stab, na.rm = TRUE) /
        mean(ipw.never.stab, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time) %>%
    mutate(
      ipw_died_never = ifelse(is.na(ipw_died_never), 0, ipw_died_never),
      
      # cumulative event
      ipw_cum_event_never =
        Reduce(
          function(prev, cur) prev + (1 - prev) * cur,
          ipw_died_never,
          accumulate = TRUE
        ),
      
      # survival
      ipw_surv_never = 1 - ipw_cum_event_never
    )
  
  # -------------------------
  # store simulation results
  # -------------------------
  
  return(
    c(ipw_results_early$ipw_surv_early,
      ipw_results_wait$ipw_surv_wait,
      ipw_results_always$ipw_surv_always,
      ipw_results_never$ipw_surv_never)
  )
  
}

boot.out.ipw <- boot(
  data = id.data,
  statistic = ipw.formula,
  R = 250
)

col.names <- c(paste0("surv.time", 1:5, ".early"),
               paste0("surv.time", 1:5, ".wait"),
               paste0("surv.time", 1:5, ".always"),
               paste0("surv.time", 1:5, ".never"),
               paste0("se.time", 1:5, ".early"),
               paste0("se.time", 1:5, ".wait"),
               paste0("se.time", 1:5, ".always"),
               paste0("se.time", 1:5, ".never"))


boot.se   <- apply(boot.out.ipw$t, 2, sd)

summary.ipw <- as.data.frame(
  t(c(boot.out.ipw$t0, boot.se))
)

colnames(summary.ipw) <- col.names
