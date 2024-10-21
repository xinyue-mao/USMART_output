
## This file contains all functions used for U-SMART

expit = function(x) exp(x) / (1 + exp(x))

gen.population.data = function (seed=5, 
                                n=1000, 
                                theta1_s1 = c(0,0,0), 
                                theta2_s1 = c(-1,1,0,0,0,0), 
                                theta1_s2 = c(0,0,0, 0,0,0, 0,0), 
                                theta2_s2 = c(-1,1, 0,0,0,0,0,0, 0, 0, 0), 
                                ut1, ut2) {
  
  set.seed(seed)
  
  ## YR(1) first stage response 
  
  R_s1_X_trt1 = c(1,0,0)   #response_stage1_Xmatrix_treatment1
  R_s1_X_trt2 = c(0,1,0)
  R_s1_X_trt3 = c(0,0,1)
  
  P_resp_s1_trt1 = expit( R_s1_X_trt1 %*% theta1_s1)  # Probability of response at stage1 with trt1
  P_resp_s1_trt2 = expit( R_s1_X_trt2 %*% theta1_s1)
  P_resp_s1_trt3 = expit( R_s1_X_trt3 %*% theta1_s1)
  
  YR_s1_trt1 = rbinom(n, 1, P_resp_s1_trt1)  # response_stage1_trt1
  YR_s1_trt2 = rbinom(n, 1, P_resp_s1_trt2)
  YR_s1_trt3 = rbinom(n, 1, P_resp_s1_trt3)
  
  ## YT(1) first stage toxicity
  
  # treatment 1
  T0_s1_X_trt1 = cbind(matrix(rep(c(1,0,1,0,0),n), nrow=n, byrow = T) ,YR_s1_trt1)   # tox0_stage1_Xmatrix_trt1
  T1_s1_X_trt1 = cbind(matrix(rep(c(0,1,1,0,0),n), nrow=n, byrow = T) ,YR_s1_trt1) 
  
  P_tox_le0_s1_trt1 = expit( T0_s1_X_trt1 %*% theta2_s1) 
  P_tox_le1_s1_trt1 = expit( T1_s1_X_trt1 %*% theta2_s1) 
  
  P_tox0_s1_trt1 = P_tox_le0_s1_trt1
  P_tox1_s1_trt1 = P_tox_le1_s1_trt1 - P_tox_le0_s1_trt1
  P_tox2_s1_trt1 = 1 - P_tox_le1_s1_trt1
  P_tox_s1_trt1 = cbind(P_tox0_s1_trt1, P_tox1_s1_trt1, P_tox2_s1_trt1)
  
  # treatment 2 
  T0_s1_X_trt2 = cbind(matrix(rep(c(1,0,0,1,0),n), nrow=n, byrow = T) ,YR_s1_trt2)
  T1_s1_X_trt2 = cbind(matrix(rep(c(0,1,0,1,0),n), nrow=n, byrow = T) ,YR_s1_trt2)
  
  P_tox_le0_s1_trt2 = expit( T0_s1_X_trt2 %*% theta2_s1)
  P_tox_le1_s1_trt2 = expit( T1_s1_X_trt2 %*% theta2_s1)
  
  P_tox0_s1_trt2 = P_tox_le0_s1_trt2
  P_tox1_s1_trt2 = P_tox_le1_s1_trt2 - P_tox_le0_s1_trt2
  P_tox2_s1_trt2 = 1 - P_tox_le1_s1_trt2
  P_tox_s1_trt2 = cbind(P_tox0_s1_trt2, P_tox1_s1_trt2, P_tox2_s1_trt2)
  
  # treatment 3 
  T0_s1_X_trt3 = cbind(matrix(rep(c(1,0,0,0,1),n), nrow=n, byrow = T) ,YR_s1_trt3)
  T1_s1_X_trt3 = cbind(matrix(rep(c(0,1,0,0,1),n), nrow=n, byrow = T) ,YR_s1_trt3)
  
  P_tox_le0_s1_trt3 = expit( T0_s1_X_trt3 %*% theta2_s1)
  P_tox_le1_s1_trt3 = expit( T1_s1_X_trt3 %*% theta2_s1)
  
  P_tox0_s1_trt3 = P_tox_le0_s1_trt3
  P_tox1_s1_trt3 = P_tox_le1_s1_trt3 - P_tox_le0_s1_trt3
  P_tox2_s1_trt3 = 1 - P_tox_le1_s1_trt3
  P_tox_s1_trt3 = cbind(P_tox0_s1_trt3, P_tox1_s1_trt3, P_tox2_s1_trt3)
  
  
  YT_s1_trt1 = YT_s1_trt2 = YT_s1_trt3 = rep(NA, n)
  for (i in 1:n) {
    YT_s1_trt1[i] = which(rmultinom(1,1,prob=P_tox_s1_trt1[i,])==1) - 1
    YT_s1_trt2[i] = which(rmultinom(1,1,prob=P_tox_s1_trt2[i,])==1) - 1
    YT_s1_trt3[i] = which(rmultinom(1,1,prob=P_tox_s1_trt3[i,])==1) - 1
  }
  
  ## YR(2) 2nd stage response
  
  R_s2_X_trt12 = cbind(matrix(rep(c(1,0,0,0,1,0),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3])  ## no tox is the reference
  R_s2_X_trt13 = cbind(matrix(rep(c(1,0,0,0,0,1),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3])
  R_s2_X_trt21 = cbind(matrix(rep(c(0,1,0,1,0,0),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3])
  R_s2_X_trt23 = cbind(matrix(rep(c(0,1,0,0,0,1),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3])
  R_s2_X_trt31 = cbind(matrix(rep(c(0,0,1,1,0,0),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3])
  R_s2_X_trt32 = cbind(matrix(rep(c(0,0,1,0,1,0),n), nrow=n, byrow = T),
                       model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3])
  
  P_resp_s2_trt12 = expit( R_s2_X_trt12 %*% theta1_s2 )
  P_resp_s2_trt13 = expit( R_s2_X_trt13 %*% theta1_s2 )
  P_resp_s2_trt21 = expit( R_s2_X_trt21 %*% theta1_s2 )
  P_resp_s2_trt23 = expit( R_s2_X_trt23 %*% theta1_s2 )
  P_resp_s2_trt31 = expit( R_s2_X_trt31 %*% theta1_s2 )
  P_resp_s2_trt32 = expit( R_s2_X_trt32 %*% theta1_s2 )
  
  summary(cbind(P_resp_s2_trt12, P_resp_s2_trt13, P_resp_s2_trt21, 
                P_resp_s2_trt23, P_resp_s2_trt31, P_resp_s2_trt32))
  
  YR_s2_trt12 = rbinom(n, 1, P_resp_s2_trt12)
  YR_s2_trt13 = rbinom(n, 1, P_resp_s2_trt13)
  YR_s2_trt21 = rbinom(n, 1, P_resp_s2_trt21)
  YR_s2_trt23 = rbinom(n, 1, P_resp_s2_trt23)
  YR_s2_trt31 = rbinom(n, 1, P_resp_s2_trt31)
  YR_s2_trt32 = rbinom(n, 1, P_resp_s2_trt32)
  
  ## YT(2) 2nd stage toxicity
  
  # treatment 1 - 2
  
  T0_s2_X_trt12 = cbind(matrix(rep(c(1,0,1,0,0,0,1,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3], YR_s2_trt12)
  T1_s2_X_trt12 = cbind(matrix(rep(c(0,1,1,0,0,0,1,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3], YR_s2_trt12)
  P_tox_le0_s2_trt12 = expit( T0_s2_X_trt12 %*% theta2_s2 )
  P_tox_le1_s2_trt12 = expit( T1_s2_X_trt12 %*% theta2_s2 )
  
  P_tox0_s2_trt12 = P_tox_le0_s2_trt12
  P_tox1_s2_trt12 = P_tox_le1_s2_trt12 - P_tox_le0_s2_trt12
  P_tox2_s2_trt12 = 1 - P_tox_le1_s2_trt12
  P_tox_s2_trt12 = cbind(P_tox0_s2_trt12, P_tox1_s2_trt12, P_tox2_s2_trt12)
  
  
  # treatment 1 - 3
  
  T0_s2_X_trt13 = cbind(matrix(rep(c(1,0,1,0,0,0,0,1),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3], YR_s2_trt13)
  T1_s2_X_trt13 = cbind(matrix(rep(c(0,1,1,0,0,0,0,1),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt1) - 1)[,2:3], YR_s2_trt13)
  P_tox_le0_s2_trt13 = expit( T0_s2_X_trt13 %*% theta2_s2 )
  P_tox_le1_s2_trt13 = expit( T1_s2_X_trt13 %*% theta2_s2 )
  
  P_tox0_s2_trt13 = P_tox_le0_s2_trt13
  P_tox1_s2_trt13 = P_tox_le1_s2_trt13 - P_tox_le0_s2_trt13
  P_tox2_s2_trt13 = 1 - P_tox_le1_s2_trt13
  P_tox_s2_trt13 = cbind(P_tox0_s2_trt13, P_tox1_s2_trt13, P_tox2_s2_trt13)
  
  # treatment 2 - 1 
  
  T0_s2_X_trt21 = cbind(matrix(rep(c(1,0,0,1,0,1,0,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3], YR_s2_trt21)
  T1_s2_X_trt21 = cbind(matrix(rep(c(0,1,0,1,0,1,0,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3], YR_s2_trt21)
  P_tox_le0_s2_trt21 = expit( T0_s2_X_trt21 %*% theta2_s2 )
  P_tox_le1_s2_trt21 = expit( T1_s2_X_trt21 %*% theta2_s2 )
  
  P_tox0_s2_trt21 = P_tox_le0_s2_trt21
  P_tox1_s2_trt21 = P_tox_le1_s2_trt21 - P_tox_le0_s2_trt21
  P_tox2_s2_trt21 = 1 - P_tox_le1_s2_trt21
  P_tox_s2_trt21 = cbind(P_tox0_s2_trt21, P_tox1_s2_trt21, P_tox2_s2_trt21)
  
  # treatment 2 - 3 
  
  T0_s2_X_trt23 = cbind(matrix(rep(c(1,0,0,1,0,0,0,1),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3], YR_s2_trt23)
  T1_s2_X_trt23 = cbind(matrix(rep(c(0,1,0,1,0,0,0,1),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt2) - 1)[,2:3], YR_s2_trt23)
  P_tox_le0_s2_trt23 = expit( T0_s2_X_trt23 %*% theta2_s2 )
  P_tox_le1_s2_trt23 = expit( T1_s2_X_trt23 %*% theta2_s2 )
  
  P_tox0_s2_trt23 = P_tox_le0_s2_trt23
  P_tox1_s2_trt23 = P_tox_le1_s2_trt23 - P_tox_le0_s2_trt23
  P_tox2_s2_trt23 = 1 - P_tox_le1_s2_trt23
  P_tox_s2_trt23 = cbind(P_tox0_s2_trt23, P_tox1_s2_trt23, P_tox2_s2_trt23)
  
  # treatment 3 - 1 
  
  T0_s2_X_trt31 = cbind(matrix(rep(c(1,0,0,0,1,1,0,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3], YR_s2_trt31)
  T1_s2_X_trt31 = cbind(matrix(rep(c(0,1,0,0,1,1,0,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3], YR_s2_trt31)
  P_tox_le0_s2_trt31 = expit( T0_s2_X_trt31 %*% theta2_s2 )
  P_tox_le1_s2_trt31 = expit( T1_s2_X_trt31 %*% theta2_s2 )
  
  P_tox0_s2_trt31 = P_tox_le0_s2_trt31
  P_tox1_s2_trt31 = P_tox_le1_s2_trt31 - P_tox_le0_s2_trt31
  P_tox2_s2_trt31 = 1 - P_tox_le1_s2_trt31
  P_tox_s2_trt31 = cbind(P_tox0_s2_trt31, P_tox1_s2_trt31, P_tox2_s2_trt31)
  
  # treatment 3 - 2 
  
  T0_s2_X_trt32 = cbind(matrix(rep(c(1,0,0,0,1,0,1,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3], YR_s2_trt32)
  T1_s2_X_trt32 = cbind(matrix(rep(c(0,1,0,0,1,0,1,0),n), nrow=n, byrow = T),
                        model.matrix(~ factor(YT_s1_trt3) - 1)[,2:3], YR_s2_trt32)
  P_tox_le0_s2_trt32 = expit( T0_s2_X_trt32 %*% theta2_s2 )
  P_tox_le1_s2_trt32 = expit( T1_s2_X_trt32 %*% theta2_s2 )
  
  P_tox0_s2_trt32 = P_tox_le0_s2_trt32
  P_tox1_s2_trt32 = P_tox_le1_s2_trt32 - P_tox_le0_s2_trt32
  P_tox2_s2_trt32 = 1 - P_tox_le1_s2_trt32
  P_tox_s2_trt32 = cbind(P_tox0_s2_trt32, P_tox1_s2_trt32, P_tox2_s2_trt32)
  
  YT_s2_trt12 = YT_s2_trt13 = YT_s2_trt21 =
    YT_s2_trt23 = YT_s2_trt31 = YT_s2_trt32 = rep(NA, n)
  for (i in 1:n) {
    YT_s2_trt12[i] = which(rmultinom(1,1,prob=P_tox_s2_trt12[i,])==1) - 1
    YT_s2_trt13[i] = which(rmultinom(1,1,prob=P_tox_s2_trt13[i,])==1) - 1
    YT_s2_trt21[i] = which(rmultinom(1,1,prob=P_tox_s2_trt21[i,])==1) - 1
    YT_s2_trt23[i] = which(rmultinom(1,1,prob=P_tox_s2_trt23[i,])==1) - 1
    YT_s2_trt31[i] = which(rmultinom(1,1,prob=P_tox_s2_trt31[i,])==1) - 1
    YT_s2_trt32[i] = which(rmultinom(1,1,prob=P_tox_s2_trt32[i,])==1) - 1
  }
  
  U_s1_trt1 = ut1[cbind(2-YR_s1_trt1, YT_s1_trt1+1)]
  U_s1_trt2 = ut1[cbind(2-YR_s1_trt2, YT_s1_trt2+1)]
  U_s1_trt3 = ut1[cbind(2-YR_s1_trt3, YT_s1_trt3+1)]
  
  U_s2_trt12 = ut2[cbind(2-YR_s2_trt12, YT_s2_trt12+1)]
  U_s2_trt13 = ut2[cbind(2-YR_s2_trt13, YT_s2_trt13+1)]
  U_s2_trt21 = ut2[cbind(2-YR_s2_trt21, YT_s2_trt21+1)]
  U_s2_trt23 = ut2[cbind(2-YR_s2_trt23, YT_s2_trt23+1)]
  U_s2_trt31 = ut2[cbind(2-YR_s2_trt31, YT_s2_trt31+1)]
  U_s2_trt32 = ut2[cbind(2-YR_s2_trt32, YT_s2_trt32+1)]
  
  raw_sim = as.data.frame(cbind(
    YR_s1_trt1, YT_s1_trt1, U_s1_trt1,
    YR_s1_trt2, YT_s1_trt2, U_s1_trt2,
    YR_s1_trt3, YT_s1_trt3, U_s1_trt3,
    YR_s2_trt12, YT_s2_trt12, U_s2_trt12,
    YR_s2_trt13, YT_s2_trt13, U_s2_trt13,
    YR_s2_trt21, YT_s2_trt21, U_s2_trt21,
    YR_s2_trt23, YT_s2_trt23, U_s2_trt23,
    YR_s2_trt31, YT_s2_trt31, U_s2_trt31,
    YR_s2_trt32, YT_s2_trt32, U_s2_trt32))
  
  sim = raw_sim
  
  responders_trt1_s1 = which(YR_s1_trt1==1)
  responders_trt2_s1 = which(YR_s1_trt2==1)
  responders_trt3_s1 = which(YR_s1_trt3==1)
  
  sim$YR_s2_trt12[responders_trt1_s1] = sim$YR_s2_trt13[responders_trt1_s1] = 
    sim$YT_s2_trt12[responders_trt1_s1] = sim$YT_s2_trt13[responders_trt1_s1] = 
    sim$U_s2_trt12[responders_trt1_s1] = sim$U_s2_trt13[responders_trt1_s1] = NA
  sim$YR_s2_trt21[responders_trt2_s1] = sim$YR_s2_trt23[responders_trt2_s1] = 
    sim$YT_s2_trt21[responders_trt2_s1] = sim$YT_s2_trt23[responders_trt2_s1] = 
    sim$U_s2_trt21[responders_trt2_s1] = sim$U_s2_trt23[responders_trt2_s1] = NA
  sim$YR_s2_trt31[responders_trt3_s1] = sim$YR_s2_trt32[responders_trt3_s1] = 
    sim$YT_s2_trt31[responders_trt3_s1] = sim$YT_s2_trt32[responders_trt3_s1] = 
    sim$U_s2_trt31[responders_trt3_s1] = sim$U_s2_trt32[responders_trt3_s1] = NA
  
  temp_U_s2_trt12 = sim$U_s2_trt12
  temp_U_s2_trt13 = sim$U_s2_trt13
  temp_U_s2_trt21 = sim$U_s2_trt21
  temp_U_s2_trt23 = sim$U_s2_trt23
  temp_U_s2_trt31 = sim$U_s2_trt31
  temp_U_s2_trt32 = sim$U_s2_trt32
  
  temp_U_s2_trt12[sim$YR_s1_trt1==1]=100
  temp_U_s2_trt13[sim$YR_s1_trt1==1]=100
  temp_U_s2_trt21[sim$YR_s1_trt2==1]=100
  temp_U_s2_trt23[sim$YR_s1_trt2==1]=100
  temp_U_s2_trt31[sim$YR_s1_trt3==1]=100
  temp_U_s2_trt32[sim$YR_s1_trt3==1]=100
  
  sim$U_trt12 = sim$U_s1_trt1 * temp_U_s2_trt12
  sim$U_trt13 = sim$U_s1_trt1 * temp_U_s2_trt13
  sim$U_trt21 = sim$U_s1_trt2 * temp_U_s2_trt21
  sim$U_trt23 = sim$U_s1_trt2 * temp_U_s2_trt23
  sim$U_trt31 = sim$U_s1_trt3 * temp_U_s2_trt31
  sim$U_trt32 = sim$U_s1_trt3 * temp_U_s2_trt32
  
  return(sim)
  
}

bound.prob = function(vec, lower_bd=0.1) {
  
  ## vec: a vector of any length summing to 1
  ## lower_bd: each component of the output vector will be larger than lower_bd
  ## and smaller than 1-lower_bd
  
  # Check if the sum of input vector is 1
  if (abs(sum(vec) - 1) > 1e-9) {
    stop("Input vector must sum to 1")
  }
  
  if (length(vec) == 2) {
    if (vec[1] <= lower_bd) {vec = c(lower_bd, 1-lower_bd)}
    if (vec[2] <= lower_bd) {vec = c(1-lower_bd, lower_bd)}
  }
  
  if (length(vec) > 2) {
    # Identify elements below the threshold and replace them with the threshold
    below_threshold = (vec <= lower_bd)
    vec[below_threshold] = lower_bd
    
    remaining = 1 - sum(vec[below_threshold])
    vec[!below_threshold] = vec[!below_threshold] / sum(vec[!below_threshold]) * remaining
  }
  ## check it one more time, in case redistribution cause probability less than 1
  if (abs(sum(vec) - 1) > 1e-9) {
    stop("Input vector must sum to 1")
  }
  
  if (length(vec) == 2) {
    if (vec[1] <= lower_bd) {vec = c(lower_bd, 1-lower_bd)}
    if (vec[2] <= lower_bd) {vec = c(1-lower_bd, lower_bd)}
  }
  
  if (length(vec) > 2) {
    # Identify elements below the threshold and replace them with the threshold
    below_threshold = (vec <= lower_bd)
    vec[below_threshold] = lower_bd
    
    remaining = 1 - sum(vec[below_threshold])
    vec[!below_threshold] = vec[!below_threshold] / sum(vec[!below_threshold]) * remaining
  }
  
  return(vec)
}

boots.max.s1 = function(dat, nboot) {
  
  max_r = rep(NA, nboot)
  n = nrow(dat)
  for (i in 1:nboot) {
    boots_id = sample(c(1:n), size=n, replace = TRUE)
    boots_sample = dat[boots_id,]
    
    U1 = mean(boots_sample$s1_U[which(boots_sample$s1_trt==1)])
    U2 = mean(boots_sample$s1_U[which(boots_sample$s1_trt==2)])
    U3 = mean(boots_sample$s1_U[which(boots_sample$s1_trt==3)])
    
    max_r[i] = which.max(c(U1, U2, U3))
  }
  return(max_r)
  ## return a vector with number of components
  ## each one component is the treatment yielding the max utility
}


boots.max.s2 = function (dat, nboot, new) {
  
  ## dat is the existing dataset, where we will bootstrap from soon
  ## nboot = number of bootstrap
  ## new is the coming dataset, we already has this patient's first stage info
  ## and we want to find the best 2nd stage trt for him/her
  
  max_r = rep(NA, nboot)
  trt = c(1,2,3)
  n = nrow(dat)
  i = 1
  while (i < nboot) {
    
    tryCatch({
      # Resample with replacement
      boots_id = sample(c(1:n), size=n, replace = TRUE)
      boots_sample = dat[boots_id,]
      
        
        data = subset(boots_sample, 
                      s1_trt == new$s1_trt &
                        s1_YR == new$s1_YR &
                        s1_YT == new$s1_YT)
        
        ## if the bootstrapped dataset does not have information for two available treatments
        ## we randomize the patient with equal prob to gather more info
        
        if (length(unique(data$s2_trt)) <2) {
          
          if (new$s1_trt == 1) {adp = c(0, 0.5, 0.5)}
          if (new$s1_trt == 2) {adp = c(0.5, 0, 0.5)}
          if (new$s1_trt == 3) {adp = c(0.5, 0.5, 0)}
          max_r[i] = 0
        }
        
        else{
          if (new$s1_trt == 1) {
            U11 = 0
            U12 = mean(subset(data, s2_trt == 2)$s2_U)
            U13 = mean(subset(data, s2_trt == 3)$s2_U)
            max_r[i] = which.max(c(U11, U12, U13))
          }
          if (new$s1_trt == 2) {
            U22 = 0
            U21 = mean(subset(data, s2_trt == 1)$s2_U)
            U23 = mean(subset(data, s2_trt == 3)$s2_U)
            max_r[i] = which.max(c(U21, U22, U23))
          }
          if (new$s1_trt == 3) {
            U33 = 0
            U31 = mean(subset(data, s2_trt == 1)$s2_U)
            U32 = mean(subset(data, s2_trt == 2)$s2_U)
            max_r[i] = which.max(c(U31, U32, U33))
          }
        }
        
     
    }, error = function(e) {
      cat("Skipping iteration", i, ": encountered an error ->", conditionMessage(e), "\n")
      # Continue to the next iteration (skip this bootstrap sample)
    })
    i = i + 1
  }
  return(max_r)
  ## return a vector with number of components
  ## each one component is the treatment yielding the max utility
}




est = function (data) {
  
  
  ## G-estimator
  
  data_s1_r = subset(data, s1_YR==1)
  data_s1_nr = subset(data, s1_YR==0)
  
  # ghatU_trt1_YT0_YR1 = mean(subset(data_s1_r, s1_trt==1 & s1_YT==0)$U)
  # ghatU_trt1_YT1_YR1 = mean(subset(data_s1_r, s1_trt==1 & s1_YT==1)$U)
  # ghatU_trt1_YT2_YR1 = mean(subset(data_s1_r, s1_trt==1 & s1_YT==2)$U)
  # 
  # ghatU_trt2_YT0_YR1 = mean(subset(data_s1_r, s1_trt==2 & s1_YT==0)$U)
  # ghatU_trt2_YT1_YR1 = mean(subset(data_s1_r, s1_trt==2 & s1_YT==1)$U)
  # ghatU_trt2_YT2_YR1 = mean(subset(data_s1_r, s1_trt==2 & s1_YT==2)$U)
  # 
  # ghatU_trt3_YT0_YR1 = mean(subset(data_s1_r, s1_trt==3 & s1_YT==0)$U)
  # ghatU_trt3_YT1_YR1 = mean(subset(data_s1_r, s1_trt==3 & s1_YT==1)$U)
  # ghatU_trt3_YT2_YR1 = mean(subset(data_s1_r, s1_trt==3 & s1_YT==2)$U)
  
  ghatU_trt1_YT0_YR1 = ghatU_trt2_YT0_YR1 = ghatU_trt3_YT0_YR1 = 10000
  ghatU_trt1_YT1_YR1 = ghatU_trt2_YT1_YR1 = ghatU_trt3_YT1_YR1 = 8000
  ghatU_trt1_YT2_YR1 = ghatU_trt2_YT2_YR1 = ghatU_trt3_YT2_YR1 = 4500
  
  P_YR1_YT0_gv_trt1 = nrow(subset(data_s1_r, s1_trt==1 & s1_YT==0)) / nrow(subset(data, s1_trt==1))
  P_YR1_YT1_gv_trt1 = nrow(subset(data_s1_r, s1_trt==1 & s1_YT==1)) / nrow(subset(data, s1_trt==1))
  P_YR1_YT2_gv_trt1 = nrow(subset(data_s1_r, s1_trt==1 & s1_YT==2)) / nrow(subset(data, s1_trt==1))
  
  P_YR1_YT0_gv_trt2 = nrow(subset(data_s1_r, s1_trt==2 & s1_YT==0)) / nrow(subset(data, s1_trt==2))
  P_YR1_YT1_gv_trt2 = nrow(subset(data_s1_r, s1_trt==2 & s1_YT==1)) / nrow(subset(data, s1_trt==2))
  P_YR1_YT2_gv_trt2 = nrow(subset(data_s1_r, s1_trt==2 & s1_YT==2)) / nrow(subset(data, s1_trt==2))
  
  P_YR1_YT0_gv_trt3 = nrow(subset(data_s1_r, s1_trt==3 & s1_YT==0)) / nrow(subset(data, s1_trt==3))
  P_YR1_YT1_gv_trt3 = nrow(subset(data_s1_r, s1_trt==3 & s1_YT==1)) / nrow(subset(data, s1_trt==3))
  P_YR1_YT2_gv_trt3 = nrow(subset(data_s1_r, s1_trt==3 & s1_YT==2)) / nrow(subset(data, s1_trt==3))
  
  
  ghatU_trt12_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==2 & s1_YT==0)$U)
  ghatU_trt12_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==2 & s1_YT==1)$U)
  ghatU_trt12_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==2 & s1_YT==2)$U)
  
  ghatU_trt13_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==3 & s1_YT==0)$U)
  ghatU_trt13_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==3 & s1_YT==1)$U)
  ghatU_trt13_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==1 & s2_trt==3 & s1_YT==2)$U)
  
  ghatU_trt21_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==1 & s1_YT==0)$U)
  ghatU_trt21_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==1 & s1_YT==1)$U)
  ghatU_trt21_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==1 & s1_YT==2)$U)
  
  ghatU_trt23_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==3 & s1_YT==0)$U)
  ghatU_trt23_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==3 & s1_YT==1)$U)
  ghatU_trt23_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==2 & s2_trt==3 & s1_YT==2)$U)
  
  ghatU_trt31_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==1 & s1_YT==0)$U)
  ghatU_trt31_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==1 & s1_YT==1)$U)
  ghatU_trt31_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==1 & s1_YT==2)$U)
  
  ghatU_trt32_YT0_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==2 & s1_YT==0)$U)
  ghatU_trt32_YT1_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==2 & s1_YT==1)$U)
  ghatU_trt32_YT2_YR0 = mean(subset(data_s1_nr, s1_trt==3 & s2_trt==2 & s1_YT==2)$U)
  
  
  
  P_YR0_YT0_gv_trt1 = nrow(subset(data_s1_nr, s1_trt==1 & s1_YT==0)) / nrow(subset(data, s1_trt==1))
  P_YR0_YT1_gv_trt1 = nrow(subset(data_s1_nr, s1_trt==1 & s1_YT==1)) / nrow(subset(data, s1_trt==1))
  P_YR0_YT2_gv_trt1 = nrow(subset(data_s1_nr, s1_trt==1 & s1_YT==2)) / nrow(subset(data, s1_trt==1))
  
  P_YR0_YT0_gv_trt2 = nrow(subset(data_s1_nr, s1_trt==2 & s1_YT==0)) / nrow(subset(data, s1_trt==2))
  P_YR0_YT1_gv_trt2 = nrow(subset(data_s1_nr, s1_trt==2 & s1_YT==1)) / nrow(subset(data, s1_trt==2))
  P_YR0_YT2_gv_trt2 = nrow(subset(data_s1_nr, s1_trt==2 & s1_YT==2)) / nrow(subset(data, s1_trt==2))
  
  P_YR0_YT0_gv_trt3 = nrow(subset(data_s1_nr, s1_trt==3 & s1_YT==0)) / nrow(subset(data, s1_trt==3))
  P_YR0_YT1_gv_trt3 = nrow(subset(data_s1_nr, s1_trt==3 & s1_YT==1)) / nrow(subset(data, s1_trt==3))
  P_YR0_YT2_gv_trt3 = nrow(subset(data_s1_nr, s1_trt==3 & s1_YT==2)) / nrow(subset(data, s1_trt==3))
  
  
  G12 = ghatU_trt1_YT0_YR1 * P_YR1_YT0_gv_trt1 +
        ghatU_trt1_YT1_YR1 * P_YR1_YT1_gv_trt1 +
        ghatU_trt1_YT2_YR1 * P_YR1_YT2_gv_trt1 +
        ghatU_trt12_YT0_YR0 * P_YR0_YT0_gv_trt1 +
        ghatU_trt12_YT1_YR0 * P_YR0_YT1_gv_trt1 +
        ghatU_trt12_YT2_YR0 * P_YR0_YT2_gv_trt1 
  
  G13 = ghatU_trt1_YT0_YR1 * P_YR1_YT0_gv_trt1 +
        ghatU_trt1_YT1_YR1 * P_YR1_YT1_gv_trt1 +
        ghatU_trt1_YT2_YR1 * P_YR1_YT2_gv_trt1 +
        ghatU_trt13_YT0_YR0 * P_YR0_YT0_gv_trt1 +
        ghatU_trt13_YT1_YR0 * P_YR0_YT1_gv_trt1 +
        ghatU_trt13_YT2_YR0 * P_YR0_YT2_gv_trt1 
  
  G21 = ghatU_trt2_YT0_YR1 * P_YR1_YT0_gv_trt2 +
        ghatU_trt2_YT1_YR1 * P_YR1_YT1_gv_trt2 +
        ghatU_trt2_YT2_YR1 * P_YR1_YT2_gv_trt2 +
        ghatU_trt21_YT0_YR0 * P_YR0_YT0_gv_trt2 +
        ghatU_trt21_YT1_YR0 * P_YR0_YT1_gv_trt2 +
        ghatU_trt21_YT2_YR0 * P_YR0_YT2_gv_trt2 
  
  G23 = ghatU_trt2_YT0_YR1 * P_YR1_YT0_gv_trt2 +
        ghatU_trt2_YT1_YR1 * P_YR1_YT1_gv_trt2 +
        ghatU_trt2_YT2_YR1 * P_YR1_YT2_gv_trt2 +
        ghatU_trt23_YT0_YR0 * P_YR0_YT0_gv_trt2 +
        ghatU_trt23_YT1_YR0 * P_YR0_YT1_gv_trt2 +
        ghatU_trt23_YT2_YR0 * P_YR0_YT2_gv_trt2 
  
  G31 = ghatU_trt3_YT0_YR1 * P_YR1_YT0_gv_trt3 +
        ghatU_trt3_YT1_YR1 * P_YR1_YT1_gv_trt3 +
        ghatU_trt3_YT2_YR1 * P_YR1_YT2_gv_trt3 +
        ghatU_trt31_YT0_YR0 * P_YR0_YT0_gv_trt3 +
        ghatU_trt31_YT1_YR0 * P_YR0_YT1_gv_trt3 +
        ghatU_trt31_YT2_YR0 * P_YR0_YT2_gv_trt3 
  
  G32 = ghatU_trt3_YT0_YR1 * P_YR1_YT0_gv_trt3 +
        ghatU_trt3_YT1_YR1 * P_YR1_YT1_gv_trt3 +
        ghatU_trt3_YT2_YR1 * P_YR1_YT2_gv_trt3 +
        ghatU_trt32_YT0_YR0 * P_YR0_YT0_gv_trt3 +
        ghatU_trt32_YT1_YR0 * P_YR0_YT1_gv_trt3 +
        ghatU_trt32_YT2_YR0 * P_YR0_YT2_gv_trt3 
  G = c(G12, G13, G21, G23, G31, G32)
  G
  ## NIPWE
  n = nrow(data)
  P_s1_trt1 = sum(data$s1_trt == 1) / n
  P_s1_trt2 = sum(data$s1_trt == 2) / n
  P_s1_trt3 = sum(data$s1_trt == 3) / n
  
  P_s2_trt12 = nrow(subset(data, s1_trt==1 & s2_trt==2)) / nrow(subset(data, s1_trt==1 & s1_YR==0))
  P_s2_trt13 = nrow(subset(data, s1_trt==1 & s2_trt==3)) / nrow(subset(data, s1_trt==1 & s1_YR==0))
  P_s2_trt21 = nrow(subset(data, s1_trt==2 & s2_trt==1)) / nrow(subset(data, s1_trt==2 & s1_YR==0))
  P_s2_trt23 = nrow(subset(data, s1_trt==2 & s2_trt==3)) / nrow(subset(data, s1_trt==2 & s1_YR==0))
  P_s2_trt31 = nrow(subset(data, s1_trt==3 & s2_trt==1)) / nrow(subset(data, s1_trt==3 & s1_YR==0))
  P_s2_trt32 = nrow(subset(data, s1_trt==3 & s2_trt==2)) / nrow(subset(data, s1_trt==3 & s1_YR==0))
  
  data$s2_trt[is.na(data$s2_trt)] = 0
  
  IPWE12 = data$s1_YR * data$U * (data$s1_trt==1) / P_s1_trt1 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==1) * (data$s2_trt==2) /(P_s1_trt1*P_s2_trt12)
  IPWE13 = data$s1_YR * data$U * (data$s1_trt==1) / P_s1_trt1 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==1) * (data$s2_trt==3) /(P_s1_trt1*P_s2_trt13)
  IPWE21 = data$s1_YR * data$U * (data$s1_trt==2) / P_s1_trt2 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==2) * (data$s2_trt==1) /(P_s1_trt2*P_s2_trt21)
  IPWE23 = data$s1_YR * data$U * (data$s1_trt==2) / P_s1_trt2 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==2) * (data$s2_trt==3) /(P_s1_trt2*P_s2_trt23)
  IPWE31 = data$s1_YR * data$U * (data$s1_trt==3) / P_s1_trt3 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==3) * (data$s2_trt==1) /(P_s1_trt3*P_s2_trt31)
  IPWE32 = data$s1_YR * data$U * (data$s1_trt==3) / P_s1_trt3 +
        (1 - data$s1_YR) * data$U * (data$s1_trt==3) * (data$s2_trt==2) /(P_s1_trt3*P_s2_trt32)
  
  IPWE = colMeans(cbind(IPWE12, IPWE13, IPWE21, IPWE23, IPWE31, IPWE32))
  IPWE
  
  W12 = data$s1_YR * (data$s1_trt==1) / P_s1_trt1 +
    (1 - data$s1_YR) * (data$s1_trt==1) * (data$s2_trt==2) /(P_s1_trt1*P_s2_trt12)
  W13 = data$s1_YR * (data$s1_trt==1) / P_s1_trt1 +
    (1 - data$s1_YR) * (data$s1_trt==1) * (data$s2_trt==3) /(P_s1_trt1*P_s2_trt13)
  W21 = data$s1_YR * (data$s1_trt==2) / P_s1_trt2 +
    (1 - data$s1_YR) * (data$s1_trt==2) * (data$s2_trt==1) /(P_s1_trt2*P_s2_trt21)
  W23 = data$s1_YR * (data$s1_trt==2) / P_s1_trt2 +
    (1 - data$s1_YR) * (data$s1_trt==2) * (data$s2_trt==3) /(P_s1_trt2*P_s2_trt23)
  W31 = data$s1_YR * (data$s1_trt==3) / P_s1_trt3 +
    (1 - data$s1_YR) * (data$s1_trt==3) * (data$s2_trt==1) /(P_s1_trt3*P_s2_trt31)
  W32 = data$s1_YR * (data$s1_trt==3) / P_s1_trt3 +
    (1 - data$s1_YR) * (data$s1_trt==3) * (data$s2_trt==2) /(P_s1_trt3*P_s2_trt32)
  
  NIPWE12 = W12 %*% data$U / sum(W12)
  NIPWE13 = W13 %*% data$U / sum(W13)
  NIPWE21 = W21 %*% data$U / sum(W21)
  NIPWE23 = W23 %*% data$U / sum(W23)
  NIPWE31 = W31 %*% data$U / sum(W31)
  NIPWE32 = W32 %*% data$U / sum(W32)
  
  NIPWE = c(NIPWE12, NIPWE13, NIPWE21, NIPWE23, NIPWE31, NIPWE32)
  NIPWE
  
  # NIPWE12 = IPWE12 %*% data$U / sum(IPWE12)
  # NIPWE13 = IPWE13 %*% data$U / sum(IPWE13)
  # NIPWE21 = IPWE21 %*% data$U / sum(IPWE21)
  # NIPWE23 = IPWE23 %*% data$U / sum(IPWE23)
  # NIPWE31 = IPWE31 %*% data$U / sum(IPWE31)
  # NIPWE32 = IPWE32 %*% data$U / sum(IPWE32)
 
  
  ## Double Robust 
  
  ghatU_trt12 = mean(subset(data, s1_trt==1 & s2_trt==2)$U)
  ghatU_trt13 = mean(subset(data, s1_trt==1 & s2_trt==3)$U)
  ghatU_trt21 = mean(subset(data, s1_trt==2 & s2_trt==1)$U)
  ghatU_trt23 = mean(subset(data, s1_trt==2 & s2_trt==3)$U)
  ghatU_trt31 = mean(subset(data, s1_trt==3 & s2_trt==1)$U)
  ghatU_trt32 = mean(subset(data, s1_trt==3 & s2_trt==2)$U)
  
  
  DRE12 = mean( data$s1_YR * data$U * (data$s1_trt==1) / P_s1_trt1 -
          data$s1_YR * ghatU_trt12 * ((data$s1_trt==1) - P_s1_trt1) / P_s1_trt1 +
          (1 - data$s1_YR) * data$U * (data$s1_trt==1) * (data$s2_trt==2) /(P_s1_trt1 * P_s2_trt12) -
          (1 - data$s1_YR) * ghatU_trt12 * ((data$s1_trt==1)*(data$s2_trt==2) - P_s1_trt1 * P_s2_trt12) / (P_s1_trt1 * P_s2_trt12) )

  DRE13 = mean( data$s1_YR * data$U * (data$s1_trt==1) / P_s1_trt1 -
    data$s1_YR * ghatU_trt13 * ((data$s1_trt==1) - P_s1_trt1) / P_s1_trt1 +
    (1 - data$s1_YR) * data$U * (data$s1_trt==1) * (data$s2_trt==3) /(P_s1_trt1 * P_s2_trt13) -
    (1 - data$s1_YR) * ghatU_trt13 * ((data$s1_trt==1)*(data$s2_trt==3) - P_s1_trt1 * P_s2_trt13) / (P_s1_trt1 * P_s2_trt13) )

  DRE21 = mean( data$s1_YR * data$U * (data$s1_trt==2) / P_s1_trt2 -
    data$s1_YR * ghatU_trt21 * ((data$s1_trt==2) - P_s1_trt2) / P_s1_trt2 +
    (1 - data$s1_YR) * data$U * (data$s1_trt==2) * (data$s2_trt==1) /(P_s1_trt2 * P_s2_trt21) -
    (1 - data$s1_YR) * ghatU_trt21 * ((data$s1_trt==2)*(data$s2_trt==1) - P_s1_trt2 * P_s2_trt21) / (P_s1_trt2 * P_s2_trt21) )

  DRE23 = mean( data$s1_YR * data$U * (data$s1_trt==2) / P_s1_trt2 -
    data$s1_YR * ghatU_trt23 * ((data$s1_trt==2) - P_s1_trt2) / P_s1_trt2 +
    (1 - data$s1_YR) * data$U * (data$s1_trt==2) * (data$s2_trt==3) /(P_s1_trt2 * P_s2_trt23) -
    (1 - data$s1_YR) * ghatU_trt23 * ((data$s1_trt==2)*(data$s2_trt==3) - P_s1_trt2 * P_s2_trt23) / (P_s1_trt2 * P_s2_trt23) )

  DRE31 = mean( data$s1_YR * data$U * (data$s1_trt==3) / P_s1_trt3 -
    data$s1_YR * ghatU_trt31 * ((data$s1_trt==3) - P_s1_trt3) / P_s1_trt3 +
    (1 - data$s1_YR) * data$U * (data$s1_trt==3) * (data$s2_trt==1) /(P_s1_trt3 * P_s2_trt31) -
    (1 - data$s1_YR) * ghatU_trt31 * ((data$s1_trt==3)*(data$s2_trt==1) - P_s1_trt3 * P_s2_trt31) / (P_s1_trt3 * P_s2_trt31) )

  DRE32 = mean( data$s1_YR * data$U * (data$s1_trt==3) / P_s1_trt3 -
    data$s1_YR * ghatU_trt32 * ((data$s1_trt==3) - P_s1_trt3) / P_s1_trt3 +
    (1 - data$s1_YR) * data$U * (data$s1_trt==3) * (data$s2_trt==2) /(P_s1_trt3 * P_s2_trt32) -
    (1 - data$s1_YR) * ghatU_trt32 * ((data$s1_trt==3)*(data$s2_trt==2) - P_s1_trt3 * P_s2_trt32) / (P_s1_trt3 * P_s2_trt32) )
  
  DRE = c(DRE12, DRE13, DRE21, DRE23, DRE31, DRE32)
  
  est = list(G=G, IPWE=IPWE, NIPWE=NIPWE, DRE=DRE)
  
  return(est)
}


eva = function (data, est) {
  
  # proportion of patients ever treated with the best treatment
  prop_best_trt = (sum(data$s1_trt==3)+sum(na.omit(data)$s2_trt==3)) / nrow(data)
  
  # sum of overall utility
  overall_ut_sample = sum(data$U)
  
  # response rate
  resp_rate = (sum(data$s1_YR==1)+sum(na.omit(data)$s2_YR==1)) / nrow(data)
  
  reg = c(12, 13, 21, 23, 31, 32)
  
  # best regime trt j-k by G-estimation
  best_reg_G = reg[which.max(est$G)]
  
  # best regime trt j-k by NIPWE
  best_reg_IPWE = reg[which.max(est$IPWE)]
  
  # best regime trt j-k by NIPWE
  best_reg_NIPWE = reg[which.max(est$NIPWE)]
  
  # best regime trt j-k by Double Robust
  best_reg_DR = reg[which.max(est$DRE)]
  
  ev = c("prop_best_trt" = prop_best_trt, 
         "overall_ut_sample" = overall_ut_sample, 
         "resp_rate" = resp_rate,
         "best_reg_G" = best_reg_G, 
         "best_reg_IPWE" = best_reg_IPWE,
         "best_reg_NIPWE" = best_reg_NIPWE, 
         "best_reg_DR" = best_reg_DR)
  
  return(ev)
  
}




simul = function(sim1, seed, n1=50, n2=70, c=1, nboot=200) {
  
  # n1 = 50       # burn in sample size at stage 1
  # n2 = 70       # burn in sample size at stage 2
  # c = 1         # tuning parameter
  # nboot = 200   # number of bootstrap in algorithm3
  
  n = nrow(sim1)
  set.seed(seed)
  
  
  trt = c(1,2,3)  # all available treatments
  

  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                                Algorithm 1                                  #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  ar1 = as.data.frame(matrix(NA, nrow=n, ncol=8))
  colnames(ar1) = c("s1_trt", "s1_YR", "s1_YT", "s1_U",
                    "s2_trt", "s2_YR", "s2_YT", "s2_U")
  p.ar1 = as.data.frame(matrix(NA, nrow=n, ncol=6))
  colnames(p.ar1) = c("s1_trt1", "s1_trt2", "s1_trt3",
                      "s2_trt1", "s2_trt2", "s2_trt3")
  
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  ar1[1:n1,1] = sample(trt, n1, replace = TRUE)
  for (i in 1:n1) {
    ar1[i,2:4] = sim1[i,(3*(ar1[i,1]-1)+1):(3*ar1[i,1])]
  }
  p.ar1[1:n1, 1:3] = 1/3
  
  u1 = 40  ## define U1>u1 as a good set
  
  # assign treatments with adaptative prob 
  for (i in (n1+1):n) {
    
    P_s1_good_trt1 = nrow(subset(ar1, s1_trt==1 & s1_U>=u1)) / nrow(subset(ar1, s1_trt==1))
    P_s1_good_trt2 = nrow(subset(ar1, s1_trt==2 & s1_U>=u1)) / nrow(subset(ar1, s1_trt==2))
    P_s1_good_trt3 = nrow(subset(ar1, s1_trt==3 & s1_U>=u1)) / nrow(subset(ar1, s1_trt==3))
    
    adaptProb = c(P_s1_good_trt1^c, P_s1_good_trt2^c, P_s1_good_trt3^c) / 
      sum(c(P_s1_good_trt1^c, P_s1_good_trt2^c, P_s1_good_trt3^c))
    
    p.ar1[i,1:3] = bound.prob(adaptProb)
    
    trt_i_s1 = sample(trt, size = 1, prob = adaptProb)
    ar1[i,1:4] = c(trt_i_s1, sim1[i,(3*(trt_i_s1-1)+1):(3*trt_i_s1)])
  }
  
  
  # for non-responders only
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  # at the second stage, treat patient with treatment different from 1st stage
  for (i in 1:n2) {
    
    prob_burn2 = c(0.5,0.5,0.5)
    prob_burn2[ar1$s1_trt[i]] = 0
    
    ar1$s2_trt[i] = sample(trt, 1, prob = prob_burn2)
    
    if (ar1$s1_trt[i]==1 & ar1$s2_trt[i]==2) {ar1[i,6:8] = sim1[i, 10:12]}
    if (ar1$s1_trt[i]==1 & ar1$s2_trt[i]==3) {ar1[i,6:8] = sim1[i, 13:15]}
    if (ar1$s1_trt[i]==2 & ar1$s2_trt[i]==1) {ar1[i,6:8] = sim1[i, 16:18]}
    if (ar1$s1_trt[i]==2 & ar1$s2_trt[i]==3) {ar1[i,6:8] = sim1[i, 19:21]}
    if (ar1$s1_trt[i]==3 & ar1$s2_trt[i]==1) {ar1[i,6:8] = sim1[i, 22:24]}
    if (ar1$s1_trt[i]==3 & ar1$s2_trt[i]==2) {ar1[i,6:8] = sim1[i, 25:27]}
  }
  
  ar1$s2_trt[ar1$s1_YR==1] = NA
  p.ar1[ar1$s1_YR==1,4:6] = NA
  p.ar1[1:n2, 4:6] = 1/2
  for (i in 1:n2) {p.ar1[i, c(ar1$s1_trt[i]+3)] = 0}
  
  
  u2 = 50
  
  # assign treatments with adaptative prob 
  for (i in (n2+1):n) {
    
    if (ar1$s1_YR[i]==1) {ar1$s2_trt[i] = NA}
    else {
      
      data = subset(ar1[1:(i-1),], 
                    s1_trt == ar1$s1_trt[i] &
                      s1_YR == ar1$s1_YR[i] &
                      s1_YT == ar1$s1_YT[i])
      
      if (length(data) <= 1 || length(unique(data$s2_trt)) <2) {
        
        if (ar1$s1_trt[i] == 1) {adp = c(0, 0.5, 0.5)}
        if (ar1$s1_trt[i] == 2) {adp = c(0.5, 0, 0.5)}
        if (ar1$s1_trt[i] == 3) {adp = c(0.5, 0.5, 0)}
        
      }
      else {
        if (ar1$s1_trt[i] == 1) {
          P_s2_good_trt12 = nrow(subset(data, s2_trt == 2 & s2_U>=u2)) / (sum(data$s2_trt == 2)+0.0001)
          P_s2_good_trt13 = nrow(subset(data, s2_trt == 3 & s2_U>=u2)) / (sum(data$s2_trt == 3)+0.0001)
          
          if (P_s2_good_trt12==0 & P_s2_good_trt13==0) {adp = c(0,0.5,0.5)}
          else {
            adaptProb = c(P_s2_good_trt12, P_s2_good_trt13) / (P_s2_good_trt12 + P_s2_good_trt13)    
            adp = c(0, 0, 0); adp[c(2,3)] = bound.prob(adaptProb)
          }
        }
        
        if (ar1$s1_trt[i] == 2) {
          P_s2_good_trt21 = nrow(subset(data, s2_trt == 1 & s2_U>=u2)) / (sum(data$s2_trt == 1)+0.0001)
          P_s2_good_trt23 = nrow(subset(data, s2_trt == 3 & s2_U>=u2)) / (sum(data$s2_trt == 3)+0.0001)
          
          if (P_s2_good_trt21==0 & P_s2_good_trt23==0) {adp = c(0.5,0,0.5)}
          else {
            adaptProb = c(P_s2_good_trt21, P_s2_good_trt23) / (P_s2_good_trt21 + P_s2_good_trt23)      
            adp = c(0, 0, 0); adp[c(1,3)] = bound.prob(adaptProb)
          }
        }
        
        if (ar1$s1_trt[i] == 3) {
          P_s2_good_trt31 = nrow(subset(data, s2_trt == 1 & s2_U>=u2)) / (sum(data$s2_trt == 1)+0.0001)
          P_s2_good_trt32 = nrow(subset(data, s2_trt == 2 & s2_U>=u2)) / (sum(data$s2_trt == 2)+0.0001)
          
          if (P_s2_good_trt31==0 & P_s2_good_trt32==0) {adp=c(0.5,0.5,0)}
          else{adaptProb = c(P_s2_good_trt31, P_s2_good_trt32) / (P_s2_good_trt31 + P_s2_good_trt32)   
          adp = c(0, 0, 0); adp[c(1,2)] = bound.prob(adaptProb)
          }
        }
        
      }
      
      p.ar1[i,4:6] = adp
      
      ar1$s2_trt[i] = sample(trt, size = 1, prob = adp)
      if (ar1$s1_trt[i]==1 & ar1$s2_trt[i]==2) {ar1[i,6:8] = sim1[i, 10:12]}
      if (ar1$s1_trt[i]==1 & ar1$s2_trt[i]==3) {ar1[i,6:8] = sim1[i, 13:15]}
      if (ar1$s1_trt[i]==2 & ar1$s2_trt[i]==1) {ar1[i,6:8] = sim1[i, 16:18]}
      if (ar1$s1_trt[i]==2 & ar1$s2_trt[i]==3) {ar1[i,6:8] = sim1[i, 19:21]}
      if (ar1$s1_trt[i]==3 & ar1$s2_trt[i]==1) {ar1[i,6:8] = sim1[i, 22:24]}
      if (ar1$s1_trt[i]==3 & ar1$s2_trt[i]==2) {ar1[i,6:8] = sim1[i, 25:27]}
    }
  }
  ar1$s2_U[which(ar1$s1_YR == 1)] = 100
  ar1$U = ar1$s1_U * ar1$s2_U
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                                Algorithm 2                                  #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  ar2 = as.data.frame(matrix(NA, nrow=n, ncol=8))
  colnames(ar2) = c("s1_trt", "s1_YR", "s1_YT", "s1_U",
                    "s2_trt", "s2_YR", "s2_YT", "s2_U")
  p.ar2 = as.data.frame(matrix(NA, nrow=n, ncol=6))
  colnames(p.ar2) = c("s1_trt1", "s1_trt2", "s1_trt3",
                      "s2_trt1", "s2_trt2", "s2_trt3")
  
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  ar2[1:n1,1] = sample(trt, n1, replace = TRUE)
  for (i in 1:n1) {
    ar2[i,2:4] = sim1[i,(3*(ar2[i,1]-1)+1):(3*ar2[i,1])]
  }
  p.ar2[1:n1, 1:3] = 1/3
  
  # assign treatments with adaptative prob 
  for (i in (n1+1):n) {
    
    E_U_s1_trt1 = mean(na.omit(ar2$s1_U[ar2$s1_trt == 1]))
    E_U_s1_trt2 = mean(na.omit(ar2$s1_U[ar2$s1_trt == 2]))
    E_U_s1_trt3 = mean(na.omit(ar2$s1_U[ar2$s1_trt == 3]))
    
    adaptProb = c(E_U_s1_trt1^c, E_U_s1_trt2^c, E_U_s1_trt3^c) / 
      sum(c(E_U_s1_trt1^c, E_U_s1_trt2^c, E_U_s1_trt3^c))
    
    p.ar2[i,1:3] = bound.prob(adaptProb)
    
    trt_i_s1 = sample(trt, size = 1, prob = adaptProb)
    ar2[i,1:4] = c(trt_i_s1, sim1[i,(3*(trt_i_s1-1)+1):(3*trt_i_s1)])
    
  }
  
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  for (i in 1:n2) {
    prob_burn2 = c(0.5,0.5,0.5)
    prob_burn2[ar2$s1_trt[i]] = 0
    
    ar2$s2_trt[i] = sample(trt, 1, prob = prob_burn2)
    
    if (ar2$s1_trt[i]==1 & ar2$s2_trt[i]==2) {ar2[i,6:8] = sim1[i, 10:12]}
    if (ar2$s1_trt[i]==1 & ar2$s2_trt[i]==3) {ar2[i,6:8] = sim1[i, 13:15]}
    if (ar2$s1_trt[i]==2 & ar2$s2_trt[i]==1) {ar2[i,6:8] = sim1[i, 16:18]}
    if (ar2$s1_trt[i]==2 & ar2$s2_trt[i]==3) {ar2[i,6:8] = sim1[i, 19:21]}
    if (ar2$s1_trt[i]==3 & ar2$s2_trt[i]==1) {ar2[i,6:8] = sim1[i, 22:24]}
    if (ar2$s1_trt[i]==3 & ar2$s2_trt[i]==2) {ar2[i,6:8] = sim1[i, 25:27]}
  }
  
  ar2$s2_trt[ar2$s1_YR==1] = NA
  p.ar2[1:n2, 4:6] = 1/2
  for (i in 1:n2) {p.ar2[i, c(ar2$s1_trt[i]+3)] = 0}
  
  # for non-responders only
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  # at the second stage, treat patient with treatment different from 1st stage
  for (i in (n2+1):n) {
    
    if (ar2$s1_YR[i]==1) {ar2$s2_trt[i] = NA}
    else {
      
      data = subset(ar2[1:(i-1),], 
                    s1_trt == ar2$s1_trt[i] &
                      s1_YR == ar2$s1_YR[i] &
                      s1_YT == ar2$s1_YT[i])
      
      if (length(data) <= 1 || length(unique(data$s2_trt)) <2) {
        
        if (ar2$s1_trt[i] == 1) {adp = c(0, 0.5, 0.5)}
        if (ar2$s1_trt[i] == 2) {adp = c(0.5, 0, 0.5)}
        if (ar2$s1_trt[i] == 3) {adp = c(0.5, 0.5, 0)}
        
      }
      else {
        if (ar2$s1_trt[i] == 1) {
          EU_s2_trt12 = mean(subset(data, s2_trt == 2)$s2_U)
          EU_s2_trt13 = mean(subset(data, s2_trt == 3)$s2_U)
          adaptProb = c(EU_s2_trt12, EU_s2_trt13) / (EU_s2_trt12 + EU_s2_trt13)    
          adp = c(0, 0, 0); adp[c(2,3)] = bound.prob(adaptProb)
        }
        
        if (ar2$s1_trt[i] == 2) {
          EU_s2_trt21 = mean(subset(data, s2_trt == 1)$s2_U)
          EU_s2_trt23 = mean(subset(data, s2_trt == 3)$s2_U)
          adaptProb = c(EU_s2_trt21, EU_s2_trt23) / (EU_s2_trt21 + EU_s2_trt23)    
          adp = c(0, 0, 0); adp[c(1,3)] = bound.prob(adaptProb)
        }
        
        if (ar2$s1_trt[i] == 3) {
          EU_s2_trt31 = mean(subset(data, s2_trt == 1)$s2_U)
          EU_s2_trt32 = mean(subset(data, s2_trt == 2)$s2_U)
          adaptProb = c(EU_s2_trt31, EU_s2_trt32) / (EU_s2_trt31 + EU_s2_trt32)    
          adp = c(0, 0, 0); adp[c(1,2)] = bound.prob(adaptProb)
        }
        
      }
      
      p.ar2[i,4:6] = adp
      
      ar2$s2_trt[i] = sample(trt, size = 1, prob = adp)
      if (ar2$s1_trt[i]==1 & ar2$s2_trt[i]==2) {ar2[i,6:8] = sim1[i, 10:12]}
      if (ar2$s1_trt[i]==1 & ar2$s2_trt[i]==3) {ar2[i,6:8] = sim1[i, 13:15]}
      if (ar2$s1_trt[i]==2 & ar2$s2_trt[i]==1) {ar2[i,6:8] = sim1[i, 16:18]}
      if (ar2$s1_trt[i]==2 & ar2$s2_trt[i]==3) {ar2[i,6:8] = sim1[i, 19:21]}
      if (ar2$s1_trt[i]==3 & ar2$s2_trt[i]==1) {ar2[i,6:8] = sim1[i, 22:24]}
      if (ar2$s1_trt[i]==3 & ar2$s2_trt[i]==2) {ar2[i,6:8] = sim1[i, 25:27]}
    }
  }
  
  ar2$s2_U[which(ar2$s1_YR == 1)] = 100
  ar2$U = ar2$s1_U * ar2$s2_U
  
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                     Algorithm 3 : bootstrap                                 #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  ar3 = as.data.frame(matrix(NA, nrow=n, ncol=8))
  colnames(ar3) = c("s1_trt", "s1_YR", "s1_YT", "s1_U",
                    "s2_trt", "s2_YR", "s2_YT", "s2_U")
  p.ar3 = as.data.frame(matrix(NA, nrow=n, ncol=6))
  colnames(p.ar3) = c("s1_trt1", "s1_trt2", "s1_trt3",
                      "s2_trt1", "s2_trt2", "s2_trt3")
  
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  ar3[1:n1,1] = sample(trt, n1, replace = TRUE)
  for (i in 1:n1) {
    ar3[i,2:4] = sim1[i,(3*(ar3[i,1]-1)+1):(3*ar3[i,1])]
  }
  p.ar3[1:n1, 1:3] = 1/3
  
  # assign treatments with adaptative prob 
  for (i in (n1+1):n) {
    
    boot_i = boots.max.s1(ar3[1:i-1,c(1,4)],nboot)
    adaptProb = c(sum(boot_i==1)^c, sum(boot_i==2)^c, sum(boot_i==3)^c) / 
      (sum(boot_i==1)^c + sum(boot_i==2)^c +sum(boot_i==3)^c)
    
    p.ar3[i,1:3] = bound.prob(adaptProb)
    
    trt_i_s1 = sample(trt, size = 1, prob = adaptProb)
    ar3[i,1:4] = c(trt_i_s1, sim1[i,(3*(trt_i_s1-1)+1):(3*trt_i_s1)])
    
  }
  
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  for (i in 1:n2) {
    prob_burn2 = c(0.5,0.5,0.5)
    prob_burn2[ar3$s1_trt[i]] = 0
    
    ar3$s2_trt[i] = sample(trt, 1, prob = prob_burn2)
    
    if (ar3$s1_trt[i]==1 & ar3$s2_trt[i]==2) {ar3[i,6:8] = sim1[i, 10:12]}
    if (ar3$s1_trt[i]==1 & ar3$s2_trt[i]==3) {ar3[i,6:8] = sim1[i, 13:15]}
    if (ar3$s1_trt[i]==2 & ar3$s2_trt[i]==1) {ar3[i,6:8] = sim1[i, 16:18]}
    if (ar3$s1_trt[i]==2 & ar3$s2_trt[i]==3) {ar3[i,6:8] = sim1[i, 19:21]}
    if (ar3$s1_trt[i]==3 & ar3$s2_trt[i]==1) {ar3[i,6:8] = sim1[i, 22:24]}
    if (ar3$s1_trt[i]==3 & ar3$s2_trt[i]==2) {ar3[i,6:8] = sim1[i, 25:27]}
    
  }
  
  ar3$s2_trt[ar3$s1_YR==1] = NA
  p.ar3[1:n2, 4:6] = 1/2
  
  for (i in 1:n2) {p.ar3[i, c(ar3$s1_trt[i]+3)] = 0}
  
  # for non-responders only
  # assign treatments with equal prob in the burn-in period and pull their YR, YT from counterfactual
  # at the second stage, treat patient with treatment different from 1st stage
  for (i in (n2+1):n) {
    
    if (ar3$s1_YR[i]==1) {
      ar3$s2_trt[i] = NA
    } else {
      boot_i = boots.max.s2(dat = ar3[1:(i-1),], nboot=nboot, new = ar3[i,])
      boot_i = na.omit(boot_i)
      
      if (sum(boot_i==0) > 0.7*length(boot_i)) {
        if (ar3$s1_trt[i] == 1) {adp = c(0, 0.5, 0.5)}
        if (ar3$s1_trt[i] == 2) {adp = c(0.5, 0, 0.5)}
        if (ar3$s1_trt[i] == 3) {adp = c(0.5, 0.5, 0)}
      } else {
        boot_i = boot_i[which(boot_i!=0)]
        if (ar3$s1_trt[i] == 1) {
          adaptProb = c( sum(boot_i==2)^c, sum(boot_i==3)^c) / (sum(boot_i==2)^c +sum(boot_i==3)^c)
          adp = c(0, 0, 0); adp[c(2,3)] = bound.prob(adaptProb)
        }
        
        if (ar3$s1_trt[i] == 2) {
          adaptProb = c( sum(boot_i==1)^c, sum(boot_i==3)^c) / (sum(boot_i==1)^c +sum(boot_i==3)^c)
          adp = c(0, 0, 0); adp[c(1,3)] = bound.prob(adaptProb)
        }
        
        if (ar3$s1_trt[i] == 3) {
          adaptProb = c( sum(boot_i==1)^c, sum(boot_i==2)^c) / (sum(boot_i==1)^c +sum(boot_i==2)^c)
          adp = c(0, 0, 0); adp[c(1,2)] = bound.prob(adaptProb)
        }
        
      }
      p.ar3[i,4:6] = adp
      ar3$s2_trt[i] = sample(trt, size = 1, prob = adp)
      
      if (ar3$s1_trt[i]==1 & ar3$s2_trt[i]==2) {ar3[i,6:8] = sim1[i, 10:12]}
      if (ar3$s1_trt[i]==1 & ar3$s2_trt[i]==3) {ar3[i,6:8] = sim1[i, 13:15]}
      if (ar3$s1_trt[i]==2 & ar3$s2_trt[i]==1) {ar3[i,6:8] = sim1[i, 16:18]}
      if (ar3$s1_trt[i]==2 & ar3$s2_trt[i]==3) {ar3[i,6:8] = sim1[i, 19:21]}
      if (ar3$s1_trt[i]==3 & ar3$s2_trt[i]==1) {ar3[i,6:8] = sim1[i, 22:24]}
      if (ar3$s1_trt[i]==3 & ar3$s2_trt[i]==2) {ar3[i,6:8] = sim1[i, 25:27]}
    }
  }
  
  ar3$s2_U[which(ar3$s1_YR == 1)] = 100
  ar3$U = ar3$s1_U * ar3$s2_U
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                          SMART w/o adaptation                               #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  rt = as.data.frame(matrix(NA, nrow=n, ncol=8))
  colnames(rt) = c("s1_trt", "s1_YR", "s1_YT", "s1_U",
                   "s2_trt", "s2_YR", "s2_YT", "s2_U")
  
  rt[1:n,1] = sample(trt, n, replace = TRUE)
  
  for (i in 1:n) {
    rt[i,2:4] = sim1[i,(3*(rt[i,1]-1)+1):(3*rt[i,1])]
  }
  
  p.rt = as.data.frame(matrix(NA, nrow=n, ncol=6))
  colnames(p.rt) = c("s1_trt1", "s1_trt2", "s1_trt3",
                     "s2_trt1", "s2_trt2", "s2_trt3")
  p.rt[1:n, 1:3] = 1/3
  
  for (i in 1:n) {
    prob_burn2 = c(0.5,0.5,0.5)
    prob_burn2[rt$s1_trt[i]] = 0
    p.rt[i,4:6] = prob_burn2
    
    rt$s2_trt[i] = sample(trt, 1, prob = prob_burn2)
    
    if (rt$s1_trt[i]==1 & rt$s2_trt[i]==2) {rt[i,6:8] = sim1[i, 10:12]}
    if (rt$s1_trt[i]==1 & rt$s2_trt[i]==3) {rt[i,6:8] = sim1[i, 13:15]}
    if (rt$s1_trt[i]==2 & rt$s2_trt[i]==1) {rt[i,6:8] = sim1[i, 16:18]}
    if (rt$s1_trt[i]==2 & rt$s2_trt[i]==3) {rt[i,6:8] = sim1[i, 19:21]}
    if (rt$s1_trt[i]==3 & rt$s2_trt[i]==1) {rt[i,6:8] = sim1[i, 22:24]}
    if (rt$s1_trt[i]==3 & rt$s2_trt[i]==2) {rt[i,6:8] = sim1[i, 25:27]}
  }
  
  
  rt$s2_trt[rt$s1_YR==1] = NA
  
  rt$s2_U[which(rt$s1_YR == 1)] = 100
  rt$U = rt$s1_U * rt$s2_U
  
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                                 Estimation                                  #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  est.ar1 = est(ar1)
  est.ar2 = est(ar2)
  est.ar3 = est(ar3)
  est.rt = est(rt)
  
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                         Evaluation Criteria                                 #
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  
  eva.ar1 = eva(ar1, est.ar1)
  eva.ar2 = eva(ar2, est.ar2)
  eva.ar3 = eva(ar3, est.ar3)
  eva.rt = eva(rt, est.rt)
  best_ind = which.max(colMeans(sim1[,28:33]))
  reg = c("12","13","21","23","31","32")
  truth = c(NA,NA,NA, rep(reg[best_ind],3))
  
  ev = list(eva.ar1=eva.ar1, eva.ar2=eva.ar2, eva.ar3=eva.ar3, eva.rt=eva.rt, truth=truth)
  
  return(result = list(ev, sample=sim1, 
                       ar1=ar1, p.ar1=p.ar1, 
                       ar2=ar2, p.ar2=p.ar2,
                       ar3=ar3, p.ar3=p.ar3,
                       smart=rt, p.smart=p.rt,
                       est.ar1=est.ar1, 
                       est.ar2=est.ar2,
                       est.ar3=est.ar3,
                       est.smart=est.rt)) 
  
}


set.seed(100000)

# parameters
pop_size = 100000

theta1_s1 = c(-0.7, 0, 0.35)
theta2_s1 = c(-0.4,2, -0.2,-0.5,-0.2, -0.2)
theta1_s2 = c(-0.2,0,0.1, -0.7,0,0.35, -0.2,-0.4)
theta2_s2 = c(-0.4,2, -0.05,-0.1,-0.05, -0.2,-0.5,-0.2, -0.2, -0.3, -0.2)

ut1 = rbind( c(100, 80, 45), c(40, 20, 5)) # since U=U1*U2, I make the ut1 worst at 5 and ut2 best at 95
ut2 = rbind( c(95, 85, 50), c(25, 15, 0))

# generate a population data
population.data = gen.population.data(seed=3, n=pop_size, theta1_s1, theta2_s1, theta1_s2, theta2_s2, ut1, ut2)





#seed = 10000+SLURM_ID
# seed = 20
# set.seed(seed)
# 
# sample_size = 1000
# patient_id = sample(1:pop_size, sample_size, replace = FALSE)
# sim1 = population.data[patient_id,]
# 
# results = simul(sim1, seed)

