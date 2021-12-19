## Figures 1 and 2 are from AddHealth so they cannot be recreated

## Figure 3

dataGenerator_binary_lik <- function(){
  set.seed(909)
  numGroup = 3
  # Priors
  prior_beta_mu   <- rep(0, numGroup)
  prior_beta_var  <- diag(1,numGroup)
  
  prior_beta_mu_indep   <- 0
  prior_beta_var_indep  <- 1
  
  prior_alpha_mu  <-0
  prior_alpha_var <- 1
  
  priorSet <- list(prior_beta_mu = prior_beta_mu, prior_beta_var = prior_beta_var, 
                   prior_beta_mu_indep = prior_beta_mu_indep, prior_beta_var_indep = prior_beta_var_indep,
                   prior_alpha_mu = prior_alpha_mu, prior_alpha_var = prior_alpha_var)
  
  # Initials
  start_beta_c <- start_beta_r <- matrix(rep(1,3), nrow = 1, ncol = numGroup)
  
  
  start_beta_c_indep <- start_beta_r_indep <-c(2,2)
  
  start_alpha <- -6
  
  initialSet <- list(start_beta_c = start_beta_c, start_beta_r = start_beta_r, start_beta_c_indep=
                       start_beta_c_indep, start_beta_r_indep = start_beta_r_indep, start_alpha = start_alpha)
  
  set.seed(51)
  n <-150
  ones <- as.matrix(rep(1,n), ncol =1 )
  U <- matrix(c(rep(c(1,0,0), 50), rep(c(0,1,0),50), rep(c(0,0,1),50)),
              ncol = 3,nrow = 150, byrow = TRUE)
  
  
  
  dataAll1 <- rnorm(150, 0, 1)
  dataAll2 <- rnorm(150, 1, 1)
  X_row1 <- matrix(rep(dataAll1, 150), nrow = 150, byrow = FALSE)
  X_col1 <- t(X_row1)
  X_row2 <- matrix(rep(dataAll2, 150), nrow = 150, byrow = FALSE)
  X_col2 <- t(X_row2)
  
  beta_r2 <- c(1,0,-1)
  beta_c2 <- c(0,-2,2)
  
  beta_r1 <- c(1,1,1)
  beta_c1 <- c(2, 2,2)
  
  lambda <- 15*matrix(c(.4,.15,.15,
                        .15,.4,.15,
                        .15,.15,.4), nrow = 3, byrow = TRUE)
  
  UbT_r1 <- U %*% beta_r1
  ## lets look at beta U^T
  dUbT_r1 <- diag(as.vector(U %*% beta_r1))
  UbT_r2 <- U %*% beta_r2
  ## lets look at beta U^T
  dUbT_r2 <- diag(as.vector(U %*% beta_r2))
  UbT_c1 <- U %*% beta_c1
  ## lets look at beta U^T
  dUbT_c1 <- diag(as.vector(U %*% beta_c1))
  
  UbT_c2 <- U %*% beta_c2
  ## lets look at beta U^T
  dUbT_c2 <- diag(as.vector(U %*% beta_c2))
  
  ## multiply by X_row to see if we get what we want
  rowXB1 <- dUbT_r1 %*%X_row1
  colXB1 <- X_col1 %*%dUbT_c1
  rowXB2 <- dUbT_r2 %*%X_row2
  colXB2 <- X_col2 %*%dUbT_c2
  
  a <- rnorm(n, 0, 1)
  b <- rnorm(n, 0, 1)
  
  pairs <- amen::rmvnorm(n, c(0,0), Sigma = matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2))
  
  
  ab_couple <-matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      ab_couple[i,j] <- pairs[i,1] + pairs[j,2]
    }
  }
  
  
  # make epsilon matrix
  EPS <- matrix(0, nrow = n, ncol = n)
  
  for(i in 1:n){
    
    for(j in 1:n){
      epsilonPair <- amen::rmvnorm(1, c(0,0), Sigma = matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2))
      
      EPS[i,j] <- epsilonPair[1]
      EPS[j,i] <- epsilonPair[2]
      
    }
  }
  
  #######
  Z = colXB1 + rowXB1 + colXB2 + rowXB2 + EPS + (U %*% lambda %*% t(U)) + ab_couple
  #######
  Y = ifelse(Z>6, 1,0 )

  
  X_r <- array(dim = c(n, n, 2))
  X_c <- array(dim = c(n, n, 2))
  X_r[,,1] <- matrix(rep(dataAll1, n), nrow = n, byrow = FALSE)
  X_c[,,1] <- t(matrix(rep(dataAll1, n), nrow = n, byrow = FALSE))
  X_r[,,2] <-matrix(rep(dataAll2, n), nrow = n, byrow = FALSE)
  X_c[,,2]<- t(matrix(rep(dataAll2, n), nrow = n, byrow = FALSE))
  
  
  
  diag(Y) <- 0
  return(list(Y = Y, X_r = X_r, X_c = X_c, trueBeta_r1 = beta_r1, trueBeta_c1 = beta_c1, trueBeta_c2 = beta_c2, trueBeta_r2 = beta_r2,
              df = cbind(dataAll1, dataAll2), initialSet = initialSet, priorSet = priorSet))
}

data_dep_zero    <- dataGenerator_binary_lik()
data_prior_0     <- data_dep_zero$priorSet
data_initial_0   <- data_dep_zero$initialSet

set.seed(112)

amen_dep_03_goodinit <- amen_master(data_dep_zero$X_r, 
                                    data_dep_zero$X_c, 
                                    NULL, 
                                    data_dep_zero$Y,
                                    iter = 150000,numGroup=3,
                                    prior_beta_mu= data_prior_0$prior_beta_mu,
                                    prior_beta_var = data_prior_0$prior_beta_var,
                                    prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                    prior_alpha_var= data_prior_0$prior_alpha_var,
                                    start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                    start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                    start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                    start_alpha= data_initial_0$start_alpha,FALSE,
                                    TRUE, keep.UV = TRUE, dcor =TRUE,
                                    symmetric = FALSE,
                                    model = "bin", odmax = 0,indepBeta=1,"all",
                                    odens = 10, 5000, FALSE)

amen_dep_03_goodinit_some_indep <- amen_master(data_dep_zero$X_r, 
                                    data_dep_zero$X_c, 
                                    NULL, 
                                    data_dep_zero$Y,
                                    iter =150000,numGroup=3,
                                    prior_beta_mu= data_prior_0$prior_beta_mu,
                                    prior_beta_var = data_prior_0$prior_beta_var,
                                    prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                    prior_alpha_var= data_prior_0$prior_alpha_var,
                                    start_beta_c= matrix(rep(1,6), nrow = 1, ncol = 3),
                                    start_beta_r = matrix(rep(1,3), nrow = 1, ncol = 3),
                                    start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                    start_alpha= data_initial_0$start_alpha,FALSE,
                                    TRUE, keep.UV = TRUE, dcor =TRUE,
                                    symmetric = FALSE,
                                    model = "bin", odmax = 0,indepBeta=1,"some",
                                    odens = 10, 5000, FALSE)

# bad initialization
amen_dep_03_badinit <- amen_master(data_dep_zero$X_r, 
                                   data_dep_zero$X_c, 
                                   NULL, 
                                   data_dep_zero$Y,
                                   iter = 150000,numGroup=3,
                                   prior_beta_mu= data_prior_0$prior_beta_mu,
                                   prior_beta_var = data_prior_0$prior_beta_var,
                                   prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                   prior_alpha_var= data_prior_0$prior_alpha_var,
                                   start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                   start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                   start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                   start_alpha= data_initial_0$start_alpha,FALSE,
                                   TRUE, keep.UV = TRUE, dcor =TRUE,
                                   symmetric = FALSE,
                                   model = "bin", odmax = 0,indepBeta=1,"all",
                                   odens = 10, 5000, TRUE)




# stay at initial est of U

amen_dep_03_initial_U <-amen_master_give_it_initialU (data_dep_zero$X_r, 
                                                          data_dep_zero$X_c, 
                                                          NULL, 
                                                          data_dep_zero$Y,
                                                          iter = 150000,numGroup=3,
                                                          prior_beta_mu= data_prior_0$prior_beta_mu,
                                                          prior_beta_var = data_prior_0$prior_beta_var,
                                                          prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                                          prior_alpha_var= data_prior_0$prior_alpha_var,
                                                          start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                                          start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                                          start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                                          start_alpha= data_initial_0$start_alpha,FALSE,
                                                          TRUE, keep.UV = TRUE, dcor =TRUE,
                                                          symmetric = FALSE,
                                                          model = "bin", odmax = 0,indepBeta=1,"all",
                                                          odens = 10, 5000, FALSE)


# stay at true U
amen_dep_03_true_U <-amen_master_giveTRUTH(data_dep_zero$X_r, 
                                                   data_dep_zero$X_c, 
                                                   NULL, 
                                                   data_dep_zero$Y,
                                                   iter = 150000,numGroup=3,
                                                   prior_beta_mu= data_prior_0$prior_beta_mu,
                                                   prior_beta_var = data_prior_0$prior_beta_var,
                                                   prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                                   prior_alpha_var= data_prior_0$prior_alpha_var,
                                                   start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                                   start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                                   start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                                   start_alpha= data_initial_0$start_alpha,FALSE,
                                                   TRUE, keep.UV = TRUE, dcor =TRUE,
                                                   symmetric = FALSE,
                                                   model = "bin", odmax = 0,indepBeta=1,"all",
                                                   odens = 10, 5000, FALSE)




# sep amen model for each inital community estimate
data<- dataGenerator_binary_lik()
set.seed(1312)
iter <- 150000
burnin <- 5000
numGroup <- 3
UV_start <- specClust(data$Y, tauAvg = TRUE,  3 ,FALSE)
U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
label_partition <- apply(U, 1, function(x) which(x ==1))
group1 <- which(label_partition==1)
group2 <- which(label_partition==2)
group3 <- which(label_partition==3)


X_r_1 <- data$X_r[group1, group1,]
X_r_2 <- data$X_r[group2, group2,]
X_r_3 <- data$X_r[group3, group3,]

X_c_1 <- data$X_c[group1, group1,]
X_c_2 <- data$X_c[group2, group2,]
X_c_3 <- data$X_c[group3, group3,]

Y_1 <- data$Y[group1, group1]
Y_2 <- data$Y[group2, group2]
Y_3 <- data$Y[group3, group3]

amen_2_R_0_1_a <-  ame(Y_1, Xdyad = NULL, Xrow =data$df[group1,] , 
                       Xcol = data$df[group1,], dcor = TRUE,
                       R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                       plot = FALSE)

amen_2_R_0_2_a <-  ame(Y_2, Xdyad = NULL, Xrow =data$df[group2,], Xcol = data$df[group2,], dcor = TRUE,
                       R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                       plot = FALSE)

amen_2_R_0_3_a <-  ame(Y_3, Xdyad = NULL, Xrow =data$df[group3,] , Xcol = data$df[group3,], dcor = TRUE,
                       R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                       plot = FALSE)
# comparing amen models
colnames(amen_2_R_0_1_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
colnames(amen_2_R_0_2_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
colnames(amen_2_R_0_3_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")




int1 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_1_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
int2 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_2_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
int3 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_3_a$BETA), prob_outer = .95, point_est = "mean")[,-1]


int <- as.data.frame(rbind(int1, int2, int3))

# basic AMEN models
sample_list <- list()
for(i in 1:5){
  samp <-fullSimulation_amen_master(data=data, anyIndep="all", family="bin", odmax=0, burnin=5000, iter=150000, 
                                    clusterNum = i,
                                    dataName = "indepDepALL")
  sample_list[[i]] <- samp
  
}



### Code for plot

#Binary Plotting

binary_indep_dep_some <-amen_dep_03_goodinit_some_indep
beta_r2 <- c(1,0,-1)
beta_c2 <- c(0,-2,2)
beta_r1 <- c(1,1,1)
beta_c1 <- c(2, 2,2)
goodinit_ub <- ubeta_posterCI_cluster_may2020(amen_dep_03_goodinit$results,
                                              beta_c1, beta_c2, beta_r1, beta_r2,
                                              150, burnin =5000, indepDep = FALSE)

somegoodinit_ub <- ubeta_posterCI_cluster_some(binary_indep_dep_some,
                                               beta_c1, beta_c2, beta_r1, beta_r2,
                                               150, burnin =5000, indepDep = FALSE)


some_indep <- plotBeta_giveAll_uni_INDEP(binary_indep_dep_some$results,
                                         beta_r1[1], 
                                         beta_c1[1], burnin = 5000)




badinit_ub <- ubeta_posterCI_cluster_may2020(amen_dep_03_badinit$results,
                                             beta_c1, beta_c2, beta_r1, beta_r2,
                                             150,burnin  = 5000,  indepDep = FALSE)

ub_bad_df <- getUB_df(amen_dep_03_badinit$results,
                      beta_c1, beta_c2, beta_r1, beta_r2,
                      150,burnin  = 0,  indepDep = FALSE)

ub_good_df <- getUB_df(amen_dep_03_goodinit$results,
                       beta_c1, beta_c2, beta_r1, beta_r2,
                       150,burnin  = 0,  indepDep = FALSE)





# collect amen models
simulIDB <- c()
for(i in 2:5){
  simulIDB<- rbind(simulIDB,   sample_list[[i]] )
  
}




amenModelsF <- rbind(simulIDB[1,]$int[-1,], simulIDB[2,]$int[-1,], simulIDB[3,]$int[-1,])
amenModelsF$model <- c(rep("R3", 4), rep("R0", 4), rep("MCMCLAMBDAU", 4))
amenModelsF$coef <- rep(c("Indep Covar (Col)","Indep Covar (Row)", "Dep Covar (Col)","Dep Covar (Row)"),3)
sepAMENMODELS <- simulIDB[4,]$int[-c(1,6,11),]
amenPlots <-plotAMEN(amenModelsF, sepAMENMODELS, beta_c1, beta_r1, beta_c2, beta_r2,FALSE)
amenPlots$amenPlot
amenPlots$sepPlot
amenPlots$dummyMelt_sep

dep_dontmoveTRUTH <- amen_dep_03_true_U
dep_dontmove      <- amen_dep_03_initial_U
inital_sep        <- int
sep_initial <- inital_sep[-c(1,6,11),]
amenPlots2 <- plotAMEN(amenModelsF, sep_initial , beta_c1, beta_r1, beta_c2, beta_r2,FALSE)


dontmove <- plot_2Dep(dep_dontmove$results$beta_cMat1,
                      dep_dontmove$results$beta_cMat2,
                      dep_dontmove$results$beta_rMat1,
                      dep_dontmove$results$beta_rMat2,
                      beta_c1, beta_c2,beta_r1,  beta_r2,indepDep = FALSE)

dontmoveTRUTH <- plot_2Dep(dep_dontmoveTRUTH$results$beta_cMat1,
                           dep_dontmoveTRUTH$results$beta_cMat2,
                           dep_dontmoveTRUTH$results$beta_rMat1,
                           dep_dontmoveTRUTH$results$beta_rMat2,
                           beta_c1, beta_c2,beta_r1,  beta_r2,indepDep = FALSE)

## first get all data for dependent rows: 4,8,12
dr1 <- cbind(as.data.frame(amenPlots$allAMEN[c(4,8,12),c(4,6,8:9)]), rep("No Comm"))
dr2 <- cbind(as.data.frame(goodinit_ub[10:12, 1:3]), rep("goodinit", 3), goodinit_ub[10:12, 4])
dr3 <- cbind(somegoodinit_ub$allInt[4:6,1:3 ], rep("goodinittell", 3), somegoodinit_ub$allInt[4:6,4])
dr4 <- cbind(as.data.frame(dontmove[c(10,11,12), c(5,7,9)]), rep("dontupdateu", 3),
             dontmove[10:12, 1])
dr5 <- cbind(amenPlots$sepAMEN[c(4,8,12), c(4,6,8)], rep("truesep", 3), 
             amenPlots$sepAMEN[c(4,8,12), 9])
dr6 <- cbind(amenPlots2$sepAMEN[c(8,4,12), c(4,6,8)], rep("dontmovesep", 3), 
             amenPlots2$sepAMEN[c(4,8,12), 9])
dr7 <- cbind(as.data.frame(dontmoveTRUTH[c(10,11,12), c(5,7,9)]), rep("dontupdateuTRUTH", 3),
             dontmoveTRUTH[10:12, 1])
dr8 <- cbind(as.data.frame(badinit_ub[10:12, 1:3]), rep("badinit", 3), badinit_ub[10:12, 4])

# dep col
dc1 <- cbind(as.data.frame(amenPlots$allAMEN[c(3,7,11),c(4,6,8:9)]), rep("No Comm"))
dc2 <- cbind(as.data.frame(goodinit_ub[7:9, 1:3]), rep("goodinit", 3), goodinit_ub[7:9, 4])
dc3 <- cbind(as.data.frame(somegoodinit_ub$allInt[1:3,1:3]), rep("goodinittell", 3), somegoodinit_ub$allInt[1:3,4])
dc4 <- cbind(as.data.frame(dontmove[c(7,8,9), c(5,7,9)]), rep("dontupdateu", 3),
             dontmove[7:9, 1])
dc5 <- cbind(amenPlots$sepAMEN[c(3,7,11), c(4,6,8)], rep("truesep", 3), 
             amenPlots$sepAMEN[c(3,7,11), 9])
dc6 <- cbind(amenPlots2$sepAMEN[c(7,3,11), c(4,6,8)], rep("dontmovesep", 3), 
             amenPlots2$sepAMEN[c(3,7,11), 9])
dc7 <- cbind(as.data.frame(dontmoveTRUTH[c(7,8,9), c(5,7,9)]), rep("dontupdateuTRUTH", 3),
             dontmoveTRUTH[7:9, 1])
dc8 <- cbind(as.data.frame(badinit_ub[7:9, 1:3]), rep("badinit", 3), badinit_ub[7:9, 4])

# indep covar
i1 <- rbind(cbind(as.data.frame(amenPlots$allAMEN[c(1,5,9),c(4,6,8,9)] ), rep("No Comm", 3)))# col
i2 <- cbind(as.data.frame(goodinit_ub[1:3, 1:3]), rep("goodinit", 3),goodinit_ub[1:3, 4])
i3 <- cbind(some_indep[1,c(5,7,9)], c("goodinittell"), c("No Comm"))
i4 <- cbind(as.data.frame(dontmove[c(1,2,3), c(5,7,9)]), rep("dontupdateu", 3),
            dontmove[1:3, 1])
i5 <- cbind(amenPlots$sepAMEN[c(1,5,9), c(4,6,8)], rep("truesep", 3), 
            amenPlots$sepAMEN[c(1,5,9), 9])
i6 <- cbind(amenPlots2$sepAMEN[c(5,1,9), c(4,6,8)], rep("dontmovesep", 3), 
            amenPlots2$sepAMEN[c(1,5,9), 9])
i7 <- cbind(as.data.frame(dontmoveTRUTH[c(1,2,3), c(5,7,9)]), rep("dontupdateuTRUTH", 3),
            dontmoveTRUTH[1:3, 1])
i8 <- cbind(as.data.frame(badinit_ub[1:3, 1:3]), rep("badinit", 3),badinit_ub[1:3, 4])


i11 <- cbind(as.data.frame(amenPlots$allAMEN[c(2,6,10),c(4,6,8,9)] ), rep("No Comm", 3))# row
i12 <- cbind(as.data.frame(goodinit_ub[4:6, 1:3]), rep("goodinit", 3),goodinit_ub[4:6, 4])
i13 <- cbind(some_indep[2,c(5,7,9)], c("goodinittell"), c("No Comm"))
i14 <- cbind(as.data.frame(dontmove[c(4,5,6), c(5,7,9)]), rep("dontupdateu", 3),
             dontmove[4:6, 1])
i15 <- cbind(amenPlots$sepAMEN[c(2,6,10), c(4,6,8)], rep("truesep", 3), 
             amenPlots$sepAMEN[c(2,6,10), 9])
i16 <- cbind(amenPlots2$sepAMEN[c(6,2,10), c(4,6,8)], rep("dontmovesep", 3), 
             amenPlots2$sepAMEN[c(2,6,10), 9])
i17 <- cbind(as.data.frame(dontmoveTRUTH[4:6, c(5,7,9)]), rep("dontupdateuTRUTH", 3),
             dontmoveTRUTH[4:6, 1])
i18 <- cbind(as.data.frame(badinit_ub[4:6, 1:3]), rep("badinit", 3),badinit_ub[4:6, 4])


colnames(i1) <- colnames(i2) <- colnames(i3) <- colnames(i11) <-colnames(i4)<-
  colnames(i12) <- colnames(i13) <- colnames(i14) <- colnames(i5) <- colnames(i15)<-colnames(i18) <-
  colnames(i6) <- colnames(i16) <- colnames(i17) <-colnames(i7) <-colnames(i8) <-   c("ll", "m", "hh", "model", "Community")

colnames(dc1) <- colnames(dc2) <- colnames(dc3) <-colnames(dc4) <-
  colnames(dr1) <- colnames(dr2) <- colnames(dr3) <- colnames(dr4) <-
  colnames(dr5) <- colnames(dr6) <- colnames(dc6) <- colnames(dc5) <-colnames(dc7) <-
  colnames(dr7) <- colnames(dc8) <-colnames(dr8) <- c("ll", "m", "hh", "model", "Community")

dep_col    <- as.data.frame(rbind(dc1, dc2, dc3, dc4, dc5, dc6, dc7,dc8))
dep_row    <- as.data.frame(rbind(dr1, dr2, dr3, dr4, dr5, dr6, dr7,dr8))

indep_col <- as.data.frame(rbind(i1, i2, i3, i4, i5, i6, i7, i8))
indep_row <- as.data.frame(rbind(i11, i12, i13,i14, i15, i16, i17, i18))

dep_col$trueBeta <- c(0,0,0, 0, -2,2, 0, -2, 2, 0, -2,2, 0, -2, 2, 0 , -2, 2, 0, -2, 2, 0, -2,2)
dep_row$trueBeta <- c(0,0,0, 1, 0, -1, 1,0, -1, 1, 0, -1,1, 0,-1, 1, 0, -1, 1, 0, -1,1,0,-1)
indep_col$trueBeta  <- c(2,2, 2, 2, 2, 2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
indep_row$trueBeta <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

## models key:
# a: R = 0 (R0)
# b. R = 3 (R3)
# c. AMEN with binary U (MCMCLAMBDAU)
# d. Oracle with known community structure (dontupdateutruth)
# e. INitialize U and dont change it (our model) (dontupdateU)
# f. run sep amen based on initial U (dontmovesep)
# g. oracel community but fit indep within each community (truesep)
# h. oracle appraoch (we know which are dep or not dep) goodinittell
# i. our approach where we let it fit comm. for all (goodinit)
dep_col$model <- as.factor(dep_col$model)
dep_row$model <- as.factor(dep_row$model)
indep_col$model <- as.factor(indep_col$model)
indep_row$model <- as.factor(indep_row$model)

levels(dep_col$model) <- c("bad", "d", "f", "g" , "i", "h", "c", "a", "b", "e")
dep_col$model <- as.character(dep_col$model)


levels(dep_row$model) <-c("bad", "d", "f", "g", "i", "h", "c", "a", "b", "e") #c("bad", "f", "e", "d", "i", "h", "c", "a", "b", "g")
dep_row$model <- as.character(dep_row$model)



levels(indep_col$model) <- c("bad", "d", "f", "g", "i", "h", "c", "a", "b", "e")
indep_col$model <- as.character(indep_col$model)



levels(indep_row$model) <- c("bad", "d", "f", "g", "i", "h", "c", "a", "b", "e")
indep_row$model <- as.character(indep_row$model)

dep_col <- subset(dep_col, model != "bad")
dep_row <- subset(dep_row, model != "bad")
indep_col <- subset(indep_col, model != "bad")
indep_row  <- subset(indep_row, model != "bad")

####### DEPENDENT COLUMN
dep_col_gg <- ggplot(subset(dep_col, Community != "No Comm"),
                     aes(x= model, y= m, col =Community))+ theme_bw()+
  ylim(-6.1, 4) + 
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ 
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T$ "))+
  #ggtitle("95% CI for Dependent Column Covariate")+
  geom_hline(size = 2,   aes(yintercept = trueBeta, col = Community)) +
  guides(col = FALSE)+  scale_color_manual(values = c("black", "dimgray", "gray"))+
  theme_bw()+ theme(axis.title.x =  element_text(size = 16), axis.text.x=element_text(size = 24,family="Courier"), 
                    # axis.text.x=element_text(size = 24), 
                    axis.title.y = element_text(size = 24),
                    axis.text.y = element_text(size = 24), 
                    title = element_text(size = 14),
                    strip.text = element_text(size = 16))



dep_col_gg1 <- ggplot(subset(dep_col, Community == "No Comm"),
                      aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ 
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\beta$"))+
  # ggtitle("95% CI for Dependent Column Covariate")+
  geom_hline(size = 2,  linetype = "dashed", aes(yintercept = trueBeta,
                                                 col = Community)) +   ylim(-6.1, 4) + 
  guides(col = FALSE)+  scale_color_manual(values = c("black", "dimgray", "gray"))+theme_bw()+ 
  theme(axis.title.x =  element_text(size = 16), axis.text.x=element_text(size = 24,family="Courier"), 
        #axis.text.x=element_text(size = 24), 
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 24), 
        title = element_text(size = 14),
        strip.text = element_text(size = 16))







grid.arrange(dep_col_gg1, dep_col_gg, nrow =1,top = "95% CI Dep. Column Covariate")

grid.arrange(dep_col_gg1, dep_col_gg, nrow =1,top = "")



dep_row_gg <- ggplot(subset(dep_row, Community != "No Comm"),
                     aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+  ylim(c(-5.5, 2.5))+
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T$ "))+
  # ggtitle("95% CI for Dependent Row Covariate")+
  geom_hline(size = 2, 
             aes(yintercept = trueBeta,
                 col = Community)) +
  guides(col = FALSE)+   scale_color_manual(values = c("black", "dimgray", "gray"))+theme_bw()+ 
  theme(  axis.title.x =  element_text(size = 16), 
          axis.text.x=element_text(size = 24,family="Courier"), 
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 24), 
          title = element_text(size = 14),
          strip.text = element_text(size = 16))


dep_row_gg1 <- ggplot(subset(dep_row, Community == "No Comm"),
                      aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ 
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\beta$ "))+
  # ggtitle("95% CI for Dependent Column Covariate")+
  geom_hline(size = 2, linetype = "dashed",
             aes(yintercept = trueBeta,
                 col = Community)) + ylim(c(-5,2.5))+
  guides(col = FALSE)+  scale_color_manual(values = c("black", "dimgray", "gray"))+theme_bw()+ 
  theme(  axis.title.x =  element_text(size = 16),   
          axis.text.x=element_text(size = 24, family="Courier"), 
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 24), 
          title = element_text(size = 14),
          strip.text = element_text(size = 16))




grid.arrange(dep_row_gg1,dep_row_gg,  nrow = 1, top = "95% CI Dep. Row Covariate")
grid.arrange(dep_row_gg1,dep_row_gg,  nrow = 1, top = "")

## Figure
#################### indep covariates
indep_col_gg1 <- ggplot(subset(indep_col, model %in% c('a','b','c','h')),
                        aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ ylim(c(-1, 3))+
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+
  ylab(TeX("$\\beta$ "))+
  #ggtitle("95% CI for Indep. Column Covariate")+
  geom_hline(size = 2, 
             aes(yintercept = trueBeta)) +
  guides(col = FALSE)+   scale_color_manual(values = c("black", "dimgray", "gray", "gray9"))+
  theme_bw()+ theme(axis.title.x =  element_text(size = 16), 
                    axis.text.x=element_text(size = 24,family = "Courier"), 
                    axis.title.y = element_text(size = 24),
                    axis.text.y = element_text(size = 24), 
                    title = element_text(size = 14),
                    strip.text = element_text(size = 16))


indep_col_gg2 <- ggplot(subset(indep_col, !(model %in% c('a','b','c','h'))),
                        aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ ylim(c(-1.2, 3))+
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+
  ylab(TeX("$\\tilde{\\beta} U^T$ "))+
  #ggtitle("95% CI for Indep. Column Covariate")+
  geom_hline(size = 2, 
             aes(yintercept = trueBeta)) +
  guides(col = FALSE)+   scale_color_manual(values = c("black", "dimgray", "gray", "gray9"))+theme_bw()+ theme(  axis.title.x =  element_text(size = 16), 
                                                                                                                 axis.text.x=element_text(size = 24, family = "Courier"), 
                                                                                                                 axis.title.y = element_text(size = 24),
                                                                                                                 axis.text.y = element_text(size = 24), 
                                                                                                                 title = element_text(size = 14),
                                                                                                                 strip.text = element_text(size = 16))


grid.arrange(indep_col_gg1, indep_col_gg2, nrow = 1)





indep_row_gg1 <- ggplot(subset(indep_row, model %in% c('a','b','c','h')),
                        aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ ylim(c(-1.5, 3))+
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+
  ylab(TeX("$\\beta$ "))+
  #ggtitle("95% CI for Indep. Column Covariate")+
  geom_hline(size = 2, 
             aes(yintercept = trueBeta)) +
  guides(col = FALSE)+   scale_color_manual(values = c("black", "dimgray", "gray", "gray9"))+theme_bw()+ theme(axis.title.x =  element_text(size = 16), 
                                                                                                               axis.text.x=element_text(size = 24, family = "Courier"), 
                                                                                                               axis.title.y = element_text(size = 24),
                                                                                                               axis.text.y = element_text(size = 24), 
                                                                                                               title = element_text(size = 14),
                                                                                                               strip.text = element_text(size = 16))


indep_row_gg2 <- ggplot(subset(indep_row, !(model %in% c('a','b','c','h'))),
                        aes(x= model, y= m, col =Community))+ theme_bw()+
  #+ facet_grid(Community~covar,scales = "free")+
  geom_errorbar(width = .8, size = 2.5, 
                aes(ymin= ll, ymax = hh), 
                position = position_dodge(width = .5))+ ylim(c(-1.5, 3))+
  geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+
  ylab(TeX("$\\tilde{\\beta} U^T$ "))+
  #ggtitle("95% CI for Indep. Column Covariate")+
  geom_hline(size = 2, 
             aes(yintercept = trueBeta)) +
  guides(col = FALSE)+   scale_color_manual(values = c("black", "dimgray", "gray", "gray9"))+theme_bw()+ theme(  axis.title.x =  element_text(size = 16), 
                                                                                                                 axis.text.x=element_text(size = 24, family = "Courier"), 
                                                                                                                 axis.title.y = element_text(size = 24),
                                                                                                                 axis.text.y = element_text(size = 24), 
                                                                                                                 title = element_text(size = 14),
                                                                                                                 strip.text = element_text(size = 16))


grid.arrange(indep_row_gg1, indep_row_gg2, nrow = 1)



## Figure 5
set.seed(1312)
iter <- 150000
burnin <- 5000
numGroup <- 3
UV_start <- specClust(data$Y, tauAvg = TRUE,  3 ,FALSE)
U <- V <- matrix(UV_start$U_hat, ncol = numGroup)

lambda <- diag(numGroup)#15*matrix(c(.4,.15,.15,
                      #.15,.4,.15,
                      #.15,.15,.4), nrow = 3, byrow = TRUE)
# spec initial
f1 <- matrixPlotter(U%*%t(U), "black", "white")

# Figure 6
# probit initial
prob_U <- probitUEstimation(c(data$Y), c(data$X_r[,,1]), c(data$X_c[,,1]),
                  c(data$X_r[,,2]),c((data$X_c[,,2])), lambda)
f2 <- matrixPlotter(prob_U, "black", "white")

# Figure 7 (true ULU)
U <- matrix(c(rep(c(1,0,0), 50), rep(c(0,1,0),50), rep(c(0,0,1),50)),
            ncol = 3,nrow = 150, byrow = TRUE)

f3 <- matrixPlotter(U %*%lambda %*%t(U), "black", "white")
grid.arrange(f3, f1, f2, nrow = 1)
# Figure 8 (timer plots)
rm(list = ls())
## simulation progression for timer plots

library(coda);library(pbivnorm)
library(bayesplot)
numGroup = 3

dataGenerator2_asymm_amen2 <- function(){
  set.seed(909)
  numGroup = 3
  # Priors
  prior_beta_mu   <- rep(0, numGroup)
  prior_beta_var  <- diag(1,numGroup)
  
  prior_beta_mu_indep   <- 0
  prior_beta_var_indep  <- 1
  
  prior_alpha_mu  <-0
  prior_alpha_var <- 1
  
  priorSet <- list(prior_beta_mu = prior_beta_mu, prior_beta_var = prior_beta_var, 
                   prior_beta_mu_indep = prior_beta_mu_indep, prior_beta_var_indep = prior_beta_var_indep,
                   prior_alpha_mu = prior_alpha_mu, prior_alpha_var = prior_alpha_var)
  
  # Initials
  start_beta_c <- start_beta_r <- matrix(rep(1,3), nrow = 1, ncol = numGroup)
  
  
  start_beta_c_indep <- start_beta_r_indep <-c(2,2)
  
  start_alpha <- -6
  
  initialSet <- list(start_beta_c = start_beta_c, start_beta_r = start_beta_r, start_beta_c_indep=
                       start_beta_c_indep, start_beta_r_indep = start_beta_r_indep, start_alpha = start_alpha)
  
  set.seed(51)
  n <-150
  ones <- as.matrix(rep(1,n), ncol =1 )
  ## Kronecker extravagenzaaaaa
  
  ## create x and row covar for one beta
  U <- matrix(c(rep(c(1,0,0), 50), rep(c(0,1,0),50), rep(c(0,0,1),50)),
              ncol = 3,nrow = 150, byrow = TRUE)
  
  
  
  dataAll1 <- rnorm(150, 0, 1)
  
  
  
  
  dataAll2 <- rnorm(150, 0, 1)
  
  
  
  
  
  
  
  
  X_row1 <- matrix(rep(dataAll1, 150), nrow = 150, byrow = FALSE)
  X_col1 <- t(X_row1)
  
  
  X_row2 <- matrix(rep(dataAll2, 150), nrow = 150, byrow = FALSE)
  X_col2 <- t(X_row2)
  
  ### generate a dep. beta 
  # beta_r <-c(1,5,10)
  # beta_c <-c(2,6,12)
  
  beta_r2 <- c(1,1,1)
  beta_c2 <- rep(1.5, 3)
  
  beta_r1 <- c(1,1,1)
  beta_c1 <- rep(.5,3)
  
  
  
  library(amen)
  epsilon <- amen::rmvnorm(n, rep(0,n), diag(.5,n))
  
  # makeSymm <- function(m) {
  #   m[lower.tri(m)] <- t(m)[lower.tri(m)]
  #   return(m)
  # }
  # 
  # epsilon <- makeSymm(epsilon)
  
  lambda <- 15*matrix(c(.4,.15,.15,
                        .15,.4,.15,
                        .15,.15,.4), nrow = 3, byrow = TRUE)
  
  UbT_r1 <- U %*% beta_r1
  ## lets look at beta U^T
  dUbT_r1 <- diag(as.vector(U %*% beta_r1))
  UbT_r2 <- U %*% beta_r2
  ## lets look at beta U^T
  dUbT_r2 <- diag(as.vector(U %*% beta_r2))
  
  
  
  UbT_c1 <- U %*% beta_c1
  ## lets look at beta U^T
  dUbT_c1 <- diag(as.vector(U %*% beta_c1))
  
  UbT_c2 <- U %*% beta_c2
  ## lets look at beta U^T
  dUbT_c2 <- diag(as.vector(U %*% beta_c2))
  
  ## multiply by X_row to see if we get what we want
  rowXB1 <- dUbT_r1 %*%X_row1
  colXB1 <- X_col1 %*%dUbT_c1
  rowXB2 <- dUbT_r2 %*%X_row2
  colXB2 <- X_col2 %*%dUbT_c2
  
  XB_check <- rowXB1 + colXB1+rowXB2 + colXB2
  
  ## okay so lets play with it
  basevecr <- c(rowXB2)
  basevecc <- c(colXB2)
  
  rewrite <- (ones %x% ((diag(1,n) * X_row2)%*%U))%*%c(beta_r2)
  rewrite2 <- (((diag(1,n) * X_col2)%*%U)%x%ones)%*%c(beta_c2)
  
  check <- cbind(basevecr, rewrite,basevecc, rewrite2)
  colnames(check) <- c("True Row", "Rewrite Row", "True Col", "Rewrite Col")
  
  a <- rnorm(n, 0, 1)
  b <- rnorm(n, 0, 1)
  
  pairs <- amen::rmvnorm(n, c(0,0), Sigma = matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2))
  
  
  ab_couple <-matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      ab_couple[i,j] <- pairs[i,1] + pairs[j,2]
    }
  }
  
  
  # make epsilon matrix
  EPS <- matrix(0, nrow = n, ncol = n)
  
  for(i in 1:n){
    
    for(j in 1:n){
      epsilonPair <- amen::rmvnorm(1, c(0,0), Sigma = matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2))
      
      EPS[i,j] <- epsilonPair[1]
      EPS[j,i] <- epsilonPair[2]
      
    }
  }
  
  #######
  Z = colXB1 + rowXB1 + colXB2 + rowXB2 + EPS + (U %*% lambda %*% t(U)) + ab_couple
  
  EZ = colXB1 + rowXB1 + colXB2 + rowXB2 + (U %*% lambda %*% t(U)) + ab_couple
  Zv <- c(Z)
  EZv <- c(EZ)
  
  epsV <- Zv-EZv
  decorrZ <- .9*Z + .9*t(Z)
  decorrEZ <- .9*EZ + .9*t(EZ)
  decorrZv <- c(decorrZ)
  decorrEZv <- c(decorrEZ)
  
  decorrEps <- decorrZv-decorrEZv
  Phi = dnorm(Zv-EZv, mean=0, sd = 1, log = TRUE)
  PhiD = dnorm(decorrZv, mean=decorrEZv)
  
  ## simualte what happens with a switch
  
  
  #######
  Y = ifelse(Z>6, 1,0 )
  sum(Y)/150^2
  
  X_r <- array(dim = c(n, n, 2))
  X_c <- array(dim = c(n, n, 2))
  X_r[,,1] <- matrix(rep(dataAll1, n), nrow = n, byrow = FALSE)
  X_c[,,1] <- t(matrix(rep(dataAll1, n), nrow = n, byrow = FALSE))
  X_r[,,2] <-matrix(rep(dataAll2, n), nrow = n, byrow = FALSE)
  X_c[,,2]<- t(matrix(rep(dataAll2, n), nrow = n, byrow = FALSE))
  
  
  
  diag(Y) <- 0
  return(list(Y = Y, X_r = X_r, X_c = X_c, trueBeta_r1 = beta_r1, trueBeta_c1 = beta_c1, trueBeta_c2 = beta_c2, trueBeta_r2 = beta_r2,
              df = cbind(dataAll1, dataAll2), initialSet = initialSet, priorSet = priorSet))
}



library(amen)
library(tictoc)
library(bayesplot)
# load data
# data set 2

#equal var and mean
data_setAMEN <- dataGenerator2_asymm_amen2()

### RUN PARALLEL JOBS 1-5000


slurmKeep <- # make job number  as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

  
  i = slurmKeep %% 10

if(i == 0){i = 10}


set.seed(slurmKeep*54 + 121210)


## R = 3

#####################################################################
es.mat <- as.data.frame(matrix(0, nrow = 1, ncol = 23))
colnames(es.mat) <- c(paste0("es", seq_len(5)), "Iter","Time", "Model", 
                      paste0("betas", seq_len(5)),paste0("betas_l", seq_len(5)),paste0("betas_u", seq_len(5)))
#rownames(es.mat) <- seq(1,20, by = 1)
if(slurmKeep <= 1000){
  iter      <- i*2000
  tic()
  amenTemp  <- ame(data_setAMEN$Y, Xdyad = NULL, Xrow = data_setAMEN$df ,Xcol = data_setAMEN$df, dcor = TRUE,
                   R=3, cvar = TRUE,rvar = TRUE,  symmetric = FALSE, model = "bin",
                   burn = 1000,odens = 10, nscan = iter, plot = FALSE, seed =slurmKeep*54 + 121210)
  
  tempTime <- toc()
  tempTime2 <- tempTime$toc - tempTime$tic
  attr(tempTime2, "elapsed") <- NULL
  es.mat[1, "Time"] <- as.numeric(tempTime2)
  
  es.mat[1, paste0("es", seq_len(5))] <- effectiveSize(as.mcmc(amenTemp$BETA))
  
  
  es.mat[1, "Iter"] <- iter
  es.mat[1, "Model"] <- "R_3"
  
  ## get posterior mean
  es.mat[1, paste0("betas", seq_len(5))] <- colMeans(amenTemp$BETA)
  es.mat[1, paste0("betas_l", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$ll
  es.mat[1, paste0("betas_u", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$hh
  
  colnames(es.mat)[9:13] <- c( "intercept", "beta_r1", "beta_r2", "beta_c1", "beta_c2")
}



### R = 0 amen with dyadic covariate as UTU
#####################################################################
UV_start    <- specClust(data_setAMEN$Y, tauAvg = TRUE,  numGroup,FALSE)
cluster_symm     <-  matrix(UV_start$U_hat, ncol = numGroup)
utu <- cluster_symm %*% t(cluster_symm)
if(slurmKeep >1000 && slurmKeep <=2000){
  iter      <- i*2000
  tic()
  
  amenTemp  <-ame(data_setAMEN$Y, Xdyad = utu, Xrow =data_setAMEN$df , Xcol = data_setAMEN$df, dcor = TRUE,
                  R=0, rvar = TRUE, cvar = TRUE, nscan = iter, 
                  symmetric = FALSE, model = "bin",burn = 1000, odens = 10, plot = FALSE, seed =slurmKeep*54 + 121210)
  
  
  
  
  tempTime <- toc()
  tempTime2 <- tempTime$toc - tempTime$tic
  attr(tempTime2, "elapsed") <- NULL
  es.mat[1, "Time"] <- as.numeric(tempTime2)
  
  es.mat[1, paste0("es", seq_len(5))] <- effectiveSize(as.mcmc(amenTemp$BETA))[-6]
  
  
  es.mat[1, "Iter"] <- iter
  es.mat[1, "Model"] <- "R_Dyad"
  
  ## get posterior mean
  es.mat[1, paste0("betas", seq_len(5))] <- colMeans(amenTemp$BETA)[-6]
  es.mat[1, paste0("betas_l", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$ll[-1]
  es.mat[1, paste0("betas_u", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$hh[-1]
  
  colnames(es.mat)[9:13] <- c( "intercept", "beta_r1", "beta_r2", "beta_c1", "beta_c2")
  
}


#####################################
######## our amen with lambda 

#####################################################################
if(slurmKeep >2000 && slurmKeep <=3000){
  
  iter      <- i*2000
  tic()
  amenTemp  <-  myAME_estimateLambda(data_setAMEN$Y, Xdyad = NULL, Xrow = data_setAMEN$df, Xcol = data_setAMEN$df, dcor = TRUE,
                                     R=3, rvar = TRUE, cvar = TRUE, nscan =iter, symmetric = FALSE,
                                     model = "bin",burn = 1000, initialU = cluster_symm,
                                     initialV = cluster_symm,odens = 10, plot = FALSE, seed =slurmKeep*54 + 121210)
  
  tempTime <- toc()
  tempTime2 <- tempTime$toc - tempTime$tic
  attr(tempTime2, "elapsed") <- NULL
  es.mat[1, "Time"] <- as.numeric(tempTime2)
  es.mat[1, paste0("es", seq_len(5))] <- effectiveSize(as.mcmc(amenTemp$BETA))
  
  
  es.mat[1, "Iter"] <- iter
  es.mat[1, "Model"] <- "R_MCMC_Lambda"
  
  ## get posterior mean
  es.mat[1, paste0("betas", seq_len(5))] <- colMeans(amenTemp$BETA)
  es.mat[1, paste0("betas_l", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$ll
  es.mat[1, paste0("betas_u", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$hh
  
  colnames(es.mat)[9:13] <- c( "intercept", "beta_r1", "beta_r2", "beta_c1", "beta_c2")
  
}


###########################################################
## standard amen with R = 0

#####################################################################
if(slurmKeep >3000 && slurmKeep <=4000){
  
  iter      <- i*2000
  tic()
  amenTemp  <-ame(data_setAMEN$Y, Xdyad = NULL, Xrow =data_setAMEN$df , Xcol = data_setAMEN$df, dcor = TRUE,
                  R=0, rvar = TRUE, cvar = TRUE, nscan = iter, 
                  symmetric = FALSE, model = "bin",burn = 1000, odens = 10, plot = FALSE,seed =slurmKeep*54 + 121210)
  
  
  
  tempTime <- toc()
  tempTime2 <- tempTime$toc - tempTime$tic
  attr(tempTime2, "elapsed") <- NULL
  es.mat[1, "Time"] <- as.numeric(tempTime2)
  es.mat[1, paste0("es", seq_len(5))] <- effectiveSize(as.mcmc(amenTemp$BETA))
  
  
  es.mat[1, "Iter"] <- iter
  es.mat[1, "Model"] <- "R_0"
  
  ## get posterior mean
  es.mat[1, paste0("betas", seq_len(5))] <- colMeans(amenTemp$BETA)
  es.mat[1, paste0("betas_l", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$ll
  es.mat[1, paste0("betas_u", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$hh
  
  colnames(es.mat)[9:13] <- c( "intercept", "beta_r1", "beta_r2", "beta_c1", "beta_c2")
  
}

#### Our amen no lambda
#####################################################################

if(slurmKeep >4000 && slurmKeep <=5000){
  
  iter      <- i*2000
  tic()
  amenTemp  <-myAME(data_setAMEN$Y, Xdyad = NULL, Xrow = data_setAMEN$df, Xcol = data_setAMEN$df, dcor = TRUE,
                    R=3, rvar = TRUE, cvar = TRUE, nscan =iter, symmetric = FALSE,
                    model = "bin",burn = 1000, initialU = cluster_symm,
                    initialV = cluster_symm,odens = 10, plot = FALSE,seed =slurmKeep*54 + 121210)
  
  
  
  tempTime <- toc()
  tempTime2 <- tempTime$toc - tempTime$tic
  attr(tempTime2, "elapsed") <- NULL
  es.mat[1, "Time"] <- as.numeric(tempTime2)
  es.mat[1, paste0("es", seq_len(5))] <- effectiveSize(as.mcmc(amenTemp$BETA))
  
  
  es.mat[1, "Iter"] <- iter
  es.mat[1, "Model"] <- "R_MCMC_noLamb"
  
  ## get posterior mean
  es.mat[1, paste0("betas", seq_len(5))] <- colMeans(amenTemp$BETA)
  es.mat[1, paste0("betas_l", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$ll
  es.mat[1, paste0("betas_u", seq_len(5))] <-  mcmc_intervals_data(as.mcmc(amenTemp$BETA), prob_outer = .95, point_est = "mean")$hh
  
  colnames(es.mat)[9:13] <- c( "intercept", "beta_r1", "beta_r2", "beta_c1", "beta_c2")
  
}







iterationNum = i

es.mat <- as.data.frame(es.mat)

Model <- es.mat$Model
saveRDS(es.mat, sprintf("%s%s_%s_%s.rds","~dir/ess_samen",
                        Model, iterationNum,slurmKeep))

# File to collect and plot timer ess/mse for showing computational advantages
library(reshape)
library(dplyr)

# Collect Data
es.finalMERGE <- c()


for(i in 1:10){
  for(j in 1:1000){
    if (file.exists(sprintf("%s%s_%s.rds","~dir/ess_samenR_3_", i,j))){
      temp <-readRDS(sprintf("%s%s_%s.rds","~dir/ess_samenR_3_", i,j))
    }else{
      temp <- NULL
    }
    es.finalMERGE<- rbind(es.finalMERGE, temp)
  }
  
}


for(i in 1:10){
  for(j in 1001:2000){
    if (file.exists(sprintf("%s%s_%s.rds","~dir/ess_samenR_DYAD_", i,j))){
      temp <-readRDS(sprintf("%s%s_%s.rds","~dir/ess_samenR_DYAD_", i,j))
    }else{
      temp <- NULL
    }
    es.finalMERGE<- rbind(es.finalMERGE, temp)
  }
  
}

# 

for(i in 1:10){
  for(j in 2001:3000){
    if (file.exists(sprintf("%s%s_%s.rds","~dir/ess_samenR_MCMC_Lambda_", i,j))){
      temp <-readRDS(sprintf("%s%s_%s.rds","~dir/ess_samenR_MCMC_Lambda_", i,j))
    }else{
      temp <- NULL
    }
    es.finalMERGE<- rbind(es.finalMERGE, temp)
  }
  
}

for(i in 1:10){
  for(j in 3001:4000){
    if (file.exists(sprintf("%s%s_%s.rds","~dir/ess_samenR_0_", i,j))){
      temp <-readRDS(sprintf("%s%s_%s.rds","~dir/ess_samenR_0_", i,j))
    }else{
      temp <- NULL
    }
    es.finalMERGE<- rbind(es.finalMERGE, temp)
  }
  
}
# 



for(i in 1:10){
  for(j in 4001:5000){
    if (file.exists(sprintf("%s%s_%s.rds","~dir/ess_samenR_MCMC_noLamb_", i,j))){
      temp <-readRDS(sprintf("%s%s_%s.rds","~dir/ess_samenR_MCMC_noLamb_", i,j))
    }else{
      temp <- NULL
    }
    es.finalMERGE<- rbind(es.finalMERGE, temp)
  }
  
}


es.finalMERGE<- as.data.frame( es.finalMERGE)

# calc. error term log(ess/(posterior mean - true value)^2)
esSummary <- es.finalMERGE %>% mutate(errorC1 = log(es4/(beta_c1 -.5)^2),
                                      errorR1 = log(es2/(beta_r1-1)^2) ,
                                      errorC2 = log(es5/(beta_c2 -.5)^2),
                                      errorR2 = log(es3/(beta_r2 -1)^2))


# get averagee error per model and number of iterations

esALL <-esSummary %>% group_by(Model, Iter) %>%
  summarise_all(funs(mean(.)))


# make df
esALL <- as.data.frame(esALL)


# melt for plottinng, get rid of unnecessary columns
meltedES <- melt(esALL, id = c("Model", "Iter"))
meltedES <- melt(esALL[,c(1,2, 24:27)], id = c("Model", "Iter"))



# create labels for coef. with tex font
esLabeller2 <- as_labeller(c(
  errorR1 = paste(expression(beta["1,r"])),
  errorR2 = paste(expression(beta["2,r"])),
  errorC1 = paste(expression(beta["1,c"])),
  errorC2 = paste(expression(beta["2,c"]))),label_parsed)



#meltMerge <- melt(esSummaryTIME, id = c("Iter", "Model", "Time"))
meltMerge <- meltedES
meltMerge$Model <- as.factor(meltMerge$Model)
meltMerge_noDyad <- subset(meltMerge, Model !="R_Dyad" )

## Create plot

ggplot(subset(meltMerge_noDyad, Model != "R_MCMC_noLamb"), aes(Iter, value, col = Model, shape = Model)) +
  #geom_smooth(se=FALSE, size =3, aes(linetype = Model)) +
  geom_point(size = 3)+ 
  facet_wrap(~variable, nrow = 2,labeller = esLabeller2)+ 
  scale_shape_manual("",values = c(15, 16, 17),labels = unname(c("AMEN: R = 0", "AMEN: R = 3",   #"AMEN: R = 0, Dyadic Covar", 
                                                                 TeX("Our Model")))) +
  theme_bw()+
  ylab(TeX("$log(\\frac{ESS}{\\beta_{err}})$"))+
  #scale_color_manual("",  values = c("midnightblue", "orange2", "gray50", "cyan3", "magenta"), 
  scale_color_manual("",  values = c("black", "gray27", 
                                     # "gray35",
                                     "gray58", "gray70"), 
                     
                     labels = unname(c("AMEN: R = 0", "AMEN: R = 3", 
                                       #"AMEN: R = 0, Dyadic Covar", 
                                       TeX("Our Model"))))+
  theme(
    axis.title.x =  element_text(size = 22),
    axis.text.x= element_text(size = 16),
    axis.title.y = element_text(
      size = 22),
    legend.text = element_text(size = 16),
    title = element_text(size = 16),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 22),
    legend.position = "bottom")

## Figure 9

amen_dep_03_badinit_progressionExperiment <- amen_master(data_dep_zero$X_r, 
                                                         data_dep_zero$X_c, 
                                                         NULL, 
                                                         data_dep_zero$Y,
                                                         iter = 15000,numGroup=3,
                                                         prior_beta_mu= data_prior_0$prior_beta_mu,
                                                         prior_beta_var = data_prior_0$prior_beta_var,
                                                         prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                                         prior_alpha_var= data_prior_0$prior_alpha_var,
                                                         start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                                         start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                                         start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                                         start_alpha= data_initial_0$start_alpha,FALSE,
                                                         TRUE, keep.UV = TRUE, dcor =TRUE,
                                                         symmetric = FALSE,
                                                         model = "bin", odmax = 0,indepBeta=1,"all",
                                                         odens = 10, 0, TRUE)



amen_dep_03_goodinit_progressionExperiment <- amen_master(data_dep_zero$X_r, 
                                                          data_dep_zero$X_c, 
                                                          NULL, 
                                                          data_dep_zero$Y,
                                                          iter = 500,numGroup=3,
                                                          prior_beta_mu= data_prior_0$prior_beta_mu,
                                                          prior_beta_var = data_prior_0$prior_beta_var,
                                                          prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                                          prior_alpha_var= data_prior_0$prior_alpha_var,
                                                          start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                                                          start_beta_r = matrix(rep(1,3), nrow = 2, ncol = 3),
                                                          start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                                          start_alpha= data_initial_0$start_alpha,FALSE,
                                                          TRUE, keep.UV = TRUE, dcor =TRUE,
                                                          symmetric = FALSE,
                                                          model = "bin", odmax = 0,indepBeta=1,"all",
                                                          odens = 10, 0, FALSE)

round1 <-matrix(amen_dep_03_badinit_progressionExperiment$results$U[1,], nrow = 150, ncol = 3)
m1<-matrixPlotter(round1 %*%t(round1), title = "Random Initialization",col2 = "grey",col1 = "black")

round1g <-matrix(amen_dep_03_goodinit_progressionExperiment$results$U[1,], nrow = 150, ncol = 3)
m1g<-matrixPlotter(round1g %*%t(round1g), title = "Smart Initialization", col2 = "grey",col1 = "black")

grid.arrange(m1, m1g, nrow = 1)

## Figure 10
plot(ub_bad_df$V1)
plot(ub_good_df$V1)
acf_bad <- acf(ub_bad_df$V1, plot = FALSE, lag = 100)
acf_good <- acf(ub_good_df$V1, plot = FALSE, lag = 100)

par(mfrow = c(1,2))
plot(acf_bad, main = "Random Initialization ACF")
plot(acf_good, main = "Smart Initialization ACF")


# Figure 11

dataGenerator_censored_binary_lik <-function(){
  
  set.seed(199)
  numGroup <- 3
  
  # priors for beta and alpha
  prior_beta_mu   <- rep(0, numGroup)
  prior_beta_var  <- diag(1,numGroup)
  
  prior_beta_mu_indep   <- 0
  prior_beta_var_indep  <- 1
  
  prior_alpha_mu  <-0
  prior_alpha_var <- 1
  
  # initialize
  start_beta_c <- start_beta_r <- matrix(rep(1,6), nrow = 2, ncol = numGroup)
  start_beta_c <- start_beta_r <- matrix(rep(-1,6), nrow = 2, ncol = numGroup)
  
  start_beta_c_indep <- start_beta_r_indep <-c(2,2)
  
  start_alpha <- -6
  
  
  priorSet <- list(prior_beta_mu = prior_beta_mu, prior_beta_var = prior_beta_var, 
                   prior_beta_mu_indep = prior_beta_mu_indep, prior_beta_var_indep = prior_beta_var_indep,
                   prior_alpha_mu = prior_alpha_mu, prior_alpha_var = prior_alpha_var)
  
  # Initials
  
  initialSet <- list(start_beta_c = start_beta_c, start_beta_r = start_beta_r, start_beta_c_indep=
                       start_beta_c_indep, start_beta_r_indep = start_beta_r_indep, start_alpha = start_alpha)
  p <- c(1/3,1/3,1/3)
  
  maintainAllOrder <- TRUE
  maxDeg <-10
  set.seed(1951)
  n <- 150
  ones <- as.matrix(rep(1,n), ncol =1 )
  ## create x and row covar for one beta
  U <- matrix(c(rep(c(1,0,0), n*p[1]), rep(c(0,1,0),n*p[2]), rep(c(0,0,1),n*p[3])),
              ncol = 3,nrow = n, byrow = TRUE)
  
  
  ## make 0 centered
  
  data1  <-rnorm(n*p[1], 0, 1)
  data2  <-rnorm(n*p[2], 0,1)
  data3  <-rnorm(n*p[3], 0, 1)
  dataAll1 <- c(data1,data2,data3)
  data1  <-rnorm(n*p[1],0, 1)
  data2  <-rnorm(n*p[2], 0,1)
  data3  <-rnorm(n*p[3], 0,1)
  dataAll2 <- c(data1,data2,data3)
  X_row1 <- matrix(rep(dataAll1, n), nrow = n, byrow = FALSE)
  X_col1 <- t(X_row1)
  X_row2 <- matrix(rep(dataAll2, n), nrow = n, byrow = FALSE)
  X_col2 <- t(X_row2)
  # set beta
  beta_r2 <- c(1,0,-1)/2
  beta_c2 <- c(0,-1,1)*1.5
  beta_r1 <- c(1,1,1)*-1
  beta_c1 <- c(1,1,1)
  lambda <- 20*matrix(c(.2,.1,.1,
                        .1,.2,.1,
                        .1,.1,.2), nrow = 3, byrow = TRUE)
  
  UbT_r1 <- U %*% beta_r1
  ## lets look at beta U^T
  dUbT_r1 <- diag(as.vector(U %*% beta_r1))
  UbT_r2 <- U %*% beta_r2
  ## lets look at beta U^T
  dUbT_r2 <- diag(as.vector(U %*% beta_r2))
  UbT_c1 <- U %*% beta_c1
  ## lets look at beta U^T
  dUbT_c1 <- diag(as.vector(U %*% beta_c1))
  UbT_c2 <- U %*% beta_c2
  ## lets look at beta U^T
  dUbT_c2 <- diag(as.vector(U %*% beta_c2))
  ## multiply by X_row to see if we get what we want
  rowXB1 <- dUbT_r1 %*%X_row1
  colXB1 <- X_col1 %*%dUbT_c1
  rowXB2 <- dUbT_r2 %*%X_row2
  colXB2 <- X_col2 %*%dUbT_c2
  XB_check <- rowXB1 + colXB1+rowXB2 + colXB2

  a <- rnorm(n, 0, 1)
  b <- rnorm(n, 0, 1)
  
  pairs <- amen::rmvnorm(n, c(0,0), Sigma = matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2))
  
  ab_couple <-matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      ab_couple[i,j] <- pairs[i,1] + pairs[j,2]
    }
  }
  
  
  # make epsilon matrix
  EPS <- matrix(0, nrow = n, ncol = n)
  
  for(i in 1:n){
    
    for(j in 1:n){
      epsilonPair <- amen::rmvnorm(1, c(0,0), Sigma = matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2))
      
      EPS[i,j] <- epsilonPair[1]
      EPS[j,i] <- epsilonPair[2]
      
    }
  }
  
  #######
  Z = colXB1 + rowXB1 + colXB2 + rowXB2 + EPS + (U %*% lambda %*% t(U)) + ab_couple

 
  X_r <- array(dim = c(n, n, 2))
  X_c <- array(dim = c(n, n, 2))
  X_r[,,1] <- matrix(rep(dataAll1, n), nrow = n, byrow = FALSE)
  X_c[,,1] <- t(matrix(rep(dataAll1, n), nrow = n, byrow = FALSE))
  X_r[,,2] <-matrix(rep(dataAll2, n), nrow = n, byrow = FALSE)
  X_c[,,2]<- t(matrix(rep(dataAll2, n), nrow = n, byrow = FALSE))
  
  
  odz = 5
  qz<-quantile(Z,1-odz/n) 
  odc <- 15
  Z<-Z-qz 

  Y <- Z*0
  for(i in 1:n){
    currRow <- Z[i,]
    orderCurr <- order(currRow, decreasing = TRUE)[1:odc]
    orderCurr_pos <- orderCurr[which(Z[i,orderCurr]>0)]
    Y[i, orderCurr_pos] <- 1
    
  }
  
  
  # top 
  trueRho <- .9
  trueAlpha <- -qz
  trueU <- U
  trueSab <-matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2)
  trueLambda <- lambda
  trueA <- pairs[,1]
  trueB <- pairs[,2]
  return(list(Y = Y, X_r = X_r, X_c = X_c, trueBeta_r1 = beta_r1, trueBeta_c1 = beta_c1, trueBeta_c2 = beta_c2, trueBeta_r2 = beta_r2,
              df = cbind(dataAll1, dataAll2), priorSet = priorSet, initialSet = initialSet, 
              trueAlpha = trueAlpha, trueU = trueU, trueLambda = trueLambda,
              trueA = trueA, trueB = trueB, trueRho = trueRho, trueSab = trueSab, trueZ = Z))
}

data_dep_zero <- dataGenerator_censored_binary_lik()

Ycbin <- 1*(data_dep_zero$Y >0)
data_prior_0     <- data_dep_zero$priorSet
data_initial_0   <- data_dep_zero$initialSet
set.seed(18983)

amen_dep_03_cbin<- amen_master(data_dep_zero$X_r,
                                data_dep_zero$X_c,
                                NULL,
                                Ycbin,
                                iter = 350000,numGroup=3,
                                prior_beta_mu= data_prior_0$prior_beta_mu,
                                prior_beta_var = data_prior_0$prior_beta_var,
                                prior_alpha_mu= data_prior_0$prior_alpha_mu,
                                prior_alpha_var= data_prior_0$prior_alpha_var,
                                start_beta_c= matrix(rep(0,6), nrow = 2, ncol = 3),#was 1 now 0
                                start_beta_r = matrix(rep(0,3), nrow = 2, ncol = 3),
                                start_beta_dyad =  matrix(rep(0,3), nrow = 2, ncol = 3),
                                start_alpha= 0,#data_initial_0$start_alpha,# changed
                                FALSE,
                                TRUE, keep.UV = TRUE, dcor =TRUE,
                                symmetric = FALSE,
                                model = "cbin", odmax = max(rowSums(Ycbin)),indepBeta=1,"all",
                                odens = 10, 
                                5000,
                                FALSE, rep(1/3,3))






ubeta_posterCI_cluster_intermediate_saves(amen_dep_03_cbin$results, beta_c1, beta_c2, 
                                          beta_r1, beta_r2,n=150,
                                          burnin=10, indepDep = FALSE,
                                          prob_choice = .95)
