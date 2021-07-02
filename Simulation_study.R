rm(list=ls())
library(fdapace)
library(ggplot2)
library(plot3D)
source('sim_functions.R')
set.seed(11235)
#Number of observation points 
ntimes     <- 101
#Define observation points in the interval [0,1]
timepoints <- seq(0,1,length.out = ntimes)

#Pre-allocate space for storing functional principal components 
#For both dimensions
dim3_fpc <- dim2_fpc <- dim1_fpc <- matrix(NA, nrow = ntimes, ncol = 3)

#First PC are constant for both dimensions 
dim2_fpc[,3] <- rep(0.1, ntimes)
dim1_fpc[,3] <- rep(0.1, ntimes)
dim3_fpc[,3] <- rep(0.1, ntimes)

#Trig basis functions for the first functional variable 
dim1_fpc[,1] <- unlist(lapply(timepoints,function(x){sqrt(2) * sin(.75* pi * x)}))
dim1_fpc[,2] <- unlist(lapply(timepoints,function(x){sqrt(2) * cos(.75* pi * x)}))
#Trig basis for the second functional variable 
dim2_fpc[,1] <- unlist(lapply(timepoints,function(x){sqrt(2) * sin(1.25 * pi * x) }))
dim2_fpc[,2] <- unlist(lapply(timepoints,function(x){sqrt(2) * -cos(1.25 * pi * x)}))
#Trig basis for the second functional variable 
dim3_fpc[,1] <- unlist(lapply(timepoints,function(x){sqrt(2) * cos(pi * x)}))
dim3_fpc[,2] <- unlist(lapply(timepoints,function(x){sqrt(2) * sin(pi * x)}))


dim1_fpc <- apply(dim1_fpc, 2, function(x){x / sqrt(sum(x^2))})
dim2_fpc <- apply(dim2_fpc, 2, function(x){x / sqrt(sum(x^2))})
dim3_fpc <- apply(dim3_fpc, 2, function(x){x / sqrt(sum(x^2))})

#Define the eigenvalues for K-L representation
dim1_eigenvalues <- c(0.75,0.25,0.05)
dim2_eigenvalues <- c(0.7,0.2,0.1)
dim3_eigenvalues <- c(0.6,0.3,0.1)

#Using tensor product to get the covariance function for X1
C_11 <- list()
for (ii in 1:3){
  C_11[[ii]] <- dim1_eigenvalues[ii] * dim1_fpc[,ii] %*% t(dim1_fpc[,ii])
}
Cov_11 <- Reduce('+',C_11)

#Using tensor product to get the covariance function for X2
C_22 <- list()
for (ii in 1:3){
  C_22[[ii]] <- dim2_eigenvalues[ii] * dim2_fpc[,ii] %*% t(dim2_fpc[,ii])
}
Cov_22 <- Reduce('+',C_22)


#Using tensor product to get the covariance function for X3
C_33 <- list()
for (ii in 1:3){
  C_33[[ii]] <- dim3_eigenvalues[ii] * dim3_fpc[,ii] %*% t(dim3_fpc[,ii])
}
Cov_33 <- Reduce('+',C_33)

#Obtain correlation surface based on the covairance function s
corr_11 <- corr_22 <- corr_33 <- matrix(NA, nrow = ntimes, ncol = ntimes)
for (ii in 1:ntimes){
  for (jj in 1:ntimes){
    corr_11[ii,jj] <- Cov_11[ii,jj] / (sqrt(Cov_11[ii,ii]) * sqrt(Cov_11[jj,jj]))
    corr_22[ii,jj] <- Cov_22[ii,jj] / (sqrt(Cov_22[ii,ii]) * sqrt(Cov_22[jj,jj]))
    corr_33[ii,jj] <- Cov_33[ii,jj] / (sqrt(Cov_33[ii,ii]) * sqrt(Cov_33[jj,jj]))
  }
}

#Get eigenvectors and values of the correlation function
corr_eigenbasis_1 <- eigen(corr_11)$vectors[,1:3]
corr_eigenvalue_1 <- eigen(corr_11)$values[1:3]
corr_eigenbasis_2 <- eigen(corr_22)$vectors[,1:3]
corr_eigenvalue_2 <- eigen(corr_22)$values[1:3]
corr_eigenbasis_3 <- eigen(corr_33)$vectors[,1:3]
corr_eigenvalue_3 <- eigen(corr_33)$values[1:3]
#Fill-in the off diagonal elements of the correlation function
corr_12 <- corr_13 <- corr_21 <- corr_23 <- corr_31<-corr_32<- matrix(NA, nrow = ntimes, ncol = ntimes)
mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_1[ii],corr_eigenvalue_2[ii])) * corr_eigenbasis_1[,ii] %*% t(corr_eigenbasis_2[,ii])
}
corr_12 <- Reduce('+',mat.list)
mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_1[ii],corr_eigenvalue_3[ii])) * corr_eigenbasis_1[,ii] %*% t(corr_eigenbasis_3[,ii])
}
corr_13 <- Reduce('+',mat.list)

mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_2[ii],corr_eigenvalue_1[ii])) * corr_eigenbasis_2[,ii] %*% t(corr_eigenbasis_1[,ii])
}
corr_21 <- Reduce('+',mat.list)
mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_2[ii],corr_eigenvalue_3[ii])) * corr_eigenbasis_2[,ii] %*% t(corr_eigenbasis_3[,ii])
}
corr_23 <- Reduce('+',mat.list)

mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_3[ii],corr_eigenvalue_1[ii])) * corr_eigenbasis_3[,ii] %*% t(corr_eigenbasis_1[,ii])
}
corr_31 <- Reduce('+',mat.list)
mat.list <- list()
for (ii in 1:3){
  mat.list[[ii]] <-  mean(c(corr_eigenvalue_3[ii],corr_eigenvalue_2[ii])) * corr_eigenbasis_3[,ii] %*% t(corr_eigenbasis_2[,ii])
}
corr_32 <- Reduce('+',mat.list)




Corr_mat <- rbind(cbind(corr_11,corr_12,corr_13),
                  cbind(corr_21,corr_22,corr_23),
                  cbind(corr_31,corr_32,corr_33))
#Obtain the eigen functiosn of the entire correlation surface 
combined.eigenbasis <- eigen(Corr_mat)$vectors[,1:3]
combined.eigenvalue <- sum(eigen(Corr_mat)$values[1:3]) * c(0.65,0.3,0.05)

par(mfrow=c(1,1))
image2D(Corr_mat)


################
#Simulate observations
################

#Preallocate space for storing simulation result 
nrep <- 200
sim_errorvar    <- matrix(NA, nrow = nrep, ncol = 3)
sim_v           <- list()
sim_v[[1]] <- sim_v[[2]] <- sim_v[[3]] <- matrix(NA, nrow = nrep, ncol = ntimes)
sim.eigenvalue  <- matrix(NA, nrow = nrep,ncol=3)
sim.eigenbasis <- list()
sim.eigenbasis[[1]] <- sim.eigenbasis[[2]]<-sim.eigenbasis[[3]] <- matrix(NA, nrow = nrep, ncol = ntimes*3)

for (sim_index in 1:nrep){
  nsim = 500
  fpc.scores <- matrix(NA, nrow = nsim, ncol = 3)
  for (pc in 1:3){
    fpc.scores[,pc] <- rnorm(nsim, mean = 0 , sd = sqrt(combined.eigenvalue[pc]))
  }
  
  combined.observation <- matrix(NA, nrow = ntimes * 3, ncol = nsim)
  for (ii in 1:nsim){
    combined.observation[,ii] <- as.vector(fpc.scores[ii,] %*% t(combined.eigenbasis))
  }
  
  #checked
  #Define v_{k}(t)
  v_mat <- matrix(NA, nrow = ntimes, ncol = 3)
  
  #X(t)
  v_mat[,1] <- unlist(lapply(timepoints, function(x){
    0.1+(1-x)^2
    
  }))
  v_mat[,2] <- unlist(lapply(timepoints, function(x){
    0.1 + .5 * (x-1)^4
  }))
  #W(t)
  v_mat[,3] <- unlist(lapply(timepoints, function(x){
    0.2+ 0.2 *cos(2 *pi *x)
  }))

  v_vec <- c(v_mat[,1],v_mat[,2],v_mat[,3])
  #Get unscaled observations
  unscale.observation <- apply(combined.observation,2,function(x){ x * sqrt(v_vec)})
  
  #Seperate combined observations into selves
  observ <- list(matrix(NA, nrow = ntimes, ncol = nsim),
                 matrix(NA, nrow = ntimes, ncol = nsim),
                 matrix(NA, nrow = ntimes, ncol = nsim))
  observ[[1]] <- unscale.observation[1:ntimes,]
  observ[[2]] <- unscale.observation[(ntimes+1):(ntimes*2),]
  observ[[3]] <- unscale.observation[(ntimes*2 + 1):(ntimes*3),]
  
  #Add measurement noises 
  error_var <- c(0.09,0.04,0.01)
  #Euclidean variable 
  observ[[1]] <- apply(observ[[1]],2,function(x){ x + rnorm(ntimes, mean = 0 , sd = sqrt(error_var[1]))})
  #Difference in radians 
  observ[[2]] <- apply(observ[[2]],2,function(x){ x + rnorm(ntimes, mean = 0 , sd = sqrt(error_var[2]))})
  #Euclidean variable 
  observ[[3]] <- apply(observ[[3]],2,function(x){ x + rnorm(ntimes, mean = 0 , sd = sqrt(error_var[3]))})
  
  ###########
  #The mean curve on S-1  is defined to be the constant 90 degree
  #Doing exponential mapping to project tangent vectors back to the manifold 
  mean_angle = unlist(lapply(timepoints,function(x){90 + 100 * (sin(2 *pi * x))^2}))
  #Preallocate space for storing exponential mapped data
  expmap_data = matrix(NA, nrow = nrow(observ[[3]]), ncol = ncol(observ[[3]]))
  
  #Exponential mapping
  for (ii in 1:ntimes){
    mean_ang_rad  <- mean_angle[ii] * pi / 180
    #Get points on the circle
    mean_coord  <-  c( cos(mean_angle[ii]), sin(mean_angle[ii]))
    #unit tangent vector towards anti-clockwise
    tangvec_pos <-  c(-sin(mean_angle[ii]), cos(mean_angle[ii]))
    #unit tangent vector towards clockwise 
    tangvec_neg <- -c(-sin(mean_angle[ii]), cos(mean_angle[ii]))
    
    #length of tangent vector v , sqrt(sum(v^2)) = distance(p,q)
    for (jj in 1:nsim){
      expmap_data[ii,jj] = mean_ang_rad + observ[[3]][ii,jj] 
    }
  }
  
  for (ii in 1:nrow(expmap_data)){
    for (jj in 1:nsim){
      if(expmap_data[ii,jj] < 0){
        expmap_data[ii,jj] <- 2*pi + expmap_data[ii,jj]
      }
      if(expmap_data[ii,jj] > 2*pi){
        expmap_data[ii,jj] <-  expmap_data[ii,jj] - 2 * pi
      }
    }
  }
  
  
  #simdata <- list(observ[[1]],expmap_data)
  
  #End of simulation of data 
  #------------------------------------------------------------------
  #1. Estimating frechet mean
  sim_mean_angle <- rep(NA,ntimes)
  for (tt in 1:ntimes){
    data_vec <- expmap_data[tt,]* 180 / pi
    sim_mean_angle[tt] = optimize(distance.fun, data_vec = data_vec,c(0,360))$minimum
  }

  #par(mfrow=c(1,1))
  #plot(timepoints, mean_angle,type='l',col='blue',lwd=2)
  #lines(timepoints, sim_mean_angle,col='purple',lwd=2)
  #legend('topleft',lty = c('solid','solid'),col=c('blue','purple'),legend=c('True','Est'))
  
  #2. Apply logarithm mapping of data on manifold (1-sphere)
  #sim_mean_angle = mean_angle
  logmap_data <- matrix(NA , nrow = nrow(expmap_data), ncol = ncol(expmap_data))
  rotate.angle <- rep(NA,101)
  for (tt in 1:ntimes){
    #vector of radians at current time tt
    data_vec <- expmap_data[tt,]
    nonna   <- which(!is.na(data_vec))
    #Change mean angle from degree to radian 
    mean_rand  <- sim_mean_angle[tt] * pi / 180
    #Get cartesian coordinates of mean radian
    mean_vec <- c(cos(mean_rand), sin(mean_rand))
    #Apply log mapping to each radian 
    for (index in 1:length(data_vec)){
      q_vec <- c(cos(data_vec[index]),sin(data_vec[index]))
      proj_p_q <- as.numeric(mean_vec %*% q_vec) * mean_vec
      u = q_vec - proj_p_q
      tang_vec <- u / sqrt(sum(u^2)) *  c(acos(mean_vec %*% q_vec))
      if (mean_rand > 0 && mean_rand < pi){
        if ( (data_vec[index] - mean_rand) >=0 && (data_vec[index] - mean_rand) <= pi){
          logmap_data[tt,index] <- sqrt(sum(tang_vec^2))
        }else{
          logmap_data[tt,index] <- -sqrt(sum(tang_vec^2))
        }
      }
      if( mean_rand ==0){
        if (data_vec[index] >=0 && data_vec[index] <= pi){
          logmap_data[tt,index] <- sqrt(sum(tang_vec^2))
        }else{
          logmap_data[tt,index] <- -sqrt(sum(tang_vec^2))
        }
      }
      if( mean_rand == (pi)){
        if (data_vec[index] >=0 && data_vec[index] <= pi){
          logmap_data[tt,index] <- -sqrt(sum(tang_vec^2))
        }else{
          logmap_data[tt,index] <- sqrt(sum(tang_vec^2))
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( (data_vec[index] - mean_rand) <=0 && (data_vec[index] - mean_rand) >= -pi){
          logmap_data[tt,index] <- -sqrt(sum(tang_vec^2))
        }else{
          logmap_data[tt,index] <- sqrt(sum(tang_vec^2))
        }
      }
      
    }
  }
  
  #perform multivariate FPCA with a normalization approach 
  combined.data <- rbind(observ[[1]],observ[[2]],logmap_data)
  
  #Doing FPCA for each dimension to find diagonal elements of the covariance functions
  #and the variance of meansurement errors
  var.dim1 <- MakeFPCAInputs(IDs = rep(1:nsim, each = ntimes),
                             tVec = rep(seq(1,101,length.out = 101), nsim),observ[[1]])
  dim1.fpca <- FPCA(var.dim1$Ly, var.dim1$Lt,
                    opt=list(userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
                             FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
  
  var.dim2 <- MakeFPCAInputs(IDs = rep(1:nsim, each = ntimes),
                             tVec = rep(seq(1,101,length.out = 101), nsim),observ[[2]])
  dim2.fpca <- FPCA(var.dim2$Ly, var.dim2$Lt, opt =list(
    userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
    FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
  
  var.dim3 <- MakeFPCAInputs(IDs = rep(1:nsim, each = ntimes),
                             tVec = rep(seq(1,101,length.out = 101), nsim),logmap_data)
  dim3.fpca <- FPCA(var.dim3$Ly, var.dim3$Lt, opt =list(
    userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
    FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
  
  sim_errorvar[sim_index,1] <- dim1.fpca$sigma2
  sim_errorvar[sim_index,2] <- dim2.fpca$sigma2
  sim_errorvar[sim_index,3] <- dim3.fpca$sigma2
  #Perform normalization 
  par(mfrow=c(1,3))
  sim_v[[1]][sim_index,] <- sample_v_1 <- diag(dim1.fpca$smoothedCov)
  plot(timepoints, v_mat[,1],type='l',col='blue');lines(timepoints, sample_v_1,col='purple')
  sim_v[[2]][sim_index,] <- sample_v_2 <- diag(dim2.fpca$smoothedCov)
  plot(timepoints, v_mat[,2],type='l',col='blue');lines(timepoints, sample_v_2,col='purple')
  sim_v[[3]][sim_index,] <- sample_v_3 <- diag(dim3.fpca$smoothedCov)
  plot(timepoints, v_mat[,3],type='l',col='blue');lines(timepoints, sample_v_3,col='purple')
  
  normalized.data <- list()
  normalized.data[[1]] <- apply(observ[[1]],2,function(x){x / sqrt(diag(dim1.fpca$smoothedCov))})
  normalized.data[[2]] <- apply(observ[[2]],2,function(x){x / sqrt(diag(dim2.fpca$smoothedCov))})
  normalized.data[[3]] <- apply(logmap_data,2,function(x){x / sqrt(diag(dim3.fpca$smoothedCov))})
  normalized.data <- Reduce('rbind',normalized.data);dim(normalized.data)
  
  var.combine <- MakeFPCAInputs(IDs = rep(1:nsim, each = ntimes*3),
                                tVec = rep(seq(0,302,length.out = 303), nsim), normalized.data)
  combine.fpca <- FPCA(var.combine$Ly, var.combine$Lt,
                       opt=list(userMu = list(t = seq(0,302,length.out = 303),mu = rep(0,303)),
                                FVEthreshold = 0.99,methodMuCovEst ='smooth'))
  sim.eigenvalue[sim_index,] <- combine.fpca$lambda
  for (fpc.index in 1:3){
    sim.eigenbasis[[fpc.index]][sim_index,] <- combine.fpca$phi[,fpc.index]
  }
}

#save.image('sim_temp.Rdata')
#load('sim_temp.Rdata')

library(ggplot2)
library(cowplot)

#check component-wise covariance and error variance estimation 
v1.df  <- data.frame(time = timepoints, v = v_mat[,1]) 
v2.df  <- data.frame(time = timepoints, v = v_mat[,2])
v3.df  <- data.frame(time = timepoints, v = v_mat[,3])


plot.1 <- ggplot(data = v1.df, aes(x = timepoints, y = v),lwd=2) + geom_line(lwd=1.2) +
          geom_line(aes(x = timepoints, y = apply(sim_v[[1]],2,mean)),linetype='dashed', size = 1.05) + 
          geom_ribbon(aes(ymin = apply(sim_v[[1]],2,function(x){quantile(x,prob =0.025)}), ymax = apply(sim_v[[1]],2,function(x){quantile(x,prob =0.975)})),
                      alpha = 0.5,fill='steelblue')+
          theme_classic() + 
          ylab( expression(paste( G[11]))) + xlab('') + 
          theme(axis.text=element_text(size=12,face='bold'),
          axis.title.y = element_text(angle = 0,size=20,face="bold",vjust =0.5))


plot.2 <- ggplot(data = v2.df, aes(x = timepoints, y = v),lwd=2) + geom_line(lwd=1.2) +
          geom_line(aes(x = timepoints, y = apply(sim_v[[2]],2,mean)),linetype='dashed', size = 1.05) + 
          geom_ribbon(aes(ymin = apply(sim_v[[2]],2,function(x){quantile(x,prob =0.025)}), ymax = apply(sim_v[[2]],2,function(x){quantile(x,prob =0.975)})),
                      alpha = 0.5,fill='steelblue')+
          theme_classic() + 
          ylab( expression(paste( G[22]))) + xlab('') + 
          theme(axis.text=element_text(size=12,face='bold'),
                axis.title.y = element_text(angle = 0,size=20,face="bold",vjust =0.5))


plot.3 <- ggplot(data = v3.df, aes(x = timepoints, y = v),lwd=2) + geom_line(lwd=1.2) +
  geom_line(aes(x = timepoints, y = apply(sim_v[[3]],2,mean)),linetype='dashed', size = 1.05) + 
  geom_ribbon(aes(ymin = apply(sim_v[[3]],2,function(x){quantile(x,prob =0.025)}), ymax = apply(sim_v[[3]],2,function(x){quantile(x,prob =0.975)})),
              alpha = 0.5,fill='steelblue')+
  theme_classic() + 
  ylab( expression(paste( G[33]))) + xlab('')+ 
  theme(axis.text=element_text(size=12,face='bold'),
  axis.title.y = element_text(angle = 0,size=20,face="bold",vjust =0.5))



plot_grid(plot.1, plot.2, plot.3, nrow= 1)


err1.df <- data.frame(y = sim_errorvar[,1])
err2.df <- data.frame(y = sim_errorvar[,2])
err3.df <- data.frame(y = sim_errorvar[,3])
plot.1 <- ggplot(data = err1.df,aes( y= y)) + geom_boxplot(lwd=1.5) + 
  geom_hline(yintercept = error_var[1], linetype='dashed',lwd=1.2) + theme_classic() + xlab('') + 
  ylab(expression(paste(sigma[1]^2))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size=12,face='bold'),
        axis.title.y = element_text(size = 20, face = 'bold',angle = 0 , vjust = 0.5))
plot.2 <- ggplot(data = err2.df,aes( y= y)) + geom_boxplot(lwd=1.5) + 
  geom_hline(yintercept = error_var[2], linetype='dashed',lwd=1.2) + theme_classic() + xlab('') + 
  ylab(expression(paste(sigma[2]^2))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size=12,face='bold'),
        axis.title.y = element_text(size = 20, face = 'bold',angle = 0 , vjust = 0.5))
plot.3 <- ggplot(data = err3.df,aes( y= y)) + geom_boxplot(lwd=1.5) + 
  geom_hline(yintercept = error_var[3], linetype='dashed',lwd=1.2) + theme_classic() + xlab('') + 
  ylab(expression(paste(sigma[3]^2))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size=12,face='bold'),
        axis.title.y = element_text(size = 20, face = 'bold',angle = 0 , vjust = 0.5))
plot_grid(plot.1, plot.2, plot.3, nrow = 1)





#check eigenvalues 
dim.fpc <- list()
dim.fpc[[1]] <- dim.fpc[[2]] <- dim.fpc[[3]] <- matrix(NA, nrow = ntimes, ncol =3)
#Plot of eigen functions for simulating data 
dim.fpc[[1]][,1] <- combined.eigenbasis[1:101,1]
dim.fpc[[1]][,2] <- combined.eigenbasis[1:101,2]
dim.fpc[[1]][,3] <- combined.eigenbasis[1:101,3]
dim.fpc[[2]][,1] <- combined.eigenbasis[102:202,1]
dim.fpc[[2]][,2] <- combined.eigenbasis[102:202,2]
dim.fpc[[2]][,3] <- combined.eigenbasis[102:202,3]
dim.fpc[[3]][,1] <- combined.eigenbasis[203:303,1]
dim.fpc[[3]][,2] <- combined.eigenbasis[203:303,2]
dim.fpc[[3]][,3] <- combined.eigenbasis[203:303,3]

sim.fpc <- list()
sim.fpc[[1]] <- sim.fpc[[2]] <- sim.fpc[[3]] <- list()
sim.fpc[[1]][[1]] <- sim.eigenbasis[[1]][,1:101]
for (jj in 1:nrep){
  if (sim.fpc[[1]][[1]][jj,1] <0){
    sim.fpc[[1]][[1]][jj,] <- -sim.fpc[[1]][[1]][jj,]
  }
}

sim.fpc[[1]][[2]] <- sim.eigenbasis[[2]][,1:101]
sim.fpc[[1]][[3]] <- sim.eigenbasis[[3]][,1:101]
for (jj in 1:nrep){
  if (sim.fpc[[1]][[3]][jj,1] <0){
    sim.fpc[[1]][[3]][jj,] <- -sim.fpc[[1]][[3]][jj,]
  }
}

sim.fpc[[2]][[1]] <- sim.eigenbasis[[1]][,102:202]
for (jj in 1:nrep){
  if (sim.fpc[[2]][[1]][jj,1] <0){
    sim.fpc[[2]][[1]][jj,] <- -sim.fpc[[2]][[1]][jj,]
  }
}
sim.fpc[[2]][[2]] <- sim.eigenbasis[[2]][,102:202]


sim.fpc[[2]][[3]] <- sim.eigenbasis[[3]][,102:202]
for (jj in 1:nrep){
  if (sim.fpc[[2]][[3]][jj,1] <0){
    sim.fpc[[2]][[3]][jj,] <- -sim.fpc[[2]][[3]][jj,]
  }
}
sim.fpc[[3]][[1]] <- sim.eigenbasis[[1]][,203:303]
for (jj in 1:nrep){
  if (sim.fpc[[3]][[1]][jj,1] >0){
    sim.fpc[[3]][[1]][jj,] <- -sim.fpc[[3]][[1]][jj,]
  }
}
sim.fpc[[3]][[2]] <- sim.eigenbasis[[2]][,203:303]
sim.fpc[[3]][[3]] <- sim.eigenbasis[[3]][,203:303]
for (jj in 1:nrep){
  if (sim.fpc[[3]][[3]][jj,1] >0){
    sim.fpc[[3]][[3]][jj,] <- -sim.fpc[[3]][[3]][jj,]
  }
}
plot(0,0,xlim=c(0,1),ylim=c(-0.2,0.2),col='white')
for (jj in 1:nrep){
  lines(timepoints, sim.fpc[[3]][[3]][jj,])
}


df.1 <- data.frame(time = timepoints, fpc1 = dim.fpc[[1]][,1])
df.2 <- data.frame(time = timepoints, fpc2 = dim.fpc[[1]][,2])
df.3 <- data.frame(time = timepoints, fpc3 = dim.fpc[[1]][,3])
df.4 <- data.frame(time = timepoints, fpc1 = dim.fpc[[2]][,1])
df.5 <- data.frame(time = timepoints, fpc2 = dim.fpc[[2]][,2])
df.6 <- data.frame(time = timepoints, fpc3 = dim.fpc[[2]][,3])
df.7 <- data.frame(time = timepoints, fpc1 = dim.fpc[[3]][,1])
df.8 <- data.frame(time = timepoints, fpc2 = dim.fpc[[3]][,2])
df.9 <- data.frame(time = timepoints, fpc3 = -dim.fpc[[3]][,3])


#geom_ribbon(aes(ymin =  apply(sim.fpc[[1]][[1]],2,function(x){quantile(x,prob =0.025)}),
#                ymax =  apply(sim.fpc[[1]][[1]],2,function(x){quantile(x,prob =0.975)})), alpha = 0.5) +
  
plot.1 <- ggplot(data = df.1, aes(x = timepoints, y = fpc1)) + 
          geom_line(aes(x = timepoints, y = fpc1),col='blue',size = 1.5,alpha = 0.5) + 
          geom_line(aes(x = timepoints, y = apply(sim.fpc[[1]][[1]],2,mean)),size = 1.5, linetype='dashed') +
          theme_classic() + xlab('') + ylab(expression(paste(phi[11]))) + 
          theme(axis.text.x = element_blank(),
                axis.text.y = element_text(face='bold',size = 11),
                axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.2 <-  ggplot(data = df.2, aes(x = timepoints, y = fpc2)) +
geom_line(aes(x = timepoints, y = fpc2),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[1]][[2]],2,mean)),size = 1.5, linetype='dashed') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[12]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.3 <-  ggplot(data = df.3, aes(x = timepoints, y = fpc3)) +
  geom_line(aes(x = timepoints, y = fpc3),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[1]][[3]],2,mean)),size = 1.5, linetype='dashed') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[13]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.4 <- ggplot(data = df.4, aes(x = timepoints, y = fpc1))+
  geom_line(aes(x = timepoints, y = fpc1),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[2]][[1]],2,mean)),size = 1.5, linetype='dotted') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[21]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.5 <- ggplot(data = df.5, aes(x = timepoints, y = fpc2))+
  geom_line(aes(x = timepoints, y = fpc2),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[2]][[2]],2,mean)),size = 1.5, linetype='dotted') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[22]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.6 <- ggplot(data = df.6, aes(x = timepoints, y = fpc3))+
  geom_line(aes(x = timepoints, y = fpc3),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[2]][[3]],2,mean)),size = 1.5, linetype='dotted') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[23]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.7 <- ggplot(data = df.7, aes(x = timepoints, y = fpc1))+
  geom_line(aes(x = timepoints, y = fpc1),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[3]][[1]],2,mean)),size = 1.5, linetype='dotdash') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[31]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.8 <- ggplot(data = df.8, aes(x = timepoints, y = fpc2))+
  geom_line(aes(x = timepoints, y = fpc2),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[3]][[2]],2,mean)),size = 1.5, linetype='dotdash') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[32]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5)) 
plot.9 <- ggplot(data = df.9, aes(x = timepoints, y = fpc3))+
  geom_line(aes(x = timepoints, y = fpc3),col='blue',size = 1.5,alpha = 0.5) + 
  geom_line(aes(x = timepoints, y = apply(sim.fpc[[3]][[3]],2,mean)),size = 1.5, linetype='dotdash') +
  theme_classic() + xlab('') + ylab(expression(paste(phi[33]))) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size = 11),
        axis.title.y  = element_text(face='bold',size = 14,angle = 0, vjust = 0.5))  

plot_grid(plot.1,plot.2,plot.3,
          plot.4,plot.5,plot.6,
          plot.7,plot.8,plot.9, nrow = 3)


