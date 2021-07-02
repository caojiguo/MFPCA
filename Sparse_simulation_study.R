rm(list=ls())
library(fdapace)
library(ggplot2)
library(plot3D)
library(scales)
source('sim_functions.R')
set.seed(11235)
#Number of observation points 
ntimes     <- 101
#Define observation points in the interval [0,1]
timepoints <- seq(0,1,length.out = ntimes)

#Pre-allocate space for storing functional principal components 
#For both dimensions
dim3_fpc <- dim2_fpc <- dim1_fpc <- matrix(NA, nrow = ntimes, ncol = 3)

#Thrid PC are constant for both dimensions 
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

par(mfrow=c(3,1))
image2D(Cov_11);image2D(corr_11)
image2D(Cov_22);image2D(corr_22)
image2D(Cov_33):image2D(corr_33)

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

nrep = 100

#Preallocate space for storing simulation result 
sparsity_vec <- seq(1,10) / 10

sim_errorvar    <- lapply(seq(1,10), function(x){matrix(NA, nrow = nrep, ncol = 3)})

sim_v           <- list()
sim_v[[1]] <- sim_v[[2]] <- sim_v[[3]] <- lapply(seq(1,10), function(x){matrix(NA, nrow = nrep, ncol = ntimes)})

sim.eigenvalue  <- lapply(seq(1,10), function(x){matrix(NA, nrow = nrep,ncol=3)})
sim.eigenbasis <- list()

sim.eigenbasis[[1]] <- sim.eigenbasis[[2]]<-sim.eigenbasis[[3]] <- lapply(seq(1,10),function(x){matrix(NA, nrow = 200, ncol = ntimes*3)})

for (sparse_index in 1:10){
  for (sim_index in 1:nrep){
    print(paste('Current simulation',sim_index,'on Sparsity',sparsity_vec[sparse_index]))
    #number of simulationf
    nsim = 500
    fpc.scores <- matrix(NA, nrow = nsim, ncol = 3)
    for (pc in 1:3){
      fpc.scores[,pc] <- rnorm(nsim, mean = 0 , sd = sqrt(combined.eigenvalue[pc]))
    }
    apply(fpc.scores,2,var);combined.eigenvalue
    
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
    
    #-------------------------
    #Sparsify each data 
    sparsity = sparsity_vec[sparse_index]
    dim1.sparsedf   <- Sparsify(t(observ[[1]]), pts = timepoints, sparsity = floor(sparsity * ntimes))
    dim2.sparsedf <- Sparsify(t(observ[[2]]), pts = timepoints, sparsity = floor(sparsity * ntimes))
    logmap.sparsedf   <- Sparsify(t(logmap_data), pts = timepoints, sparsity = floor(sparsity * ntimes))
    par(mfrow=c(1,3))
    #Doing FPCA for each dimension to find diagonal elements of the covariance functions
    #and the variance of meansurement errors
    dim1.fpca <- FPCA(dim1.sparsedf$Ly, dim1.sparsedf$Lt,
                      opt=list(userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
                               FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid = 101))
    plot(dim1.fpca$workGrid,diag(dim1.fpca$smoothedCov),type='l');lines(timepoints,v_mat[,1])
    dim2.fpca <- FPCA(dim2.sparsedf$Ly, dim2.sparsedf$Lt, opt =list(
      userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
      FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid = 101))
    plot(dim2.fpca$workGrid,diag(dim2.fpca$smoothedCov),type='l');lines(timepoints,v_mat[,2])
    dim3.fpca <- FPCA(logmap.sparsedf$Ly, logmap.sparsedf$Lt, opt =list(
      userMu = list(t = seq(0,101,length.out = 101),mu = rep(0,101)),
      FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid = 101))
    plot(dim3.fpca$workGrid,diag(dim3.fpca$smoothedCov),type='l');lines(timepoints,v_mat[,3])
    par(mfrow=c(1,1))
    
    sim_errorvar[[sparse_index]][sim_index,1] <- dim1.fpca$sigma2
    sim_errorvar[[sparse_index]][sim_index,2] <- dim2.fpca$sigma2
    sim_errorvar[[sparse_index]][sim_index,3] <- dim3.fpca$sigma2
    
    
    sim_v[[1]][[sparse_index]][sim_index,] <- sample_v_1 <- diag(dim1.fpca$smoothedCov)
    plot(timepoints, v_mat[,1],type='l',col='blue');lines(timepoints, sample_v_1,col='purple')
    sim_v[[2]][[sparse_index]][sim_index,] <- sample_v_2 <- diag(dim2.fpca$smoothedCov)
    plot(timepoints, v_mat[,2],type='l',col='blue');lines(timepoints, sample_v_2,col='purple')
    sim_v[[3]][[sparse_index]][sim_index,] <- sample_v_3 <- diag(dim3.fpca$smoothedCov)
    plot(timepoints, v_mat[,3],type='l',col='blue');lines(timepoints, sample_v_3,col='purple')
    
    #Normalization
    normalized.df    <- list()
    normalized.df$Lt <- list()
    normalized.df$Ly <- list()
    for (jj in 1:nsim){

      #Concatenating timepoints
      normalized.df$Lt[[jj]] <- c(dim1.sparsedf$Lt[[jj]], dim2.sparsedf$Lt[[jj]] + 1.01, logmap.sparsedf$Lt[[jj]] + 2.02) * 100

      observ.index <- dim1.fpca$workGrid %in%  dim1.sparsedf$Lt[[jj]]
      normvec <- sqrt(diag(dim1.fpca$smoothedCov))[observ.index]
      dim1.normed <- dim1.sparsedf$Ly[[jj]] / normvec

      observ.index <- dim2.fpca$workGrid  %in%  dim2.sparsedf$Lt[[jj]]
      normvec <- sqrt(diag(dim2.fpca$smoothedCov))[observ.index]
      dim2.normed  <- dim2.sparsedf$Ly[[jj]] / normvec

      observ.index <- dim3.fpca$workGrid %in%  logmap.sparsedf$Lt[[jj]]
      normvec <- sqrt(diag(dim3.fpca$smoothedCov))[observ.index]
      log.normed <- logmap.sparsedf$Ly[[jj]] / normvec

      #Concatenating observations
      normalized.df$Ly[[jj]] <- c(dim1.normed, dim2.normed, log.normed)

    }

    combine.fpca <- FPCA(normalized.df$Ly, normalized.df$Lt,
                         opt=list(userMu = list(t = seq(0,302,length.out = 303),mu = rep(0,303)),
                                  FVEthreshold = 0.90,methodMuCovEst ='smooth',nRegGrid = 303, verbose=T))

    for (fpc.index in 1:2){
      sim.eigenbasis[[fpc.index]][[sparse_index]][sim_index,] <- combine.fpca$phi[,fpc.index]
    }

    #save.image('sparse_simulation.Rdata')
  }
}
#load saved data for plots 
load('sparse_simulation.Rdata')

#examine result
true.eigenbasis <- list()
true.eigenbasis[[1]] <- true.eigenbasis[[2]] <- true.eigenbasis[[3]] <- matrix(NA, nrow = 101, ncol = 3)
for (jj in 1:3){
  true.eigenbasis[[jj]][,1] <- combined.eigenbasis[1:101,jj]
}
for (jj in 1:3){
  true.eigenbasis[[jj]][,2] <- combined.eigenbasis[102:202,jj]
}
for (jj in 1:3){
  true.eigenbasis[[jj]][,3] <- combined.eigenbasis[203:303,jj]
}

#check v
par(mfrow=c(1,3))
plot.index = 23
plot(timepoints, v_mat[,1],type='l',lwd=2)
for (jj in 1:9){
  lines(timepoints,sim_v[[1]][[jj]][plot.index,],col='blue')
}
lines(timepoints, sim_v[[1]][[10]][plot.index,], col='red',lwd=1.5)

plot(timepoints, v_mat[,2],type='l',lwd=2)
for (jj in 1:9){
  lines(timepoints,sim_v[[2]][[jj]][plot.index,],col='blue')
}
lines(timepoints, sim_v[[2]][[10]][plot.index,], col='red',lwd=1.5)

plot(timepoints, v_mat[,3],type='l',lwd=2)
for (jj in 1:9){
  lines(timepoints,sim_v[[3]][[jj]][plot.index,],col='blue')
}
lines(timepoints, sim_v[[3]][[10]][plot.index,], col='red',lwd=1.5)


nquadpts     <- length(timepoints)
quadwts      <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
quadwts      <- c(1,rep(c(4,2),(nquadpts-1)/2))
quadwts[nquadpts] <- 1
quadwts      <- quadwts/ (3 * ntimes)

#mean integrated squared errors
v.MISE <- list()
v.plot <- list()
for (dim in 1:3){
  v.MISE[[dim]] <- matrix(NA, nrow = 50, ncol = 10)
  for (jj in 1:10){
    for (ii in 1:50){
      diff <- sim_v[[dim]][[jj]][ii,] - v_mat[,1]
      v.MISE[[dim]][ii,jj] <-  sum((sqrt(quadwts) * diff)^2)
    }
  }
  
  v1.error  <- as.vector(v.MISE[[dim]])
  v1.group  <- factor(rep(c('10%','20%','30%','40%','50%',
                            '60%','70%','80%','90%','100%' ),each = 50),
                      levels = c('100%','90%','80%','70%','60%',
                                 '50%','40%','30%','20%','10%' ))
  v1.df <- data.frame(group = v1.group, error = v1.error)
  
  v.plot[[dim]] <- ggplot(data =v1.df,aes(x = group, y = error)) + geom_boxplot(width=0.5) + theme_classic() + 
    scale_y_continuous(labels = format_format(big.mark = " ", decimal.mark = ".", scientific = TRUE)) + 
    theme(axis.text = element_text(size = 14,face='bold')) + xlab('') + ylab('')
  
}
library(cowplot)
plot_grid(v.plot[[1]],v.plot[[2]],v.plot[[3]],nrow=3)



#summarized in table 
apply(v.MISE[[1]],2,mean)
sqrt(apply(v.MISE[[1]],2,var))
apply(v.MISE[[2]],2,mean)
sqrt(apply(v.MISE[[2]],2,var))
apply(v.MISE[[3]],2,mean)
sqrt(apply(v.MISE[[3]],2,var))



#ISE for 1st fpc
for (ii in 1:50){
  for (jj in 1:10){
    if (sim.eigenbasis[[1]][[jj]][ii,1] <0){
      sim.eigenbasis[[1]][[jj]][ii,] <- -sim.eigenbasis[[1]][[jj]][ii,]
    }
  }
}
fpc.1.ISE <- list()
fpc1.plot <- list()
for (dim in 1:3){
  fpc.1.ISE[[dim]] <- matrix(NA, nrow = 50, ncol = 10)
  for (jj in 1:10){
    for (ii in 1:50){
        if (dim == 1){
          diff <- sim.eigenbasis[[1]][[jj]][ii,1:101] - true.eigenbasis[[1]][,1]
          fpc.1.ISE[[dim]][ii,jj]  <- sum((sqrt(quadwts) * diff)^2)
        }
        if (dim == 2){
          diff <- sim.eigenbasis[[1]][[jj]][ii,102:202] - true.eigenbasis[[1]][,2]
          fpc.1.ISE[[dim]][ii,jj]  <- sum((sqrt(quadwts) * diff)^2)
        }
        if (dim == 3){
          diff <- sim.eigenbasis[[1]][[jj]][ii,203:303] - true.eigenbasis[[1]][,3]
          fpc.1.ISE[[dim]][ii,jj]  <- sum((sqrt(quadwts) * diff)^2)
        }
    }
  }

  
  v1.error  <- as.vector(fpc.1.ISE[[dim]])
  v1.group  <- factor(rep(c('10%','20%','30%','40%','50%',
                            '60%','70%','80%','90%','100%' ),each = 50),
                      levels = c('100%','90%','80%','70%','60%',
                                 '50%','40%','30%','20%','10%' ))
  v1.df <- data.frame(group = v1.group, error = v1.error)
  yupper = c(3e-4,9e-4,1e-3)
  fpc1.plot[[dim]] <- ggplot(data =v1.df,aes(x = group, y = error)) + geom_boxplot(width=0.5) + theme_classic() + 
    scale_y_continuous(labels = format_format(big.mark = " ", decimal.mark = ".", scientific = TRUE),limits = c(0,yupper[dim])) + 
    theme(axis.text = element_text(size = 14,face='bold')) + xlab('') + ylab('')
}
plot_grid(fpc1.plot[[1]],fpc1.plot[[2]],fpc1.plot[[3]],nrow =3)



