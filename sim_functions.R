distance.fun <- function(mean.angle,data_vec){
  x_vec <- cos(data_vec * pi / 180)
  y_vec <- sin(data_vec * pi / 180)

  x_mean <- cos(mean.angle * pi / 180)
  y_mean <- sin(mean.angle * pi / 180)
  meanvec <- c(x_mean, y_mean)

  datamat <- cbind(x_vec, y_vec)

  dotprod <- rep(NA, length(data_vec))

  loc <- which(!is.na(data_vec))
  
  for (jj in 1:length(loc)){
    index = loc[jj]
    dotprod[index] <- datamat[index,] %*% meanvec
  }
  sqr.dist <- mean((acos(dotprod))^2,na.rm=TRUE)
  return(sqr.dist)
}
getdata <- function(x){

  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:252){
    filename <- paste('~/sfuvault/Wind Data/File ',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori/sfuvault/Wind Data/File ',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
  datamat   <- matrix(NA, nrow = numofrep, ncol = 24)
  rownames  <- rep(NA, numofrep)
  year      <- seq(1998,2018,1)
  day       <- seq(1,31,1)
  index     <- 1
  indicator <- matrix(NA, nrow = numofrep, ncol = 2)
  for (jj in 1:21){
    for (kk in 1:31){
      tempmat           <- dec.data[dec.data$Year == year[jj],]
      tempmat           <- tempmat[tempmat$Day  == day[kk],]
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }

  van.data <- list(spd.data = spd.data, ang.data = datamat)
  #----------------------------------------------------
  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:21){
    filename <- paste('~/sfuvault/toronto wind data/Toronto File',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori/sfuvault/toronto wind data/Toronto File',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
  datamat   <- matrix(NA, nrow = numofrep, ncol = 24)
  rownames  <- rep(NA, numofrep)
  year      <- seq(1998,2018,1)
  day       <- seq(1,31,1)
  index     <- 1
  indicator <- matrix(NA, nrow = numofrep, ncol = 2)
  for (jj in 1:21){
    for (kk in 1:31){
      tempmat           <- dec.data[dec.data$Year == year[jj],]
      tempmat           <- tempmat[tempmat$Day  == day[kk],]
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }
  ton.data <- list(spd.data = spd.data, ang.data = datamat)

  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:21){
    filename <- paste('~/sfuvault/halifax wind data/Halifax File',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori//sfuvault/halifax wind data/Halifax File',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
  datamat   <- matrix(NA, nrow = numofrep, ncol = 24)
  rownames  <- rep(NA, numofrep)
  year      <- seq(1998,2018,1)
  day       <- seq(1,31,1)
  index     <- 1
  indicator <- matrix(NA, nrow = numofrep, ncol = 2)
  for (jj in 1:21){
    for (kk in 1:31){
      tempmat           <- dec.data[dec.data$Year == year[jj],]
      tempmat           <- tempmat[tempmat$Day  == day[kk],]
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }
  hal.data <- list(spd.data = spd.data, ang.data = datamat)
  #-------------------------------------------------------------

  return(list(van.data = van.data,
              ton.data = ton.data,
              hal.data = hal.data))
}

llr.mean <- function(time,data,bw){
  s0_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s0_vec[ii] <- get_s_0(ii,h_bw =bw)
  }
  
  s1_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s1_vec[ii] <- get_s_1(ii,h_bw =bw)
  }
  
  s2_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s2_vec[ii] <- get_s_2(ii,h_bw =bw)
  }
  
  r0_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    r0_vec[ii] <- get_r_0(ii,obs =data,h_bw =bw)
  }
  
  r1_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    r1_vec[ii] <- get_r_1(ii,obs = data,h_bw =bw)
  }
  
  return((r0_vec * s2_vec - r1_vec * s1_vec)/(s0_vec*s2_vec - (s1_vec)^2))
}
get_s_0 <- function(t,h_bw){
  s_0_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_0_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)
    }
  }
  return(sum(s_0_mat,na.rm=TRUE)/nrow(s_0_mat)/ncol(s_0_mat))
}
get_s_1 <- function(t,h_bw){
  s_1_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_1_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*((timemat[ii,jj] - t)/h_bw)
    }
  }
  return(sum(s_1_mat,na.rm=TRUE)/nrow(s_1_mat)/ncol(s_1_mat))
}

get_s_2 <- function(t,h_bw){
  s_2_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_2_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*(((timemat[ii,jj] - t)/h_bw)^2)
    }
  }
  return(sum(s_2_mat,na.rm=TRUE)/nrow(s_2_mat)/ncol(s_2_mat))
}

get_r_0 <- function(t,obs,h_bw){
  r_0_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      r_0_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*obs[ii,jj]
    }
  }
  return(sum(r_0_mat,na.rm=TRUE)/nrow(r_0_mat)/ncol(r_0_mat))
}

get_r_1 <- function(t,obs,h_bw){
  r_1_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      r_1_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*((timemat[ii,jj] - t)/h_bw)*obs[ii,jj]
    }
  }
  return(sum(r_1_mat,na.rm=TRUE)/nrow(r_1_mat)/ncol(r_1_mat))
}

epane  <- function(x){
  if (abs(x) <= 1){
    return(3/4*(1-x^2))
  }else{
    return(0)
  }
  
}



# for (tt in 1:ntimes){
#   data_vec  = allcity.data[[ii]]$ang.data[,tt] * 10
#   mean_rand = ang.meanmat[ii,tt] * pi / 180
#   #Get cartesian coordinates of mean radian
#   mean_vec <- c(cos(mean_rand), sin(mean_rand))
#   for (index in 1:length(data_vec)){
#     data_rand <- data_vec[index] * pi / 180
#     #Vector from original to data point on circle
#     q_vec <- c(cos(data_rand),sin(data_rand))
#     #Projection of q onto p
#     proj_p_q <- as.numeric(mean_vec %*% q_vec) * mean_vec
#     #Tangent vector from proj_p_q to q
#     u = q_vec - proj_p_q
#     #normalized u and scale it with great circle distance (geodesic distance): acos(p,q)
#     tang_vec <- u / sqrt(sum(u^2)) * c( acos(mean_vec %*% q_vec))
#
#     if(is.na(data_vec[index])){
#       if (mean_rand > 0 && mean_rand < pi){
#         if ( angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if( mean_rand ==0){
#         if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if( mean_rand == (pi)){
#         if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
#           angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if (mean_rand > pi && mean_rand < (2 * pi)){
#         if ( angle.temp[index,tt] <=0 && angle.temp[index,tt] >= -pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#
#
#     }else{
#       if (mean_rand > 0 && mean_rand < pi){
#         if ( (data_rand - mean_rand) >=0 && (data_rand - mean_rand) <= pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if( mean_rand ==0){
#         if (data_rand >=0 && data_rand <= pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if( mean_rand == (pi)){
#         if (data_rand >=0 && data_rand <= pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#       if (mean_rand > pi && mean_rand < (2 * pi)){
#         if ( (data_rand - mean_rand) <=0 && (data_rand - mean_rand) >= -pi){
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }else{
#           angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
#         }
#       }
#
#     }
#
#   }
#
# }


S_pq <- function(p,q,s,t,h_R,timemat){
  tempmat <- matrix(NA, nrow = nrow(timemat), ncol = ncol(timemat))
  sum_vec <- rep(NA, nrow = nrow(timemat))
  for (ii in 1:nrow(timemat)){
    index = 1
    tempvec <- rep(NA, length(timemat[ii,])*( length(timemat[ii,])-1))
    for (kk in 1:ncol(timemat)){
      for (jj in 1:ncol(timemat)){
        if( jj != kk){
          tempvec[index] <- (((timemat[ii,jj] - s)/h_R)^p)*(((timemat[ii,kk] - t)/h_R)^q)*
            (1/h_R)*epane((timemat[ii,jj] - s)/h_R)*
            (1/h_R)*epane((timemat[ii,kk] - t)/h_R)
          index = index + 1
        }else{
          
        }
      }
    }
    sum_vec[ii] <- sum(tempvec)
  }
  return(sum(sum_vec) / length(sum_vec) / nrow(timemat))
}


R_pq <- function(p,q,s,t,h_R,timemat,obs){
  tempmat <- matrix(NA, nrow = nrow(timemat), ncol = ncol(timemat))
  sum_vec <- rep(NA, nrow = nrow(timemat))
  for (ii in 1:nrow(timemat)){
    index = 1
    tempvec <- rep(NA, length(timemat[ii,])*( length(timemat[ii,])-1))
    for (kk in 1:ncol(timemat)){
      for (jj in 1:ncol(timemat)){
        if( jj != kk){
          tempvec[index] <- obs[ii,jj] * obs[ii,kk]*(((timemat[ii,jj] - s)/h_R)^p)*(((timemat[ii,kk] - t)/h_R)^q)*
            (1/h_R)*epane((timemat[ii,jj] - s)/h_R)*
            (1/h_R)*epane((timemat[ii,kk] - t)/h_R)
          index = index + 1
        }else{
          
        }
      }
    }
    sum_vec[ii] <- sum(tempvec,na.rm= TRUE)
  }
  return(sum(sum_vec) / length(sum_vec) / nrow(timemat))
}


getcov <- function(s,t,mean_vec,h_R,obs, timemat){
  S_00  <- S_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_01  <- S_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_10  <- S_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_11  <- S_pq(p = 1 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_02  <- S_pq(p = 0 , q = 2 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_20  <- S_pq(p = 2 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  R_00  <- R_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)
  R_01  <- R_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)
  R_10  <- R_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)
  
  A_1   <- S_20 * S_02 - (S_11)^2
  A_2   <- S_10 * S_02 - S_01 * S_11
  A_3   <- S_01 * S_20 - S_10 * S_11
  B     <- A_1 * S_00 - A_2 * S_10 - A_3 * S_01
  
  cov.val <- (1/B) * (A_1 * R_00 - A_2 * R_10 - A_3 * R_01) - mean_vec[s+1] * mean_vec[t+1]
  return(cov.val)
}

getcrosscov <- function(s, t, h_R,obs_1,obs_2, timemat){
  S_00  <- S_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_01  <- S_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_10  <- S_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_11  <- S_pq(p = 1 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_02  <- S_pq(p = 0 , q = 2 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_20  <- S_pq(p = 2 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  Rc_00  <- Rc_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs_1 = obs_1, obs_2 = obs_2)
  Rc_01  <- Rc_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat, obs_1 = obs_1, obs_2 = obs_2)
  Rc_10  <- Rc_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat,  obs_1 = obs_1, obs_2 = obs_2)
  
  A_1   <- S_20 * S_02 - (S_11)^2
  A_2   <- S_10 * S_02 - S_01 * S_11
  A_3   <- S_01 * S_20 - S_10 * S_11
  B     <- A_1 * S_00 - A_2 * S_10 - A_3 * S_01
  
  cov.val <- (1/B) * (A_1 * Rc_00 - A_2 * Rc_10 - A_3 * Rc_01)
  return(cov.val)
}


