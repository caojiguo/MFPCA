#clear workspace 
rm(list=ls())
set.seed(1234)
#load dependencies 
library(ggplot2)
library(cowplot)
library(fdapace)
library(plot3D)
library(RColorBrewer)
#some functions for runing FPCA_MVDM
source('functions.R')

###########################################
#Data has been saved to a Rdata file, 'winddata.rds', from the original csv files which are 
#available from climate.weather.gov.ca
###########################################
#get wind data from Vancouver, Toronto, and Halifax
data <- readRDS('winddata.rds')

#calculating missing entries in wind angle observations
missing.df <- apply(data$van.data$ang.data,1,function(x){length(which(is.na(x)))}) / 24
#Converting into proportions rather than frequencies 

#plot of the distribution of portions of missing entries
ggplot(data = data.frame(numofmissing = missing.df),aes(x = numofmissing)) + 
  geom_histogram(breaks=seq(0,1,length.out = 11)[-1], 
                 col="black",  
                 fill="grey", 
                 alpha=.7) + theme_classic() + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 14, face = 'bold')) +
  scale_x_continuous(breaks=seq(0,1,length.out = 11)[-1],
                     labels = scales::percent(x = seq(0,1,length.out = 11)[-1],accuracy = 1)) #+ 
  #scale_y_continuous(labels = scales::percent(x= c(0,0.02,0.04,0.06,0.08),accuracy = 1))


#wind rose plot
#have been createed through another python file , see 'windrose_plot.py'
plot.windrose(spd = as.vector(data$van.data$spd.data),
              dir = as.vector(data$van.data$ang.data * 10),
              spdseq = seq(0,80,10),spdmin = 0, spdmax = 82)



##########################################################
#1a. Estimate mean funcitons of wind direction, Riemannian variable, and 
#    speed, Euclidean variable.
#1b. Perform log-mapping on direction variable.
#2. Estimate covariance functions of both wind direction and speed.
#3. Normalization
#4. FPCA_MVDM
#5. Repeat 1-4 for Vancouver, Toronto and Halifax
##########################################################

#Total curve display for cities
tempdata_ang <- data$van.data$ang.data * 10
tempdata_spd <- data$van.data$spd.data


plot_size = 2
plot_sample = sample(sort(missing.df, index.return=TRUE,decreasing=TRUE)$ix[1:100], size = plot_size)
#plot_sample = sort(missing.df, index.return=TRUE,decreasing=TRUE)$ix[1:plot_size]


ang.df <- data.frame(ID = rep(1:plot_size, each = 24),
                     Time = rep(1:24, plot_size),
                     Direction = c(tempdata_ang[plot_sample,]))
ang.df$ID = as.character(ang.df$ID)
plot_1 = ggplot(data = ang.df, aes(x = Time, y = Direction,color = ID))  + geom_line(size=1.2) + 
      theme_minimal() + theme(axis.text=element_text(size=12,face='bold'),
                              axis.title=element_text(size=14,face="bold") ,
                              axis.title.x = element_blank(),
                              axis.text.x  = element_blank(),
                              legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) + 
  scale_y_continuous(breaks = seq(0,360,60))
  

spd.df <- data.frame(ID = rep(1:plot_size, each = 24),
                     Time = rep(1:24, plot_size),
                     Speed = c(tempdata_spd[plot_sample,]))
spd.df$ID = as.character(ang.df$ID)
plot_2 = ggplot(data = spd.df, aes(x = Time, y = Speed,color = ID))  + geom_line(size=1.2) + 
  theme_minimal() + theme(axis.text=element_text(size=12,face='bold'),
                          axis.title=element_text(size=14,face="bold") ,
                          legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) 
plot_grid(plot_1,plot_2,nrow = 2, align = 'v')

###
#1a. Vancouver
###
#Estimate the mean functions of vancouer wind speed and direction 
#Mean function of wind direction (angle)
ntimes = 24
van_mean_angle <- rep(NA,ntimes)
for (tt in 1:ntimes){
  #Conert radians to degree
  data_vec      <- (data$van.data$ang.data * 10)[,tt]
  van_mean_angle[tt] <- optimize(distance.fun, data_vec = data_vec,c(0,360))$minimum
}

#mean function of wind speed of vancouver 
timemat <- matrix(rep(0:23,651),byrow=TRUE, nrow = 651)
van_mean_speed <- llr.mean(time = seq(0,23,1),data = data$van.data$spd.data, bw = 3)
plot(seq(0,23,1), van_mean_speed,type='l')

###
#1b. Vancouver
###
van.ang.logmap <- matrix(NA , nrow = nrow(data$van.data$ang.data), ncol = ncol(data$van.data$ang.data))
rotate.angle <- rep(NA,24)
for (ii in 1:24){
  print(paste('Hour',ii))
  data_vec          <- data$van.data$ang.data[,ii] * 10
  
  mean_coord     <- c(cos(van_mean_angle[ii] * pi / 180), sin(van_mean_angle[ii] * pi / 180))
  
  nonna   <- which(!is.na(data_vec))
  tempmat <- matrix(NA,nrow = length(data_vec),ncol =2 )
  for (index in 1:length(nonna)){
    x_vec <- c(cos(data_vec[nonna[index]] * pi / 180),sin(data_vec[nonna[index]] * pi / 180))
    vec_u <- x_vec - as.vector(mean_coord %*% x_vec) * mean_coord
    v_vec <- (vec_u / sqrt(as.vector(vec_u %*% vec_u))) * acos(as.vector(mean_coord %*% x_vec))
    tempmat[nonna[index],] <- v_vec
  }
  
  tempmean <- apply(tempmat,2,function(x) return(mean(x,na.rm=TRUE)))
  
  
  if (coef(lm(tempmat[,2] ~ tempmat[,1]))[2] > 0){
    angle <- 2*pi - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }else{
    angle <- - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }
  
  rotate.angle[ii] <- angle
  rotation.mat <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
  
  rotated.points <- matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
  for (index in 1:length(nonna)){
    prep_vec <- tempmat[nonna[index],]
    prep_vec <- prep_vec - tempmean
    rotated.points[nonna[index],] <-  rotation.mat   %*% prep_vec
  }
  tempvec <- rep(NA, length(data_vec))
  for (jj in 1:length(nonna)){
    tempvec[nonna[jj]] <- rotated.points[nonna[jj],1]
  }
  van.ang.logmap[,ii] <- tempvec
}
#check mean of logmapped data = 0 
apply(van.ang.logmap,2,function(x){mean(x, na.rm=TRUE)})
#subtract the mean function from speed 
van.spd.demean <- t(apply(data$van.data$spd.data, 1, function(x){ x - van_mean_speed}))

###
#2. Vancouver
###
#Covariance of wind direction 
van.ang.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(van.ang.logmap))
van.ang.FPCA <- FPCA(van.ang.L$Ly, van.ang.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid=24))
#G(s,s)

v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(van.ang.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_classic() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=14, face = 'bold'),
        axis.title=element_text(size=14,face="bold"))
# timegrid <- expand.grid(seq(0,23,length.out = 101),seq(0,23,length.out = 101))[c(2,1)]
# van.ang.cov <- as.data.frame(cbind(timegrid, as.vector(van.ang.FPCA$fittedCov)))
# colnames(van.ang.cov) <- c('x','y','cov')
# ggplot(van.ang.cov, aes(x = x, y = y, z = cov)) + geom_tile(aes(fill = cov))
#   stat_contour(geom = 'polygon', aes(fill = ..level..)) + 
#   geom_tile(aes(fill = cov)) +
#   stat_contour(bins = 15) 


#Covariance of wind speed
van.spd.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(van.spd.demean))
van.spd.FPCA <- FPCA(van.spd.L$Ly, van.spd.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid = 24))

v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(van.spd.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_classic() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=14, face = 'bold'),
        axis.title=element_text(size=14,face="bold"))

library(tidyr)
library(ggplot2)
covplotdata =  data.frame(cbind(expand.grid(seq(0,24,length.out=201),seq(0,24,length.out=201)),as.vector(van.ang.FPCA$fittedCov)))
colnames(covplotdata) = c('x','y','Cov')
b <- seq(-0.3,1,0.2)
pal = c(hcl(0,100, seq(20,100, length.out=25)), hcl(240,100, seq(100,20, length.out=25)))
ggplot(data = covplotdata, aes(x=x,y=y,z=Cov))+geom_raster(aes(fill=Cov)) + 
 ylab('') + xlab('') +  scale_fill_gradientn(colours =rev(brewer.pal(11,"RdYlBu")), breaks = b, labels = format(b))+
  theme(axis.line=element_blank())+theme_bw()+
  theme(axis.ticks.length=unit(.2, "cm"),axis.ticks = element_line(size=1.5))+
  theme(axis.text=element_text(size=20, face = 'bold'),
        axis.title=element_text(size=20,face="bold")) + 
  theme(axis.text.x = element_text(vjust = 1), axis.text.y = element_text(hjust = 0.8))+
  scale_x_continuous(breaks=seq(0,24,4), expand = c(0, 0)) + scale_y_continuous(breaks=seq(0,24,4), expand = c(0, 0)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + theme(legend.title=element_blank())+
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20,draw.ulim = FALSE, draw.llim = FALSE))+
  theme(legend.text = element_text(size=16,face="bold"))

b = seq(20,70,10)
covplotdata =  data.frame(cbind(expand.grid(seq(0,24,length.out=201),seq(0,24,length.out=201)),as.vector(van.spd.FPCA$fittedCov)))
colnames(covplotdata) = c('x','y','Cov')

ggplot(data = covplotdata, aes(x=x,y=y,z=Cov))+geom_raster(aes(fill=Cov)) + 
  ylab('') + xlab('') +  scale_fill_gradientn(colours =rev(brewer.pal(11,"RdYlBu")), breaks = b, labels = format(b))+
  theme(axis.line=element_blank())+theme_bw()+
  theme(axis.ticks.length=unit(.2, "cm"),axis.ticks = element_line(size=1.5))+
  theme(axis.text=element_text(size=20, face = 'bold'),
        axis.title=element_text(size=20,face="bold")) + 
  theme(axis.text.x = element_text(vjust = 1), axis.text.y = element_text(hjust = 0.8))+
  scale_x_continuous(breaks=seq(0,24,4), expand = c(0, 0)) + scale_y_continuous(breaks=seq(0,24,4), expand = c(0, 0)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + theme(legend.title=element_blank())+
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20,draw.ulim = FALSE, draw.llim = FALSE))+
  theme(legend.text = element_text(size=16,face="bold"))


#3. Vancouver
###
van.spd.nlized  <- t(apply(van.spd.demean,1,function(x){ x / sqrt(diag(van.spd.FPCA$fittedCov))}))
van.ang.nlized  <- t(apply(van.ang.logmap,1,function(x){ x / sqrt(diag(van.ang.FPCA$fittedCov))}))
nlized.data <- cbind(van.ang.nlized,van.spd.nlized)


###
#4. Vancouver
###
van.nlized.L   <- MakeFPCAInputs(IDs = rep(1:651, each=48), 
                             tVec=rep(seq(0,47,1),651), t(nlized.data))
van.nlized.FPCA <- FPCA(van.nlized.L$Ly, van.nlized.L$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth',nRegGrid = 48))
van.nlized.FPCA$cumFVE
#Taking leading 4 FPCS
van.spd.fpc <- van.ang.fpc <- matrix(NA, ncol = 4, nrow = 24)
van.ang.fpc.plots <- list()
van.spd.fpc.plots <- list()
for (jj in 1:4){
  van.ang.fpc[,jj] <- van.nlized.FPCA$phi[1:24,jj]
  van.spd.fpc[,jj] <- van.nlized.FPCA$phi[25:48,jj]
}





#################################
#Do mixed FPCA for Toronto and Halifax
#Estimate the mean functions of Toronto wind speed and direction 
#Mean function of wind direction (angle)
ntimes = 24
ton_mean_angle <- rep(NA,ntimes)
for (tt in 1:ntimes){
  #Conert radians to degree
  data_vec      <- (data$ton.data$ang.data * 10)[,tt]
  ton_mean_angle[tt] <-  optimize(distance.fun, data_vec = data_vec,c(0,360))$minimum
}

#mean function of wind direction of Toronto 
timemat <- matrix(rep(0:23,651),byrow=TRUE, nrow = 651)
ton_mean_speed <- llr.mean(time = seq(0,23,1),data = data$ton.data$spd.data, bw = 3)
plot(seq(0,23,1), ton_mean_speed,type='l')
###############
#Perform log-mappping to direction data 
ton.ang.logmap <- matrix(NA , nrow = nrow(data$ton.data$ang.data), ncol = ncol(data$ton.data$ang.data))
rotate.angle <- rep(NA,24)
for (ii in 1:24){
  print(paste('Hour',ii))
  data_vec          <- data$ton.data$ang.data[,ii] * 10
  
  mean_coord     <- c(cos(ton_mean_angle[ii] * pi / 180), sin(ton_mean_angle[ii] * pi / 180))
  
  nonna   <- which(!is.na(data_vec))
  tempmat <- matrix(NA,nrow = length(data_vec),ncol =2 )
  for (index in 1:length(nonna)){
    x_vec <- c(cos(data_vec[nonna[index]] * pi / 180),sin(data_vec[nonna[index]] * pi / 180))
    vec_u <- x_vec - as.vector(mean_coord %*% x_vec) * mean_coord
    v_vec <- (vec_u / sqrt(as.vector(vec_u %*% vec_u))) * acos(as.vector(mean_coord %*% x_vec))
    tempmat[nonna[index],] <- v_vec
  }
  
  tempmean <- apply(tempmat,2,function(x) return(mean(x,na.rm=TRUE)))
  
  
  if (coef(lm(tempmat[,2] ~ tempmat[,1]))[2] > 0){
    angle <- 2*pi - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }else{
    angle <- - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }
  
  rotate.angle[ii] <- angle
  rotation.mat <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
  
  rotated.points <- matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
  for (index in 1:length(nonna)){
    prep_vec <- tempmat[nonna[index],]
    prep_vec <- prep_vec - tempmean
    rotated.points[nonna[index],] <-  rotation.mat   %*% prep_vec
  }
  tempvec <- rep(NA, length(data_vec))
  for (jj in 1:length(nonna)){
    tempvec[nonna[jj]] <- rotated.points[nonna[jj],1]
  }
  ton.ang.logmap[,ii] <- tempvec
}
#check mean of logmapped data = 0 
apply(ton.ang.logmap,2,function(x){mean(x, na.rm=TRUE)})
#subtract the mean function from both wind direction and speed 
ton.spd.demean <- t(apply(data$ton.data$spd.data, 1, function(x){ x - ton_mean_speed}))



#Covariance of wind direction 
ton.ang.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(ton.ang.logmap))
ton.ang.FPCA <- FPCA(ton.ang.L$Ly, ton.ang.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid=24))
#G(s,s)
par(mar = c(1, 1, 1, 1))
image2D(ton.ang.FPCA$fittedCov,colkey = FALSE,xlab='',ylab='',
        xaxt='n',yaxt='n')
v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(ton.ang.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_minimal() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
# timegrid <- expand.grid(seq(0,23,length.out = 101),seq(0,23,length.out = 101))[c(2,1)]
# ton.ang.cov <- as.data.frame(cbind(timegrid, as.vector(ton.ang.FPCA$fittedCov)))
# colnames(ton.ang.cov) <- c('x','y','cov')
# ggplot(ton.ang.cov, aes(x = x, y = y, z = cov)) + geom_tile(aes(fill = cov))
#   stat_contour(geom = 'polygon', aes(fill = ..level..)) + 
#   geom_tile(aes(fill = cov)) +
#   stat_contour(bins = 15) 


#Covariance of wind speed
ton.spd.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(ton.spd.demean))
ton.spd.FPCA <- FPCA(ton.spd.L$Ly, ton.spd.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid= 24))

image2D(ton.spd.FPCA$fittedCov,colkey = FALSE,xlab='',ylab='',
        xaxt='n',yaxt='n')
v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(ton.spd.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_minimal() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

#############
#perform normalization 
ton.ang.FPCA <- FPCA(ton.ang.L$Ly, ton.ang.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
ton.spd.FPCA <- FPCA(ton.spd.L$Ly, ton.spd.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth'))

ton.spd.nlized  <- t(apply(ton.spd.demean,1,function(x){ x / sqrt(diag(ton.spd.FPCA$fittedCov))}))
ton.ang.nlized  <- t(apply(ton.ang.logmap,1,function(x){ x / sqrt(diag(ton.ang.FPCA$fittedCov))}))
nlized.data <- cbind(ton.ang.nlized,ton.spd.nlized)
ton.nlized.L   <- MakeFPCAInputs(IDs = rep(1:651, each=48), 
                                 tVec=rep(seq(0,47,1),651), t(nlized.data))
ton.nlized.FPCA <- FPCA(ton.nlized.L$Ly, ton.nlized.L$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth',nRegGrid = 48))
ton.nlized.FPCA$cumFVE
#Taking leading 6 FPCS
ton.spd.fpc <- ton.ang.fpc <- matrix(NA, ncol = 4, nrow = 24)
ton.ang.fpc.plots <- list()
ton.spd.fpc.plots <- list()
for (jj in 1:4){
  ton.ang.fpc[,jj] <- ton.nlized.FPCA$phi[1:24,jj]
  ton.spd.fpc[,jj] <- ton.nlized.FPCA$phi[25:48,jj]
}




#Estimate the mean functions of Halifax wind speed and direction 
#Mean function of wind direction (angle)
ntimes = 24
hal_mean_angle <- rep(NA,ntimes)
for (tt in 1:ntimes){
  #Conert radians to degree
  data_vec      <- (data$hal.data$ang.data * 10)[,tt]
  hal_mean_angle[tt] <-  optimize(distance.fun, data_vec = data_vec,c(0,360))$minimum
}

#mean function of wind speed of Halifax 
timemat <- matrix(rep(0:23,651),byrow=TRUE, nrow = 651)
hal_mean_speed <- llr.mean(time = seq(0,23,1),data = data$hal.data$spd.data, bw = 3)
plot(seq(0,23,1), hal_mean_speed,type='l')
###############
#Perform log-mappping to direction data 
hal.ang.logmap <- matrix(NA , nrow = nrow(data$hal.data$ang.data), ncol = ncol(data$hal.data$ang.data))
rotate.angle <- rep(NA,24)
for (ii in 1:24){
  print(paste('Hour',ii))
  data_vec          <- data$hal.data$ang.data[,ii] * 10
  
  mean_coord     <- c(cos(hal_mean_angle[ii] * pi / 180), sin(hal_mean_angle[ii] * pi / 180))
  
  nonna   <- which(!is.na(data_vec))
  tempmat <- matrix(NA,nrow = length(data_vec),ncol =2 )
  for (index in 1:length(nonna)){
    x_vec <- c(cos(data_vec[nonna[index]] * pi / 180),sin(data_vec[nonna[index]] * pi / 180))
    vec_u <- x_vec - as.vector(mean_coord %*% x_vec) * mean_coord
    v_vec <- (vec_u / sqrt(as.vector(vec_u %*% vec_u))) * acos(as.vector(mean_coord %*% x_vec))
    tempmat[nonna[index],] <- v_vec
  }
  
  tempmean <- apply(tempmat,2,function(x) return(mean(x,na.rm=TRUE)))
  
  
  if (coef(lm(tempmat[,2] ~ tempmat[,1]))[2] > 0){
    angle <- 2*pi - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }else{
    angle <- - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }
  
  rotate.angle[ii] <- angle
  rotation.mat <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
  
  rotated.points <- matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
  for (index in 1:length(nonna)){
    prep_vec <- tempmat[nonna[index],]
    prep_vec <- prep_vec - tempmean
    rotated.points[nonna[index],] <-  rotation.mat   %*% prep_vec
  }
  tempvec <- rep(NA, length(data_vec))
  for (jj in 1:length(nonna)){
    tempvec[nonna[jj]] <- rotated.points[nonna[jj],1]
  }
  hal.ang.logmap[,ii] <- tempvec
}
#check mean of logmapped data = 0 
apply(hal.ang.logmap,2,function(x){mean(x, na.rm=TRUE)})
#subtract the mean function from both wind direction and speed 
hal.spd.demean <- t(apply(data$hal.data$spd.data, 1, function(x){ x - hal_mean_speed}))


#Covariance of wind direction 
hal.ang.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(hal.ang.logmap))
hal.ang.FPCA <- FPCA(hal.ang.L$Ly, hal.ang.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid=24))
#G(s,s)
par(mar = c(1, 1, 1, 1))
image2D(hal.ang.FPCA$fittedCov,colkey = FALSE,xlab='',ylab='',
        xaxt='n',yaxt='n')
v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(hal.ang.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_minimal() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
# timegrid <- expand.grid(seq(0,23,length.out = 101),seq(0,23,length.out = 101))[c(2,1)]
# hal.ang.cov <- as.data.frame(cbind(timegrid, as.vector(hal.ang.FPCA$fittedCov)))
# colnames(hal.ang.cov) <- c('x','y','cov')
# ggplot(hal.ang.cov, aes(x = x, y = y, z = cov)) + geom_tile(aes(fill = cov))
#   stat_contour(geom = 'polygon', aes(fill = ..level..)) + 
#   geom_tile(aes(fill = cov)) +
#   stat_contour(bins = 15) 


#Covariance of wind speed
hal.spd.L    <- MakeFPCAInputs(IDs = rep(1:651, each=24), 
                               tVec=rep(seq(0,23,1),651), t(hal.spd.demean))
hal.spd.FPCA <- FPCA(hal.spd.L$Ly, hal.spd.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth',nRegGrid = 24))

image2D(hal.spd.FPCA$fittedCov,colkey = FALSE,xlab='',ylab='',
        xaxt='n',yaxt='n')
v.df <- data.frame(time = seq(0,23,length.out = 24), v = diag(hal.spd.FPCA$fittedCov))
ggplot(data = v.df, aes(x = time, y = v),lwd=2) + geom_line(lwd=1.2) + theme_minimal() + 
  ylab('') + xlab('') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

#############
#perform normalization 
hal.ang.FPCA <- FPCA(hal.ang.L$Ly, hal.ang.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
hal.spd.FPCA <- FPCA(hal.spd.L$Ly, hal.spd.L$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth'))
hal.spd.nlized  <- t(apply(hal.spd.demean,1,function(x){ x / sqrt(diag(hal.spd.FPCA$fittedCov))}))
hal.ang.nlized  <- t(apply(hal.ang.logmap,1,function(x){ x / sqrt(diag(hal.ang.FPCA$fittedCov))}))
nlized.data <- cbind(hal.ang.nlized,hal.spd.nlized)
hal.nlized.L   <- MakeFPCAInputs(IDs = rep(1:651, each=48), 
                                 tVec=rep(seq(0,47,1),651), t(nlized.data))
hal.nlized.FPCA <- FPCA(hal.nlized.L$Ly, hal.nlized.L$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth',nRegGrid = 48))
hal.nlized.FPCA$cumFVE
#Taking leading 6 FPCS
hal.spd.fpc <- hal.ang.fpc <- matrix(NA, ncol = 4, nrow = 24)
hal.ang.fpc.plots <- list()
hal.spd.fpc.plots <- list()
for (jj in 1:4){
  hal.ang.fpc[,jj] <- hal.nlized.FPCA$phi[1:24,jj]
  hal.spd.fpc[,jj] <- hal.nlized.FPCA$phi[25:48,jj]
}

####################
#FPCA_MVDM finished, the subsequent codes are producing plots
######
van.nlized.FPCA$cumFVE
ton.nlized.FPCA$cumFVE
hal.nlized.FPCA$cumFVE

#fpc plots of all three cities

ang.fpcplots <- list()
for (jj in 1:4){
  plot.df <- data.frame(time = rep(seq(0,23,length.out = 24),3),
                   fpc  = c(van.ang.fpc[,jj], ton.ang.fpc[,jj], hal.ang.fpc[,jj]),
                   city = rep(c('Vancouver','Toronto','Halifax'),each = 24))
  ang.fpcplots[[jj]] <- ggplot(data = plot.df,aes(x = time, y = fpc,group = city)) + 
                geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
                geom_hline(yintercept = 0,linetype = 'dashed') + 
                theme_classic() + xlab('') + ylab('') +
                theme(axis.text.y  = element_text(size = 16, face = 'bold')) + 
                scale_x_continuous(breaks=seq(0,24,4),
                                   labels=seq(0,24,4)) + theme(axis.text.x = element_text(size = 16, face='bold'))
}

spd.fpcplots <- list()
for (jj in 1:4){
  plot.df <- data.frame(time = rep(seq(0,23,length.out = 24),3),
                        fpc  = c(van.spd.fpc[,jj], ton.spd.fpc[,jj], hal.spd.fpc[,jj]),
                        city = rep(c('Vancouver','Toronto','Halifax'),each = 24))
  spd.fpcplots[[jj]] <- ggplot(data = plot.df,aes(x = time, y = fpc,group = city)) + 
    geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
    geom_hline(yintercept = 0,linetype='dashed') + 
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.y  = element_text(size = 16, face = 'bold')) + 
    scale_x_continuous(breaks=seq(0,24,4),
                       labels=seq(0,24,4)) + theme(axis.text.x = element_text(size = 16, face='bold'))
}

#Figure 8
plot_grid(ang.fpcplots[[1]],ang.fpcplots[[2]],ang.fpcplots[[3]],ang.fpcplots[[4]],
          spd.fpcplots[[1]],spd.fpcplots[[2]],spd.fpcplots[[3]],spd.fpcplots[[4]],
          nrow = 2)


#Reconstruction 
recon_Z_van = van.nlized.FPCA$xiEst %*% t(van.nlized.FPCA$phi)
recon_Z_ton = ton.nlized.FPCA$xiEst %*% t(ton.nlized.FPCA$phi)
recon_Z_hal = hal.nlized.FPCA$xiEst %*% t(hal.nlized.FPCA$phi)

recon_V_van_ang = t(apply(recon_Z_van[,1:24],1,function(x){ x * sqrt(diag(van.ang.FPCA$fittedCov))}))
recon_V_van_spd = t(apply(recon_Z_van[,25:48],1,function(x){ x * sqrt(diag(van.spd.FPCA$fittedCov))}))
recon_V_ton_ang = t(apply(recon_Z_ton[,1:24],1,function(x){ x * sqrt(diag(ton.ang.FPCA$fittedCov))}))
recon_V_ton_spd = t(apply(recon_Z_ton[,25:48],1,function(x){ x * sqrt(diag(ton.spd.FPCA$fittedCov))}))
recon_V_hal_ang = t(apply(recon_Z_hal[,1:24],1,function(x){ x * sqrt(diag(hal.ang.FPCA$fittedCov))}))
recon_V_hal_spd = t(apply(recon_Z_hal[,25:48],1,function(x){ x * sqrt(diag(hal.spd.FPCA$fittedCov))}))
recon_X_van = recon_V_van_spd + van_mean_speed
recon_X_ton = recon_V_ton_spd + ton_mean_speed
recon_X_hal = recon_V_hal_spd + hal_mean_speed


angle.est  <- matrix(NA, nrow = 651, ncol = 24)
angle.temp <- recon_V_van_ang
for (tt in 1:24){
  data_vec  = (data$van.data$ang.data * 10)[,tt]
  mean_rand = van_mean_angle[tt]
  #Get cartesian coordinates of mean radian
  mean_vec <- c(cos(mean_rand), sin(mean_rand))
  for (index in 1:length(data_vec)){
    data_rand <- data_vec[index] * pi / 180
    #Vector from original to data point on circle
    q_vec <- c(cos(data_rand),sin(data_rand))
    #Projection of q onto p
    proj_p_q <- as.numeric(mean_vec %*% q_vec) * mean_vec
    #Tangent vector from proj_p_q to q
    u = q_vec - proj_p_q
    #normalized u and scale it with great circle distance (geodesic distance): acos(p,q)
    tang_vec <- u / sqrt(sum(u^2)) * c( acos(mean_vec %*% q_vec))
    
    if(is.na(data_vec[index])){
      if (mean_rand > 0 && mean_rand < pi){
        if ( angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( angle.temp[index,tt] <=0 && angle.temp[index,tt] >= -pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
      
    }else{
      if (mean_rand > 0 && mean_rand < pi){
        if ( (data_rand - mean_rand) >=0 && (data_rand - mean_rand) <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( (data_rand - mean_rand) <=0 && (data_rand - mean_rand) >= -pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
    }
    
  }
  
}
angle.est <- t(apply(angle.temp,1, function(x) {x + van_mean_angle * pi / 180}))
recon_W_van = angle.est


van_ori_data_ang    = data.frame(x = seq(0,23,1), y = data$van.data$ang.data[402,] * 10)
van_recon_data_ang  = data.frame(x = seq(0,23,1), y = recon_W_van[402,] * 180 / pi)

van_ori_data_spd    = data.frame(x = seq(0,23,length.out =24), y = data$van.data$spd.data[23,])
van_recon_data_spd  = data.frame(x = seq(0,23,length.out =24),y = recon_X_van[23,])


#Toronto
angle.est  <- matrix(NA, nrow = 651, ncol = 24)
angle.temp <- recon_V_ton_ang
for (tt in 1:24){
  data_vec  = (data$ton.data$ang.data * 10)[,tt]
  mean_rand = ton_mean_angle[tt]
  #Get cartesian coordinates of mean radian
  mean_vec <- c(cos(mean_rand), sin(mean_rand))
  for (index in 1:length(data_vec)){
    data_rand <- data_vec[index] * pi / 180
    #Vector from original to data point on circle
    q_vec <- c(cos(data_rand),sin(data_rand))
    #Projection of q onto p
    proj_p_q <- as.numeric(mean_vec %*% q_vec) * mean_vec
    #Tangent vector from proj_p_q to q
    u = q_vec - proj_p_q
    #normalized u and scale it with great circle distance (geodesic distance): acos(p,q)
    tang_vec <- u / sqrt(sum(u^2)) * c( acos(mean_vec %*% q_vec))
    
    if(is.na(data_vec[index])){
      if (mean_rand > 0 && mean_rand < pi){
        if ( angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( angle.temp[index,tt] <=0 && angle.temp[index,tt] >= -pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
      
    }else{
      if (mean_rand > 0 && mean_rand < pi){
        if ( (data_rand - mean_rand) >=0 && (data_rand - mean_rand) <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( (data_rand - mean_rand) <=0 && (data_rand - mean_rand) >= -pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
    }
    
  }
  
}
angle.est <- t(apply(angle.temp,1, function(x) {x + ton_mean_angle * pi / 180}))
recon_W_ton = angle.est

missing.df <- apply(data$ton.data$ang.data,1,function(x){length(which(is.na(x)))}) / 24
sort(missing.df, index.return=TRUE)

ton_ori_data_ang    = data.frame(x = seq(0,23,1), y = data$ton.data$ang.data[286,] * 10)
ton_recon_data_ang  = data.frame(x = seq(0,23,1), y = recon_W_ton[286,] * 180 / pi)


#Halifax
angle.est  <- matrix(NA, nrow = 651, ncol = 24)
angle.temp <- recon_V_hal_ang
for (tt in 1:24){
  data_vec  = (data$hal.data$ang.data * 10)[,tt]
  mean_rand = hal_mean_angle[tt]
  #Get cartesian coordinates of mean radian
  mean_vec <- c(cos(mean_rand), sin(mean_rand))
  for (index in 1:length(data_vec)){
    data_rand <- data_vec[index] * pi / 180
    #Vector from original to data point on circle
    q_vec <- c(cos(data_rand),sin(data_rand))
    #Projection of q onto p
    proj_p_q <- as.numeric(mean_vec %*% q_vec) * mean_vec
    #Tangent vector from proj_p_q to q
    u = q_vec - proj_p_q
    #normalized u and scale it with great circle distance (geodesic distance): acos(p,q)
    tang_vec <- u / sqrt(sum(u^2)) * c( acos(mean_vec %*% q_vec))
    
    if(is.na(data_vec[index])){
      if (mean_rand > 0 && mean_rand < pi){
        if ( angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (angle.temp[index,tt] >=0 && angle.temp[index,tt] <= pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( angle.temp[index,tt] <=0 && angle.temp[index,tt] >= -pi){
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
      
    }else{
      if (mean_rand > 0 && mean_rand < pi){
        if ( (data_rand - mean_rand) >=0 && (data_rand - mean_rand) <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(-angle.temp[index,tt]) * mean_vec + sin(-angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand ==0){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if( mean_rand == (pi)){
        if (data_rand >=0 && data_rand <= pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      if (mean_rand > pi && mean_rand < (2 * pi)){
        if ( (data_rand - mean_rand) <=0 && (data_rand - mean_rand) >= -pi){
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }else{
          angle.est[index,tt] <- acos((cos(angle.temp[index,tt]) * mean_vec + sin(angle.temp[index,tt]) * c(-sin(mean_rand),cos(mean_rand)))[1])
        }
      }
      
    }
    
  }
  
}
angle.est <- t(apply(angle.temp,1, function(x) {x + hal_mean_angle * pi / 180}))
recon_W_hal = angle.est

missing.df <- apply(data$hal.data$ang.data,1,function(x){length(which(is.na(x)))}) / 24
sort(missing.df, index.return=TRUE)

hal_ori_data_ang    = data.frame(x = seq(0,23,1), y = data$hal.data$ang.data[130,] * 10)
hal_recon_data_ang  = data.frame(x = seq(0,23,1), y = recon_W_hal[130,] * 180 / pi)



van_plot_recon_ang = ggplot(van_ori_data_ang, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#619CFF') + geom_line(data = van_recon_data_ang,color='#619CFF', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold") ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) + 
  scale_y_continuous(breaks = seq(0,360,60))

van_plot_recon_spd = ggplot(van_ori_data_spd, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#619CFF') + geom_line(data = van_recon_data_spd,color='#619CFF', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold") ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4))



ton_plot_recon_ang = ggplot(ton_ori_data_ang, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#F8766D') + geom_line(data = ton_recon_data_ang,color='#F8766D', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold") ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) + 
  scale_y_continuous(breaks = seq(0,360,30))


ton_ori_data_spd    = data.frame(x = seq(0,23,length.out =24), y = data$ton.data$spd.data[286,])
ton_recon_data_spd  = data.frame(x = seq(0,23,length.out =24),y = recon_X_ton[286,])
ton_plot_recon_spd = ggplot(ton_ori_data_spd, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#F8766D') + geom_line(data = ton_recon_data_spd,color='#F8766D', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold") ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4))


hal_plot_recon_ang = ggplot(hal_ori_data_ang, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#00BA38') + geom_line(data = hal_recon_data_ang,color='#00BA38', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold") ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') + annotate('point',x = 7, y = 20, color='#00BA38', shape = 16, size = 8)+ 
        annotate('text', x = 7 , y = 45, label='A',size = 7) + 
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) + 
  scale_y_continuous(breaks = seq(0,360,60),labels=seq(0,360,60))

hal_ori_data_spd    = data.frame(x = seq(0,23,length.out =24), y = data$hal.data$spd.data[286,])
hal_recon_data_spd  = data.frame(x = seq(0,23,length.out =24),y = recon_X_hal[286,])
hal_plot_recon_spd = ggplot(hal_ori_data_spd, aes(x = x, y = y)) + geom_point(shape=0,size = 2,color='#00BA38') + geom_line(data = hal_recon_data_spd,color='#00BA38', aes(x = x, y = y),size = 1.2)+ 
  theme(axis.text = element_blank()) +theme_classic() + 
  theme(axis.text=element_text(size=12,face='bold'),
        axis.title=element_text(size=14,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks=seq(0,24,4),labels=seq(0,24,4)) +
  scale_y_continuous(breaks = seq(0,25,5),labels=seq(0,25,5))


plot_grid(hal_plot_recon_ang,van_plot_recon_ang, ton_plot_recon_ang,
          hal_plot_recon_spd,van_plot_recon_spd, ton_plot_recon_spd,nrow =2,align='v')


#
#-------------------
#Finite sample performance of mixed FPCA
#Number of simulated data in each simulation 
nsimdata   <- 500
#Number of finite sample (simulation) to perform 
nsim       <- 100 
#Number of components to take from the real data analysis 
ncomponent <- 7
#Preallocate space for storing estimated FPC of simulated data
van.simfpc <- list()
ton.simfpc <- list()
hal.simfpc <- list()
for (jj in 1:ncomponent){
  van.simfpc[[jj]] <- ton.simfpc[[jj]] <- hal.simfpc[[jj]] <- matrix(NA, nrow = nsim, ncol = 202)
}


for (iter in 1:nsim){
  print(iter)
  #Simulate scores for each city based on the eigenvalues 
  #Preallocate space 
  van.score  <- matrix(NA, nrow = nsimdata , ncol = ncomponent)
  ton.score  <- matrix(NA, nrow = nsimdata , ncol = ncomponent)
  hal.score  <- matrix(NA, nrow = nsimdata , ncol = ncomponent)
  #Simulate scores, scores ~ N(0, eigenvalues)
  for (ii in 1:ncomponent){
    van.score[,ii] <- rnorm(n = nsimdata, mean = 0 , sd = sqrt(van.nlized.FPCA$lambda[ii]))
  }
  for (ii in 1:ncomponent){
    ton.score[,ii] <- rnorm(n = nsimdata, mean = 0 , sd = sqrt(ton.nlized.FPCA$lambda[ii]))
  }
  for (ii in 1:ncomponent){
    hal.score[,ii] <- rnorm(n = nsimdata, mean = 0 , sd = sqrt(hal.nlized.FPCA$lambda[ii]))
  }
  
  van.logsimdata <- van.score %*%  t(apply(van.nlized.FPCA$phi,2, function(x) {x / sqrt(sum(x^2))})[,1:ncomponent])
  ton.logsimdata <- ton.score %*%  t(apply(ton.nlized.FPCA$phi,2, function(x) {x / sqrt(sum(x^2))})[,1:ncomponent])
  hal.logsimdata <- hal.score %*%  t(apply(hal.nlized.FPCA$phi,2, function(x) {x / sqrt(sum(x^2))})[,1:ncomponent])
  
  #Perform mixed FPCA on simulated data 
  L1      <- MakeFPCAInputs(IDs = rep(1:nsimdata, each=202), tVec=rep(seq(0,201,1),nsimdata), t(van.logsimdata)) 
  L2      <- MakeFPCAInputs(IDs = rep(1:nsimdata, each=202), tVec=rep(seq(0,201,1),nsimdata), t(ton.logsimdata)) 
  L3      <- MakeFPCAInputs(IDs = rep(1:nsimdata, each=202), tVec=rep(seq(0,201,1),nsimdata), t(hal.logsimdata)) 
  
  (van.simpca <- FPCA(L1$Ly, L1$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth')))
  (ton.simpca <- FPCA(L2$Ly, L2$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth')))
  (hal.simpca <- FPCA(L3$Ly, L3$Lt, opt=list(FVEthreshold = 0.99,methodMuCovEst = 'smooth')))
  
  for (jj in 1:6){
    van.simfpc[[jj]][iter,] <- van.simpca$phi[,jj] * sqrt(sum(van.nlized.FPCA$phi[,jj]^2))
    ton.simfpc[[jj]][iter,] <- ton.simpca$phi[,jj] * sqrt(sum(ton.nlized.FPCA$phi[,jj]^2))
    hal.simfpc[[jj]][iter,] <- hal.simpca$phi[,jj] * sqrt(sum(hal.nlized.FPCA$phi[,jj]^2))
  }
}

#Plotting finite sample performance 
#Vancouver 
#Vancouver wind direction
van.sim.fpc.angle <- list()

for (jj in 1:6){
  van.sim.fpc.angle[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  van.sim.fpc.angle[[jj]] <- van.simfpc[[jj]][,1:101]
}

for (ii in 1:nsim){
  if ( van.sim.fpc.angle[[1]][ii,1] < 0 ){
    van.sim.fpc.angle[[1]][ii,] <- -van.sim.fpc.angle[[1]][ii,]
  }
  
}

for (ii in 1:nsim){
  if ( van.sim.fpc.angle[[2]][ii,1] > 0 ){
    van.sim.fpc.angle[[2]][ii,] <- -van.sim.fpc.angle[[2]][ii,]
  }
  
}


#Vancouver wind speed
van.sim.fpc.speed <- list()
for (jj in 1:6){
  van.sim.fpc.speed[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  van.sim.fpc.speed[[jj]] <- van.simfpc[[jj]][,102:202]
}


for (ii in 1:nsim){
  if (van.sim.fpc.speed[[1]][ii,1] < 0){
    van.sim.fpc.speed[[1]][ii,] <- -van.sim.fpc.speed[[1]][ii,]
  }
}

for (ii in 1:nsim){
  if (van.sim.fpc.speed[[2]][ii,1] > 0){
    van.sim.fpc.speed[[2]][ii,] <- -van.sim.fpc.speed[[2]][ii,]
  }
}



van.direction.plots <- list()
van.speed.plots     <- list()
for (jj in 1:6){
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = van.nlized.FPCA$phi[1:101,jj],
                        lower = apply(van.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.025)}), upper = apply(van.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.975)}))
  van.direction.plots[[jj]] <- ggplot(data = plot.df, aes(x = time, y = fpc)) + 
                                geom_line(aes(x = time, y = fpc),col='blue',size = 1.5) + 
                                geom_ribbon(aes(ymin=lower, ymax=upper), 
                                            alpha=0.3,)+
                                theme_classic() + xlab('') + ylab('') + 
                                theme(axis.text.x = element_blank(),
                                      axis.text.y = element_text(face='bold',size = 12),
                                      axis.title  = element_text(face='bold',size = 14)) 
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = van.nlized.FPCA$phi[102:202,jj],
                        lower = apply(van.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.025)}),upper = apply(van.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.975)}))
  van.speed.plots[[jj]] <-  ggplot(data = plot.df, aes(x = time, y = fpc)) + 
    geom_line(aes(x = time, y = fpc),col='blue',size = 1.5) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), 
                alpha=0.3,)+
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold',size = 12),
          axis.title  = element_text(face='bold',size = 14)) 
  
  
}

#Toronto
#Toronto wind direction
ton.sim.fpc.angle <- list()
for (jj in 1:6){
  ton.sim.fpc.angle[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  ton.sim.fpc.angle[[jj]] <- ton.simfpc[[jj]][,1:101]
}

for (ii in 1:nsim){
  if ( ton.sim.fpc.angle[[1]][ii,1] > 0 ){
    ton.sim.fpc.angle[[1]][ii,] <- -ton.sim.fpc.angle[[1]][ii,]
  }
  
}

for (ii in 1:nsim){
  if ( ton.sim.fpc.angle[[2]][ii,1] > 0 ){
    ton.sim.fpc.angle[[2]][ii,] <- -ton.sim.fpc.angle[[2]][ii,]
  }
  
}


#Toronto wind speed
ton.sim.fpc.speed <- list()
for (jj in 1:6){
  ton.sim.fpc.speed[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  ton.sim.fpc.speed[[jj]] <- ton.simfpc[[jj]][,102:202]
}


for (ii in 1:nsim){
  if (ton.sim.fpc.speed[[1]][ii,1] > 0){
    ton.sim.fpc.speed[[1]][ii,] <- -ton.sim.fpc.speed[[1]][ii,]
  }
}

for (ii in 1:nsim){
  if (ton.sim.fpc.speed[[2]][ii,101] < 0){
    ton.sim.fpc.speed[[2]][ii,] <- -ton.sim.fpc.speed[[2]][ii,]
  }
}

ton.direction.plots <- list()
ton.speed.plots     <- list()
for (jj in 1:6){
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = ton.nlized.FPCA$phi[1:101,jj],
                        lower = apply(ton.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.025)}), upper = apply(ton.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.975)}))
  ton.direction.plots[[jj]] <- ggplot(data = plot.df, aes(x = time, y = fpc)) + 
    geom_line(aes(x = time, y = fpc),col='green',size = 1.5) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), 
                alpha=0.3,)+
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold',size = 12),
          axis.title  = element_text(face='bold',size = 14)) 
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = ton.nlized.FPCA$phi[102:202,jj],
                        lower = apply(ton.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.025)}),upper = apply(ton.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.975)}))
  ton.speed.plots[[jj]] <- ggplot(data = plot.df, aes(x = time, y = fpc)) + 
    geom_line(aes(x = time, y = fpc),col='green',size = 1.5) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), 
                alpha=0.3,)+
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold',size = 12),
          axis.title  = element_text(face='bold',size = 14)) 
  
}
#Halifax 
#Halifax wind direction
hal.sim.fpc.angle <- list()
for (jj in 1:6){
  hal.sim.fpc.angle[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  hal.sim.fpc.angle[[jj]] <- hal.simfpc[[jj]][,1:101]
}

for (ii in 1:nsim){
  if ( hal.sim.fpc.angle[[1]][ii,40] > 0 ){
    hal.sim.fpc.angle[[1]][ii,] <- -hal.sim.fpc.angle[[1]][ii,]
  }
  
}


for (ii in 1:nsim){
  if ( hal.sim.fpc.angle[[2]][ii,40] > 0 ){
    hal.sim.fpc.angle[[2]][ii,] <- -hal.sim.fpc.angle[[2]][ii,]
  }
  
}


#Halifax wind speed
hal.sim.fpc.speed <- list()
for (jj in 1:6){
  hal.sim.fpc.speed[[jj]] <- matrix(NA, nrow = nsim, ncol = 101)
  hal.sim.fpc.speed[[jj]] <- hal.simfpc[[jj]][,102:202]
}


for (ii in 1:nsim){
  if (hal.sim.fpc.speed[[1]][ii,40] < 0){
    hal.sim.fpc.speed[[1]][ii,] <- -hal.sim.fpc.speed[[1]][ii,]
  }
}

for (ii in 1:nsim){
  if (hal.sim.fpc.speed[[2]][ii,40] > 0){
    hal.sim.fpc.speed[[2]][ii,] <- -hal.sim.fpc.speed[[2]][ii,]
  }
}

hal.direction.plots <- list()
hal.speed.plots     <- list()
for (jj in 1:6){
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = hal.nlized.FPCA$phi[1:101,jj],
                        lower = apply(hal.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.025)}), upper = apply(hal.sim.fpc.angle[[jj]],2,function(x){quantile(x,prob =0.975)}))
  hal.direction.plots[[jj]] <- ggplot(data = plot.df, aes(x = time, y = fpc)) + 
    geom_line(aes(x = time, y = fpc),col='red',size = 1.5) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), 
                alpha=0.3,)+
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold',size = 12),
          axis.title  = element_text(face='bold',size = 14)) 
  
  plot.df <- data.frame(time = seq(0,1,length.out = 101), fpc = hal.nlized.FPCA$phi[102:202,jj],
                        lower = apply(hal.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.025)}),upper = apply(hal.sim.fpc.speed[[jj]],2,function(x){quantile(x,prob =0.975)}))
  hal.speed.plots[[jj]] <- ggplot(data = plot.df, aes(x = time, y = fpc)) + 
    geom_line(aes(x = time, y = fpc),col='red',size = 1.5) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), 
                alpha=0.3,)+
    theme_classic() + xlab('') + ylab('') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold',size = 12),
          axis.title  = element_text(face='bold',size = 14)) 
  
}


plot_grid(van.direction.plots[[1]],van.direction.plots[[2]],
          ton.direction.plots[[1]],ton.direction.plots[[2]],
          hal.direction.plots[[1]],hal.direction.plots[[2]],
          nrow = 3)
plot_grid(van.speed.plots[[1]],van.speed.plots[[2]],
          ton.speed.plots[[1]],ton.speed.plots[[2]],
          hal.speed.plots[[1]],hal.speed.plots[[2]],
          nrow = 3)


#merged plot

fpc1.ang.df <- data.frame(time = rep(seq(0,23,length.out = 101),3 ),
                          fpc = c(van.nlized.FPCA$phi[1:101,1], ton.nlized.FPCA$phi[1:101,1], hal.nlized.FPCA$phi[1:101,1]),
                          lower = c(apply(van.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(ton.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(hal.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.025)})),
                          upper = c(apply(van.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(ton.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(hal.sim.fpc.angle[[1]],2,function(x){quantile(x,prob =0.975)})),
                          city = rep(c('Vancouver','Toronto','Halifax'), each = 101))
fpc1.spd.df <- data.frame(time = rep(seq(0,23,length.out = 101),3 ),
                          fpc = c(van.nlized.FPCA$phi[102:202,1], ton.nlized.FPCA$phi[102:202,1], hal.nlized.FPCA$phi[102:202,1]),
                          lower = c(apply(van.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(ton.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(hal.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.025)})),
                          upper = c(apply(van.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(ton.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(hal.sim.fpc.speed[[1]],2,function(x){quantile(x,prob =0.975)})),
                          city = rep(c('Vancouver','Toronto','Halifax'), each = 101))
fpc2.ang.df <- data.frame(time = rep(seq(0,23,length.out = 101),3 ),
                          fpc = c(van.nlized.FPCA$phi[1:101,2], ton.nlized.FPCA$phi[1:101,2], hal.nlized.FPCA$phi[1:101,2]),
                          lower = c(apply(van.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(ton.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(hal.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.025)})),
                          upper = c(apply(van.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(ton.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(hal.sim.fpc.angle[[2]],2,function(x){quantile(x,prob =0.975)})),
                          city = rep(c('Vancouver','Toronto','Halifax'), each = 101))
fpc2.spd.df <- data.frame(time = rep(seq(0,23,length.out = 101),3 ),
                          fpc = c(van.nlized.FPCA$phi[102:202,2], ton.nlized.FPCA$phi[102:202,2], hal.nlized.FPCA$phi[102:202,2]),
                          lower = c(apply(van.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(ton.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.025)}),
                                    apply(hal.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.025)})),
                          upper = c(apply(van.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(ton.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.975)}),
                                    apply(hal.sim.fpc.speed[[2]],2,function(x){quantile(x,prob =0.975)})),
                          city = rep(c('Vancouver','Toronto','Halifax'), each = 101))

plot1 <- ggplot(data = fpc1.ang.df,aes(x = time, y = fpc,group = city)) + 
  geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
  geom_hline(yintercept = 0,linetype='dashed') + 
  theme_classic() + xlab('') + ylab('') + 
  theme(axis.text.y  = element_text(size = 12, face = 'bold'),
        legend.position = 'none') + 
  geom_ribbon(aes(ymin=lower, ymax=upper,fill=city,col=city), alpha=0.3,)+
  scale_x_continuous(breaks=seq(0,23,4),
                     labels=seq(0,23,4)) + theme(axis.text.x = element_text(size = 11, face='bold'))

plot2 <- ggplot(data = fpc1.spd.df,aes(x = time, y = fpc,group = city)) + 
  geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
  geom_hline(yintercept = 0,linetype='dashed') + 
  theme_classic() + xlab('') + ylab('') + 
  theme(axis.text.y  = element_text(size = 12, face = 'bold'),
        legend.position = 'none') + 
  geom_ribbon(aes(ymin=lower, ymax=upper,fill=city,col=city), alpha=0.3,)+
  scale_x_continuous(breaks=seq(0,23,4),
                     labels=seq(0,23,4)) + theme(axis.text.x = element_text(size = 11, face='bold'))

plot3 <- ggplot(data = fpc2.ang.df,aes(x = time, y = fpc,group = city)) + 
  geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
  geom_hline(yintercept = 0,linetype='dashed') + 
  theme_classic() + xlab('') + ylab('') + 
  theme(axis.text.y  = element_text(size = 12, face = 'bold'),
        legend.position = 'none') + 
  geom_ribbon(aes(ymin=lower, ymax=upper,fill=city,col=city), alpha=0.3,)+
  scale_x_continuous(breaks=seq(0,23,4),
                     labels=seq(0,23,4)) + theme(axis.text.x = element_text(size = 11, face='bold'))

plot4 <- ggplot(data = fpc2.spd.df,aes(x = time, y = fpc,group = city)) + 
  geom_line(aes(col=city),size = 1.2, show.legend = FALSE) + 
  geom_hline(yintercept = 0,linetype='dashed') + 
  theme_classic() + xlab('') + ylab('') + 
  theme(axis.text.y  = element_text(size = 12, face = 'bold'),
        legend.position = 'none') + 
  geom_ribbon(aes(ymin=lower, ymax=upper,fill=city,col=city), alpha=0.3,)+
  scale_x_continuous(breaks=seq(0,23,4),
                     labels=seq(0,23,4)) + theme(axis.text.x = element_text(size = 11, face='bold'))
plot_grid(plot1,plot3,
          plot2,plot4, nrow = 2)
