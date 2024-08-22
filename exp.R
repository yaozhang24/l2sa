source('utils.R')





#Load the data
load("nhanes.fish.rda")
Z <- as.numeric(nhanes.fish$fish.level == "high")
X_0 <- nhanes.fish[, c("gender", "income","income.missing", "age",  "race", "education", "smoking.ever","smoking.now")]
X_0$race <- factor(X_0$race)
X_0 <-model.matrix(~ . - 1, X_0)
X <- X_0[, -c(10)] 
Y <- log2(nhanes.fish$o.LBXTHG)


num_treated <- sum(Z==1)
num_control <- sum(Z==0)
n <- num_control + num_treated 
prop_0 <- num_control/n
prop_1 <- 1- prop_0 



#Number of folds in cross-fitting

K <- 9





#MSM

#Sensitivity parameters
Gamma_list <- seq(1, 11, by=0.5)
mean_l <- c() 



#Generate the sensitivity curve
for (Gamma in Gamma_list){
    
    #compute the EIFs
    out_lower_1 <- psi_cf_real_data(X,Y,Z,Gamma,K,"treated","lower")
    out_upper_0 <- psi_cf_real_data(X,Y,Z,Gamma,K,"control","upper")
  
    #average the EIFs
    mean_l <- c(mean_l, mean(out_lower_1-out_upper_0))
    
}

#Plot the sensitivity curve
png("real_3.png",width=7,height=7, units = 'in', res = 500)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(x=Gamma_list,y=mean_l,type="l", bty="L", xlab=expression(Gamma),ylab="ATE",main = "MSM", cex.lab = 2,cex.main = 2,
     xaxt = 'n',,yaxt = 'n', xaxs = "i", xlim = c(1,10), ylim = c(-1,2 ))
axis(1, at = seq(1,11, by = 2), las=1,cex.axis=2.0)
axis(2, at = seq(-1, 4, by = 1), las=1,cex.axis=2.0)
lines(Gamma_list,mean_l,col="paleturquoise3",lwd =5)
abline(h = 0, col="black",lwd =2)


dev.off()





#L2


#Sensitivity parameters
lambda_list <- c(seq(0,4, by=0.1))
m_l = c()
mean_l = c()


#Generate the sensitivity curve
for (lambda in lambda_list){

  #compute the EIFs
  out_lower_1_ <- psi_12_cf_real_data(X,Y,Z,lambda,K,"treated","lower")
  out_upper_0_ <- psi_12_cf_real_data(X,Y,Z,lambda,K,"control","upper")
  
  #average the EIFs
  m_l <- c(m_l,  mean(prop_1*out_lower_1_[[1]] + prop_0*out_upper_0_[[1]]))
  mean_l <- c(mean_l,  mean(out_lower_1_[[2]] - out_upper_0_[[2]]))
  

}


#Plot the sensitivity curve
png("real_4.png",width=7,height=7, units = 'in', res = 500)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(x=m_l,y=mean_l,type="l", bty="L",  xlab=expression(italic(Sigma)),ylab="ATE",
     main = expression(bold(italic(L^2)*"-analysis")), 
     cex.lab = 2,cex.main = 2,
     xaxt = 'n',,yaxt = 'n',xaxs = "i",xlim = c(1,2), ylim = c(-1,2 ))
axis(1, at = seq(1, 2, by = 0.2), las=1,cex.axis=2.0)
axis(2, at = seq(-1, 4, by = 1), las=1,cex.axis=2.0)


abline(h = 0, col="black",lwd =2)
lines(m_l,mean_l,col="tomato",lwd =5)


dev.off()


