library(glmnet)








psi_cf_real_data <-function(X, Y, Z,Gamma,K,opt1,opt2){
  
  ###MSM###

  set.seed(1234)
  
  if(opt1=="control"){
    Z <- 1-Z
  }
  
  alpha <- Gamma/(1+Gamma)
  
  if (opt2=="upper"){
    
    sign_change <- 1
    
  } else if(opt2=="lower"){
    
    sign_change <- -1
    Y <- Y*sign_change
    
  } else{
    stop("incorrect opt")
  }
  
  
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  
  folds <- split(sample(1:n),1:K)
  phi_list <- c()
  
  for (k in 1:K){ 
    i<- folds[[k]]
    
    #Split the data
    X_ave <- X[i,]
    X_eta <- X[-i,]
    
    Y_ave <- Y[i]
    Y_eta <- Y[-i]
    
    Z_ave <- Z[i]
    Z_eta <- Z[-i]
    
    q_hat <- c()
    mu_plus_hat <- c()
    mu_minus_hat <- c()
    
    
    X_eta_1 <- X_eta[Z_eta ==1,]
    Y_eta_1 <- Y_eta[Z_eta ==1]
    
    
    
    #Fitting a outcome model
    Y_model <- lm(Y~., data.frame(Y = Y_eta_1, X = X_eta_1)) 
    mean_Y <- predict(Y_model, data.frame(Y = Y_eta_1, X = X_eta_1)) 
    residual_hat <- drop (Y_eta_1) - mean_Y
    mu_X <- predict(Y_model,data.frame(Y = Y_ave, X = X_ave))
    
    
    #Fitting a treatment model
    e_model <- glmnet( X_eta,Z_eta, family= "binomial", alpha=0, lambda=0.01) 
    e_hat <- predict(e_model, X_ave, type="response")
    e_hat_2 = e_hat
    e_hat_2[e_hat_2>0.9] = 0.9
    e_hat_2[e_hat_2<0.1] = 0.1
    e_hat = e_hat_2
    
    
    mu_hat <- c()
    for (i in 1:nrow(X_ave)) {
      
      x <- X_ave[i,]
      
      #Kernel estimator
      Kx <- kernel(X_eta_1,x) 
      Wx <- Kx /sum(Kx)
      
      mean_vector <- mu_X[i] + residual_hat
      
      if (Gamma==1){
        mu_hat <- c(mu_hat, Wx %*% mean_vector)
      }
      
      #Quantile and Truncated expectations
      q_hat_x <- nw_q_xi(mean_vector,Wx,alpha,1,"gamma")
      mu_plus_hat_x <- nw(q_hat_x,mean_vector,Wx,1,"1+")
      mu_minus_hat_x <- nw(q_hat_x,mean_vector,Wx,1,"1-")
      
      
      q_hat <- c(q_hat,q_hat_x)
      mu_plus_hat <- c(mu_plus_hat,mu_plus_hat_x)
      mu_minus_hat <- c(mu_minus_hat,mu_minus_hat_x)
      
    }
    

    #EIF formula
    phi <- function(Q,E,MU_1,MU_2){
      
      
      a <- ((1-alpha)*Q + (Y_ave-Q)*(Y_ave>Q) - MU_1)*Z_ave*((1-Gamma) + Gamma/E)
      b <- (alpha*Q + (Y_ave-Q)*(Y_ave<Q) - MU_2)*Z_ave*((1-1/Gamma) + (1/Gamma)/E)
      c <- ((1-Gamma)*Z_ave+Gamma)* MU_1
      d <- ((1-1/Gamma)*Z_ave+1/Gamma)*MU_2
      a + c + b + d
    }
    
    
    if (Gamma==1){
      #Use AIPW under no unconfounding
      phi_k <- (Y_ave - mu_hat)*Z_ave/e_hat + mu_hat
      phi_k <- phi_k*sign_change
      
    } else{
      #Use EIF under MSM
      phi_k <- phi(q_hat,e_hat,mu_plus_hat,mu_minus_hat)*sign_change
    }
    

    phi_list <- c(phi_list,phi_k)
    
    
  }
  
  return(phi_list)
  
}















psi_12_cf_real_data <- function(X, Y, Z,lambda,K,opt1,opt2){
  
  ###L2-Lagrangian formulation###
  
  set.seed(1234)
  
  if(opt1=="control"){
    Z <- 1-Z
  }
  
  if (opt2=="upper"){
    
    sign_change <- -1
    
    Y <- Y*sign_change
    
  } else if(opt2=="lower"){
    
    sign_change <- 1
    
  } else{
    stop("incorrect opt")
  }
  
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  
  folds <- split(sample(1:n),1:K)
  
  phi_1_list <- c()
  phi_2_list <- c()
  
  for (k in 1:K){ 
    i<- folds[[k]]
    
    #Split the data
    X_ave <- X[i,]
    X_eta <- X[-i,]
    
    Y_ave <- Y[i]
    Y_eta <- Y[-i]
    
    Z_ave <- Z[i]
    Z_eta <- Z[-i]
    
    
    xi_hat <- c()
    y_p_hat <- c()
    y_1_hat <- c()
    E_h_1_hat <- c()
    E_h_2_hat <- c()
    mu_hat <- c()
    
    X_eta_1 <- X_eta[Z_eta ==1,]
    Y_eta_1 <- Y_eta[Z_eta ==1]
    
    
    #Fitting a outcome model
    Y_model <- lm(Y~. , data.frame(Y = Y_eta_1, X = X_eta_1)) 
    mean_Y <- predict(Y_model, data.frame(Y = Y_eta_1, X = X_eta_1)) 
    residual_hat <- drop (Y_eta_1) - mean_Y
    mu_X <- predict(Y_model,data.frame(Y = Y_ave, X = X_ave))
    
    
    #Fitting a treatment model
    e_model <- glmnet(X_eta, Z_eta, family= "binomial", alpha=0,lambda=0.01)
    e_hat <- predict(e_model, X_ave, type="response") 
    e_hat_2 = e_hat
    e_hat_2[e_hat_2>0.9] = 0.9
    e_hat_2[e_hat_2<0.1] = 0.1
    e_hat = e_hat_2
    

    if (lambda != 0){
      
      
      for (i in 1:nrow(X_ave)) {
        
        x <- X_ave[i,]
        
        #Kernel estimator
        Kx <- kernel(X_eta_1,x) 
        Wx <- Kx /sum(Kx)
        mean_vector <- mu_X[i] + residual_hat
        mu_x <- Wx %*% mean_vector
        
        
        #Threshold and Truncated expectations
        e_x <- e_hat[i]
        xi_hat_x <- nw_q_xi(mean_vector,Wx,lambda,1,"lambda",e_x)
        y_p <-nw(xi_hat_x,mean_vector,Wx,1,"p")
        y_1 <-nw(xi_hat_x,mean_vector,Wx,1,"1-")
        y_2 <-nw(xi_hat_x,mean_vector,Wx,1,"2-")
        E_h_1_x <- e_x*mu_x + lambda*(xi_hat_x *y_1 -y_2)
        E_h_2_x <- e_x*(2-e_x) + lambda^2*(xi_hat_x^2*y_p - 2*xi_hat_x*y_1 + y_2) 
        
        
        mu_hat <- c(mu_hat,mu_x)
        xi_hat <- c(xi_hat,xi_hat_x)
        y_p_hat <- c(y_p_hat,y_p)
        y_1_hat <- c(y_1_hat,y_1)
        E_h_1_hat <- c(E_h_1_hat,E_h_1_x)
        E_h_2_hat <- c(E_h_2_hat,E_h_2_x)
      }
      
      h_hat <- e_hat + lambda*(xi_hat-Y_ave)*(Y_ave<=xi_hat)
      
      
      
    } else{
      
      for (i in 1:nrow(X_ave)) {
        
        x <- X_ave[i,]
        
        Kx <- kernel(X_eta_1,x) 
        Wx <- Kx /sum(Kx)
        mean_vector <- mu_X[i] + residual_hat
        mu_x <- Wx %*% mean_vector
        
        mu_hat <- c(mu_hat, mu_x)
        
      }
      
      
    }

    
    #EIF
    phi_1 <-function(H,Y_P,E_H_2,E){  (2*(1-E)*(1-H)/Y_P+H^2- E_H_2)*Z_ave/E + 2*(1-E)*(Z_ave-E+(E-Z_ave)/Y_P) + E_H_2   }
    
    phi_2 <-function(MU,H,Y_P,Y_1,E_H_1,E){ ((1-H)/Y_P*Y_1+H*Y_ave - E_H_1)*Z_ave/E + (E-Z_ave)/Y_P*Y_1 + (Z_ave-E)*MU + E_H_1  }
    
    
    #one-step estimation
    
    if (lambda==0){
      #Use AIPW under no unconfounding
      phi_1_k <- rep(1,length(mu_hat))
      phi_2_k <- (Y_ave - mu_hat)*Z_ave/e_hat + mu_hat
      phi_2_k <- phi_2_k*sign_change
      
    } else{
      #Use EIF under L^2 model
      phi_1_k <-phi_1(h_hat,y_p_hat,E_h_2_hat,e_hat)
      phi_2_k <-phi_2(mu_hat,h_hat,y_p_hat,y_1_hat,E_h_1_hat,e_hat)*sign_change
      
    }
    

    phi_1_list <- c(phi_1_list,phi_1_k)
    phi_2_list <- c(phi_2_list,phi_2_k)
    
  }
  
  
  return(list(phi_1_list, phi_2_list))
  
  
}





kernel <-function(X,x,h_x=NULL){
  
  #Nadaraya Watson kernel estimator
  
  if (is.null(dim(X))){
    if (is.null(h_x)){
      rep(1,length(X))
    } else{
      dnorm(x,mean=X, sd= h_x)
    }
  } else{
    if (is.null(h_x)){
      rep(1,dim(X)[1])
    } else{
      apply(X,1, function(Xi) prod(dnorm(x,mean=Xi,sd=h_x)) )
    }
  }
}




nw <- function(y,Y,Wx,h,opt){
  
  #Kernel estimators of expectations
  
  if (opt == "d"){
    Ky <-  dnorm(y, mean=Y, sd=h)
  }  else if(opt == "p"){
    Ky <-  pnorm(y, mean=Y, sd=h)
  }  else if(opt =="1-"){
    Ky <- truncnorm_1_minus(y, mean=Y, sd=h)
  }  else if(opt =="1+"){
    Ky <- truncnorm_1_plus(y, mean=Y, sd=h)
  }  else if(opt == "2-"){
    Ky <- truncnorm_2_minus(y, mean=Y, sd=h)
  }  else if(opt == "2+"){
    Ky <- truncnorm_2_plus(y, mean=Y, sd=h)
  }  else{
    stop("invalid option.")
  }
  
  return(Wx %*% drop(Ky))
}



nw_q_xi <- function(Y,Wx,para,h,index,e=NULL){
  
  
  #Compute the quantile or threshold
  
  f_para <- function(xi){
    
    if (index=="gamma"){
      
      tt <- (nw(xi,Y,Wx,h,"p")-para)
      
      return(tt)
      
    } else if(index=="lambda"){
      y_p <- nw(xi,Y,Wx,h,"p")
      y_1 <- nw(xi,Y,Wx,h,"1-")
      tt <- (xi*y_p - y_1 - (1-e)/para)
      
      return(tt)
      
    } else{
      stop("wrong index")
    }
  }
  
  
  
  a <- -5000000
  b <- 5000000
  c <- uniroot(f_para,interval = c(a,b))$root 
  
  
  return(c)
  
}





#Truncated expectations for Gaussian distributions


truncnorm_1_minus <- function(y,mean,sd){
  beta <- (y-mean)/sd
  return(mean*pnorm(beta)-sd*dnorm(beta))
}   

truncnorm_1_plus <- function(y,mean,sd){
  beta <- (y-mean)/sd
  return(mean*(1-pnorm(beta))+sd*dnorm(beta))
}   

truncnorm_2_minus <- function(y,mean,sd){
  beta <- (y-mean)/sd
  sum <- (mean^2+sd^2)*pnorm(beta) - 2*mean*sd*dnorm(beta)-sd^2*beta*dnorm(beta)
  return(sum)
}  

truncnorm_2_plus <- function(y,mean,sd){
  beta <- (y-mean)/sd
  sum <- (mean^2+sd^2)*(1-pnorm(beta))+2*mean*sd*dnorm(beta)+sd^2*beta*dnorm(beta)
  return(sum)
}  


