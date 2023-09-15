likehood <- function(y, x, theta){
   K = ncol(y)-1
   n = nrow(y)
   P = ncol(x)
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   likehood_value <- 0
   eta <- 0.5 + numeric(length = K)
   for(needle_t in 2:n){
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      cum_eta <- cumsum(eta)
      likehood_value <- likehood_value  - sum(y[needle_t, 2:(K+1)]*cum_eta)+log(1+sum(exp(cum_eta)))
   }
   likehood_value 
}


score <- function(y, x, theta){
   K = ncol(y)-1
   n = nrow(y)
   P = ncol(x)
   number_of_parameters <- length(theta)
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   #partial_eta_initial <-  numeric(length = number_of_parameters)
   matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
   eta <- 0.5 + numeric(length = K)
   score_value <- numeric(number_of_parameters)
   for (needle_t in 2:n) {
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      for (needle_k in 1:K) {
         matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k]) 
         matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
         matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
         matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
      }
      inverse_cum_y <- rev(cumsum(rev(y[needle_t, 2:(K+1)])))
      cumulate_eta <- cumsum(eta)
      err <- inverse_cum_y - rev(cumsum(rev(cumulate_eta)))
      score_value <- score_value + apply(-err * t(matrix_partial_eta), 2, sum)
   }
 score_value
}


x_path <- function(n,P) matrix(data = rnorm(n*P), nrow = n, ncol = P)


x_path2 <- function(n,P) matrix(data = rcauchy(n*P), nrow = n, ncol = P)

paths <- function(n,P,n_cat,theta, eta_0 = 0.5){
   K <- n_cat-1
   x <- x_path(n,P)
   y <- matrix(data = 0, ncol = n_cat, nrow = n)
   y[1,1] <- 1
   ## vérifier la taille de theta
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   
   eta <- eta_0 + numeric(length = K)
   
   for(needle_t in 2:n){
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      for (needle_k in 1:K) {
         numerator <- exp(cumsum(eta))
         pi <- numerator  /(1+sum( numerator) )
      }
      y[needle_t, ] <- rmultinom(1, 1, c(max(0,1-sum(pi)), pi))[,1]
   }
   return(list(y = y, x = x))
}


paths2 <- function(n,P,n_cat,theta, eta_0 = 0.5){
   K <- n_cat-1
   x <- x_path(n,P)
   y <- matrix(data = 0, ncol = n_cat, nrow = n)
   y[1,1] <- 1

   ## vérifier la taille de theta
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   
   eta <- eta_0 + numeric(length = K)
   
   for(needle_t in 2:n){
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + sum(beta * eta)

for (needle_k in 1:K) {
   numerator <- exp(cumsum(eta))
   pi <- numerator  /(1+sum( numerator) )
}
 y[needle_t, ] <- rmultinom(1, 1, c(max(0,1-sum(pi)), pi))[,1]
   }
   return(list(y = y, x = x))
}

bipaths <- function(n,P,n_cat,theta, eta_0 = 0.5){
   K <- n_cat-1
   x1 <- x_path(n,P)
   x2 <-  x_path(n,P)
   y1 <- matrix(data = 0, ncol = n_cat, nrow = n)
   y1[1,1] <- 1
   y2 <- matrix(data = 0, ncol = n_cat, nrow = n)
   y2[1,1] <- 1
   ## vérifier la taille de theta
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   
   eta1 <- eta_0 + numeric(length = K)
   eta2 <- eta_0 + numeric(length = K)
   
   for(needle_t in 2:n){
      eta1  <- omega  + sum(gamma * x1[needle_t-1, ]) + sum(alpha * y1[needle_t-1, 2:(K+1)]) + beta * eta1
      eta2  <- omega  + sum(gamma * x2[needle_t-1, ]) + sum(alpha * y2[needle_t-1, 2:(K+1)]) + beta * eta2
      
      for (needle_k in 1:K) {
         numerator1 <- exp(cumsum(eta1))
         pi1 <- numerator1  /(1+sum( numerator1) )
         numerator2 <- exp(cumsum(eta2))
         pi2 <- numerator2  /(1+sum( numerator2) )
      }
      u <- runif(1)
      y1[needle_t, ] <- (which.max(u-cumsum(c(1-sum(pi1), pi1))<0)==1:n_cat) * 1
      y2[needle_t, ] <- (which.max((1-u)-cumsum(c(1-sum(pi2), pi2))<0)==1:n_cat) * 1
   }
   return(list(y1 = y1, x1 = x1, y2 = y2, x2 = x2))
}

bipaths_concorde <- function(n,P,n_cat,theta, theta_other,  eta_0 = 0.5){
   K <- n_cat-1
   x1 <- x_path(n,P)
   x2 <-  x_path(n,P)
   y1 <- matrix(data = 0, ncol = n_cat, nrow = n)
   y1[1,1] <- 1
   y2 <- matrix(data = 0, ncol = n_cat, nrow = n)
   y2[1,1] <- 1
   ## vérifier la taille de theta
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   
   omega_other <- theta_other[1:K]
   gamma_other <- theta_other[(K+1):(K+P)]
   alpha_other <- theta_other[(K+P+1):(2*K + P)]
   beta_other <- theta_other[(2*K + P+1):(3*K + P)]
   
   eta1 <- eta_0 + numeric(length = K)
   eta2 <- eta_0 + numeric(length = K)
   
   for(needle_t in 2:n){
      eta1  <- omega  + sum(gamma * x1[needle_t-1, ]) + sum(alpha * y1[needle_t-1, 2:(K+1)]) + beta * eta1
      eta2  <- omega_other  + sum(gamma_other * x2[needle_t-1, ]) + sum(alpha_other * y2[needle_t-1, 2:(K+1)]) + beta_other * eta2
      
      for (needle_k in 1:K) {
         numerator1 <- exp(cumsum(eta1))
         pi1 <- numerator1  /(1+sum( numerator1) )
         numerator2 <- exp(cumsum(eta2))
         pi2 <- numerator2  /(1+sum( numerator2) )
      }
      u <- runif(1)
      y1[needle_t, ] <- (which.max(u-cumsum(c(1-sum(pi1), pi1))<0)==1:n_cat) * 1
      y2[needle_t, ] <- (which.max(u-cumsum(c(1-sum(pi2), pi2))<0)==1:n_cat) * 1
   }
   return(list(y1 = y1, x1 = x1, y2 = y2, x2 = x2))
}



sequence_score <- information_matrix <- function(y, x, theta){
   K = ncol(y)-1
   n = nrow(y)
   P = ncol(x)
   number_of_parameters <- length(theta)
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   score <- matrix(0, nrow = n-10, ncol = number_of_parameters)
   matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
   eta <- 0.5 + numeric(length = K)
   for (needle_t in 2:n) {
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      for (needle_k in 1:K) {
         matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k]) 
         matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
         matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
         matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
      }
      inverse_cum_y <- rev(cumsum(rev(y[needle_t, 2:(K+1)])))
      cumulate_eta <- cumsum(eta)
      err <- inverse_cum_y - rev(cumsum(rev(exp(cumulate_eta))))/(1+sum(exp(cumulate_eta)))
      score_t <- apply(-err * t(matrix_partial_eta), 2, sum)
      if (needle_t > 10){
         score[needle_t-10, ] <-  t(score_t)
      }
   }
   score
}


information_matrix <- function(y, x, theta){
      K = ncol(y)-1
      n = nrow(y)
      P = ncol(x)
      number_of_parameters <- length(theta)
      omega <- theta[1:K]
      gamma <- theta[(K+1):(K+P)]
      alpha <- theta[(K+P+1):(2*K + P)]
      beta <- theta[(2*K + P+1):(3*K + P)]
      #partial_eta_initial <-  numeric(length = number_of_parameters)
      information <- matrix(0, ncol = number_of_parameters, nrow = number_of_parameters)
      matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
      eta <- 0.5 + numeric(length = K)
      for (needle_t in 2:n) {
         eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
         for (needle_k in 1:K) {
            matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k]) 
            matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
            matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
            matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
         }
         inverse_cum_y <- rev(cumsum(rev(y[needle_t, 2:(K+1)])))
         cumulate_eta <- cumsum(eta)
         err <- inverse_cum_y - rev(cumsum(rev(exp(cumulate_eta))))/(1+sum(exp(cumulate_eta)))
         score_t <- apply(-err * t(matrix_partial_eta), 2, sum)
         if (needle_t > 10){
            information <- information + score_t %*% t(score_t)
         }
      }
      information/(n-10)
   }


sensibility_matrix <- function(y, x, theta){
   K = ncol(y)-1
   n = nrow(y)
   P = ncol(x)
   number_of_parameters <- length(theta)
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
   #partial_eta_initial <-  numeric(length = number_of_parameters)
   matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
   sensibility <- matrix(0, ncol = number_of_parameters, nrow = number_of_parameters)
   eta <- 0.5 + numeric(length = K)
   for (needle_t in 2:n) {
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      for (needle_k in 1:K) {
         matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k]) 
         matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
         matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
         matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
      }
      cumulate_eta <- cumsum(eta)
      
      sensibility_t <- matrix(0, ncol = number_of_parameters, nrow = number_of_parameters)
      

      if (needle_t > 10){
      
      #### calcul de h_un
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 1:K){
            right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
         }
      h_un <-  (1/((1+sum(exp(cumulate_eta)))^2))* right_none  %*% t(matrix_partial_eta[,1])
      
     
      #### calcul de h_deux
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 2:K){
         right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
      }
      h_deux <- (1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[2:K]) * matrix_partial_eta[,1]  + (1 + exp(cumulate_eta)[1]) * right_none) %*% t(matrix_partial_eta[,2])
      
      sensibility_t <- sensibility_t + h_un + h_deux
      if(K>2){
      for (needle_k in 3:K) {
         left_none <- matrix_partial_eta[,1] 
         for(needle_ell in 2:(needle_k-1)){
            left_none <- left_none + (1+sum(exp(cumulate_eta)[1:(needle_ell-1)])) * matrix_partial_eta[,needle_ell]
         }
         right_none <- numeric(number_of_parameters)
         for(needle_ell in needle_k:K){
         right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
         }
         sensibility_t <- sensibility_t + (1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[needle_k:K]) * left_none + (1 + sum(exp(cumulate_eta)[1:(needle_k -1)])) * right_none) %*% t(matrix_partial_eta[,needle_k])
      }
      }
      }
      sensibility <- sensibility + sensibility_t
      
   }
   sensibility/(n-10)
}



portemanteau_test <- function(y, x, theta){
   K = ncol(y)-1
   n = nrow(y)
   P = ncol(x)
   number_of_parameters <- length(theta)
   omega <- theta[1:K]
   gamma <- theta[(K+1):(K+P)]
   alpha <- theta[(K+P+1):(2*K + P)]
   beta <- theta[(2*K + P+1):(3*K + P)]
  
   matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
   eta <- 0.5 + numeric(length = K)
   Sensibility <-  sensibility_matrix(y, x,  theta = theta)
   Information <- information_matrix(y, x,  theta = theta)
   c_10 <- numeric(number_of_parameters)
   c_20 <- numeric(number_of_parameters)
   c_3plus0 <- matrix(0,nrow = number_of_parameters, ncol = K-2)
   overline_mu2 <- matrix(0, nrow = K, ncol = K)
   G <- matrix(0,ncol = K, nrow  = number_of_parameters)
   rho <- numeric(K)
   err_moins1 <- rev(cumsum(rev(exp(cumsum(rep(1/2, K))))))/(1+sum(exp(cumsum(rep(1/2, K)))))
   for (needle_t in 2:n) {
      eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
      for (needle_k in 1:K) {
         matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k]) 
         matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
         matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
         matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
      }
      
      inverse_cum_y <- rev(cumsum(rev(y[needle_t, 2:(K+1)])))
      cumulate_eta <- cumsum(eta)
      err <- inverse_cum_y - rev(cumsum(rev(exp(cumulate_eta))))/(1+sum(exp(cumulate_eta)))
      
      if(needle_t>10){
         
         for(needle_k in 1:K){
           G[,needle_k] <- G[,needle_k] + err_moins1[needle_k] * err[needle_k] * apply(err* t(matrix_partial_eta) %*% t(solve(Sensibility)), 2, sum)
         }
         
         #### calcul de h_un
         right_none <- numeric(number_of_parameters)
         for(needle_ell in 1:K){
            right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
         }
         c_10 <- c_10 -  err[1]*(1/((1+sum(exp(cumulate_eta)))^2))* right_none   
         
         
         #### calcul de h_deux
         right_none <- numeric(number_of_parameters)
         for(needle_ell in 2:K){
            right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
         }
         c_20 <- c_20 -  err[2]*(1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[2:K]) * matrix_partial_eta[,1]  + (1 + exp(cumulate_eta)[1]) * right_none) 
         
         if(K>2){
            for (needle_k in 3:K) {
               left_none <- matrix_partial_eta[,1] 
               for(needle_ell in 2:(needle_k-1)){
                  left_none <- left_none + (1+sum(exp(cumulate_eta)[1:(needle_ell-1)])) * matrix_partial_eta[,needle_ell]
               }
               right_none <- numeric(number_of_parameters)
               for(needle_ell in needle_k:K){
                  right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
               }
               c_3plus0[, needle_k-2] <- c_3plus0[, needle_k-2] - err[needle_k]*(1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[needle_k:K]) * left_none + (1 + sum(exp(cumulate_eta)[1:(needle_k -1)])) * right_none) 
            }
         }
         
         rho_cur <- err_moins1 * err
         rho <- rho + rho_cur
         overline_mu2 <- overline_mu2 + rho_cur %*% t(rho_cur)
      }
      
      err_moins1 <- err  
   }
   matrix_C <- cbind(c_10, c_20, c_3plus0)/(n-10)
   G <- G/(n-10)
   rho <- rho/(n-10)
   overline_mu2 <- overline_mu2/(n-10)

   matrix_W <- overline_mu2 + t(matrix_C) %*% solve(Sensibility)%*%Information%*%t(solve(Sensibility)) %*% matrix_C + t(matrix_C) %*% G + t(G) %*% matrix_C
   statistic <- (n-10) * rho %*% solve(matrix_W) %*% rho
   p_value <- pchisq(statistic, df = K, lower.tail = FALSE)
   list(statistic = statistic, p_value = p_value)
}



poly_autoregression <- function(y, x){
   likehood_data <- function(theta) likehood(y=y, x = x, theta = theta)
   repeat_opt <- 10
   K = ncol(y)-1
   P = ncol(x)
   
   list_values <- numeric(repeat_opt)*NA
   mat_param <- matrix(NA, ncol = 3*K+P , nrow = repeat_opt)
   
   for(needle_rep in 1:repeat_opt){
      print(needle_rep)
      parauto <- rexp(2*K)
      theta_init <- c(rnorm(K+P), parauto/(2*sum(parauto)))
      
      res_opt <- tryCatch({
         optim(theta_init, fn = likehood_data, control = list(maxit = 1e6), method = "L-BFGS-B", 
              lower = c(rep(-Inf, 2*K+P),  rep(-0.999999, K)),
              upper = c(rep(Inf, 2*K+P),  rep(0.999999, K)))
      }, error = function(e) {
         return(NA)
      })
      
      if(length(res_opt)>1){
      list_values[needle_rep] <- res_opt$value
      mat_param[needle_rep, ] <- res_opt$par
      }
   }
   parameters <- mat_param[which.min(list_values),]
   loss_value <- min(list_values, na.rm = TRUE)
   Information <-  information_matrix(y, x, parameters)
   Sensibility <-  sensibility_matrix(y, x, parameters)
   parameters_covariance <- solve(Sensibility)%*%Information%*%t(solve(Sensibility))/(nrow(y)-1)
   standard_deviation <- sqrt(diag(parameters_covariance))
   AIC <- 2*loss_value + 2 * (3 * K + P)
   portemanteau = portemanteau_test(y = y, x = x, theta = parameters)
   list(parameters = parameters, standard_deviation = standard_deviation, AIC =  AIC,
        loss_value =  loss_value, Information = Information, Sensibility = Sensibility,
        parameters_covariance = parameters_covariance, nsample = nrow(y), 
        number_of_parameters = (3 * K + P), portemanteau = portemanteau, y=y, x = x) #
}

comparison_test <- function(model_poly_autoregressive1, model_poly_autoregressive2){
 nsample <- model_poly_autoregressive1$nsample  
 hat_theta_2 <- model_poly_autoregressive2$parameters
 matrix_v <- model_poly_autoregressive1$nsample * model_poly_autoregressive1$parameters_covariance +  model_poly_autoregressive2$nsample * model_poly_autoregressive2$parameters_covariance  
 global_statistic <- (nsample-1) * (hat_theta_1 - hat_theta_2) %*% solve(matrix_v) %*% (hat_theta_1 - hat_theta_2)
 p_value <- pchisq(global_statistic, df = length(hat_theta_2), lower.tail = FALSE)
 list(global_statistic = global_statistic, p_value = p_value)
}

comparison_test_dependent <- function(model_poly_autoregressive1, model_poly_autoregressive2){
   nsample <- model_poly_autoregressive1$nsample 
   hat_theta_1 <- model_poly_autoregressive1$parameters
   hat_theta_2 <- model_poly_autoregressive2$parameters
   inv_sensibility_1 <- solve(model_poly_autoregressive1$Sensibility)
   inv_sensibility_2 <- solve(model_poly_autoregressive2$Sensibility)
   score_1 <- sequence_score(model_poly_autoregressive1$y, model_poly_autoregressive1$x, model_poly_autoregressive1$parameters)
   score_2 <- sequence_score(model_poly_autoregressive2$y, model_poly_autoregressive2$x, model_poly_autoregressive2$parameters)
   S <- matrix(0, ncol = length(hat_theta_2), nrow = length(hat_theta_2))
   for (needle_t in 1:nrow(score_1)) {
      S <- S + score_1[needle_t,] %*% t(score_2[needle_t,])
   }
   S <- S/nrow(score_1)
   matrix_v <- (nsample-1)* model_poly_autoregressive1$parameters_covariance +  (nsample-1)* model_poly_autoregressive2$parameters_covariance  +  inv_sensibility_1 %*% S %*% t(inv_sensibility_2) + inv_sensibility_2 %*% t(S) %*% t(inv_sensibility_1)
   global_statistic <- (nsample-1) * (hat_theta_1 - hat_theta_2) %*% solve(matrix_v) %*% (hat_theta_1 - hat_theta_2)
   p_value <- pchisq(global_statistic, df = length(hat_theta_2), lower.tail = FALSE)
   list(global_statistic = global_statistic, p_value = p_value, test_variable =abs(sqrt(nsample-1) * (hat_theta_1 - hat_theta_2) / sqrt(diag(matrix_v))))
}



