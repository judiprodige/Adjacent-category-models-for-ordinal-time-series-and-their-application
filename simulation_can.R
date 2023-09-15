source("model_v5.R")
n_length <- c(70, 100, 300, 500, 1000)
P_value <-  5
n_cat_val <-  4

parameters_list <- vector(mode = "list", length = length(n_length))
sd_list <- vector(mode = "list", length = length(n_length))

theta_value <- c(c(1.2, 0.7, 0.5), c(-0.8, 1.5, -1.5, 2, 2), c(0.3, -0.3, 0.5), c(0.8, -0.2, 0.3))
times_simulations <- 499

for (needle_paths_length in seq_along(n_length)) {
   matrix_parameters_simulation <- matrix(nrow = times_simulations, ncol = length(theta_value))
   matrix_sd_simulation <- matrix(nrow = times_simulations, ncol = length(theta_value))
   for(needle_simulation in 1:times_simulations){
      if(needle_simulation %% 50 == 1) {print(paste(needle_paths_length, needle_simulation))}
      data_simulated <- paths(n = n_length[needle_paths_length], P = P_value, n_cat = n_cat_val, theta = theta_value)
      y_sim <- data_simulated$y
      x_sim <- data_simulated$x
      matrix_parameters_simulation[needle_simulation,] <- poly_autoregression(y = y_sim, x = x_sim)$parameters
      I <-  information_matrix(y_sim, x_sim, matrix_parameters_simulation[needle_simulation,])
      J <-  sensibility_matrix(y_sim, x_sim, matrix_parameters_simulation[needle_simulation,])
      matrix_sd_simulation[needle_simulation,] <- sqrt(diag(solve(J)%*%I%*%t(solve(J))/n_length[needle_paths_length]))
   }
   parameters_list[[needle_paths_length]] <- matrix_parameters_simulation
   save(parameters_list, file = "simulation_parameters.R")
   save(sd_list, file = "simulation_sd.R")
}



