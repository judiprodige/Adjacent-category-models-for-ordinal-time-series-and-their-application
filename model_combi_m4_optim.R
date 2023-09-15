source("model_v5.R")

library(doParallel)

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-2) 
registerDoParallel(cluster)


load("climBerth_f.RData")
load("climTemis_f.RData")
load("temis_cat20_95.RData")
load("berth_cat20_95.RData")


#numero des colonne contenant les covariables
v1 <- c(13:20, 29:32)
v2 <- c(21:32)
v3 <- c(33:44)
v4 <- c(41:52)

#nombre de covariables
n_cova1 <- length(v1)
n_cova2 <- length(v2)
n_cova3 <- length(v3)
n_cova4 <- length(v4)

#liste de toutes les combinaisons possibles de covariables
list_comb1 <- lapply(1:length(v1),  function(i) combn(v1, i))
list_comb2 <- lapply(1:length(v2),  function(i) combn(v2, i))
list_comb3 <- lapply(1:length(v3),  function(i) combn(v3, i))
list_comb4 <- lapply(1:length(v4),  function(i) combn(v4, i))

combi_list <- list(list_comb1, list_comb2, list_comb3, list_comb4)



#liste de parametres pour les modeles
#la longueur = nombre de combinaisons possibles 

param <- list()

#list_model_temis <- list()
#list_model_berth <- list()

tot_comb <- 2^n_cova1-1

for( k in 1:length(combi_list)){
  
  compt <- tot_comb*(k-1)
  compt2 <- k*tot_comb
  
  param[[compt2]] <- combi_list[[k]][[n_cova1]][,1]

  
  for(i in 1:(n_cova1-1)){
    for(j in 1:length(list_comb1[[i]][1,])){
      
      compt <- compt +1
      
      param[[compt]] <- combi_list[[k]][[i]][,j]

      }
  }
  print(compt2)
}

save(param, file = paste("m4b2/param",".RData",sep = ""))


list_model_berth <- vector("list",length =length(param))

list_model_temis <- vector("list",length =length(param))

clusterExport(cluster, c("list_model_berth", "list_model_temis", "param", "berth_cat20_95", "temis_cat20_95", "climBerth_saison", "climTemis_saison"))


foreach(t=1:length(param)) %dopar% {
  
  an <- 4:75
  list_model_berth[[t]] <- assign(paste("model_berth", t,sep = "_"), poly_autoregression(berth_cat20_95[an,],as.matrix(climBerth_saison[an,param[[t]]])))
  
  list_model_temis[[t]] <- assign(paste("model_Temis", t,sep = "_"), poly_autoregression(temis_cat20_95[an,],as.matrix(climTemis_saison[an,param[[t]]])))
  
  
  mod_berth <-list_model_berth[[t]]
  mod_temis <-list_model_temis[[t]]
    
  save(mod_berth, file = paste("m4b2/model_berth_",t,".RData",sep = ""))
  save(mod_temis, file = paste("m4t2/model_temis_",t,".RData",sep = ""))
  
}



#Stop cluster
stopCluster(cluster)

