

load("climBerth_f.RData")
load("climTemis_f.RData")
load("m4b2/param.RData")


m_berth <- list.files(path = 'm4b2/', pattern="model_berth" , full.names = T)

m_temis <- list.files(path = 'm4t2/', pattern="model_temis", full.names = T)


res_bert <- data.frame(ID=integer(),
                       numemo_model=integer(),
                       nb_cova=integer(),
                       #covariables=character(),
                       portemanteau=integer(),
                       pval_signif=logical(),
                       AIC=integer(),
                       #parametres=double(),
                       n_param_sign=integer(),
                       stringsAsFactors=FALSE)

res_temis <- data.frame(ID=integer(),
                        numemo_model=integer(),
                       nb_cova=integer(),
                       #covariables=character(),
                       portemanteau=integer(),
                       pval_signif=logical(),
                       AIC=integer(),
                       #parametres=double(),
                       n_param_sign=integer(),
                       stringsAsFactors=FALSE)

IDb <- 0

for(i in 1:length(m_berth)){
  load(m_berth[i])
  IDb <- IDb+1
  num_modb <- as.integer(sub(".*berth_", "", sub(".RData.*", "", m_berth[i])))
  n_varb <- length(param[[num_modb]])
  #cova <- colnames(climBerth_saison)[param[[num_mod]]]
  t_statb <- mod_berth$portemanteau[2]
  t_signb <- t_statb>0.05
  m_AICb <- mod_berth$AIC
  #m_param <- mod_berth$parameters
  n_param_sigb <- sum(abs(mod_berth$parameters/mod_berth$standard_deviation)>1.80)
  
  res_bert[i,] <- c(IDb, num_modb,n_varb, t_statb, t_signb, m_AICb, n_param_sigb)
  
  
  print(IDb)
}


IDt <- 0

for(j in 1:length(m_temis)){
  load(m_temis[j])
  IDt <- IDt+1
  num_modt <- as.integer(sub(".*temis_", "", sub(".RData.*", "", m_temis[j])))
  n_vart <- length(param[[num_modt]])
  #cova <- colnames(climTemis_saison)[param[[num_modt]]]
  t_statt <- mod_temis$portemanteau[2]
  t_signt <- t_statt>0.05
  m_AICt <- mod_temis$AIC
  #m_param <- mod_berth$parameters
  n_param_sigt <- sum(abs(mod_temis$parameters/mod_temis$standard_deviation)>1.80)
  
  res_temis[j,] <- c(IDt, num_modt, n_vart, t_statt, t_signt, m_AICt, n_param_sigt)
  
  
  print(IDt)
}

library(dplyr)

temis_bon <- res_temis %>% filter(pval_signif==TRUE)
berth_bon <- res_bert %>% filter(pval_signif==TRUE)


#meilleur modèle temis
load(m_temis[13244])
mod_temis$parameters/mod_temis$standard_deviation
colnames(climTemis_saison)[param[[718]]]


#meilleur modèle berth
load(m_berth[10855])
mod_berth$parameters/mod_berth$standard_deviation
colnames(climBerth_saison)[param[[5028]]]



temis_b <- temis_bon %>% filter(n_param_sign!=0)
berth_b <- berth_bon %>% filter(n_param_sign!=0)

save(temis_b, file = "bon_models_temis.RData")
save(berth_b, file = "bon_model_berth.RData")
