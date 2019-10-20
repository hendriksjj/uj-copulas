
# results_dir <- 'C:\\Users\\F5056756\\Google Drive\\Jurgens\\MCOM Financial Economics\\0. Research Methodology\\1. Paper\\results'
# setwd('C:\\Users\\F5056756\\Google Drive\\Jurgens\\MCOM Financial Economics\\0. Research Methodology\\1. Paper\\clean-data\\0. FULL_DATA')

library(dplyr)
library(MTS)
library(psych)
library(tseries)
library(rugarch)
library(VineCopula)

##### FUNCTIONS #####
# Helper Functions for data
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}
country_sector <- function(x) {
  s <- strsplit(x, " ")[[1]]
  len <- length(s)
  Country <- paste(s[1:(len - 1)], collapse = ' ')
  Sector <- s[len]
  if(Sector == 'Financials'){
    Sector <- 'Financial'
  }
  if(Sector == 'Industrials'){
    Sector <- 'Industrial'
  }
  if(Sector == 'Resources'){
    Sector <- 'Resource'
  }
  data.frame(Country_Sector = x, Country, Sector)
}

# Writes Descriptive Statistics.csv
descr_stats <- function(df_indices){
  library(psych)
  library(tseries)
  colnames(df_indices) <- colnames(df_indices) %>%
    stringr::str_replace_all('_', ' ') %>%
    lapply(FUN = simpleCap) %>%
    unlist()
  
  descr <- psych::describe(df_indices) %>%
    as.data.frame() %>%
    select(n, mean, sd, skew, kurtosis)
  
  descr$Country_Sector <- row.names(descr)
  
  jb_p_val <- apply(X = df_indices, 2 , FUN = function(x) jarque.bera.test(x)$p.value) %>%
    data.frame()
  colnames(jb_p_val) <- 'Jarque-Bera p-value'
  jb_p_val$Country_Sector <- row.names(jb_p_val)
  
  jb_test_stat <- apply(X = df_indices, 2 , FUN = function(x) jarque.bera.test(x)$statistic) %>%
    data.frame()
  colnames(jb_test_stat) <- 'Jarque-Bera test statistic'
  jb_test_stat$Country_Sector <- row.names(jb_test_stat)
  
  descr <- descr %>%
    inner_join(jb_test_stat) %>%
    inner_join(jb_p_val)
  
  descr <- descr %>%
    inner_join(descr$Country_Sector %>%
                 lapply(country_sector) %>%
                 rbind_list())
  
  df_ret <- descr %>%
    select(Country, Sector, Mean = mean, `Standard Deviation` = sd, Skewness = skew, Kurtosis = kurtosis, `Jarque-Bera test statistic`) %>% 
    mutate_if(is.numeric, round, 4)  
  df_ret$Country[duplicated(df_ret$Country)] = ''
  write.csv(df_ret, 'finalResults//Descriptive Statistics.csv')
  return(TRUE)
}

# Writes correlations.csv
corr <- function(df_full){
  df_indices <- df_full[,3:ncol(df_full)]
  cc <- colnames(df_indices) %>%
    stringr::str_replace_all('_', ' ') %>%
    lapply(FUN = simpleCap) %>%
    unlist() %>%
    lapply(country_sector) %>%
    rbind_list() %>%
    mutate(' ' = paste(Country, Sector)) %>%
    select(' ') %>%
    as.vector()
  
  colnames(df_indices) <- cc$` `
  write.csv(round(cor(df_indices), 4), 'finalResults//correlations.csv')
  return(TRUE)
}

# Provides best ARMA(p_max, q_max)-GARCH(r_max, s_max) model specification for index y
fun_garch_spec <- function(y, p_max, q_max, r_max, s_max){
  library(rugarch)
  select_train <- 1:floor(1 * length(y))
  select_test <- setdiff(1:length(y), select_train)
  y_train <- y[select_train]
  y_test <- y[select_test]
  fit_save <- list()
  garch.solvers <- c('nlminb', 'solnp', 'lbfgs', 'gosolnp', 'nloptr', 'hybrid')
  
  
  
  
  count = 1
  for(p in 0:p_max){
    for(q in 0:q_max){
      for(r in 1:r_max){
        for(s in 1:s_max){
          if(count == 1){
            model_spec <- data.frame(ModelNum = count, ModelType = 'GARCH', AR = p,MA = q, G = r, ARCH = s, Model = "sGARCH", SubModel = 0, Dist = "std", MSE_in_Sample = 0, MSE_out_Sample = 0, solver = '', stringsAsFactors = FALSE, bds_ts = 0, bds_p = 0)
          }
          
          model_spec <- rbind(model_spec, list(count,'GARCH',p,q,r,s,"sGARCH", 0, "std", 0, 0, '', 0, 0))
          count = count + 1
        }
      }
    }
  }
  
  # sample(1:nrow(model_spec), nrow(model_spec))
  for(i in 1:nrow(model_spec)){
    
    print('----------------------------------------------------------')
    print(paste('Model ',i))
    
    if(model_spec[i,'ModelType']=='GARCH'){
      spec=ugarchspec(mean.model = list(armaOrder = c(model_spec[i,'AR'], model_spec[i,'MA']), include.mean = TRUE, archm = FALSE, external.regressors = NULL), variance.model = list(model = model_spec[i,'Model'], submodel = if(model_spec[i,'SubModel']==0) NULL else model_spec[i,'SubModel'], garchOrder = c(model_spec[i,'G'], model_spec[i,'ARCH'])),distribution.model = model_spec[i,'Dist'])
      
      
      fit <- NULL
      forc <- NULL
      attempt <- 1
      while((is.null(fit) || is.null(forc)) && attempt <= length(garch.solvers)){
        try(fit <- ugarchfit(spec, data=y, out.sample=floor(0.75 * length(y)), solver = garch.solvers[attempt]),silent = TRUE)
        try(forc <- ugarchforecast(fit, n.ahead=(length(y) - floor(0.75 * length(y)))), silent = TRUE)
        
        attempt <- attempt + 1
      }
      
      model_spec[i, 'solver'] <- garch.solvers[attempt - 1]
      
      
      
      if(!is.null(forc)){
        forc_resid <- y_test - forc@forecast$seriesFor
        data_resid <- fit@fit$residuals
        if(i==1){
          resid_save = data.frame(forc_resid)
          rownames(resid_save) <- NULL
          colnames(resid_save) <- paste(i)
          
          actual_resid_save = data.frame(data_resid)
          rownames(actual_resid_save) <- NULL
          colnames(actual_resid_save) <- paste(i)
        } else {
          resid_save_temp <- data.frame(data.frame(forc_resid))
          rownames(resid_save_temp) <- NULL
          colnames(resid_save_temp) <- paste(i)
          resid_save <- cbind.data.frame(resid_save, resid_save_temp)
          
          actual_resid_save_temp <- data.frame(data.frame(data_resid))
          rownames(actual_resid_save_temp) <- NULL
          colnames(actual_resid_save_temp) <- paste(i)
          actual_resid_save <- cbind.data.frame(actual_resid_save, actual_resid_save_temp)
        }
        model_spec[i, 'solver'] <- garch.solvers[attempt - 1]
        model_spec[i,'MSE_in_Sample'] <- sum(fit@fit$residuals * fit@fit$residuals)/length(fit@fit$residuals)
        model_spec[i,'MSE_out_Sample'] <- sum(forc_resid * forc_resid)/length(forc_resid)
        sigma_square_std <- log(fit@fit$residuals ** 2)
        
        # copdata <- pobs(sigma_square_std)
        
        bds <- bds.test(sigma_square_std)
        model_spec[i,'bds_ts']  <- bds$statistic[1]
        model_spec[i,'bds_p']  <- bds$p.value[1]
        if(bds$p.value[1] > 0.05){
          i <- nrow(model_spec)
        }
      } else {
        model_spec[i,'MSE_in_Sample'] <- NULL
        model_spec[i,'MSE_out_Sample'] <- NULL
        model_spec[i,'bds_ts']  <- NULL
        model_spec[i,'bds_p']  <- NULL
      }
      
      fit_save[i] <- fit
      
    }
  }
  
  # return(model_spec[order(model_spec$MSE_out_Sample),][1,])  
  if(max(model_spec$bds_p) > 0.05){
    return(model_spec[order(model_spec$bds_p, decreasing=TRUE),])    
  } else {
    return(model_spec[order(model_spec$bds_ts, decreasing=FALSE),])    
  }
  
}

# Writes marginal_models.csv
model_output <- function(spec_save){
  country_sec <- spec_save$index %>%
    stringr::str_replace_all('_', ' ') %>%
    lapply(FUN = simpleCap) %>%
    unlist() %>%
    lapply(country_sector) %>%
    rbind_list()
  
  country_sec$index <- spec_save$index
  
  df_return <- spec_save %>%
    select(index, p = AR, q = MA, r = G, s = ARCH, `Out of sample MSE` = MSE_out_Sample, `BDS p-value`=bds_p) %>%
    mutate_if(is.numeric, round, 4) %>%
    inner_join(country_sec) %>%
    select(Country, Sector, p, q, r, s, `BDS p-value`)
  
  df_return$Country[duplicated(df_return$Country)] <- ''
  
  write.csv(df_return, 'finalResults//marginal_models.csv')
  
  return(TRUE)
}

##### CODE #####

