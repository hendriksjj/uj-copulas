options(scipen = 99)

library(VineCopula)
library(psych)
library(dplyr)
library(kdecopula)
library(tseries)

# Helper functions
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

copula_results <- function(copdata, rvm){
  copula_families <- rvm$family
  copula_par_1 <- rvm$par
  copula_par_2 <- rvm$par2
  copula_names <- rvm$names
  count = 0
  for(i in 1:(ncol(copula_families) - 1)){
    for(j in (i + 1):(nrow(copula_families))){
      left_index_name <- copula_names[i] %>%
        stringr::str_replace_all('_', ' ') %>%
        lapply(FUN = simpleCap) %>%
        unlist()
      right_index_name <- copula_names[j] %>%
        stringr::str_replace_all('_', ' ') %>%
        lapply(FUN = simpleCap) %>%
        unlist()
      
      copula_family <- df_copula_class[df_copula_class$Code == copula_families[j, i],]$Copula %>%
        as.character() %>%
        stringr::str_squish() %>%
        simpleCap()
      
      if(count == 0){
        df_ret <- data.frame(Indices = paste(left_index_name, right_index_name, sep = ' and '), 
                             copula_family,
                             Param1 = copula_par_1[j, i],
                             Param2 = copula_par_2[j, i])      
      } else {
        df_ret <- rbind(df_ret,
                        data.frame(Indices = paste(left_index_name, right_index_name, sep = ' and '), 
                                   copula_family,
                                   Param1 = copula_par_1[j, i],
                                   Param2 = copula_par_2[j, i])
        )
      }
      count = count + 1
      
    }
  }
  
  write.csv(df_ret, 'finalResults//Copula Results.csv') 
  return(TRUE)
}

# Upper and Lower dependency calculations
U_TDC <- function(xy, upper_perc){
  N <- nrow(xy)
  i_u <- N * upper_perc 
  cop_est <- sum((xy[,1] <= i_u/N) & (xy[,2] <= i_u/N))/N
  ret_val <- (1 - 2 * (i_u/N) + cop_est)/(1 - i_u/N)
  return(ret_val)
}
L_TDC <- function(xy, lower_perc){
  N <- nrow(xy)
  i_l <- N * lower_perc 
  cop_est <- sum((xy[,1] <= i_l/N) & (xy[,2] <= i_l/N))/N
  ret_val <- cop_est/(1 - i_l/N)
  return(ret_val)
}

# Simulations
simu_TDC <- function(rvm, lower_perc, upper_perc, sim_size, number_of_sims){
  
  indices <- rvm$names
  combos <- list()
  for(i in 1:(length(indices) - 1)){
    for(j in (i+1):length(indices)){
      print('--')
      print(i)
      print(j)
      name <- paste(indices[i],'and',indices[j])
      combos[[name]] <- c(indices[i], indices[j])
      print('--')
    }
  }
  
  
  # Start of main simulation
  # Create df to store values
  df_return <- data.frame(SimNu = integer(),
                          combo_name = character(),
                          left_country = character(),
                          right_country = character(),
                          U_TDC = numeric(),
                          L_TDC = numeric(),
                          stringsAsFactors = FALSE)
  
  for(i in 1:number_of_sims){
    print('-----')
    print(paste('Simulation number:',i))
    print('-----')
    copdata_sim <- RVineSim(N = sim_size, RVM=rvm, U = NULL)  
    for(combo in combos){
      
      df_append <- data.frame(SimNu = i,
                              combo_name = paste(combo, collapse = ' and '),
                              left_country = combo[1],
                              right_country = combo[2],
                              U_TDC = U_TDC(copdata_sim[,combo], upper_perc),
                              L_TDC = L_TDC(copdata_sim[,combo], lower_perc))
      
      df_return <- rbind(df_return,
                         df_append)
    }
  }
  
  
  mean_val <- aggregate(df_return[,c('U_TDC','L_TDC')], list(df_return$combo_name), mean)
  
  
  percs <- c((1 - 0.99)/2, (1 - 0.95)/2,(1 - 0.90)/2,
             1 - (1 - 0.90)/2, 1 - (1 - 0.95)/2, 1 - (1 - 0.99)/2)
  
  conf_intervals <- aggregate(df_return[,c('U_TDC','L_TDC')], list(df_return$combo_name), FUN = quantile, probs  = percs)
  
  
  conf_levels <- list(c('5%', '95%'),
                      c('2.5%', '97.5%'),
                      c('0.5%', '99.5%'))
  
  
  # U_TDC
  U_TDC_insig_10 <- (0 >= conf_intervals$U_TDC[,conf_levels[[1]]][,1]) & (0 <= conf_intervals$U_TDC[,conf_levels[[1]]][,2])
  U_TDC_insig_5 <- (0 >= conf_intervals$U_TDC[,conf_levels[[2]]][,1]) & (0 <= conf_intervals$U_TDC[,conf_levels[[2]]][,2])
  U_TDC_insig_1 <- (0 >= conf_intervals$U_TDC[,conf_levels[[3]]][,1]) & (0 <= conf_intervals$U_TDC[,conf_levels[[3]]][,2])
  
  mean_val$U_TDC_sig_10 <- !U_TDC_insig_10
  mean_val$U_TDC_sig_5 <- !U_TDC_insig_5
  mean_val$U_TDC_sig_1 <- !U_TDC_insig_1
  
  # U_TDC
  L_TDC_insig_10 <- (0 >= conf_intervals$L_TDC[,conf_levels[[1]]][,1]) & (0 <= conf_intervals$L_TDC[,conf_levels[[1]]][,2])
  L_TDC_insig_5 <- (0 >= conf_intervals$L_TDC[,conf_levels[[2]]][,1]) & (0 <= conf_intervals$L_TDC[,conf_levels[[2]]][,2])
  L_TDC_insig_1 <- (0 >= conf_intervals$L_TDC[,conf_levels[[3]]][,1]) & (0 <= conf_intervals$L_TDC[,conf_levels[[3]]][,2])
  
  mean_val$L_TDC_sig_10 <- !L_TDC_insig_10
  mean_val$L_TDC_sig_5 <- !L_TDC_insig_5
  mean_val$L_TDC_sig_1 <- !L_TDC_insig_1
  
  return(list(aggregated_results = mean_val,
              raw_simulations = df_return)
  )
}

# Results df buildup
df_return_from_agg <- function(aggregated_results){
  agg_str <- aggregated_results %>%
    mutate(`Upper Tail 1` = as.character(round(U_TDC, 4))) %>%
    mutate(`Lower Tail 1` = as.character(round(L_TDC, 4))) %>%
    mutate(`Upper Tail` = ifelse(U_TDC_sig_1, paste(`Upper Tail 1`, '***', sep = ''),
                                 ifelse(U_TDC_sig_5, paste(`Upper Tail 1`, '**', sep = ''),
                                        ifelse(U_TDC_sig_10, paste(`Upper Tail 1`, '*', sep = ''),
                                               `Upper Tail 1`)))) %>%
    mutate(`Lower Tail` = ifelse(L_TDC_sig_1, paste(`Lower Tail 1`, '***', sep = ''),
                                 ifelse(L_TDC_sig_5, paste(`Lower Tail 1`, '**', sep = ''),
                                        ifelse(L_TDC_sig_10, paste(`Lower Tail 1`, '*', sep = ''),
                                               `Lower Tail 1`)))) %>%
    select(CountryGroup = Group.1, `Upper Tail`, `Lower Tail`)
  df_ret <- data.frame(matrix(data = '',nrow = 17, ncol = 17), stringsAsFactors = FALSE)
  
  df_ret[1,3] <- 'Brazil'
  df_ret[1,6] <- 'Russia'
  df_ret[1,9] <- 'India'
  df_ret[1,12] <- 'China'
  df_ret[1,15] <- 'South Africa'
  
  df_ret[3,1] <- 'Brazil'
  df_ret[6,1] <- 'Russia'
  df_ret[9,1] <- 'India'
  df_ret[12,1] <- 'China'
  df_ret[15,1] <- 'South Africa'
  
  df_ret[2,3] <- 'Financial'
  df_ret[2,4] <- 'Industrial'
  df_ret[2,5] <- 'Resource'
  df_ret[2,6] <- 'Financial'
  df_ret[2,7] <- 'Industrial'
  df_ret[2,8] <- 'Resource'
  df_ret[2,9] <- 'Financial'
  df_ret[2,10] <- 'Industrial'
  df_ret[2,11] <- 'Resource'
  df_ret[2,12] <- 'Financial'
  df_ret[2,13] <- 'Industrial'
  df_ret[2,14] <- 'Resource'
  df_ret[2,15] <- 'Financial'
  df_ret[2,16] <- 'Industrial'
  df_ret[2,17] <- 'Resource'
  
  df_ret[3,2] <- 'Financial'
  df_ret[4,2] <- 'Industrial'
  df_ret[5,2] <- 'Resource'
  df_ret[6,2] <- 'Financial'
  df_ret[7,2] <- 'Industrial'
  df_ret[8,2] <- 'Resource'
  df_ret[9,2] <- 'Financial'
  df_ret[10,2] <- 'Industrial'
  df_ret[11,2] <- 'Resource'
  df_ret[12,2] <- 'Financial'
  df_ret[13,2] <- 'Industrial'
  df_ret[14,2] <- 'Resource'
  df_ret[15,2] <- 'Financial'
  df_ret[16,2] <- 'Industrial'
  df_ret[17,2] <- 'Resource'
  
  for( i in 3:nrow(df_ret)){
    df_ret[i,i] = '.'
  }
  
  #agg_row_count =79, 80, 81
  for(agg_row_count in 1:nrow(agg_str)){
    agg_row <- agg_str[agg_row_count,]
    x <- as.character(agg_row$CountryGroup)
    cs_cs <- strsplit(x, split = ' and ') %>%
      unlist()
    
    cs_1 <- cs_cs[1] %>%
      stringr::str_replace_all('_', ' ') %>%
      lapply(FUN = simpleCap) %>%
      unlist() %>%
      country_sector()
    
    cs_2 <- cs_cs[2] %>%
      stringr::str_replace_all('_', ' ') %>%
      lapply(FUN = simpleCap) %>%
      unlist() %>%
      country_sector()
    
    # Row
    R_Country <- which(cs_1$Country == df_ret[,1])
    R_Sector <- which(cs_1$Sector == df_ret[,2])
    row_select <- R_Sector[which(R_Sector >= R_Country & R_Sector <= (R_Country + 2))]
    
    # Column
    C_Country <- which(cs_2$Country == as.character(df_ret[1,]))
    C_Sector <- which(cs_2$Sector == as.character(df_ret[2,]))
    col_select <- C_Sector[which(C_Sector >= C_Country & C_Sector <= (C_Country + 2))]
    
    if(row_select > col_select){
      df_ret[row_select, col_select] <- agg_row$`Lower Tail`  
      df_ret[col_select, row_select] <- agg_row$`Upper Tail`
    } else {
      df_ret[row_select, col_select] <- agg_row$`Upper Tail`  
      df_ret[col_select, row_select] <- agg_row$`Lower Tail`    
    }
    
  }
  return(df_ret)
}
