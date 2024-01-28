getIndexes <- function(dfs, s_df, dbs){
  indexes <- array(integer(3 * s_df * dbs), dim = c(3,s_df,dbs))
  for(i in 1:s_df){
    for(j in 1:dbs){
      indexes[1,i,j] <- i + s_df
      indexes[2,i,j] <- nrow(dfs[[j]]) - s_df - 12 + i
      indexes[3,i,j] <- nrow(dfs[[j]]) - s_df + i
    }
  }
  return(indexes)
}


initializeStates <- function(dfs, Quantiles, dbs){
  Quantiles[5] <- Quantiles[5] + 1
  for(i in 1:dbs){
    for(n in 1:nrow(dfs[[i]])){
      if(is.na(dfs[[i]][[n, 'CGM']])){
        
      }
      else{
        dfs[[i]][n, 'State'] <- which.max(dfs[[i]][[n, 'CGM']] < Quantiles)
      }
    }
  }
  return(dfs)
}




initializeDist <- function(df, target, filter){
  formula <- as.formula(paste(target,"~",filter))
  Mean <- unlist(aggregate(formula, df, mean)[[target]])
  Variance <- unlist(aggregate(formula, df, var)[[target]])
  Scale <- Variance/Mean
  Shape <- Mean/Scale
  return(list(shape = Shape, scale = Scale))
}


getTrans <- function(df){
  trans <- matrix(integer(25), nrow = 5, byrow = TRUE)
  df <- df[!is.na(df$time_spend),]
  for(i in 1:nrow(df)){
    trans[df[i, 'coming_from'], df[i, 'going_to']] <- trans[df[i, 'coming_from'], df[i, 'going_to']] + 1
  }
  trans <- apply(trans, 2, function(x) x / rowSums(trans))
  print(trans)
  return(trans)
}



getSoj <- function(dfs, dbs){
  
  df <- data.frame(coming_from = c(), going_to = c(), time_spend = c(), patient_number = c())
  for(j in 1:dbs){
    series <- dfs[[j]]
    print(series)
    a <- series[1] #current state
    b <- c() # unique sequence of states
    c <- c() # time spend in each state
    count <- 0
    for(i in 1:length(series)){
      if(i != length(series) & is.na(series[i]) & is.na(series[i+1])){
        count <- count + 1
      }
      else if((is.na(a) & !is.na(series[i])) | (!is.na(a) & is.na(series[i])) | (a != series[i])){
        b <- append(b, a)
        c <- append(c, count)
        count <- 1
        a <- series[i]
      }
      else if((is.na(a) & is.na(series[i])) | (a == series[i])){
        count <- count + 1
      }
    }
    
    if((is.na(b[length(b)]) & !is.na(series[length(series)])) | (!is.na(b[length(b)]) & is.na(series[length(series)])) | (b[length(b)] != series[length(series)])){
      b <- append(b, series[length(series)])
      c <- append(c, 1)
    }

    for(i in 1:length(b)){
      if(i == 1 | i == length(b)){
        c[i] = NA
      }
      else if(is.na(b[i-1]) | is.na(b[i+1]) | is.na(b[i])){
        c[i] = NA
      }
      else{
        
      }
    }
    coming_from <- append(NA, b)
    going_to <- append(b, NA)
    time_spend <- append(NA, c)
    patient_number <- rep(j, length(coming_from))
    df_tmp <- data.frame(coming_from = coming_from, going_to = going_to, time_spend = time_spend, patient_number = patient_number)
    df <- rbind(df, df_tmp)
  }
  return(df)
}
