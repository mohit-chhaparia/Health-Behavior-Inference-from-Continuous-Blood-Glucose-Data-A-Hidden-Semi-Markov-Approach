interpolate_CGM <- function(df){
  timestamps <- as.POSIXct(df$`Timestamp (YYYY-MM-DDThh:mm:ss)`, format="%Y-%m-%dT%H:%M:%S")
  start <- timestamps[1]
  end <- timestamps[length(timestamps)]
  timeSeq <- seq(from = start, to = end, by = "5 min")
  
  
  # Assuming time_seq is a vector of POSIXct times
  # Apply the function to each time point in the sequence
  interpolated_values <- sapply(timeSeq, interpolate_between_nearest, df = df)
  
  # Combine the time points and their interpolated values into a new data frame
  result_df <- data.frame(
    `Timestamp (YYYY-MM-DDThh:mm:ss)` = timeSeq,
    CGM = interpolated_values
  )
  
  return(result_df)
}

interpolate_between_nearest <- function(time_point, df) {
  
  # Check if there's a value exactly at the target time point
  exact_value <- df %>%
    filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` == time_point) %>%
    select(CGM) %>%
    pull()
  
  # If there's an exact match, return it without interpolation
  if(length(exact_value) == 1) {
    return(exact_value)
  }
  
  # Find the nearest point before the target time within -5 minutes
  before_point <- df %>%
    filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` < time_point) %>%
    filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` >= (time_point - minutes(5))) %>%
    slice_max(`Timestamp (YYYY-MM-DDThh:mm:ss)`)
  
  # Find the nearest point after the target time within +5 minutes
  after_point <- df %>%
    filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` > time_point) %>%
    filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` <= (time_point + minutes(5))) %>%
    slice_min(`Timestamp (YYYY-MM-DDThh:mm:ss)`)
  
  if(nrow(before_point) == 0 && nrow(after_point) == 1){
    tp <- after_point$`Timestamp (YYYY-MM-DDThh:mm:ss)`
    before_point <- df %>%
      filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` < time_point) %>%
      filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` >= (tp - minutes(10))) %>%
      slice_max(`Timestamp (YYYY-MM-DDThh:mm:ss)`)
  }
  else if(nrow(before_point) == 1 && nrow(after_point) == 0){
    tp <- before_point$`Timestamp (YYYY-MM-DDThh:mm:ss)`
    after_point <- df %>%
      filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` > time_point) %>%
      filter(`Timestamp (YYYY-MM-DDThh:mm:ss)` <= (tp + minutes(10))) %>%
      slice_min(`Timestamp (YYYY-MM-DDThh:mm:ss)`)
  }
  
  # If either before or after point is missing, return NA
  if(nrow(before_point) == 0 | nrow(after_point) == 0) {
    return(NA)
  }
  
  # Use approx to interpolate between the two nearest points
  approximated_value <- round(approx(
    x = c(before_point$`Timestamp (YYYY-MM-DDThh:mm:ss)`, after_point$`Timestamp (YYYY-MM-DDThh:mm:ss)`),
    y = c(before_point$CGM, after_point$CGM),
    xout = time_point
  )$y, 2)
  
  return(approximated_value)
}

timeDiff <- function(dfs){
  df <- data.frame(
    id = integer(),
    uniqueValues = I(list()),
    outOfRangeValues = I(list()),
    outOfRangeIndexes = I(list())
  )
  ignore_vals <- 295:305
  for(i in 1:20){
    unique_vals <- unique(diff(as.POSIXct(dfs[[i]]$`Timestamp (YYYY-MM-DDThh:mm:ss)`, format="%Y-%m-%dT%H:%M:%S")))
    outOfRange_vals <- unique_vals[!unique_vals %in% ignore_vals]
    outOfRange_Ind <- which(diff(as.POSIXct(dfs[[i]]$`Timestamp (YYYY-MM-DDThh:mm:ss)`, format="%Y-%m-%dT%H:%M:%S")) %in% outOfRange_vals)
    
    new_records <- list(
      list(id=i,
           uniqueValues=I(list(unique_vals)),
           outOfRangeValues=I(list(outOfRange_vals)),
           outOfRangeIndexes=I(list(outOfRange_Ind)))
    )
    df <- rbind(df, as.data.frame(new_records))
  }
  return(df)
}


library(readxl)
library(writexl)


setwd("SET_WORKING_DIRECTORY")


file_path_pattern_read <- "FILE_NAME-%d.xlsx"
file_path_pattern_write <- "FILE_NAME-%d.csv"

dfs=list();
df1=lst();
for (i in 1:20) {
  print(i)
  file_name_read <- sprintf(file_path_pattern_read, 1000 + i) # creates the filename
  df <- read_excel(file_name_read)
  df$CGM <- as.integer(ifelse(df$CGM == 'Low', '40', ifelse(df$CGM == 'High', '400', df$CGM)))
  df$`Timestamp (YYYY-MM-DDThh:mm:ss)` <- as.POSIXct(df$`Timestamp (YYYY-MM-DDThh:mm:ss)`, format="%Y-%m-%dT%H:%M:%S")
  df1[[i]] <- df
  df <- df[!is.na(df$CGM),]
  df <- df[df$`Event Type` == 'EGV',]
  df <- interpolate_CGM(df)
  df$`Patient Number` <- i
  dfs[[i]] <- df # reads the file and assigns it to the variable
  file_name_write <- sprintf(file_path_pattern_write, 1000 + i)
  write.csv(df, file_name_write)
  rm(df)
}
df_time <- timeDiff(df1)
write.csv(df_time, "FILE_NAME.csv")