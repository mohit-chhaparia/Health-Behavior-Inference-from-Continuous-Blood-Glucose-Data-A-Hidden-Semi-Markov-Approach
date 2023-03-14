# LOAD LIBRARIES
library(readxl)
library(writexl)
library(mvtnorm)
library(mhsmm)
library(IMIFA)
library(ggplot2)
library(gridExtra)

setwd("SET_WORKING_DIRECTORY")

# LOAD ALL THE DATASETS
df_P1 <- read_excel("PATIENT_1_DATASET.xlsx")
df_P1 <- df_P1[is.na(df_P1$CGM) == 0,]
df_P2 <- read_excel("PATIENT_2_DATASET.xlsx")
df_P2 <- df_P2[is.na(df_P2$CGM) == 0,]
df_P3 <- read_excel("PATIENT_3_DATASET.xlsx")
df_P3 <- df_P3[is.na(df_P3$CGM) == 0,]
df_P4 <- read_excel("PATIENT_4_DATASET.xlsx")
df_P4 <- df_P4[is.na(df_P4$CGM) == 0,]
df_P5 <- read_excel("PATIENT_5_DATASET.xlsx")
df_P5 <- df_P5[is.na(df_P5$CGM) == 0,]
df_P6 <- read_excel("PATIENT_6_DATASET.xlsx")
df_P6 <- df_P6[is.na(df_P6$CGM) == 0,]
df_P7 <- read_excel("PATIENT_7_DATASET.xlsx")
df_P7 <- df_P7[is.na(df_P7$CGM) == 0,]
df_P8 <- read_excel("PATIENT_8_DATASET.xlsx")
df_P8 <- df_P8[is.na(df_P8$CGM) == 0,]
df_P9 <- read_excel("PATIENT_9_DATASET.xlsx")
df_P9 <- df_P9[is.na(df_P9$CGM) == 0,]
df_P10 <- read_excel("PATIENT_10_DATASET.xlsx")
df_P10 <- df_P10[is.na(df_P10$CGM) == 0,]
df_P11 <- read_excel("PATIENT_11_DATASET.xlsx")
df_P11 <- df_P11[is.na(df_P11$CGM) == 0,]
df_P12 <- read_excel("PATIENT_12_DATASET.xlsx")
df_P12 <- df_P12[is.na(df_P12$CGM) == 0,]
df_P13 <- read_excel("PATIENT_13_DATASET.xlsx")
df_P13 <- df_P13[is.na(df_P13$CGM) == 0,]
df_P14 <- read_excel("PATIENT_14_DATASET.xlsx")
df_P14 <- df_P14[is.na(df_P14$CGM) == 0,]
df_P15 <- read_excel("PATIENT_15_DATASET.xlsx")
df_P15 <- df_P15[is.na(df_P15$CGM) == 0,]
df_P16 <- read_excel("PATIENT_16_DATASET.xlsx")
df_P16 <- df_P16[is.na(df_P16$CGM) == 0,]
df_P17 <- read_excel("PATIENT_17_DATASET.xlsx")
df_P17 <- df_P17[is.na(df_P17$CGM) == 0,]
df_P18 <- read_excel("PATIENT_18_DATASET.xlsx")
df_P18 <- df_P18[is.na(df_P18$CGM) == 0,]
df_P19 <- read_excel("PATIENT_19_DATASET.xlsx")
df_P19 <- df_P19[is.na(df_P19$CGM) == 0,]
df_P20 <- read_excel("PATIENT_20_DATASET.xlsx")
df_P20 <- df_P20[is.na(df_P20$CGM) == 0,]


# USER-DEFINED MSTEP.GAMMA FUNCTION
mstep.gamma <- function(x, wt=NULL){
  tmp1 = apply(wt, 2, function(w) weighted.params(x,w))
  return(list(shape=tmp1[1,], scale=tmp1[2,]))
}

weighted.params <- function(x, wt){
  tmp1=weighted.mean(x,wt)
  sum1=0
  for(i in 1:length(x)){
    sum1 = sum1 + wt[i]*(x[i]-tmp1)^2
  }
  tmp2 = sum1/sum(wt)
  scale1 = tmp2/tmp1
  shape1 = tmp1/scale1
  return(c(shape1,scale1))
}


# USER-DEFINED RGAMMA FUNCTION
rgamma.hsmm <- function(j, model){
  rgamma(n=1, shape = model$parms.emission$shape[j], scale = model$parms.emission$scale[j])
}


# USER DEFINED DGAMMA FUNCTION
dgamma.hsmm <- function(x, j, model){
  dgamma(x=x, shape = model$parms.emission$shape[j], scale = model$parms.emission$scale[j])
}


# USER-DEFINED TIME IN STATE FUNCTION
time_in_state <- function(pred_list){
  last_predicted_state <- pred_list[length(pred_list)]
  count <- 1
  for(i in 1:length(pred_list)){
    if(last_predicted_state == pred_list[length(pred_list)-i]){
      count <- count + 1
    }
    else{
      return(count)
    }
  }
  return(count)
}


# USER-DEFINED PREDICT FUNCTION
predict_new <- function(model, observed_CGM_list, predict_at = 60){
  pred_from_hsmm <- predict(model, observed_CGM_list)$s
  
  
  last_state <- pred_from_hsmm[length(pred_from_hsmm)]
  original_state <- last_state
  tau <- time_in_state(pred_from_hsmm)
  cumsum_d <- matrix(integer(length(model$model$sojourn$d)), nrow = length(model$model$sojourn$d) / 5, byrow = TRUE)
  for(i in 1:5){
    cumsum_d[,i] <- cumsum(model$model$sojourn$d[,i])
  }
  runif_vals <- runif(n = 10000, min = cumsum_d[tau,last_state] , max = max(cumsum_d[,last_state]))
  state_predictions <- matrix(integer(10000 * predict_at / 5), nrow = 10000, byrow = TRUE)
  runif_records <- matrix(integer(10000 * predict_at / 5), nrow = 10000, byrow = TRUE)
  time_matrix <- matrix(integer(10000 * predict_at / 5), nrow = 10000, byrow = TRUE)
  for(i in 1:10000){
    runif_records[i,1] <- runif_vals[i]
    FCV <- 1
    flag <- 0
    if(length(unique(cumsum_d[, last_state]==runif_vals[i]))==2){
      time <- which(cumsum_d[, last_state]==runif_vals[i]) - tau
      time_matrix[i,1] <- time
    }
    else{
      time <- which.max(cumsum_d[, last_state] > runif_vals[i])-1 - tau
      time_matrix[i,1] <- time
    }
    if(time <= tau){
      truncated_time <- 0
      flag <- 1
    }
    else{
      truncated_time <- time
    }
    repeat{
      if(flag==1){
        time <- which.max(cumsum_d[, last_state] > runif(n = 1, min = 0 , max = max(cumsum_d[,last_state])))
        time_matrix[i,FCV] <- time
        if(time < 1){

        }
      }
      while(time >= 1 & FCV <= floor(predict_at/5)){
        state_predictions[i,FCV] <- last_state
        FCV <- FCV + 1
        time <- time - 1
      }
      if(FCV == floor(predict_at/5)+1){
        last_state <- original_state
        break
      }
      flag <- 1
      if(last_state == 1){
        last_state <- 2
      }
      else if(last_state == 5){
        last_state <- 4
      }
      else{
        new_runif <- runif(1,0,1)
        transit_prob <- new_runif
        runif_records[i, FCV] <- new_runif
        if(model$model$transition[last_state,last_state-1] <= model$model$transition[last_state,last_state+1]){
          if(transit_prob <= model$model$transition[last_state,last_state-1]){
            last_state <- last_state - 1
          }
          else{
            last_state <- last_state + 1
          }
        }
        else{
          if(transit_prob <= model$model$transition[last_state,last_state+1]){
            last_state <- last_state + 1
          }
          else{
            last_state <- last_state - 1
          }
        }
      }
    }
  }
  # FRACTION OF BEING IN A PARTICULAR STATE AT EACH TIME POINT
  fraction_of_state <- matrix(integer(predict_at), nrow = 5, byrow = TRUE)
  for(i in 1:(predict_at/5)){
    pred_count <- table(state_predictions[,i])
    m <- 1
    for(j in 1:5){
      pred_count_header_str <- names(pred_count)
      pred_count_header_int <- strtoi(pred_count_header_str)
      if(m > length(pred_count_header_int)){
        
      }
      else{
        if(j == pred_count_header_int[m]){
          fraction_of_state[j,i] <- pred_count[pred_count_header_str[m]]
          m <- m + 1
        }
        else{
          fraction_of_state[j,i] <- 0
        }
      }
    }
    fraction_of_state[,i] <- fraction_of_state[,i] / length(state_predictions[,i])
  }
  prob_of_state <- matrix(integer(15), nrow = 3, byrow = TRUE)
  # ROW 1 GIVES PROBABILITY OF HYPOGLYCEMIA FOR EACH STATE
  # ROW 2 GIVES PROBABILITY OF NORMAL RANGE FOR EACH STATE
  # ROW 3 GIVES PROBABILITY OF HYPERGLYCEMIA FOR EACH STATE
  for(i in 1:5){
    prob_of_state[1,i] <- pgamma(70, shape = model$model$parms.emission$shape[i], scale = model$model$parms.emission$scale[i])
    prob_of_state[3,i] <- 1 - pgamma(180, shape = model$model$parms.emission$shape[i], scale = model$model$parms.emission$scale[i])
    prob_of_state[2,i] <- 1 - prob_of_state[3,i] - prob_of_state[1,i]
  }
  
  # ROW 1 GIVES PROBABILITY OF HYPOGLYCEMIA FOR EACH TIME POINT
  # ROW 2 GIVES PROBABILITY OF NORMAL RANGE FOR EACH TIME POINT
  # ROW 3 GIVES PROBABILITY OF HYPERGLYCEMIA FOR EACH TIME POINT
  prob_at_time_points <- prob_of_state %*% fraction_of_state
  
  emis_prediction <- matrix(integer(10000 * predict_at / 5), nrow = 10000, byrow = TRUE)
  for(i in 1:nrow(state_predictions)){
    for(j in 1:ncol(state_predictions)){
      emis_prediction[i,j] <- rgamma(
        n = 1,
        shape = model$model$parms.emission$shape[state_predictions[i,j]],
        scale = model$model$parms.emission$scale[state_predictions[i,j]]
        )
    }
  }
  emis_prob_Tm_Pts_30mins <- integer(predict_at / 5)
  emis_prob_Tm_Pts_60mins <- integer(predict_at / 5)
  for(i in 1:nrow(state_predictions)){
    emis_prob_Tm_Pts_30mins[i] <- any(emis_prediction[i,1:6] < 70)
    emis_prob_Tm_Pts_60mins[i] <- any(emis_prediction[i,] < 70)
  }
  Prob_of_Hypo_30min <- sum(emis_prob_Tm_Pts_30mins) / 10000
  Prob_of_Hypo_60min <- sum(emis_prob_Tm_Pts_60mins) / 10000
  
  list_to_return <- list(
    pred_from_hsmm[length(pred_from_hsmm)],
    state_predictions,
    fraction_of_state,
    prob_of_state,
    prob_at_time_points,
    time_matrix,
    Prob_of_Hypo_30min,
    Prob_of_Hypo_60min,
    emis_prediction
    )
  
  return(list_to_return)
}


# NUMBER OF WINDOWS OR NUMBER OF PREDICTIONS OF 12 TIME POINTS EACH FOR EACH PATIENT
s_df = 288

# GENERATING DATA START AND END POINTS FOR EACH ITERATION AND EACH DATASET
indexes <- array(integer(3 * s_df * 20), dim = c(3,s_df,20))
dfs <- list(df_P1, df_P2, df_P3, df_P4, df_P5,
            df_P6, df_P7, df_P8, df_P9, df_P10,
            df_P11, df_P12, df_P13, df_P14, df_P15,
            df_P16, df_P17, df_P18, df_P19, df_P20
            )
for(i in 1:s_df){
  for(j in 1:20){
    indexes[1,i,j] <- i
    indexes[2,i,j] <- nrow(dfs[[j]]) - 300 + i
    indexes[3,i,j] <- nrow(dfs[[j]]) - 288 + i
  }
}

# ASSIGNING INITIAL STATES AS PER QUANTILES OF THE GLUCOSE VALUES OF ALL THE DATASETS COMBINED
df_quantile <- rbind(df_P1, df_P2, df_P3, df_P4, df_P5,
                     df_P6, df_P7, df_P8, df_P9, df_P10,
                     df_P11, df_P12, df_P13, df_P14, df_P15,
                     df_P16, df_P17, df_P18, df_P19, df_P20
                     )
Quantiles <- quantile(unlist(df_quantile[,'CGM']), probs = c(0.15, 0.375, 0.625, 0.85))
for(o in 1:20){
  for(n in 1:nrow(dfs[[o]])){
    if(is.na(dfs[[o]][n, 'CGM'])){
      
    }
    else if(dfs[[o]][n, 'CGM'] < Quantiles[1]){
      dfs[[o]][n, 'State'] <- 1
    }
    else if(dfs[[o]][n, 'CGM'] < Quantiles[2]){
      dfs[[o]][n, 'State'] <- 2
    }
    else if(dfs[[o]][n, 'CGM'] < Quantiles[3]){
      dfs[[o]][n, 'State'] <- 3
    }
    else if(dfs[[o]][n, 'CGM'] < Quantiles[4]){
      dfs[[o]][n, 'State'] <- 4
    }
    else{
      dfs[[o]][n, 'State'] <- 5
    }
  }
}
df_P1 <- dfs[[1]]
df_P2 <- dfs[[2]]
df_P3 <- dfs[[3]]
df_P4 <- dfs[[4]]
df_P5 <- dfs[[5]]
df_P6 <- dfs[[6]]
df_P7 <- dfs[[7]]
df_P8 <- dfs[[8]]
df_P9 <- dfs[[9]]
df_P10 <- dfs[[10]]
df_P11 <- dfs[[11]]
df_P12 <- dfs[[12]]
df_P13 <- dfs[[13]]
df_P14 <- dfs[[14]]
df_P15 <- dfs[[15]]
df_P16 <- dfs[[16]]
df_P17 <- dfs[[17]]
df_P18 <- dfs[[18]]
df_P19 <- dfs[[19]]
df_P20 <- dfs[[20]]


# INITIALIZING DATAFRAME TO STORE OUTPUT
Tm_Pts <- 12
# 20 DATASET * 1 PREDICTION PER ITERATION * 12 TIMEPOINTS PER ITERATION * 288 ITERATIONS OR PREDICTIONS
Time_Point <- integer(20 * 1 * s_df * Tm_Pts)
State <- integer(20 * 1 * s_df * Tm_Pts)
CGM <- integer(20 * 1 * s_df * Tm_Pts)
Predicted_State <- integer(20 * 1 * s_df * Tm_Pts)
Hypo_Prob <- integer(20 * 1 * s_df * Tm_Pts)
Patient_Number <- integer(20 * 1 * s_df * Tm_Pts)
Actual_State_Com <- integer(20 * 1 * s_df * Tm_Pts)
Pred_State_Com <- integer(20 * 1 * s_df * Tm_Pts)
Prob_of_Hypo_30min <- integer(20 * 1 * s_df * Tm_Pts)
Prob_of_Hypo_60min <- integer(20 * 1 * s_df * Tm_Pts)
Time_Stamp <- as.character(integer(20 * 1 * s_df * Tm_Pts))

output <- data.frame(
  Patient_Number,
  Time_Point,
  State,
  CGM,
  Predicted_State,
  Hypo_Prob,
  Actual_State_Com,
  Pred_State_Com,
  Prob_of_Hypo_30min,
  Prob_of_Hypo_60min,
  Time_Stamp
  )

output_index <- 1


# MODEL TRAINING AND PREDICTION
for(i in 1:(s_df)){
  print(i)

  df_train <- rbind(df_P1[indexes[1,i,1]:indexes[2,i,1],], df_P2[indexes[1,i,2]:indexes[2,i,2],],
                    df_P3[indexes[1,i,3]:indexes[2,i,3],], df_P4[indexes[1,i,4]:indexes[2,i,4],],
                    df_P5[indexes[1,i,5]:indexes[2,i,5],], df_P6[indexes[1,i,6]:indexes[2,i,6],],
                    df_P7[indexes[1,i,7]:indexes[2,i,7],], df_P8[indexes[1,i,8]:indexes[2,i,8],],
                    df_P9[indexes[1,i,9]:indexes[2,i,9],], df_P10[indexes[1,i,10]:indexes[2,i,10],],
                    df_P11[indexes[1,i,11]:indexes[2,i,11],], df_P12[indexes[1,i,12]:indexes[2,i,12],],
                    df_P13[indexes[1,i,13]:indexes[2,i,13],], df_P14[indexes[1,i,14]:indexes[2,i,14],],
                    df_P15[indexes[1,i,15]:indexes[2,i,15],], df_P16[indexes[1,i,16]:indexes[2,i,16],],
                    df_P17[indexes[1,i,17]:indexes[2,i,17],], df_P18[indexes[1,i,18]:indexes[2,i,18],],
                    df_P19[indexes[1,i,19]:indexes[2,i,19],], df_P20[indexes[1,i,20]:indexes[2,i,20],]
  )

  if(any(rowSums(is.na(df_train))==8)){
    df_train <- na.omit(df_train)
  }

  rownames(df_train) <- 1:nrow(df_train)

  df_test <- dfs
  
  
  
  j = 1
  Coming_From <- c(0)
  Entering <- c(0)
  Time_Spend <- c(0)
  CGM_Entering <- c(0)
  CGM_Leaving <- c(0)
  df_Sojourn <- data.frame(CGM_Entering, CGM_Leaving, Coming_From, Entering, Time_Spend)
  df_Sojourn[1,"Coming_From"] = df_train[1,"State"]
  df_Sojourn[1, "CGM_Entering"] = df_train[1,"CGM"]
  count = 1
  j = 1

  for(k in 1:nrow(df_train)){

    if(is.na(df_train[k,"State"]) | is.na(df_train[k+1,"State"])) {
      
    } else if(df_train[k,"State"]==df_train[k+1,"State"]) {
      count <- count + 1
    } else {
      df_Sojourn[j+1,"Coming_From"] = df_train[k+1,"State"]
      df_Sojourn[j,"Time_Spend"] = count
      df_Sojourn[j,"Entering"] = df_train[k+1,"State"]
      df_Sojourn[j,"CGM_Leaving"] = df_train[k, "CGM"]
      df_Sojourn[j+1, "CGM_Entering"] = df_train[k+1,"CGM"]
      count = 1
      j = j + 1
    }
  }

  df_Sojourn <- na.omit(df_Sojourn)
  rownames(df_Sojourn) <- 1:nrow(df_Sojourn)
  
  index <- c(1, 1, 1, 1, 1)
  State1 <- c(0)
  State2 <- c(0)
  State3 <- c(0)
  State4 <- c(0)
  State5 <- c(0)
  df_StateTimes <- data.frame(State1, State2, State3, State4, State5)
  for(j in 1:nrow(df_Sojourn)) {
    if(df_Sojourn[j,"Coming_From"] == 1) {
      df_StateTimes[index[1],"State1"] = df_Sojourn[j,"Time_Spend"]
      index[1] = index[1] + 1
    }
    else if(df_Sojourn[j,"Coming_From"] == 2) {
      df_StateTimes[index[2],"State2"] = df_Sojourn[j,"Time_Spend"]
      index[2] = index[2] + 1
    }
    else if(df_Sojourn[j,"Coming_From"] == 3) {
      df_StateTimes[index[3],"State3"] = df_Sojourn[j,"Time_Spend"]
      index[3] = index[3] + 1
    }
    else if(df_Sojourn[j,"Coming_From"] == 4) {
      df_StateTimes[index[4],"State4"] = df_Sojourn[j,"Time_Spend"]
      index[4] = index[4] + 1
    }
    else {
      df_StateTimes[index[5],"State5"] = df_Sojourn[j,"Time_Spend"]
      index[5] = index[5] + 1
    }
  }
  rownames(df_StateTimes) <- 1:nrow(df_StateTimes)

  J <- 5
  init <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  
  # INITIALIZING THE TRANSITION MATRIX
  count <- integer(25)
  for(j in 1:nrow(df_Sojourn)) {
    tmp <- count[(unlist(df_Sojourn[j, "Coming_From"]) - 1) * 5 + unlist(df_Sojourn[j, "Entering"])]
    count[(unlist(df_Sojourn[j, "Coming_From"]) - 1) * 5 + unlist(df_Sojourn[j, "Entering"])] <- tmp + 1
  }

  trans <- matrix(round((count / sum(count)), 4), nrow = 5, byrow = TRUE)
  for(k in 1:5) {
    for(j in 1:5){
      if(j-k==1 | j-k == -1) {
        
      } else {
        trans[k,j] <- 0.0000
      }
    }
  }

  for(j in 1:5){
    trans[j,] <- trans[j,] / sum(trans[j,])
  }

  # INITIALIZING THE EMISSION DISTRIBUTION
  Mean_Emis <- integer(5)
  Count_Emis <- integer(5)
  for(j in 1:nrow(df_train)) {
    if(is.na(df_train[j,"CGM"])){
      
    }
    else if(df_train[j,"CGM"] < Quantiles[1]){
      Mean_Emis[1] <- Mean_Emis[1] + df_train[j,"CGM"]
      Count_Emis[1] <- Count_Emis[1] + 1
    }
    else if(df_train[j,"CGM"] < Quantiles[2]) {
      Mean_Emis[2] <- Mean_Emis[2] + df_train[j,"CGM"]
      Count_Emis[2] <- Count_Emis[2] + 1
    }
    else if(df_train[j,"CGM"] < Quantiles[3]) {
      Mean_Emis[3] <- Mean_Emis[3] + df_train[j,"CGM"]
      Count_Emis[3] <- Count_Emis[3] + 1
    }
    else if(df_train[j,"CGM"] < Quantiles[4]) {
      Mean_Emis[4] <- Mean_Emis[4] + df_train[j,"CGM"]
      Count_Emis[4] <- Count_Emis[4] + 1
    }
    else {
      Mean_Emis[5] <- Mean_Emis[5] + df_train[j,"CGM"]
      Count_Emis[5] <- Count_Emis[5] + 1
    }
  }

  for(j in 1:5){
    Mean_Emis[j] <- unlist(Mean_Emis[j])/unlist(Count_Emis[j])
  }
  Mean_Emis <- unlist(Mean_Emis)

  Var_Emis <- integer(5)
  
  for(j in 1:nrow(df_train)) {
    if(is.na(df_train[j, "CGM"])){
      
    }
    else if(df_train[j,"CGM"] < Quantiles[1]) {
      Var_Emis[1] <- Var_Emis[1] + (Mean_Emis[1] - df_train[j,"CGM"])^2
    }
    else if(df_train[j,"CGM"] < Quantiles[2]) {
      Var_Emis[2] <- Var_Emis[2] + (Mean_Emis[2] - df_train[j,"CGM"])^2
    }
    else if(df_train[j,"CGM"] < Quantiles[3]) {
      Var_Emis[3] <- Var_Emis[3] + (Mean_Emis[3] - df_train[j,"CGM"])^2
    }
    else if(df_train[j,"CGM"] < Quantiles[4]) {
      Var_Emis[4] <- Var_Emis[4] + (Mean_Emis[4] - df_train[j,"CGM"])^2
    }
    else {
      Var_Emis[5] <- Var_Emis[5] + (Mean_Emis[5] - df_train[j,"CGM"])^2
    }
  }

  Alpha_Emis <- integer(5)
  Beta_Emis <- integer(5)
  for(j in 1:5) {
    Var_Emis[j] <- unlist(Var_Emis[j])/(unlist(Count_Emis[j]) - 1)
    Beta_Emis[j] <- unlist(Var_Emis[j])/unlist(Mean_Emis[j])
    Alpha_Emis[j] <- unlist(Mean_Emis[j])/unlist(Beta_Emis[j])
  }

  Var_Emis <- unlist(Var_Emis)
  emis <- list(shape = unlist(Alpha_Emis), scale = unlist(Beta_Emis))
  
  len_Sojourn <- integer(5)
  Mean_Sojourn <- integer(5)
  Var_Sojourn <- integer(5)
  Alpha_Sojourn <- integer(5)
  Beta_Sojourn <- integer(5)
  
  # INITIALIZING THE SOJOURN DISTRIBUTION
  for(j in 1:5) {
    len_Sojourn[j] <- length(na.omit(unlist(df_StateTimes[,colnames(df_StateTimes)[j]])))
    Mean_Sojourn[j] <- mean(na.omit(unlist(df_StateTimes[,colnames(df_StateTimes)[j]])))
    Var_Sojourn[j] <- var(na.omit(unlist(df_StateTimes[,colnames(df_StateTimes)[j]])))
    Beta_Sojourn[j] <- Var_Sojourn[j]/Mean_Sojourn[j]
    Alpha_Sojourn[j] <- Mean_Sojourn[j]/Beta_Sojourn[j]
  }
  
  M <- 288
  d <- cbind(
    dgamma(1:M, shape = Alpha_Sojourn[1], scale = Beta_Sojourn[1]),
    dgamma(1:M, shape = Alpha_Sojourn[2], scale = Beta_Sojourn[2]),
    dgamma(1:M, shape = Alpha_Sojourn[3], scale = Beta_Sojourn[3]),
    dgamma(1:M, shape = Alpha_Sojourn[4], scale = Beta_Sojourn[4]),
    dgamma(1:M, shape = Alpha_Sojourn[5], scale = Beta_Sojourn[5])
  )
  
  N <- length(x = df_train$CGM)
  train <- list(x = df_train$CGM, N = N)
  class(train) <- "hsmm.data"

  # MODEL TRAINING
  startval <- hsmmspec(init, trans, emis, list(d = d, type = "gamma"), dens.emission = dgamma.hsmm, rand.emission = rgamma.hsmm, mstep = mstep.gamma)
  h.activity <- hsmmfit(train, startval, mstep = mstep.gamma, maxit = 100, M = M, lock.transition = FALSE, lock.d = FALSE, graphical = TRUE)
  
  # 1 SET OF PREDICTION OF 12 TIMEPOINTS FOR EACH PATIENT
  for(j in 1:20){

    for(m in 1:1){

      if(indexes[1,i,j]==0){
        
      }
      else{
        test <- list(x = (df_test[[j]])[indexes[1,i,j]:(indexes[2,i,j]+m-1),]$CGM, N = indexes[2,i,j] - indexes[1,i,j] + m)
        class(test) <- "hsmm.data"
        
        test_pred <- predict_new(h.activity, test, predict_at =  (5 * Tm_Pts))
        for(n in 1:Tm_Pts){

          output[output_index, 'Patient_Number'] <- j
          output[output_index,'State'] <- (df_test[[j]])[indexes[2,i,j]+m-1+n,]$State
          output[output_index,'CGM'] <- (df_test[[j]])[indexes[2,i,j]+m-1+n,]$CGM
          output[output_index, 'Time_Point'] <- n
          output[output_index, 'Predicted_State'] <- which.max(test_pred[[5]][,n])
          output[output_index, 'Hypo_Prob'] <- test_pred[[5]][1,n]
          output[output_index, 'Prob_of_Hypo_30min'] <- test_pred[[7]]
          output[output_index, 'Prob_of_Hypo_60min'] <- test_pred[[8]]
          output[output_index, 'Time_Stamp'] <- (df_test[[j]])[indexes[2,i,j]+m-1+n,]$`Timestamp (YYYY-MM-DDThh:mm:ss)`
          output_index <- output_index + 1
        }
      }
    }
  }
}


# GENERATES SENSITIVITY AND SPECIFICITY DATA FOR EACH TIMEPOINT USING THRESHOLDS FROM 0 TO 1 WITH A STEP SIZE OF 0.001

TPR <- matrix(integer(1001 * Tm_Pts), nrow = 1001)
FPR <- matrix(integer(1001 * Tm_Pts), nrow = 1001)
Ratio <- matrix(integer(1001 * Tm_Pts), nrow = 1001)

for(i in 1:Tm_Pts){

  df <- output[output$Time_Point == i,]
  for(j in 1:1001){

    ratio <- (j - 1) / 1000
    confusion_matrix <- matrix(integer(4), nrow = 2)

    for(m in 1:nrow(df)){
      if(is.na(df[m,'Hypo_Prob'])){
        
      }
      else if(df[m,'Hypo_Prob'] >= ratio){
        pred_state <- 1
      }
      else{
        pred_state <- df[m,'Predicted_State']
      }
      if(is.na(df[m,'State'])){
        
      }
      else if(pred_state == df[m,'State'] & pred_state == 1){
        confusion_matrix[1,1] <- confusion_matrix[1,1] + 1
      }
      else if(pred_state != 1 & df[m,'State'] != 1){
        confusion_matrix[2,2] <- confusion_matrix[2,2] + 1
      }
      else if(df[m,'State'] == 1){
        confusion_matrix[2,1] <- confusion_matrix[2,1] + 1
      }
      else{
        confusion_matrix[1,2] <- confusion_matrix[1,2] + 1
      }
    }
    # SENSITIVITY
    TPR[j,i] <- confusion_matrix[1,1] / (confusion_matrix[1,1] + confusion_matrix[2,1])
    # (1 - SPECIFICITY)
    FPR[j,i] <- confusion_matrix[1,2] / (confusion_matrix[1,2] + confusion_matrix[2,2])
    Ratio[j,i] <- ratio
  }
}


# PLOTS ROC CURVE FOR EACH TIMEPOINT USING THE ABOVE GENERATED DATA

p1 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 1], y = TPR[, 1])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 1") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p2 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 2], y = TPR[, 2])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 2") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p3 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 3], y = TPR[, 3])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 3") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p4 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 4], y = TPR[, 4])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 4") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p5 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 5], y = TPR[, 5])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 5") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p6 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 6], y = TPR[, 6])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 6") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p7 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 7], y = TPR[, 7])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 7") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p8 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 8], y = TPR[, 8])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 8") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p9 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 9], y = TPR[, 9])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 9") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p10 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 10], y = TPR[, 10])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 10") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p11 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 11], y = TPR[, 11])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 11") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

p12 <- ggplot() +
  geom_line(mapping = aes(x = FPR[, 12], y = TPR[, 12])) +
  xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle("TimePoint 12") +
  xlim(0,1) + ylim(0,1) + coord_fixed(ratio = 1)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 3)


# PLOTS PREDICTED PROBABILITY OF HYPOGLYCEMIA AGAINST THE OBSERVED GLUCOSE VALUES FOR EACH TIMEPOINT
par(mfrow=c(3, 4))
plot(x = output[output[, 'Time_Point'] == 1, 'CGM'] , y = output[output[, 'Time_Point'] == 1, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 1',
     )
plot(x = output[output[, 'Time_Point'] == 2, 'CGM'] , y = output[output[, 'Time_Point'] == 2, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 2',
     )
plot(x = output[output[, 'Time_Point'] == 3, 'CGM'] , y = output[output[, 'Time_Point'] == 3, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 3',
     )
plot(x = output[output[, 'Time_Point'] == 4, 'CGM'] , y = output[output[, 'Time_Point'] == 4, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 4',
     )
plot(x = output[output[, 'Time_Point'] == 5, 'CGM'] , y = output[output[, 'Time_Point'] == 5, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 5',
     )
plot(x = output[output[, 'Time_Point'] == 6, 'CGM'] , y = output[output[, 'Time_Point'] == 6, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 6',
     )
plot(x = output[output[, 'Time_Point'] == 7, 'CGM'] , y = output[output[, 'Time_Point'] == 7, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 7',
)
plot(x = output[output[, 'Time_Point'] == 8, 'CGM'] , y = output[output[, 'Time_Point'] == 8, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 8',
)
plot(x = output[output[, 'Time_Point'] == 9, 'CGM'] , y = output[output[, 'Time_Point'] == 9, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 9',
)
plot(x = output[output[, 'Time_Point'] == 10, 'CGM'] , y = output[output[, 'Time_Point'] == 10, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 10',
)
plot(x = output[output[, 'Time_Point'] == 11, 'CGM'] , y = output[output[, 'Time_Point'] == 11, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 11',
)
plot(x = output[output[, 'Time_Point'] == 12, 'CGM'] , y = output[output[, 'Time_Point'] == 12, 'Hypo_Prob'],
     xlab = 'Glucose (mg/dL)' , ylab = 'Probability for Hypoglycemic Event', main = 'Timepoint 12',
)


# GENERATES CONFUSION MATRIX FOR A PREDICTION THRESHOLD OF 0.05
confusion_matrix <- array(integer(4 * Tm_Pts), dim = c(2,2,Tm_Pts))
for(i in 1:Tm_Pts){

  df <- output[output$Time_Point == i,]
  rownames(df) <- 1:nrow(df)
  for(m in 1:nrow(df)){
    
    ratio <- 0.05
    if(df[m,'Hypo_Prob'] >= ratio){
      pred_state <- 1
      
    }
    else{
      pred_state <- df[m,'Predicted_State']
    }
    if(is.na(df[m,'State'])){
      
    }
    else if(pred_state == df[m,'State'] & pred_state == 1){
      confusion_matrix[1,1,i] <- confusion_matrix[1,1,i] + 1
    }
    else if(pred_state != 1 & df[m,'State'] != 1){
      confusion_matrix[2,2,i] <- confusion_matrix[2,2,i] + 1
    }
    else if(df[m,'State'] == 1){
      confusion_matrix[2,1,i] <- confusion_matrix[2,1,i] + 1
    }
    else{
      confusion_matrix[1,2,i] <- confusion_matrix[1,2,i] + 1
    }
  }
}


# ALIGNING INITIAL DATA TO SHOW 3 STATES INSTEAD OF 5
for(i in 1:nrow(output)){
  if(is.na(output[i,'State'])){
    
  }
  else if(output[i, 'State'] == 1){
    output[i, 'Actual_State_Com'] = 1
  }
  else if(output[i, 'State'] == 5){
    output[i, 'Actual_State_Com'] = 3
  }
  else{
    output[i, 'Actual_State_Com'] = 2
  }
  if(is.na(output[i,'Hypo_Prob'])){
    
  }
  else if(output[i, 'Hypo_Prob'] >= 0.05){
    output[i, 'Pred_State_Com'] = 1
  }
  else{
    output[i, 'Pred_State_Com'] = output[i, 'Predicted_State']
  }
}


# PLOTTING INITIALIZED AND THE PREDICTED STATE
df <- output[output[, 'Time_Point'] == 1 & output[,'Patient_Number'] == 1,]
par(mfrow=c(2, 1))
plot(x = 1:nrow(df), y = df[, 'Actual_State_Com'], main = 'Actual State')
plot(x = 1:nrow(df), y = df[, 'Pred_State_Com'], main = 'Predicted State')



write.csv(output,'LOCATION/FILE_NAME.csv')
write.csv(h.activity$model$sojourn$d,'LOCATION/FILE_NAME.csv')
write.csv(TPR,'LOCATION/FILE_NAME.csv')
write.csv(FPR,'LOCATION/FILE_NAME.csv')


# GENERATING STATE SEQUENCE USING THE IN-BUILD PREDICT FUNCTION
df_StateSequence <- data.frame()
for(i in 1:20){
  InitializedState <- dfs[[i]][, 'State'][[1]]
  CGM <- dfs[[i]][, 'CGM'][[1]]
  PatientNumber <- rep(i, length(InitializedState))
  PredictedState <- predict(h.activity, CGM)$s
  df <- data.frame(PatientNumber, CGM, InitializedState, PredictedState)
  df_StateSequence <- rbind(df_StateSequence, df)
}

write.csv(df_StateSequence,'LOCATION/FILE_NAME.csv')
write.csv(h.activity$model$d,'LOCATION/FILE_NAME.csv')
write.csv(h.activity$model$D,'LOCATION/FILE_NAME.csv')


# TRAINED EMISSION DISTRIBUTION SHAPE AND SCALE PARAMETERS
print(h.activity$model$parms.emission)
# TRAINED TRANSITION MATRIX
print(h.activity$model$transition)
# TRAINED SOJOURN DISTRIBUTION SHAPE AND SCALE PARAMETERS
print(h.activity$model$sojourn$shape)
print(h.activity$model$sojourn$scale)
# LATENT STATE INITIALIZATION QUANTILES
print(Quantiles)
# MINIMUM GLUCOSE VALUE IN THE DATASET
print(min(df_P1$CGM))
# MAXIMUM GLUCOSE VALUE IN THE DATASET
print(max(df_P1$CGM))
# TRAINED INITIAL DISTRIBUTION
print(h.activity$model$init)



# INITIALIZED EMISSION DISTRIBUTION MEAN AND STANDARD DEVIATION VALUES
print(Mean_Emis)
print(sqrt(Var_Emis))
# INITIALIZED SOJOURN DISTRIBUTION MEAN AND STANDARD DEVIATION VALUES
print(Mean_Sojourn)
print(sqrt(Var_Sojourn))
# INITIALIZED TRANSITION MATRIX
print(trans)


# SAVING OUTPUT WORKSPACE
save.image('LOCATION/FILE_NAME.RData')
