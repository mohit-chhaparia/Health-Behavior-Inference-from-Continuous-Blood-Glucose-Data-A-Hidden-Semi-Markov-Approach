library(readxl)
library(writexl)
library(mhsmm)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(lubridate)
library(dipsaus)


setwd("SET_WORKING_DIRECTORY")
source("Training and Prediction Functions.R")
source("Initialization Functions.R")
setwd("SET_WORKING_DIRECTORY")


patient_number <- 1
file_path_pattern <- "FILE_NAME-%d.csv" # Creates generic filename pattern
dfs=list(); # Creates an empty list to hold datasets
file_name <- sprintf(file_path_pattern, 1000 + patient_number) # Creates the filename
df <- read.csv(file_name) # Reads the file
dfs[[1]] <- df # Assigns it to the datasets list
rm(df) # Deletes extra copy of the dataset


s_df <- 288 # Number of windows or iterations
dbs <- 1 # Number of datasets
Tm_Pts <- 12 # Timepoints to be predicted






# Quantiles to initialize all the required variables using first 1 day CGM data (288 records) from all the datasets
Quantiles <- quantile(unlist(lapply(dfs, function(x) x$CGM[1:288])), 
                      probs = c(0.15, 0.375, 0.625, 0.85, 1.00), na.rm = TRUE)

# Assign initial states for all the datasets based on Quantiles of first 1 day CGM data from all the datasets
dfs <- initializeStates(dfs, Quantiles, dbs)


############################################################
## Dataframe structure initialization to store the output ##
############################################################
# 1 DataFrame * 1 Predictions * (Tm_Pts + 1) TimePoints * s_df Subsets of each DataFrame

# Records time points from 0 to 12 (or 0 to 60 minutes)
Time_Point <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of state values from 1 to 5 as per the initialization using the quantiles approach
State <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of observed CGM (Glucose in mg/dL) values
CGM <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of predicted states from 1 to 3 without using any threshold (1 represents
# Hypoglycemic state; 2 represents Normal state; 3 represents Hyperglycemic state)
Predicted_State <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Predicted probability of hypoglycemic state using the Monte Carlo approach
Hypo_Prob <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Patient number from 1 to 20
Patient_Number <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of observed state from 1 to 3 (1 represents hypoglycemic state with Glucose
# values less than 70 mg/dL; 2 represents normal state with Glucose values from 70 mg/dL
# to less than 180 mg/dL; 3 represents hyperglycemic state with Glucose values of
# 180 mg/dL and above).
Actual_State_Com <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of predicted states from 1 to 3 using a threshold of 0.05 (1 represents
# Hypoglycemic state; 2 represents Normal state; 3 represents Hyperglycemic state)
Pred_State_Com <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of observed hypoglycemic event within 30 minutes with values 0, 1, or NA
# 1 represents atleast 1 observed glucose value exists which is less than 70 mg/dL within the 30 minutes interval
# 0 represents no observed glucose value exists which is less than 70 mg/dL within the 30 minutes interval
# NA represents the following cases
# Case 1: Value for 0 to 30 minutes interval is stored at time point 1, value for 
# 5 to 35 minutes interval is stored at time point 2 and so on. Hence, no value is
# stored at time points 0, 8, 9, 10, 11, and 12
# Case 2: All observed Glucose values are missing from the required window
Actual_Hypo_Event_30min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Predicted probability of an hypoglycemic event within the next 30 minutes
# Value of time point 1 represents the probability of hypoglycemic event from 0 to 30 minutes,
# value of time point 2 represents the probability of hypoglycemic event from 5 to 35 minutes,
# and so on. Value is NA for time points 0, 8, 9, 10, 11, and 12.
Prob_of_Hypo_30min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Record of observed hypoglycemic event within 60 minutes with values 0, 1, or NA
# 1 represents atleast 1 observed glucose value exists which is less than 70 mg/dL within the 60 minutes interval
# 0 represents no observed glucose value exists which is less than 70 mg/dL within the 60 minutes interval
# NA represents the following cases
# Case 1: Value for 0 to 60 minutes interval is stored at time point 1. Hence, the value
# at all the remaining time points is NA.
# Case 2: All observed Glucose values are missing from the required window
Actual_Hypo_Event_60min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Predicted probability of an hypoglycemic event within the next 60 minutes
# Value of time point 1 represents the probability of hypoglycemic event from 0 to 60 minutes.
# Value is NA for all the remaining time points.
Prob_of_Hypo_60min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Date and time for the observed Glucose value.
Time_Stamp <- as.character(integer(dbs * 1 * s_df * (Tm_Pts + 1)))

# Minimum observed Glucose value within a 30 minute window
# A value of NA represents the following cases:
# Case 1: All the observed glucose value within the required 30 minute window are NAs.
# Case 2: Value for 0 to 30 minutes interval is stored at time point 1,
# value for 5 to 35 minutes interval is stored at time point 2, and so on.
# Hence, value is NA for time points 0, 8, 9, 10, 11, and 12.
Min_CGM_30min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

# Minimum observed Glucose value within a 60 minute window
# A value of NA represents the following cases:
# Case 1: All the observed glucose value within the required 60 minute window are NAs.
# Case 2: Value for 0 to 60 minutes interval is stored at time point 1.
# Values for all the remaining time points are NAs.
Min_CGM_60min <- integer(dbs * 1 * s_df * (Tm_Pts + 1))

output <- data.frame(
  Patient_Number,
  Time_Point,
  State,
  CGM,
  Predicted_State,
  Hypo_Prob,
  Actual_State_Com,
  Pred_State_Com,
  Actual_Hypo_Event_30min,
  Prob_of_Hypo_30min,
  Actual_Hypo_Event_60min,
  Prob_of_Hypo_60min,
  Time_Stamp,
  Min_CGM_30min,
  Min_CGM_60min
)
output_index <- 1


############################################################
## Initializing parameters for the model training process ##
############################################################

J <- 5 # Number of States
init <- c(0.15, 0.225, 0.25, 0.225, 0.15) # Initial Distribution (Same as the quantiles used in the initialization process)
df_Sojourn <- getSoj(lapply(dfs, function(x) x$State[1:288]), dbs) # A dataframe with entering & leaving states, time spend in the state, and patient number.
trans <- getTrans(df_Sojourn) # Initialize the transition matrix
emis <- initializeDist(do.call(rbind,lapply(dfs, function(x) x[1:288,])), 'CGM', 'State') # Initialize the emission distribution
soj <- initializeDist(df_Sojourn, 'time_spend', 'coming_from') # Initialize the sojourn distribution
M <- 144 # Maximum time spend in a state at a time (12 hours or 144 records)
d <- do.call(cbind, lapply(1:J, function(m) dgamma(1:M, shape = soj$shape[m], scale = soj$scale[m])))

rm(df_Sojourn)



###################################
## Model training and prediction ##
###################################
indexes <- getIndexes(dfs, s_df, dbs) # Get indexes for the moving window approach
# dimension of indexes is 3 x 288 x 1

for(i in 1:s_df){
  print(sprintf("Window Number: %i", i))
  
  # df_train is the training data set for ith iteration
  df_train <- do.call(rbind, lapply(1:length(dfs), function(m) dfs[[m]][indexes[1,i,m]:indexes[2,i,m],]))
  
  # N is a vector indicating how many records are present from the patient's data set
  # in the training data set (df_train).
  N <- as.vector(table(df_train$`Patient.Number`))
  train <- list(x = df_train$CGM, N = N)
  class(train) <- "hsmm.data"
  
  startval <- hsmmspec(init, trans, emis, list(d = d, type = "gamma"), dens.emission = dgamma.hsmm, rand.emission = rgamma.hsmm, mstep = mstep.gamma)
  h.activity <- hsmmfit(train, startval, mstep = mstep.gamma, maxit = 100, M = M, lock.transition = FALSE, lock.d = FALSE, graphical = TRUE)
  
  for(j in 1:dbs){
    if(indexes[1,i,j]==0){
      print('Error in Indexes!!!!')
    }
    else{
      test <- list(x = dfs[[j]][indexes[1,i,j]:indexes[2,i,j],]$CGM, N = indexes[2,i,j] - indexes[1,i,j] + 1)
      class(test) <- "hsmm.data"
      
      test_pred <- predict_new(h.activity, test, predict_at =  (5 * Tm_Pts))
      for(n in 0:Tm_Pts){
        output[output_index, 'Patient_Number'] <- j
        output[output_index,'State'] <- (dfs[[j]])[indexes[2,i,j]+n,]$State
        output[output_index,'CGM'] <- (dfs[[j]])[indexes[2,i,j]+n,]$CGM
        output[output_index, 'Time_Point'] <- n
        output[output_index, 'Time_Stamp'] <- dfs[[j]][indexes[2,i,j]+n,]$`Timestamp..YYYY.MM.DDThh.mm.ss.`
        output[output_index, 'Predicted_State'] <- ifelse(n>0, which.max(test_pred[[5]][,n]), NA)
        output[output_index, 'Hypo_Prob'] <- ifelse(n>0, test_pred[[5]][1,n], NA)
        output[output_index, 'Actual_Hypo_Event_30min'] <- ifelse(n<8 & n>0,ifelse(any((dfs[[j]])[(indexes[2,i,j]+n):(indexes[2,i,j]+n+5),]$CGM < 70), 1, 0),NA)
        output[output_index, 'Prob_of_Hypo_30min'] <- ifelse(n>0, test_pred[[7]][n], NA)
        output[output_index, 'Actual_Hypo_Event_60min'] <- ifelse(n<2 & n>0,ifelse(any((dfs[[j]])[(indexes[2,i,j]+n):(indexes[2,i,j]+n+11),]$CGM < 70), 1, 0),NA)
        output[output_index, 'Prob_of_Hypo_60min'] <- ifelse(n>0, test_pred[[8]][n], NA)
        output[output_index, 'Min_CGM_30min'] <- ifelse(n<8 & n>0, min((dfs[[j]])[(indexes[2,i,j]+n):(indexes[2,i,j]+n+5),]$CGM, na.rm = TRUE), NA)
        output[output_index, 'Min_CGM_60min'] <- ifelse(n==1, min((dfs[[j]])[(indexes[2,i,j]+n):(indexes[2,i,j]+n+11),]$CGM, na.rm = TRUE), NA)
        output_index <- output_index + 1
      }
    }
  }
}


output$Actual_State_Com <- ifelse(is.na(output$CGM), NA, ifelse(output$CGM < 70, 1, ifelse(output$CGM >= 180, 3, 2)))
output$Pred_State_Com <- ifelse(is.na(output$Hypo_Prob), NA, ifelse(output$Hypo_Prob >= 0.05, 1, output$Predicted_State))


# To generate state sequence using the in-built predict function and the trained model
df_StateSequence <- data.frame()
for(i in 1:dbs){
  InitializedState <- dfs[[i]][, 'State'][[1]]
  CGM <- dfs[[i]][, 'CGM'][[1]]
  PatientNumber <- rep(i, length(InitializedState))
  PredictedState <- predict(h.activity, CGM)$s
  df <- data.frame(PatientNumber, CGM, InitializedState, PredictedState)
  df_StateSequence <- rbind(df_StateSequence, df)
}

setwd("SET_WORKING_DIRECTORY")
write.csv(output,'FILE_NAME.csv')
write.csv(h.activity$model$sojourn$d,'FILE_NAME.csv')
write.csv(df_StateSequence,'FILE_NAME.csv')
write.csv(h.activity$model$d,'FILE_NAME.csv')
write.csv(h.activity$model$D,'FILE_NAME.csv')
save.image('FILE_NAME.RData')