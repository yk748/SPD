library("ecp")
library ("e1071")
library("Matrix")
library("nnet")
library ("tree")
library("randomForest")
library("gbm")
library("glmnet")
library("ggplot2")
library("tidyr")
library("RColorBrewer")
library("patchwork")

source("library.R")

##################################################################
# Data illustration:
##################################################################
# Choose the subject:
load("subject7.Rdata") # Choose subject 7
# load("subject10.Rdata") # Choose subject 10

# Number of samples:
TT <- dim(data)[1]

##################################################################
# Running:
##################################################################
# Set up for the experiment:
output <- list()

iter <- 1; cnt <- 1
cnt_max <- TT
window_size <- 5
# step_size <- 2 # c = 2 in the manuscript
step_size <- 5 # c = 5 in the manuscript

# Main loop:
while(cnt <= cnt_max & iter <= 50){
  
  # ------------------------------------------------------------- #
  # Indicate the iteration (the number of trials in total) and count
  # ------------------------------------------------------------- #
  cat(paste0("#--------------------------------------------------# \n"))
  cat(paste0("Iteration ",iter,".\n"))
  cat(paste0("Count ",cnt,".\n"))
  cat(paste0("#--------------------------------------------------# \n"))
  hold_out <- c((5+cnt):(5+cnt+window_size-1))
  
  # Condition check for hold_out
  if (hold_out[length(hold_out)] > TT){
    cat(paste0("Break; The hold out samples exceed the last observation.\n"))
    break
  }else if(sum(is.na(match(t_star,hold_out))) < length(t_star)){
    cat(paste0("Break; The hold out samples contain the finalized change points.\n"))
    cnt <- cnt + step_size
    next
  }
  
  # ------------------------------------------------------------- #
  # Data preparation
  # ------------------------------------------------------------- #
  X_name <- c("active_time", "step_count",
              "conversation_percent","tic_voiced_time",
              "time_at_home", "sleep_duration",
              "sleep_interruption",
              "total_activity_duation","total_location_duration",
              "travel_diameter","radius_of_gyration")
  Y_name <- c("stress","arousal","valence")
  
  X_full <- data.frame(data$active_time, 
                       data$step_count,
                       data$conversation_percent,
                       data$tic_voiced_time,
                       data$time_at_home,
                       data$sleep_duration,
                       data$sleep_interruption,
                       data$travel_diameter,
                       data$total_activity_duration,
                       data$total_location_duration,
                       data$radius_of_gyration)
  colnames(X_full) <- X_name
  
  Y_full <- data.frame(data$stress,
                       data$arousal,
                       data$valence)
  colnames(Y_full) <- Y_name
  
  X_train <- X_full[-c(hold_out),]
  X_test <- X_full[c(hold_out),]
  Y_train <- Y_full[-c(hold_out),]
  Y_test <- Y_full[c(hold_out),]
  
  # ------------------------------------------------------------- #
  # Construct adjusted responses
  # ------------------------------------------------------------- #
  Y_adj <- response_adjustment(Y_train,Y_test)$series
  Y_train$adj_stress <- Y_adj$train
  Y_test$adj_stress <- Y_adj$test
  Y_name <- c(Y_name,"adj_stress")
  
  # ------------------------------------------------------------- #
  # Algorithm:
  # ------------------------------------------------------------- #
  cat(paste0("#--------------------------------------------------# \n"))
  cat(paste0("Algorithm begins: \n"))
  cat(paste0("#--------------------------------------------------# \n"))
  
  # Step 1. CPD
  Step1 <- SPD_step1(X_train,Y_train,X_name,Y_name)
  
  if (length(Step1$t_hat) == 0){
    cat("Terminated because there was no change point.\n")
    cnt <- cnt + step_size
    next
  }else{
    cat(paste0("Detected change points are ",
               length(Step1$t_hat),".\n"))
  }
  
  # Step 2: Validation of significance
  Step2 <- SPD_step2(Step1)
  
  if (length(Step2$t_hat) == 0){
    cat("Terminated because all segments are not significant.\n")
    cnt <- cnt + step_size
    next
  }else{
    for (k in 1:(length(Step1$t_hat)+1)){
      cat(paste0("p-value for the ",k,"th segment is ",
                 round(Step2$p_val[k],7),", ",Step2$decision[k],".\n"))
    }
  }
  
  # Step 3: Validation of distinction in equations
  Step3 <- SPD_step3(Step2)
  
  if (length(Step3$t_hat) == 0){
    cat("Terminated because all segments are not distinct.\n")
    cnt <- cnt + step_size
    next
  }else{
    for (k in 1:length(Step2$t_hat)){
      if ( Step3$decision[k] == "skip"){
        cat(paste0("Skip ",k," and ",k+1,
                   "th segment since there is no difference in equations.\n"))
      }else{
        cat(paste0("p-value for the ",k," and ",k+1,"th segment is "
                   ,round(Step3$p_val[k],7),", ",Step3$decision[k],".\n"))
      }
    }
  }
  
  # Step 4: Validation of distinction in responses
  Step4 <- SPD_step4(Step3)
  
  if (length(Step4$t_hat) == 0){
    cat("Terminated because all responses are not distinct.\n")
    cnt <- cnt + step_size
    next
  }else{
    for (k in 1:length(Step3$t_hat)){
      if ( Step4$decision[k] == "skip"){
        cat(paste0("Skip ",k," and ",k+1,
                   "th segment since there is no difference in distributions.\n"))
      }else{
        cat(paste0("p-value for the ",k," and ",k+1,"th segment is "
                   ,round(Step4$p_val[k],7),", ",Step4$decision[k],".\n"))
      }
    }
  }
  
  # Step 5: Prediction step
  Final_segment <- SPD_stress_ind(Step4)
  Step5 <- SPD_step5(Final_segment,Y_train,X_train,Y_test,X_test,X_name,window_size)
  
  # ------------------------------------------------------------- #
  # Store:
  # ------------------------------------------------------------- #
  output[[iter]] <- list(hold_out = hold_out,
                         result = Final_segment,
                         classification = Step5)
  
  iter <- iter + 1
  cnt <- cnt + step_size
}

# save(output,step_size,window_size,file="subj7_c2.RData") # If subject 7 with c=2
# save(output,step_size,window_size,file="subj7_c5.RData") # If subject 7 with c=5

# save(output,step_size,window_size,file="subj10_c2.RData") # If subject 10 with c=2 
# save(output,step_size,window_size,file="subj10_c5.RData") # If subject 10 with c=5


##################################################################
# Result:
##################################################################
# load("data_SPD_final_subj7_step2.RData"); t_star <- c(31,92)
# load("data_SPD_final_subj7_step5.RData"); t_star <- c(31,92)

# load("data_SPD_final_subj10_step2.RData"); t_star <- 69
# load("data_SPD_final_subj10_step5.RData"); t_star <- 69

metric_output <- SPD_metric(output,iter,step_size,t_star)
round(metric_output$table,2)

