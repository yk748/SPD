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
# Running:
##################################################################
load("subject_id_final.RData")

for (subj_no in 1:6){
  
  # ------------------------------------------------------------- #
  # Set up for the experiment:
  load(paste0(subject_id[subj_no],".RData"))
  t_star <- subject_t_star[[subj_no]]
  TT <- dim(data)[1] # Number of samples:
  
  for (trial in 1:2){
    
    output <- list()
    iter <- 1; cnt <- 1
    cnt_max <- TT
    window_size <- 5
    
    if (trial == 1){
      step_size <- 2 # w = 2 in the manuscript
    }else{
      step_size <- 5 # w = 5 in the manuscript
    }
    cat(paste0("#--------------------------------------------------# \n"))
    cat(paste0("Step size ",step_size,".\n"))
    
    # ------------------------------------------------------------- #
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
    save(output,step_size,window_size,data,
         file=paste0("result_step",step_size,"_",subject_id[subj_no],".RData"))
  }
}

##################################################################
# Result:
##################################################################
load("subject_id_final.RData")

for (subj_no in 1:6){
  
  for (trial in 1:2){
    
    # ------------------------------------------------------------- #
    # load the result file:
    if (trial == 1){
      step_size <- 2 # w = 2 in the manuscript
    }else{
      step_size <- 5 # w = 5 in the manuscript
    }
    
    file_name <- paste0("result_step",step_size,"_",subject_id[subj_no],".RData")
    load(file_name)
    
    t_star <- subject_t_star[[subj_no]]
    metric_output <- SPD_metric(output,step_size,t_star)
    
    # ------------------------------------------------------------- #
    # Print the output:
    cat(paste0("#--------------------------------------------------# \n"))
    cat(paste0("Subject ",subject_id[subj_no],"\n"))
    cat(paste0("Step size ",step_size,"\n"))
    cat(paste0("Number of iteration is ",length(output),"\n"))
    print(round(metric_output$table,2))
    cat(paste0("#--------------------------------------------------# \n"))
  }
}

