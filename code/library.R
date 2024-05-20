# ------------------------------------------------------------- #
# response_adjustment
# ------------------------------------------------------------- #
response_adjustment <- function(Y_train,Y_test){
  
  mean_arousal <- mean(Y_train$arousal)
  mean_valence <- mean(Y_train$valence)
  
  TT <- dim(Y_train)[1]
  Y_adj_train <- vector("numeric",TT)
  t_train_pls <- vector("numeric",length=TT)
  t_train_mns <- vector("numeric",length=TT)
  for (t in 1:TT){
    adj_score <- pam_grid(mean_arousal,mean_valence,
                          Y_train$arousal[t],
                          Y_train$valence[t])
    Y_adj_train[t] <- Y_train$stress[t] + adj_score
    if (adj_score > 0){
      t_train_pls[t] <- adj_score
    }else if (adj_score < 0){
      t_train_mns[t] <- adj_score
    }else{
      next
    }
  }
  
  tt <- dim(Y_test)[1]
  Y_adj_test <- vector("numeric",tt)
  t_test_pls <- vector("numeric",length=tt)
  t_test_mns <- vector("numeric",length=tt)
  for (t in 1:tt){
    adj_score <- pam_grid(mean_arousal,mean_valence,
                          Y_test$arousal[t],
                          Y_test$valence[t])
    Y_adj_test[t] <- Y_test$stress[t] + adj_score
    if (adj_score > 0){
      t_test_pls[t] <- adj_score
    }else if (adj_score < 0){
      t_test_mns[t] <- adj_score
    }else{
      next
    }
  }
  
  output <- list()
  output$series <- list(train=Y_adj_train,test=Y_adj_test)
  output$time <- list(train_add = t_train_pls,
                      train_sub = t_train_mns,
                      test_add = t_test_pls,
                      test_sub = t_test_mns)
  return(output)
}

# ------------------------------------------------------------- #
pam_grid <- function(mean_arousal,mean_valence,
                     arousal_t,valence_t){
  
  if (arousal_t > mean_arousal & 
      valence_t < mean_valence ){
    
    if ( arousal_t > ((mean_arousal+5)/2) &
         valence_t > ((mean_valence+1)/2) ){
      adj <- 1 # grid 8
    }else if( arousal_t < ((mean_arousal+5)/2) &
              valence_t > ((mean_valence+1)/2) ){
      adj <- 0.75 # grid 7
    }else if( arousal_t > ((mean_arousal+5)/2) &
              valence_t < ((mean_valence+1)/2) ){
      adj <- 0.5 # grid 6
    }else{
      adj <- 0.25 # grid 5
    }
    
  }else if(arousal_t < mean_arousal &
           valence_t > mean_valence ){
    
    if ( arousal_t > ((mean_arousal+1)/2) &
         valence_t > ((mean_valence+5)/2) ){
      adj <- -1 # grid 12
    }else if( arousal_t < ((mean_arousal+1)/2) &
              valence_t > ((mean_valence+5)/2) ){
      adj <- -0.75 # grid 11
    }else if( arousal_t > ((mean_arousal+1)/2) &
              valence_t < ((mean_valence+5)/2) ){
      adj <- -0.5 # grid 10
    }else{
      adj <- -0.25 # grid 9
    }
    
  }else{
    adj <- 0
  }
  return(adj)
}


# ------------------------------------------------------------- #
# SPD_step1
# ------------------------------------------------------------- #
SPD_step1 <- function(X_train,Y_train,X_name,Y_name){
  
  # Use Matteson and James (2014, JASA)
  output_1st_mmt <- e.divisive(X_train, sig.lvl = 0.05, R = 199, alpha = 1)
  output_2nd_mmt <- e.divisive(X_train, sig.lvl = 0.05, R = 199, alpha = 2)
  output_seg <- sort(unique(c(output_1st_mmt$estimates,
                              output_2nd_mmt$estimates)))
  
  # Construct the segments and cosegments
  if (length(output_seg[-c(1,length(output_seg))]) == 0){
    
    data_seg <- list()
    t_hat <- NULL
    
  }else{
    
    win_start <- output_seg[-length(output_seg)]
    win_end <- c(output_seg[-c(1,length(output_seg))]-1,
                 output_seg[length(output_seg)]-1)
    data_seg <- list()
    for (k in 1:length(win_start)){
      data_seg[[k]] <- data.frame(time = c(win_start[k]:win_end[k]),
                                  Y_train[c(win_start[k]:win_end[k]),],
                                  X_train[c(win_start[k]:win_end[k]),])
      colnames(data_seg[[k]]) <- c("time",Y_name,X_name)
    }
    
    t_hat <- output_seg[-c(1,length(output_seg))]
    
  }
  
  output <- list()
  output$data_seg <- data_seg
  output$t_hat <- t_hat
  output$X_name <- X_name
  return(output)
}



# ------------------------------------------------------------- #
# SPD_step2
# ------------------------------------------------------------- #
SPD_step2 <- function(Step1){
  
  num_seg <- length(Step1$t_hat)+1
  data_seg <- Step1$data_seg
  alpha <- 0.05
  X_name <- Step1$X_name
  
  # Step 2: Test of significance
  sig_idx <- NULL; insig_idx <- NULL; TT_k <- vector("numeric",num_seg)
  decision <- vector("numeric",num_seg)
  current_seg <- NULL
  pval <- vector("numeric",num_seg)
  
  for (k in 1:num_seg){
    
    TT_k[k] <- dim(data_seg[[k]])[1]
    
    if ( sum(diff(data_seg[[k]]$adj_stress) == 0) == (TT_k[k]-1)){ 
      # if adjusted stress has no fluctuation,
      decision[k] <- "No fluctuation, merged"
      pval[k] <- NA
      
    }else if( TT_k[k] <= length(X_name)){ 
      # If sample length is not enough (less than equal to 11),
      decision[k] <- "Sample size shortage, merged"
      pval[k] <- NA
      
    }
    
    if (decision[k]!=0){
      current_seg <- rbind(current_seg,data_seg[[k]])
      next
      
    }else{
      current_seg <- data_seg[[k]]
      
    }
    
    reg <- lm(adj_stress~.,current_seg[,c("adj_stress",X_name)])
    null <- lm(adj_stress~1,current_seg[,c("adj_stress",X_name)])
    anova <- anova(reg,null)
    pval[k] <- anova$`Pr(>F)`[2]
    
    if(is.nan(pval[k])){
      decision[k] <- "NaN produced, merged"
      pval[k] <- NA
      current_seg <- rbind(current_seg,data_seg[[k]])
      next
      
    }else{
      if (pval[k] < alpha){ # If pval < alpha
        sig_idx <- c(sig_idx,k)
        decision[k] <- "significant"
        
      }else{ # If pval >= alpha
        insig_idx <- c(insig_idx,k)
        decision[k] <- "insignificant"
        
      }
    }
  }
  
  # Rearragement after step 2:
  if (sum(decision=="significant") == 0){
    sig_data_seg <- list()
    sig_seg_start <- NULL
    
  }else{
    TT <- sum(TT_k)
    # Based on Step 1 result:
    win_start <- c(1,Step1$t_hat); win_end <- c(Step1$t_hat-1,TT)
    # Refinement:
    sig_data_seg <- list(); 
    sig_seg_start <- sig_seg_end <- NULL
    cnt_sig <- 0; cumm_idx <- NULL; 
    effective_idx <- c(1:length(win_start))[!is.na(match(decision,c("insignificant","significant")))]
    
    for (k in 1:length(win_start)){
      
      if (is.na(match(k,effective_idx)!=0)){
        cumm_idx <- c(cumm_idx,k)
        if (k == length(win_start)){
          
        }else{
          next 
        }
        
      }else{
        if (is.na(match(k,sig_idx))){ # if this segment is not significant,
          cumm_idx <- c(cumm_idx,k)
          
          cnt_sig <- cnt_sig + 1
          sig_data_seg[[cnt_sig]] <- do.call("rbind",lapply(X=cumm_idx,
                                                            FUN=function(X)data_seg[[X]]))
          
          sig_seg_start <- c(sig_seg_start,win_start[min(cumm_idx)])
          sig_seg_end <- c(sig_seg_end,win_end[max(cumm_idx)])
          cumm_idx <- NULL
          
        }else{ # if this segment is significant,
          cnt_sig <- cnt_sig + 1
          sig_data_seg[[cnt_sig]] <- data_seg[[k]]
          sig_seg_start <- c(sig_seg_start,win_start[k])
          sig_seg_end <- c(sig_seg_end,win_end[k])
          
          # If cummulative segments that have to be merged exist,
          if (sum(cumm_idx)>0){
            cnt_sig <- cnt_sig + 1
            sig_data_seg[[cnt_sig]] <- do.call("rbind",lapply(X=cumm_idx,
                                                              FUN=function(X)data_seg[[X]]))
            
            sig_seg_start <- c(sig_seg_start,win_start[min(cumm_idx)])
            sig_seg_end <- c(sig_seg_end,win_end[max(cumm_idx)])
            cumm_idx <- NULL
          }
        }
      }
    } 
  }
  
  output <- list()
  output$data_seg <- sig_data_seg
  output$t_hat <- sig_seg_start[-1]
  output$p_val <- pval
  output$decision <- decision
  output$X_name <- X_name
  return(output)
}


# ------------------------------------------------------------- #
# SPD_step3
# ------------------------------------------------------------- #
SPD_step3 <- function(Step2){
  
  num_seg <- length(Step2$t_hat)+1
  data_seg <- Step2$data_seg
  alpha <- 0.05
  X_name <- Step2$X_name
  
  dis_idx <- list()
  indis_idx <- list()
  F_stat <- vector("numeric",(num_seg-1))
  pval <- vector("numeric",(num_seg-1))
  decision <- vector("numeric",(num_seg-1))
  TT <- dim(data_seg[[1]])[1]
  for (k in 1:(num_seg-1)){
    
    dat1 <- data_seg[[k]][,c(X_name,"adj_stress")]
    dat2 <- data_seg[[(k+1)]][,c(X_name,"adj_stress")]
    dat_c <- rbind(dat1,dat2)
    N1 <- dim(dat1)[1]; N2 <- dim(dat2)[1]
    num_param <- dim(dat_c)[2]
    
    TT <- TT + N2
    
    reg1 <- lm(adj_stress~.,dat1)
    S1 <- sum(reg1$residuals^2)
    
    reg2 <- lm(adj_stress~.,dat2)
    S2 <- sum(reg2$residuals^2)
    
    
    reg_c <- lm(adj_stress~.,dat_c)
    S_c <- sum(reg_c$residuals^2)
    
    df1 <- num_param
    df2 <- N1+N2-2*num_param
    F_stat[k] <- ( (S_c - (S1+S2))/df1 )/( (S1+S2)/df2 )
    pval[k] <- 1- pf(F_stat[k],df1,df2)
    
    if (pval[k] < alpha){
      dis_idx[[k]] <- c(k,k+1)
      indis_idx[[k]] <- 0
      decision[k] <- "distinct"
    }else{
      dis_idx[[k]] <- 0
      indis_idx[[k]] <- c(k,k+1)
      decision[k] <- "indistinct"
    }   
  }
  
  # Rearragement after step 3:
  if (sum(decision == "distinct") == 0){
    final_data_seg <- list()
    final_seg_start <- NULL
    
  }else{
    cnt_dis <- 1
    merge_idx <- vector("numeric",length=num_seg)
    for (k in 1:(num_seg-1)){
      if (k == 1){
        
        merge_idx[k] <- cnt_dis
        if ( (num_seg-1) == 1){
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
          break
        }
        
      }else if( k == (num_seg-1)){
        
        if ( length(indis_idx[[(k-1)]])==2 & length(indis_idx[[k]])==2 ){ 
          # if (k-1)th and kth segments are not distinct and kth and (k+1)th segments are not distinct either,
          merge_idx[k] <- cnt_dis
          merge_idx[(k+1)] <- cnt_dis
        }else if( length(indis_idx[[(k-1)]])==2 & length(indis_idx[[k]])==1 ){
          # if (k-1)th and kth segments are not distinct and kth and (k+1)th segments are distinct,
          merge_idx[k] <- cnt_dis
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
        }else if( length(indis_idx[[(k-1)]])==1 & length(indis_idx[[k]])==2 ){ 
          # if (k-1)th and kth segments are distinct but kth and (k+1)th segments are not distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
          merge_idx[(k+1)] <- cnt_dis
        }else{
          # if (k-1)th and kth segments are distinct and kth and also (k+1)th segments are distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
        }
        
      }else{
        
        if ( length(dis_idx[[(k-1)]]) == 2){ 
          # if (k-1)th and kth segments are distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
        }else if ( length(dis_idx[[(k-1)]]) == 1 ){ 
          # if (k-1)th and kth segments are not distinct,
          merge_idx[k] <- cnt_dis
        }
        
      }
    }
    
    final_data_seg <- list()
    merge_list <- split(seq_along(merge_idx), merge_idx)
    win_start <- c(1,Step2$t_hat); win_end <- c(Step2$t_hat-1,TT)
    final_seg_start <- NULL ; final_seg_end <- NULL
    for (k in 1:length(unique(merge_idx))){
      final_data_seg[[k]] <- do.call("rbind",lapply(X=merge_list[[k]],
                                                    FUN=function(X)data_seg[[X]]))
      final_seg_start[k] <- win_start[min(merge_list[[k]])]
      final_seg_end[k] <- win_end[max(merge_list[[k]])]
    }
  }
  
  output <- list()
  output$data_seg <- final_data_seg
  output$t_hat <- final_seg_start[-1]
  output$p_val <- pval
  output$decision <- decision
  # output$notice <- notice
  return(output)
}


# ------------------------------------------------------------- #
# SPD_step4
# ------------------------------------------------------------- #
SPD_step4 <- function(Step3){
  
  num_seg <- length(Step3$t_hat)+1
  data_seg <- Step3$data_seg
  alpha <- 0.05
  
  data_Y_table <- array(NA,dim=c(num_seg,4))
  for (k in 1:num_seg){
    for (j in c(2:5)){
      data_Y_table[k,(j-1)] <- sum(data_seg[[k]][,"stress"] == j)
    }
  }
  dimnames(data_Y_table) <- list(segment = as.character(1:num_seg),
                                 value = c("2","3","4","5"))
  TT <- sum(data_Y_table)
  
  dis_idx <- list()
  indis_idx <- list()
  pval <- vector("numeric",(num_seg-1))
  decision <- vector("numeric",(num_seg-1))
  
  for (k in 1:(num_seg-1)){
    
    focus_table <- data_Y_table[c(k,(k+1)),]
    x1 <- focus_table[1,which(colSums(focus_table) != 0)]
    x2 <- focus_table[2,which(colSums(focus_table) != 0)]
    x <- x1 + x2 
    
    if (all(x > 5) & sum(x == 0) < 2){
      
      # all x_i > 5 and # of non-zero categories is 0 or 1
      
      pval[k] <- chisq.test(cbind(x1,x2))$p.value
      
    }else if( any(x <= 5) & sum(x == 0) < 2){
      
      # some x_i <= 5 but # of non-zero categories is 0 or 1
      x1 <- x1[-which(x <= 5)]
      x2 <- x2[-which(x <= 5)]
      pval[k] <- chisq.test(cbind(x1,x2))$p.value
      
    }else if (sum(x == 0) == 2){
      
      # only 2 categories are non-zero.
      n1 <- sum(x1); n2 <- sum(x2)
      x1 <- x1[1]; x2 <- x2[1]
      pval[k] <- prop.test(x = c(x1,x2), n = c(n1,n2))$p.value
      
    }else{ 
      # only 1 category is non-zero
      pval[k] <- NULL
      
    }
    
    if (is.null(pval[k])){
      dis_idx[[k]] <- 0
      indis_idx[[k]] <- c(k,k+1)
      decision[k] <- "skip"
    }else{
      if (pval[k] < alpha){
        dis_idx[[k]] <- c(k,k+1)
        indis_idx[[k]] <- 0
        decision[k] <- "distinct"
      }else{
        dis_idx[[k]] <- 0
        indis_idx[[k]] <- c(k,k+1)
        decision[k] <- "indistinct"
      }   
    }
  }
  
  
  # Rearragement after step 4:
  if (sum(decision == "distinct") == 0){
    final_data_seg <- list()
    final_seg_start <- NULL
    
  }else{
    cnt_dis <- 1
    merge_idx <- vector("numeric",length=num_seg)
    for (k in 1:(num_seg-1)){
      if (k == 1){
        
        merge_idx[k] <- cnt_dis
        if ( (num_seg-1) == 1){
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
          break
        }
        
      }else if( k == (num_seg-1)){
        
        if ( length(indis_idx[[(k-1)]])==2 & length(indis_idx[[k]])==2 ){ 
          # if (k-1)th and kth segments are not distinct and kth and (k+1)th segments are not distinct either,
          merge_idx[k] <- cnt_dis
          merge_idx[(k+1)] <- cnt_dis
        }else if( length(indis_idx[[(k-1)]])==2 & length(indis_idx[[k]])==1 ){
          # if (k-1)th and kth segments are not distinct and kth and (k+1)th segments are distinct,
          merge_idx[k] <- cnt_dis
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
        }else if( length(indis_idx[[(k-1)]])==1 & length(indis_idx[[k]])==2 ){ 
          # if (k-1)th and kth segments are distinct but kth and (k+1)th segments are not distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
          merge_idx[(k+1)] <- cnt_dis
        }else{
          # if (k-1)th and kth segments are distinct and kth and also (k+1)th segments are distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
          cnt_dis <- cnt_dis + 1
          merge_idx[(k+1)] <- cnt_dis
        }
        
      }else{
        
        if ( length(dis_idx[[(k-1)]]) == 2){ 
          # if (k-1)th and kth segments are distinct,
          cnt_dis <- cnt_dis + 1
          merge_idx[k] <- cnt_dis
        }else if ( length(dis_idx[[(k-1)]]) == 1 ){ 
          # if (k-1)th and kth segments are not distinct,
          merge_idx[k] <- cnt_dis
        }
        
      }
    }
    
    final_data_seg <- list()
    merge_list <- split(seq_along(merge_idx), merge_idx)
    win_start <- c(1,Step3$t_hat); win_end <- c(Step3$t_hat-1,TT)
    final_seg_start <- NULL ; final_seg_end <- NULL
    for (k in 1:length(unique(merge_idx))){
      final_data_seg[[k]] <- do.call("rbind",lapply(X=merge_list[[k]],
                                                    FUN=function(X)data_seg[[X]]))
      final_seg_start[k] <- win_start[min(merge_list[[k]])]
      final_seg_end[k] <- win_end[max(merge_list[[k]])]
    } 
  }
  
  output <- list()
  output$data_seg <- final_data_seg
  output$t_hat <- final_seg_start[-1]
  output$p_val <- pval
  output$decision <- decision
  return(output)
}


# ------------------------------------------------------------- #
# SPD_stress_ind
# ------------------------------------------------------------- #
SPD_stress_ind <- function(Step4,dicho=FALSE){
  
  num_seg <- length(Step4$t_hat)+1
  data_seg <- Step4$data_seg
  
  data_label <- list()
  if (dicho == FALSE){
    for (k in 1:num_seg){
      data_label[[k]] <- data.frame(label = rep(k,dim(data_seg[[k]])[1]),data_seg[[k]])
    }
  }else{
    for (k in 1:num_seg){
      if (k == 1){
        data_label[[1]] <- data.frame(label = rep(1,dim(data_seg[[k]])[1]),data_seg[[k]])
      }else if (k == 2){
        data_label[[2]] <- data.frame(label = rep(2,dim(data_seg[[k]])[1]),data_seg[[k]])
      }else if (k > 1 & k %% 2 == 1){
        data_label[[1]] <- rbind(data_label[[1]],
                                 data.frame(label = rep(1,dim(data_seg[[k]])[1]),data_seg[[k]]))
      }else if (k > 1 & k %% 2 == 0){
        data_label[[2]] <- rbind(data_label[[2]],
                                 data.frame(label = rep(2,dim(data_seg[[k]])[1]),data_seg[[k]]))
      }
    } 
  }
  data_tot <- do.call("rbind",data_label)
  return(data_tot)
}



# ------------------------------------------------------------- #
# SPD_step5
# ------------------------------------------------------------- #
SPD_step5 <- function(Final_segment,Y_train,X_train,Y_test,X_test,X_name,window_size){
  
  num_labels <- length(unique(Final_segment$label))
  # ------------------------------------------------------------- #
  # Benchmarks:
  Resp <- Pred <- array(NA,dim=c(num_labels,window_size))
  
  cluster_Y <- kmeans(Y_train$adj_stress,centers=num_labels)
  dat_lm <- data.frame(Y=Y_train$adj_stress,X_train)
  lm_full <- lm(Y~.,dat_lm)
  Y_pred <- predict(lm_full,X_test)
  for (k in 1:num_labels){
    Y_k_bar <- mean(Y_train$adj_stress[which(cluster_Y$cluster == k)])
    for (t in 1:window_size){
      Resp[k,t] <- (Y_test$adj_stress[t] - Y_k_bar)^2
      Pred[k,t] <- (Y_pred[t] - Y_k_bar)^2 
    }
  }
  decision_resp <- mapply(x=1:window_size,function(x)which.min(Resp[,x]))
  cat(paste0("The decision from responses is ",decision_resp,"\n"))
  decision_pred <- mapply(x=1:window_size,function(x)which.min(Pred[,x]))
  cat(paste0("The decision from prediction is ",decision_pred,"\n"))
  
  # ------------------------------------------------------------- #
  # Step 5:
  # Logistic (multinomial) regression
  tryCatch(
    {
      if (num_labels == 2){
        logistic_data <- Final_segment[,c("label",X_name)]
        logistic_data$label[which(logistic_data$label==1)] <- 1
        logistic_data$label[which(logistic_data$label!=1)] <- 0
        
        reg_fit <- glm(label~., family = binomial, data=logistic_data)
        pred_reg <- predict(reg_fit,as.data.frame(X_test),type="link")
        decision_reg <- ifelse(pred_reg >= 0.5,1,2)
        decision_reg <- as.integer(decision_reg)
        # cat(paste0("The decision from logistic regression is ",decision_reg,".\n")) 
      }else{
        mn_data <- Final_segment[,c("label",X_name)]
        reg_fit <- nnet::multinom(label~.,mn_data)
        
        pred_reg <- predict(lm_fit,X_test,type="probs")
        decision_reg <- apply(pred_reg,1,function(x) which.max(x))
        decision_reg <- as.integer(decision_reg)
        cat(paste0("The decision from multinomial regression is ",decision_reg,".\n"))
      }
    }, error = function(msg){
      cat(paste0("Regression has not been completed.\n"))
      return(NA)}
  )
  
  # ------------------------------------------------------------- #
  # Standard svm
  tryCatch(
    {
      if (num_labels == 2){
        svm_data <- Final_segment[,c("label",X_name)]
        svm_data$label[which(svm_data$label==1)] <- 1
        svm_data$label[which(svm_data$label!=1)] <- -1
        
        tune.out <- tune(svm, label~.,data=svm_data,
                         kernel = "linear",
                         ranges = list(cost=10^(-2:2)))
        
        svmfit <- svm(label~., data=svm_data, kernel="linear",
                      cost=tune.out$best.parameters)
        pred_svm <- predict(svmfit,X_test)
        decision_svm <- ifelse(sign(pred_svm) > 0, 1,2)
        decision_svm <- as.integer(decision_svm)
        cat(paste0("The decision from SVM are ",decision_svm,".\n"))
      }else{
        pred_svm <- array(NA,dim=c(window_size,length(num_labels)))
        for (k in 1:num_labels){
          svm_data <- Final_segment[,c("label",X_name)]
          svm_data$label[which(svm_data$label==k)] <- 1
          svm_data$label[which(svm_data$label!=k)] <- -1
          
          tune.out <- tune(svm, label~.,data=svm_data,
                           kernel = "linear",
                           ranges = list(cost=10^(-2:2)))
          
          svmfit <- svm(label~., data=svm_data, kernel="linear",
                        cost=tune.out$best.parameters)
          pred_svm[,k] <- predict(svmfit,X_test)
        }
        decision_svm <- apply(pred_svm,1,function(x) which.max(x))
        decision_svm <- as.integer(decision_svm)
        cat(paste0("The decision from SVM is ",decision_svm,".\n"))
      }
    }, error = function(msg){
      cat(paste0("SVM has not been completed.\n"))
      return(NA)}
  )
  
  # ------------------------------------------------------------- #
  # LDA:
  tryCatch(
    {
      TT <- dim(Final_segment)[1]
      pi_hat <- mapply(x=1:num_labels,function(x)sum(Final_segment$label==x)/TT)
      mu_hat <- t(mapply(x=1:num_labels,function(x) 
        colMeans(Final_segment[which(Final_segment$label==x),X_name])))
      Cov_hat <- t(Final_segment[,X_name]) %*% as.matrix(Final_segment[,X_name])/TT
      dis_fn_LDA <- array(NA,dim=c(window_size,num_labels))
      for (t in 1:window_size){
        for (k in 1:num_labels){
          dis_fn_LDA[t,k] <- as.matrix(X_test[t,]) %*% solve(Cov_hat) %*% mu_hat[k,] - 0.5*t(mu_hat[k,]) %*% solve(Cov_hat) %*% mu_hat[k,] + log(pi_hat[k])
        }
      }
      decision_LDA <- apply(dis_fn_LDA,1,function(x) which.max(x))
      cat(paste0("The decision from LDA is ",decision_LDA,".\n"))
    }, error = function(msg){
      cat(paste0("LDA has not been completed.\n"))
      return(NA)}
  )
  
  # ------------------------------------------------------------- #
  # RF:
  tryCatch(
    {
      rf_data <- Final_segment[,c("label",X_name)]
      rf_fit <- randomForest(factor(label)~.,rf_data , mtry=floor(11/3), 
                             ntree=500, importance=TRUE)
      pred_rf <- predict(rf_fit,X_test)
      decision_rf <- as.numeric(pred_rf)
      cat(paste0("The decision from random forest is ",decision_rf,".\n"))
    }, error = function(msg){
      cat(paste0("RF has not been completed.\n"))
      return(NA)}
  )
  
  
  # ------------------------------------------------------------- #
  # Boosting:
  tryCatch(
    {
      boost_data <- Final_segment[,c("label",X_name)]
      if (num_labels == 2){
        dist <- "bernoulli"
        boost_data$label[which(boost_data$label==1)] <- 1
        boost_data$label[which(boost_data$label!=1)] <- 0
      }else{
        dist <- "multinomial"
      }
      
      boost_fit <- gbm(label~.,data = boost_data, distribution = dist, cv.folds = 10,
                       shrinkage = 0.01, n.minobsinnode = 10, n.trees = 5000)
      pred_boost <- predict.gbm(object = boost_fit,
                                newdata = as.data.frame(X_test), 
                                n.trees = 5000, type = "response")
      
      if (num_labels == 2){
        decision_boost <- ifelse(pred_boost > 0.5,1,2)
      }else{
        decision_boost <- apply(pred_boost,1,function(x) which.max(x))
      }
      decision_boost <- as.integer(decision_boost)
      cat(paste0("The decision from boosting is ",decision_boost,"\n"))
    }, error = function(msg){
      cat(paste0("Boosting has not been completed.\n"))
      return(NA)}
  )
  
  # ------------------------------------------------------------- #
  if(!exists("decision_reg")){
    reg_fit <- NULL
    pred_glm <- NULL
    decision_reg <- rep(0,window_size)
  }
  if(!exists("decision_svm")){
    svmfit <- NULL
    pred_svm <- NULL
    decision_svm <- rep(0,window_size)
  }
  if(!exists("decision_LDA")){
    dis_fn_LDA <- NULL
    decision_LDA <- rep(0,window_size)
  }
  if(!exists("decision_rf")){
    rf_fit <- NULL
    decision_rf <- rep(0,window_size)
  }
  if(!exists("decision_boost")){
    boost_fit <- NULL
    decision_boost <- rep(0,(step_size+1))
  }
  
  # ------------------------------------------------------------- #
  output <- list()
  output$resp <- decision_resp
  output$pred <- decision_pred
  output$reg <- decision_reg
  output$svm <- decision_svm
  output$lda <- decision_LDA
  output$rf <- decision_rf
  output$boost <- decision_boost
  
  output$clust_spec <- cluster_Y
  output$lm_spec <- lm_full
  output$reg_spec <- reg_fit
  output$svm_spec <- svmfit
  output$lda_spec <- dis_fn_LDA
  output$rf_spec <- rf_fit
  output$boost_spec <- boost_fit
  return(output)
}


# ------------------------------------------------------------- #
# SPD_metric
# ------------------------------------------------------------- #
SPD_metric <- function(output,iter,step_size,t_star){
  
  num_seg <- length(t_star)+1
  num_iter <- length(output)
  true_labels <- array(NA,dim=c(num_iter,window_size))
  for (iter in 1:num_iter){
    hold_out <- output[[iter]]$hold_out
    
    for (h in 1:window_size){
      true_labels[iter,h] <- num_seg - sum(hold_out[h] < t_star)
    }
  }
  
  pred_labels <- array(NA,dim=c(num_iter,window_size,7))
  for (iter in 1:num_iter){
    pred_labels[iter,,1] <- output[[iter]]$classification$reg
    pred_labels[iter,,2] <- output[[iter]]$classification$svm
    pred_labels[iter,,3] <- output[[iter]]$classification$lda
    pred_labels[iter,,4] <- output[[iter]]$classification$rf
    pred_labels[iter,,5] <- output[[iter]]$classification$boost
    pred_labels[iter,,6] <- output[[iter]]$classification$resp
    pred_labels[iter,,7] <- output[[iter]]$classification$pred
  }
  
  rec_k <- array(NA,dim=c(num_seg,7))
  prec_k <- array(NA,dim=c(num_seg,7))
  f1_k <- array(NA,dim=c(num_seg,7))
  acc <- vector("numeric",7)
  for (l in 1:7){
    acc[l] <- sum(pred_labels[,,l] == true_labels)/(window_size*num_iter)
    for (k in 1:num_seg){
      rec_k[k,l] <-  sum(pred_labels[,,l] == k & true_labels == k)/sum(true_labels == k)
      prec_k[k,l] <-  sum(pred_labels[,,l] == k & true_labels == k)/sum(pred_labels[,,l] == k)
      f1_k[k,l] <- (2*rec_k[k,l]*prec_k[k,l])/(rec_k[k,l]+prec_k[k,l])
    }
  }
  rec_k[is.nan(rec_k)] <- 0
  prec_k[is.nan(prec_k)] <- 0
  f1_k[is.nan(f1_k)] <- 0
  rec <- colMeans(rec_k)
  prec <- colMeans(prec_k)
  
  table <- rbind(rec,prec,acc,f1_k)
  colnames(table) <- c("Reg","SVM","LDA","RF","Boost","Resp","Pred")
  name_method <- c("Recall","Precision","Accuracy")
  for (i in 1:num_seg){
    name_method <- c(name_method,paste0("F1_",i))
  }
  rownames(table) <- name_method
  
  output <- list()
  output$true_labels <- true_labels
  output$pred_labels <- pred_labels
  output$recall <- rec
  output$precision <- prec
  output$accuracy <- acc
  output$f1 <- f1_k
  output$table <- table
  return(output)
}

