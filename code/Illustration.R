library("ecp")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("tidyr")
library("RColorBrewer")
library("patchwork")

source("library.R")

##################################################################
# Figure 3 in Section 3.2 & Figure 6 in Appendix:
##################################################################
# Choose the subject: Uncommentize the following
load("sternrelief9.RData")
# load("gayevskyreds18.RData")

# Number of samples:
TT <- dim(data)[1]

# Passive sensing data:
ts_passive  <- ts.union(act_time = ts(data$active_time),
                        step_cnt = ts(data$step_count),
                        conv_pc = ts(data$conversation_percent),
                        tic_voice = ts(data$tic_voiced_time),
                        home_time = ts(data$time_at_home),
                        slp_dr = ts(data$sleep_duration),
                        slp_int = ts(data$sleep_interruption),
                        trv_dm = ts(data$travel_diameter),
                        act_dr = ts(data$total_activity_duration),
                        lc_dr = ts(data$total_location_duration),
                        gyration = ts(data$radius_of_gyration))
times <- attr(ts_passive,"tsp")
df_passive <- cbind(as.data.frame(ts_passive), time = seq(times[1], times[2], 1/times[3]))
tibble_passive <- pivot_longer(df_passive, 1:11)
tibble_passive$name <- factor(tibble_passive$name,levels=unique(tibble_passive$name))

# Active data:
Y1 <- Y2 <- data.frame(data$stress,
                       data$arousal,
                       data$valence)
colnames(Y1) <- colnames(Y2) <- c("stress","arousal","valence")
resp_adj <- response_adjustment(Y1,Y2)

ts_active  <- ts.union(arousal = ts(data$arousal),
                       valence = ts(data$valence),
                       stress = ts(data$stress),
                       adj_stress = ts(resp_adj$series$train))
df_active <- cbind(as.data.frame(ts_active), time = seq(times[1], times[2], 1/times[3]))

tibble_active <- pivot_longer(df_active, 1:4)
tibble_active$name <- factor(tibble_active$name,levels=unique(tibble_active$name))


mean_level <- data.frame(name=c("arousal","valence"),
                         vals = c(mean(data$arousal),mean(data$valence)),
                         stringsAsFactors = FALSE)
mean_level$name <- factor(mean_level$name,levels=unique(mean_level$name))
add_time <- c(1:TT)[resp_adj$time$train_add > 0]
sub_time <- c(1:TT)[resp_adj$time$train_sub < 0]
add_score <- data.frame(name=replicate(length(add_time),"adj_stress"),
                        vals = add_time,
                        stringsAsFactors = FALSE)
add_score$name <- factor(add_score$name,levels=unique(add_score$name))
sub_score <- data.frame(name=replicate(length(sub_time),"adj_stress"),
                        vals = sub_time,
                        stringsAsFactors = FALSE)
sub_score$name <- factor(sub_score$name,levels=unique(sub_score$name))

# Plot initial data:
paired <- colorRampPalette(brewer.pal(12,"Paired"))
pl_active <- ggplot(tibble_active, aes(time, value, colour=name)) +
  geom_line() +
  scale_x_continuous(name="Time",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="",breaks=c(1:6)) +
  scale_color_manual(values=c("violetred","dodgerblue","black","black")) +
  geom_vline(data=add_score,aes(xintercept=vals),col="violetred") +
  geom_vline(data=sub_score,aes(xintercept=vals),col="dodgerblue") +
  geom_hline(data=mean_level,aes(yintercept=vals),linetype="dashed") +
  facet_grid(name~., switch = "y", scales = "fixed",
             labeller=labeller(name = act.name.labs)) +
  labs(title = "Active variables") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size=10, face=2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        legend.position = "none")

pl_passive <- ggplot(tibble_passive, aes(time, value, colour=name)) +
  geom_line(alpha=1.5) +
  scale_x_continuous(name="Time",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="") +
  scale_color_manual(values=paired(11)) +
  facet_grid(name~., switch = "y", scales = "free_y",
             labeller=labeller(name = pass.name.labs)) +
  labs(title = "Passive sensing variables") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size = 10, face = 2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# Choose the subject: Uncommentize the following
pdf(file ="./initialization1.pdf", width=13, height=13)
# pdf(file ="./initialization2.pdf", width=13, height=13)
par(mar=c(0,0,0,0))
ggarrange(pl_active,pl_passive,nrow=2,align="v",heights=c(5,12))
dev.off()


##################################################################
# Figure 4 in Section 3.2 & Figure 7 in Appendix:
##################################################################
X_name <- c("active_time", "step_count",
            "conversation_percent","tic_voiced_time",
            "time_at_home", "sleep_duration",
            "sleep_interruption",
            "total_activity_duation","total_location_duration",
            "travel_diameter","radius_of_gyration")
Y_name <- c("arousal","valence","stress","adj_stress")

X_train <- data.frame(data$active_time, 
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
Y_train <- data.frame(data$arousal,
                      data$valence,
                      data$stress,
                      resp_adj$series$train)

set.seed(123)
# Step 1. CPD
Step1 <- SPD_step1(X_train,Y_train,X_name,Y_name)
# Step 2: Validation of significance
Step2 <- SPD_step2(Step1)
# Step 3: Validation of distinction in equations
Step3 <- SPD_step3(Step2)
# Step 4: Validation of distinction in responses
Step4 <- SPD_step4(Step3)


# Display results:
tibble_stress <- tibble_active[which(tibble_active$name=="adj_stress"|
                                       tibble_active$name=="stress"),]

# Result after Step 1:
pl_stress_step1 <- ggplot(tibble_stress, aes(time, value, colour=name)) +
  geom_line() +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="",breaks=c(1:6)) +
  scale_color_manual(values=c("black","black")) +
  geom_vline(data=add_score,aes(xintercept=vals),linetype="dashed",col="violetred") +
  geom_vline(data=sub_score,aes(xintercept=vals),linetype="dashed",col="dodgerblue") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  facet_grid(name~., switch = "y", scales = "fixed",
             labeller=labeller(name = act.name.labs)) +
  labs(title = "Result of Step 1") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size=8, face=2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=6),
        legend.position = "none")

pl_passive_step1 <- ggplot(tibble_passive, aes(time, value, colour=name)) +
  geom_line(alpha=1.5) +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  scale_color_manual(values=paired(11)) +
  facet_grid(name~., switch = "y", scales = "free_y",
             labeller=labeller(name = pass.name.labs)) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size = 8, face = 2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")


# Result after Step 2:
pl_stress_step2 <- ggplot(tibble_stress, aes(time, value, colour=name)) +
  geom_line() +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="",breaks=c(1:6)) +
  scale_color_manual(values=c("black","black")) +
  geom_vline(data=add_score,aes(xintercept=vals),linetype="dashed",col="violetred") +
  geom_vline(data=sub_score,aes(xintercept=vals),linetype="dashed",col="dodgerblue") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  facet_grid(name~., switch = "y", scales = "fixed",
             labeller=labeller(name = act.name.labs)) +
  labs(title = "Result of Step 2") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size=8, face=2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=6),
        legend.position = "none")

pl_passive_step2 <- ggplot(tibble_passive, aes(time, value, colour=name)) +
  geom_line(alpha=1.5) +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  scale_color_manual(values=paired(11)) +
  facet_grid(name~., switch = "y", scales = "free_y",
             labeller=labeller(name = pass.name.labs)) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size = 8, face = 2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")


# Result after Step 3:
pl_stress_step3 <- ggplot(tibble_stress, aes(time, value, colour=name)) +
  geom_line() +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="",breaks=c(1:6)) +
  scale_color_manual(values=c("black","black")) +
  geom_vline(data=add_score,aes(xintercept=vals),linetype="dashed",col="violetred") +
  geom_vline(data=sub_score,aes(xintercept=vals),linetype="dashed",col="dodgerblue") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  geom_vline(xintercept=Step3$t_hat,col="orange",linetype=2,size=2) +
  facet_grid(name~., switch = "y", scales = "fixed",
             labeller=labeller(name = act.name.labs)) +
  labs(title = "Result of Step 3") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size=8, face=2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=6),
        legend.position = "none")

pl_passive_step3 <- ggplot(tibble_passive, aes(time, value, colour=name)) +
  geom_line(alpha=1.5) +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  geom_vline(xintercept=Step3$t_hat,col="orange",linetype=2,size=2) +
  scale_color_manual(values=paired(11)) +
  facet_grid(name~., switch = "y", scales = "free_y",
             labeller=labeller(name = pass.name.labs)) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size = 8, face = 2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# Result after Step 4:
pl_stress_step4 <- ggplot(tibble_stress, aes(time, value, colour=name)) +
  geom_line() +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="",breaks=c(1:6)) +
  scale_color_manual(values=c("black","black")) +
  geom_vline(data=add_score,aes(xintercept=vals),linetype="dashed",col="violetred") +
  geom_vline(data=sub_score,aes(xintercept=vals),linetype="dashed",col="dodgerblue") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  geom_vline(xintercept=Step3$t_hat,col="orange",linetype=2,size=2) +
  geom_vline(xintercept=Step4$t_hat,col="red",linetype=1,size=2) +
  facet_grid(name~., switch = "y", scales = "fixed",
             labeller=labeller(name = act.name.labs)) +
  labs(title = "Result of Step 4") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size=8, face=2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=6),
        legend.position = "none")

pl_passive_step4 <- ggplot(tibble_passive, aes(time, value, colour=name)) +
  geom_line(alpha=1.5) +
  scale_x_continuous(name="",breaks=seq(times[1],times[2],10)) +
  scale_y_continuous(name="") +
  geom_vline(xintercept=Step1$t_hat,col="green",linetype=4,size=2) +
  geom_vline(xintercept=Step2$t_hat,col="yellow",linetype=3,size=2) +
  geom_vline(xintercept=Step3$t_hat,col="orange",linetype=2,size=2) +
  geom_vline(xintercept=Step4$t_hat,col="red",linetype=1,size=2) +
  scale_color_manual(values=paired(11)) +
  facet_grid(name~., switch = "y", scales = "free_y",
             labeller=labeller(name = pass.name.labs)) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "npc"),
        strip.text = element_text(size = 8, face = 2),
        strip.text.y.left = element_text(angle=0,hjust=1),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# Choose the subject: Uncommentize the following
pdf(file = "./result_illustration1.pdf", width=13, height=13)
# pdf(file = "./result_illustration2.pdf", width=13, height=13)
par(mar=c(0,0,0,0))
ggarrange(
  ggarrange(pl_stress_step1,pl_passive_step1,nrow=2,align="v",heights=c(3,9)),
  ggarrange(pl_stress_step2,pl_passive_step2,nrow=2,align="v",heights=c(3,9)),
  ggarrange(pl_stress_step3,pl_passive_step3,nrow=2,align="v",heights=c(3,9)),
  ggarrange(pl_stress_step4,pl_passive_step4,nrow=2,align="v",heights=c(3,9)),
  nrow=2,ncol=2)
dev.off()

