################################################################################
#
# Visual Search Main Analysis
# Author: Iskra Todorova
# Last Update: 18-09-2024
# R version : R 4.4.1
#
################################################################################
#
# Before you begin:
# - Loading /installing packages
# - Setting working directory
#
################################################################################
#
#
# Outline:
# 1) Preparing Data
#   1.1) Load data 
#   1.2) Organize
# 2) Visual Search (VS) Data
#   2.1) Fixation Duration
#     2.1.1) Model fixation duration proportions
#     2.1.2) Robust model, add dispersion variability for the predictors
#   2.2) Clean Data VS
#   2.3) Hits
#     2.3.1) Hit Model
#     2.3.2) Adding fixation proportions to the hit model
#     2.3.3) Additional Analysis for Hits, based on task relevant variables
#   2.4) Time to target
#     2.4.1) Model time to target
#     2.4.2)Robust model for reaction time
# 3) Check the ids for descriptive Statistic for the sample
#
#
################################################################################

# load packages ----

pkgs<-c("groundhog",
        "dplyr",                  
        "here",
        "ggplot2",
        "lme4",
        "lmerTest",
        "emmeans",
        "performance",
        "car",
        "glmmTMB",
        "robustlmm",
        "broom.mixed",
        "DHARMa")          

# check if required packages are installed
installed_packages = pkgs %in% rownames(installed.packages())

# install packages if not installed
if (any(installed_packages == FALSE)) {
  install.packages(pkgs[!installed_packages])
}

lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Package", pkg, "not found."))
  }
})

# run if code doesn't work due to package updates or unavailable functions
# groundhog.day <- '2024-04-20'
# groundhog.library(pkgs,groundhog.day, tolerate.R.version = '4.2.2')

# Set paths ----

# set dynamically based on where this script is saved

set_here()


# 1) Preparing Data -----
#   1.1) Load data----
# Load visual search data
load(here("all_data_preprocessed_280824.Rdata"))

rm(df)
rm(df_agg)

# load matched data = demo data T2
load(here("matched.Rdata"))

#   1.2) Organize ----

hit_summary <-df_trial %>%
  group_by(target_position) %>%
  summarize(n_hits = sum(missed_target, na.rm = T), n_total = n())
#-- -> no hits on position 735_815_900_980

# create data set only with the matched data ids
# ids structure in matched_data differ, make it as in pic variable in df_trial
matched_data <- matched_data %>%
  mutate(ID_Studie = sprintf("%03d", as.numeric(ID_Studie)))

d<-df_trial[df_trial$pic %in% matched_data$ID_Studie,]

# check unique timepoint names
table(d$timepoint)

#  select Baseline and Follow-Up Data, discard the rest
vs <- d %>% 
  filter(timepoint %in% c("T2", "K", "K_FU","FU2"))

table(vs$timepoint)

# recode timepoint values
vs <- vs %>%
  mutate(timepoint = dplyr::recode(timepoint, 
                                   "K" = "T2", 
                                   "FU2" = "FU",
                                   "K_FU" = "FU"))

# add only two levels in timepoint
vs$timepoint <- factor(vs$timepoint, levels = c("T2", "FU"))

# 2) Visual Search Data ----
  
# convert missed_target variable to hit, easier for understanding,
# missed_target is False and True, meaning False is a Hit
vs$hit <- !vs$missed_target


#   2.1) Fixation duration ---- 

# Data set with only 32 trials 
vs_32 <- vs[vs$trial %in% levels(vs$trial)[1:32], ]

# remove NAs in fixation duration
vs_fixation_dur_32 <- vs_32 %>% 
  filter(!is.na(fixation_duration))

hist(vs_fixation_dur_32$fixation_duration)

### calculate fixation duration proportions
# Calculate with a trial length  of 4 seconds: 2s fix cross + 1.5 stimulus + 0.5 jitter
vs_fixation_dur_32$fixation_proportions <- vs_fixation_dur_32$fixation_duration/4

hist(vs_fixation_dur_32$fixation_proportions)

stats_fdp_32 <- vs_fixation_dur_32 %>%
  group_by(group, timepoint) %>%
  summarise(
    mean_fixation = mean(fixation_proportions, na.rm = TRUE) * 100,  
    sd_fixation = sd(fixation_proportions, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# mean fixation durations proportions by participant, group and timepoint
vs_prop_pic_32 <- vs_fixation_dur_32 %>%
  group_by(pic, group, timepoint) %>%
  summarise(mean_fixation_proportion = mean(fixation_proportions, na.rm = TRUE), .groups = "drop")

# violin fix dur participant, group and timepoint
fix_dur_violin <- ggplot(vs_prop_pic_32, aes(x = group, y = mean_fixation_proportion, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = "red", outlier.shape = 1) + 
  facet_wrap(~timepoint) + 
  labs(#title = "Fixation Duration Mean Proportion by Trial",
       x = "Group",
       y = "Proportions") +
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        plot.title = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(fill = "white", color = NA))
fix_dur_violin 

ggsave("Fixation-proportions-id-timepoint-group.png", fix_dur_violin, width = 10, height = 8)

#     2.1.1) Model fixation duration proportions ----
# Use Beta Regression, better for proportions , make sure no values = 0 or 1
m_fd_prop_32 <- glmmTMB(fixation_proportions ~ group * timepoint + (1 | pic),
                     data = vs_fixation_dur_32,
                     family = beta_family())
summary(m_fd_prop_32)
check_model(m_fd_prop_32)

#     2.1.2) Robust model, add dispersion variability for the predictors ----
m_fd_rob <- glmmTMB(fixation_proportions ~ group * timepoint + (1 | pic),
                            dispformula = ~ group * timepoint,  
                            data = vs_fixation_dur_32,
                            family = beta_family())
summary(m_fd_rob)

performance::check_model(m_fd_rob)
performance::check_collinearity(m_fd_rob)
performance::performance(m_fd_rob)

emmeans(m_fd_rob, ~  group | timepoint, type= "response")
summary(contrast(emmeans(m_fd_rob, ~  group | timepoint, type = "response"), "pairwise"), infer = c(TRUE, TRUE))

# Plot EMMs
emm_fdp <- as.data.frame(emmeans(m_fd_rob, ~  group | timepoint, type= "response"))

fix_dur_emms <- ggplot(emm_fdp, aes(x = timepoint, y = response, group = group, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  labs(#title = "Estimated Probabilities of Fixation Duration Proportions",
       x = "Timepoint",
       y = "Probability ") +
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        #plot.title = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(fill = "white", color = NA))

fix_dur_emms

ggsave("EMMS-fix-dur.png", fix_dur_emms, width = 10, height = 8)

#   2.2) Clean Data VS  ----

# Time to target 
# exclude the cases that have values under 0.1
vs_fixation_dur_32 <- vs_fixation_dur_32 %>% 
  filter(time_to_target > 0.1)

# Hits 
# exlude all NAs in Hit
vs_fixation_dur_32 <- vs_fixation_dur_32 %>% 
  filter(!is.na(hit))

# hit duration
# check if there are values below 0.1 s
any(vs_fixation_dur_32$hit_duration<0.1)
#-->none

vs_df <- vs_fixation_dur_32

#   2.3) Hits -----
#     2.3.1) Hit Model ----
vs_df$hit <- factor(ifelse(vs_df$hit, "Hit", "Miss"), levels = c("Miss", "Hit"))

m_hit <- glmer(hit ~ group * timepoint + (1 | pic), vs_df,  family = binomial(link = "logit"))
summary(m_hit)

performance::check_model(m_hit)
sim_res <- simulateResiduals(m_hit)
plot(sim_res)
performance(m_hit)
vif(m_hit)
check_collinearity(m_hit)


emm_hit <- emmeans(m_hit, ~ group | timepoint, type = "response")
emm_hit

contrast(emm_hit, "pairwise")
summary(contrast(emm_hit, "pairwise"), infer = c(TRUE, TRUE))

# Plot EMM
emm_hits <- as.data.frame(emm_hit)

plot_hit_prob <- ggplot(emm_hits, aes(x = timepoint, y = prob, group = group, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  labs(#title = "Probability of Hit",
       x = "Timepoint",
       y = "Probability") +
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        #plot.title = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(fill = "white", color = NA))


plot_hit_prob
ggsave("EMMS-hit.png", plot_hit_prob, width = 10, height = 8)


#     2.3.2) adding fixation proportions to the hit model ----

m_hit_fd <- glmer(hit ~ group * timepoint + fixation_proportions +(1 | pic), vs_df,  family = binomial())
summary(m_hit_fd)

# model comparison
anova(m_hit, m_hit_fd)

performance::check_model(m_hit_fd)
sim_res_hfd <- simulateResiduals(m_hit_fd)
plot(sim_res_hfd)
performance(m_hit_fd)
vif(m_hit_fd)
check_collinearity(m_hit_fd)

emm_hit_p <- emmeans(m_hit_fd, ~ fixation_proportions, type="response")
emm_hit_p

#     2.3.3) Additional Analysis for Hits, based on task relevant variables ----
m_hit_add <- glmer(hit ~ target + target_position  +(1 | pic), vs_df,  family = binomial(link = "logit"))
summary(m_hit_add)

performance::check_model(m_hit_add)
sim_res_add <- simulateResiduals(m_hit_add)
plot(sim_res_add)
performance(m_hit_add)
vif(m_hit_add)
check_collinearity(m_hit_add)

# target letters
emmeans(m_hit_add, ~ target , type = "response")
contrast(emmeans(m_hit_add, ~ target , type = "response"), "pairwise",  infer = c(TRUE, TRUE))

letter_type_prop <- summary(emmeans(m_hit_add, ~ target , type = "response"))
letter_prob_df <- as.data.frame(letter_type_prop)

# letters plot
bar_plot_letters <- ggplot(letter_prob_df, aes(x = target, y = prob, fill = target)) +
  geom_bar(stat = "identity", position = "dodge", fill= "#00BFC4", color = "black", alpha = 0.5) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, position = position_dodge(0.9)) +  
  labs(x = "Target Letter", y = "Hit Probability")+
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        plot.title = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(fill = "white", color = NA))

bar_plot_letters

ggsave("Letters_prob_2.png", bar_plot_letters, width = 8, height = 6)

# target position
emmeans(m_hit_add, ~ target_position , type = "response")
contrast(emmeans(m_hit_add, ~ target_position , type = "response"), "pairwise")

# Target positions plot
all_target_positions <- data.frame(
    Position = c("985_65_1150_230", "485_565_650_730", "1235_315_1400_480",
                 "1235_565_1400_730", "735_65_900_230", "485_315_650_480", "985_815_1150_980", "735_815_900_980"),
    x1 = c(985,  485, 1235, 1235, 735, 485, 985, 735),
    y1 = c(65,  565, 315, 565, 65, 315, 815, 815),
    x2 = c(1150, 650, 1400, 1400, 900, 650, 1150, 900),
    y2 = c(230, 730, 480, 730, 230, 480, 980, 980),
    prob = c(0.872, 0.743, 0.847, 0.831, 0.874, 0.845, 0.814, NA),
    se = c(0.0306, 0.0411, 0.0324, 0.0373, 0.0301, 0.0324, 0.0396, NA),
    significant = c(FALSE,TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)  
  )
  
  # colors
  normal_color <- "#F8766D"
  significant_color <- "#00BFC4"  
  
target_pos<-  ggplot(all_target_positions, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2)) +
    #geom_rect(aes(xmin = 0, xmax = 1920, ymin = 0, ymax = 1080), fill = NA, color = "black", size = 0.5) +
    geom_rect(aes(fill = ifelse(significant, significant_color, normal_color)), alpha = 0.5) +
    geom_text(aes(x = (x1 + x2) / 2, y = (y1 + y2) / 2, 
                  label = sprintf("M = %.3f\nSE = %.3f%s", prob, se, ifelse(significant, "*", ""))),
              size = 2.5, fontface = "bold",
              hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(name = "X Coordinate", limits = c(0, 1920)) +
    scale_y_continuous(name = "Y Coordinate", limits = c(0, 1080)) +
    labs(fill = "Target Position") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
          axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
          plot.title = element_text(size = 12, face = "bold"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_blank(), 
         #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          legend.position = "none")
target_pos  

ggsave("target-positions-probabilities.png", target_pos , width = 10, height = 8)   
  
#   2.4) Time to target -----

# filter only Hits in each df

vs_rt <- vs_df %>% 
  filter(hit == "Hit")

#     2.4.1) Model time to target----
m_rt <- lmer(scale(time_to_target) ~ group * timepoint  + fixation_proportions + (1 | pic) , data = vs_rt)
anova(m_rt)

performance::check_model(m_rt)
sim_res_rt <- simulateResiduals(m_rt)
plot(sim_res_rt)


check_heteroscedasticity(m_rt)
vif(m_rt)
performance(m_rt)

summary(m_rt)

emmeans(m_rt, ~ group)
emmeans(m_rt, ~ fixation_proportions)
emmeans(m_rt, ~ group | timepoint)
emmeans(m_rt, ~ group * timepoint)

contrast(emmeans(m_rt, ~ group), infer = c(TRUE, TRUE))

contrast(emmeans(m_rt, ~ group * timepoint), infer = c(TRUE, TRUE) )
contrast(emmeans(m_rt, ~ group * timepoint))

emm_rt <- as.data.frame(emmeans(m_rt, ~ group * timepoint))

plot_rt_prob <- ggplot(emm_rt, aes(x = timepoint, y = emmean, group = group, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  labs(#title = "Probability of Hit",
    x = "Timepoint",
    y = "Mean Reaction") +
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        #plot.title = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(fill = "white", color = NA))


plot_rt_prob
ggsave("EMMS-rt.png", plot_hit_prob, width = 10, height = 8)

#     2.4.2)Robust model for reaction time ----
m_rt_rob <- rlmer(scale(time_to_target) ~ group * timepoint  + fixation_proportions +  (1 | pic) , data = vs_rt)
summary(m_rt_rob)

model <- "m_rt"
model_rob <- "m_rt_rob"

output <- as.data.frame(tidy(get(model_rob))[tidy(get(model_rob))$effect=="fixed",])
output$df <- summary(get(model))$coefficients[,3]# df Satterthwaite
output$p.rob <- 2*pt(abs(summary(get(model_rob))$coefficients[,3]), summary(get(gsub(model_rob,pattern = "_rob",replacement = "")))$coefficients[,3], lower=FALSE)
output[4:length(output)] <- round(output[4:length(output)],digits = 3)
output#[grepl(output$term,pattern = "group"),

performance::check_model(m_rt_rob)
sim_res_rt_rob <- simulateResiduals(m_rt_rob)
plot(sim_res_rt_rob)
check_heteroscedasticity(m_rt_rob)

vif(m_rt_rob)
performance(m_rt_rob)

# 3) Check the ids for descriptive Statistic for the sample ----
# vs_df used for fixation durations and hits
ASD_T2_ids <- vs_df %>%
  filter(group == "ASD", timepoint == "T2") %>%
  select(pic) %>%
  distinct()
#--> 52
ASD_FU_ids <- vs_df %>%
  filter(group == "ASD", timepoint == "FU") %>%
  select(pic) %>%
  distinct()
#--> 18
TD_T2_ids <- vs_df %>%
  filter(group == "TD", timepoint == "T2") %>%
  select(pic) %>%
  distinct()
#--> 46
TD_FU_ids <- vs_df %>%
  filter(group == "TD", timepoint == "FU") %>%
  select(pic) %>%
  distinct()
#-->15

# vs_rt used to calculatet reaction time only in hit cases
ASD_T2_ids_rt <- vs_rt %>%
  filter(group == "ASD", timepoint == "T2") %>%
  select(pic) %>%
  distinct()
#--> 50
ASD_FU_ids_rt <- vs_rt %>%
  filter(group == "ASD", timepoint == "FU") %>%
  select(pic) %>%
  distinct()
#--> 16
TD_T2_ids_rt <- vs_rt %>%
  filter(group == "TD", timepoint == "T2") %>%
  select(pic) %>%
  distinct()
#--> 45
TD_FU_ids_rt <- vs_rt %>%
  filter(group == "TD", timepoint == "FU") %>%
  select(pic) %>%
  distinct()
#-->15

trial_number_fd <- vs_df %>%
  group_by(pic,group, timepoint) %>%
  summarise(
    Trials = n(),  
    .groups = "drop"
  )

fd_mean_sd_sum <- vs_df %>%
  group_by(group, timepoint) %>%
  summarise(
    Trials = n(),  
    Mean_FD = mean(fixation_duration, na.rm = TRUE),  
    SD_FD = sd(fixation_duration, na.rm = TRUE),
    .groups = "drop"
  )

summary_vs_fd <- trial_number_fd %>% 
  group_by(group,timepoint) %>% 
  summarise(
    IDs = n_distinct(pic),
    trial_m = mean(Trials, na.rm= T),
    trial_sd = sd(Trials, na.rm= T),
    trial_median = median(Trials, na.rm=T),
    trial_min = min(Trials, na.rm=T),
    trial_max = max(Trials, na.rm=T),
    .groups ="drop"
  ) %>% 
  merge(fd_mean_sd_sum)

trial_number_rt <- vs_rt %>%
  group_by(pic,group, timepoint) %>%
  summarise(
    Trials = n(),  
    .groups = "drop"
  )

rt_mean_sd_sum <- vs_rt %>%
  group_by(group, timepoint) %>%
  summarise(
    Trials = n(),  
    Mean_FD = mean(time_to_target, na.rm = TRUE),  
    SD_FD = sd(time_to_target, na.rm = TRUE),
    .groups = "drop"
  )

summary_vs_rt <- trial_number_rt %>% 
  group_by(group,timepoint) %>% 
  summarise(
    IDs = n_distinct(pic),
    trial_m = mean(Trials, na.rm= T),
    trial_sd = sd(Trials, na.rm= T),
    trial_median = median(Trials, na.rm=T),
    trial_min = min(Trials, na.rm=T),
    trial_max = max(Trials, na.rm=T),
    .groups ="drop"
  ) %>% 
  merge(rt_mean_sd_sum)


