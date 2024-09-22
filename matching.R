################################################################################
#
# Demographic Data Analysis- Matching
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
# 1) Preparing demographic data
#   1.1) Get data 
#   1.2) Exclusions
#   1.3) Create separate data set for the matching
# 2) Data Summary Statistics
#   2.1) Data Summary
#   2.2) Tests
# 3) Matching 
#   3.1) Perform matching
#   3.2) Data Summary
#   3.3) Tests
#   3.4) Visualizations
# 4) Data set for the analysis
#   4.1) Create data set for Baseline and Follow-Up, based on matching df
#   4.2) Create data set only with Follow-Up data
#   4.3) Data Summary
#   4.4) Tests for FU df
#   4.5) Visualizations
# 5) Consort Chart
# 6) Sample description  

################################################################################

# load packages ----

pkgs<-c("groundhog",
        "dplyr",                  
        "here",
        "modelsummary",
        "summarytools",
        "ggplot2",
        "hrbrthemes",
        "viridis",
        "MatchIt",
        "gridExtra",
        "psych",
        "car")      

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


# 1) Preparing demographic data -----

#   1.1) Get data ----

load(here("preprocessed_demo.Rdata"))

#   1.2) Exclusions ----

# create a vector to exclude specific ids from the data; 
# based on preprocessing script "data_processing_visualsearch_0622"

exclude_ids<-c('036','085','925','942','965')

# Outliers check

# check if there are IQ values smaller than 30 and higher than 140, both groups
IQ_outliers <- d %>%
  filter(IQ < 30 | IQ > 140)
# results: 10 IDs on different timepoints, none in T2

# CBCL outliers in TD 
CBCL_outliers <- d %>%
  filter((CBCL == 65 | CBCL > 65) & group =="TD")
#results: 1 observation, FU

# clean the data set
d_1 <- d %>% 
  # exclusion: specific ids
  filter(!ID_Studie %in% exclude_ids)

#   1.3) Create separate data set for the matching, selecting only T2 ----
d_T2 <- d_1 %>% 
  # select only the baseline measurements for the matching, which equals T2
  filter(timepoint == "T2") %>%
  # convert ASD to 1 and TD to 0 for the matching
  mutate(group = ifelse(group == 'ASD', 1, 0)) %>% 
  # remove NAs in IQ
  filter(!is.na(IQ))

# 2) Data Summary Statistics ----

#   2.1) Data Summaries----

stats<- stby(
  data = d_T2,
  INDICES = d_T2$group, 
  FUN = descr, 
  stats = "common" 
)
print(stats)

#   2.2) Tests----

# Age
#Check Assumptions
shapiro.test(d_T2$Age[d_T2$group == "1"])
#---> not significant
shapiro.test(d_T2$Age[d_T2$group == "0"]) 
#---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(Age ~ group, data = d_T2)
#---> significant

# Test age
# check assumptions
shapiro.test(d_T2$test_age[d_T2$group == "1"])
#---< significant
shapiro.test(d_T2$test_age[d_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(test_age ~ group, data = d_T2)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(d_T2$IQ[d_T2$group == "1"])
#---< significant
shapiro.test(d_T2$IQ[d_T2$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(IQ ~ group, data = d_T2)
# ---> significant



# CBCL
# check assumptions
shapiro.test(d_T2$CBCL[d_T2$group == "1"])
#---< significant
shapiro.test(d_T2$CBCL[d_T2$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(CBCL ~ group, data = d_T2)
# ---> significant

# SRS
# check assumptions
shapiro.test(d_T2$SRS_Score[d_T2$group == "1"])
#---< significant
shapiro.test(d_T2$SRS_Score[d_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(SRS_Score ~ group, data = d_T2)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(d_T2$RBSR_Score[d_T2$group == "1"])
#---< significant
shapiro.test(d_T2$RBSR_Score[d_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(RBSR_Score ~ group, data = d_T2)
# ---> significant

# chi square gender before matching
gender_table <- table(d_T2$gender, d_T2$group)
chi_square_test <- chisq.test(gender_table)
chi_square_test
print(chi_square_test$expected)

# 3) Matching ----

#   3.1) Perform matching ----

set.seed(100)
all.match<-matchit(group~test_age+gender,
                   data=d_T2,
                   method='nearest',discard='both',
                   ratio=8, # match four controls to each ASD
                   replace=T,caliper=0.1)
all.match
summary(all.match)
# Unmatched ASD = 24, TD = 8
matched_data <- match.data(all.match)

save(matched_data, file = "matched.Rdata")

#for comparison only with test_age, drops ASD= 23, TD=9
set.seed(100)
all.match.<-matchit(group~test_age,
                   data=d_T2,
                   method='nearest',discard='both',
                   ratio=8, # match four controls to each ASD
                   replace=T,caliper=0.1)
all.match.
summary(all.match.)
matched_data <- match.data(all.match)

#   3.2) Data Summaries ----

# create data set only with the matched data ids
# ids structure in matched_data differ, make it as in pic variable in df_trial
matched_data <- matched_data %>%
  mutate(ID_Studie = sprintf("%03d", as.numeric(ID_Studie)))

d<-df_trial[df_trial$pic %in% matched_data$ID_Studie,]

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

# filter ids with eyetracking data
unique_ids_T2 <- vs %>%
  filter(timepoint == "T2") %>%
  select(pic) %>%
  distinct()

# selct from matching data only the ids that have eye-tracking data
matched_data<-matched_data[matched_data$ID_Studie %in% unique_ids_T2$pic,]

print(stats_matched_data)

#   3.3) Tests ----

# Age
#Check Assumptions
shapiro.test(matched_data$Age[matched_data$group == "1"])
#---> not significant
shapiro.test(matched_data$Age[matched_data$group == "0"]) 
#---> not significant 

# Levene Test
leveneTest(Age~factor(group), matched_data)
#--> significant

# Welch's t-test
t.test(Age~group,matched_data, var.equal = F)


# Test age
# check assumptions
shapiro.test(matched_data$test_age[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$test_age[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(test_age ~ group, data = matched_data)
# ---> not significant

# RBS-R
# check assumptions
shapiro.test(matched_data$IQ[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$IQ[matched_data$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(IQ ~ group, data = matched_data)
# ---> significant



# CBCL
# check assumptions
shapiro.test(matched_data$CBCL[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$CBCL[matched_data$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(CBCL ~ group, data = matched_data)
# ---> significant

# SRS
# check assumptions
shapiro.test(matched_data$SRS_Score[matched_data$group == "1"])
#---< not significant
shapiro.test(matched_data$SRS_Score[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(SRS_Score ~ group, data = matched_data)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(matched_data$RBSR_Score[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$RBSR_Score[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(RBSR_Score ~ group, data = matched_data)
# ---> significant

# chi-square test after matching
gender_group_table_2 <- table(matched_data$gender, matched_data$group)
chi_square_test_3 <- chisq.test(gender_group_table_2)
chi_square_test_3
print(chi_square_test_3$expected) # over 5


#   3.4) Visualization ----

# Density plots before & after matching

# create labels 
gender_labels <- c("1"="Male" , "2"="Female" )
legend_labels <- c("TD","ASD")

# convert group to factor for the plots
matched_data$group <- as.factor(matched_data$group)
d_T2$group <- as.factor(d_T2$group)

# density plots

x_limits <- c(0, 80)
y_limits <- c(0, 0.15)

p1 <- ggplot(d_T2, aes(x=test_age,  fill = group)) + 
  theme_bw() + 
  theme(legend.position="top") +
  scale_fill_viridis(discrete=TRUE, labels = legend_labels) +
  geom_density(alpha=0.5)+
  facet_wrap(~ gender, labeller = labeller(gender=gender_labels))+
  labs(title = "Before Matching",
       x = "Developmental Age",
       y = "Density",
       fill= "Group") +
  theme_ipsum() +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90), 
        plot.title = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5),
  )
print(p1)

p2 <- ggplot(matched_data, aes(x=test_age, fill = group)) + 
  theme_bw() + 
  theme(legend.position="top") +
  scale_fill_viridis(discrete=TRUE, labels = legend_labels) +
  geom_density(alpha=0.5)+
  facet_wrap(~ gender, labeller = labeller(gender=gender_labels)) +
  labs(title = "After Mathing",
       x = "Developmental Age",
       y = "Density",
       fill= "Group") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_ipsum() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5))
print(p2)

Matching_Plot <- grid.arrange( p1, p1, nrow=2)

ggsave("matching_density_plot.png", Matching_Plot, width = 10, height = 8)


# 4) Data set for the analysis ----

#   4.1) Create data set for Baseline and Follow-Up, based on matching df----

# Create a data frame based only on the ids from the matched data
d_m<-d[d$ID_Studie %in% matched_data$ID_Studie,]

# drop all unnecessary timepoints, recode FU2 to FU
df <- d_m %>% 
  filter(timepoint %in% c("T2", "FU", "FU2")) %>% 
  mutate(timepoint = ifelse(timepoint == 'FU2', 'FU', timepoint))

stats_table<-describeBy(df %>% select(Age, test_age, ADI_Toddler_ges, ADOS, CBCL, SRS_Score, RBSR_Score), 
           group = list(df$group, df$timepoint), mat=T)
print(stats_table)

#   4.2) Create data set only with Follow-Up data ----
d_FU <- subset(df, timepoint == "FU")

stat_FU<- stby(
  data = d_FU,
  INDICES = d_FU$group, 
  FUN = descr, 
  stats = "common" 
)

print(stat_FU)

duplicates_FU <- d_FU %>%
  group_by(ID_Studie, timepoint, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

print(duplicates_FU)
# results- 20 duplicates

duplicate_IDS <- d_FU %>%
  filter(ID_Studie %in% duplicates_FU$ID_Studie)

# double entries for these ids, remove the one with NA values in SRS 
# (the SRS scores are unavailable in the double entries)

df_FU <- d_FU %>%
  filter(!is.na(SRS_Score)) 

save(df_FU, file = "df_FU.Rdata")

#### Do not execute the following lines, somehow for some of the IDs that have NAs in IQ, we have eye-tracking data
#### see script visual-search-analysis for the data used for sample description of FU
# ids with NAs in test_age
#ids_na_testage <- df_FU %>%
#  filter(is.na(test_age)) %>%  
#  distinct(ID_Studie)
# 21 --> drop them, because we also don't have eye-tracking data for them, they probably discontinued participant
#df_FU <- df_FU %>% 
#  filter(!is.na(test_age))

#   4.3) Data Summary ----

# summary stats FU
stats_FU<- stby(
  data = df_FU,
  INDICES = df_FU$group, 
  FUN = descr, 
  stats = "common" 
)

print(stats_FU)

#   4.4) Tests for FU df  ----

# check the matching variables, do a t-test and chi-square to see if something change in the significance values
t.test(test_age~group, df_FU)

gender_group_table_4 <- table(df_FU$gender, df_FU$group)
chi_square_test_4 <- chisq.test(gender_group_table_4)
chi_square_test_4

# Age
#Check Assumptions
shapiro.test(d_FU$Age[d_FU$group == "ASD"])
#---> not significant
shapiro.test(d_FU$Age[d_FU$group == "TD"]) 
#---> not significant 

# Levene Test
leveneTest(Age~group, d_FU)
#--> not significant

t.test(Age~group,matched_data)


# Test age
# check assumptions
shapiro.test(matched_data$test_age[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$test_age[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(test_age ~ group, data = matched_data)
# ---> not significant

# RBS-R
# check assumptions
shapiro.test(matched_data$IQ[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$IQ[matched_data$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(IQ ~ group, data = matched_data)
# ---> significant



# CBCL
# check assumptions
shapiro.test(matched_data$CBCL[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$CBCL[matched_data$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(CBCL ~ group, data = matched_data)
# ---> significant

# SRS
# check assumptions
shapiro.test(matched_data$SRS_Score[matched_data$group == "1"])
#---< not significant
shapiro.test(matched_data$SRS_Score[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(SRS_Score ~ group, data = matched_data)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(matched_data$RBSR_Score[matched_data$group == "1"])
#---< significant
shapiro.test(matched_data$RBSR_Score[matched_data$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(RBSR_Score ~ group, data = matched_data)
# ---> significant

#   4.5) Visualizations ----

# Density plots Baseline and Follow Up

x_limit <- c(0, 120)
y_limit <- c(0, 0.15)

p3 <- ggplot(matched_data, aes(x=test_age, fill = group)) + 
  theme_bw() + 
  theme(legend.position="top") +
  scale_fill_viridis(discrete=TRUE, labels = legend_labels) +
  geom_density(alpha=0.5)+
  facet_wrap(~ gender, labeller = labeller(gender=gender_labels)) +
  labs(title = "Baseline",
       x = "Developmental Age",
       y = "Density",
       fill= "Group") +
  scale_x_continuous(limits = x_limit) +
  scale_y_continuous(limits = y_limit) +
  theme_ipsum() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5))
print(p3)

p4 <- ggplot(d_FU, aes(x=test_age,  fill = group)) + 
  theme_bw() + 
  theme(legend.position="top") +
  scale_fill_viridis(discrete=TRUE, labels = legend_labels) +
  geom_density(alpha=0.5)+
  facet_wrap(~ gender, labeller = labeller(gender=gender_labels))+
  labs(title = "Follow-Up",
       x = "Developmental Age",
       y = "Density",
       fill= "Group") +
  scale_x_continuous(limits = x_limit) +
  scale_y_continuous(limits = y_limit) +
  theme_ipsum() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),  
        axis.title.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 90),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5))
print(p4)

Matching_Plot_2 <- grid.arrange( p3, p4, nrow=2)

ggsave("matching_density_plot_baseline-fu.png", Matching_Plot_2, width = 10, height = 8)

#### 5) Data for CONSORT Flow Chart ----

# IDs at the beginning ASD+TD
unique_ids <- n_distinct(d$ID_Studie)
print(unique_ids)
# 151

# IDs devided by group and timepoint
unique_ids_group <- d %>%
  group_by(group, timepoint) %>% 
  summarize(unique_ids_group = n_distinct(ID_Studie), .groups = "drop")
# ids at baseline
# ASD = 76
# TD = 75
# ids at FU --> not for the chart
# ASD = 23+46
# TD = 30

# after excluding the exclude_ids, based on assessment problems
unique_ids_group_d_1 <- d_1 %>%
  group_by(group) %>%  
  summarize(unique_ids_group_d_1 = n_distinct(ID_Studie))
# ASD =75
# TD =73

# After filtering for NAs in IQ - dataset d_T2, see 2.1)
# ASD = 71
# TD = 72

# Matching T2, see summary matched_data 3.2)
# ASD = 63
# TD =48
###---> final sample for the visual search

# Timepoint FU
unique_ids_group_df_FU <- df_FU %>%
  group_by(group) %>%  
  summarize(unique_ids_group_df_FU = n_distinct(ID_Studie))
# ASD = 41
# TD = 17
# ---> duplicates due to demo data preprocessing ---> exclude duplicates
# ---> for 20 ids in the ASD group there is no IQ data, for some of them there is eye-tracking data
# further Consort Flow Chart data with Visual Search data


# 6) Sample Description and partially Consort Chart----
# Load visual search data
load(here("all_data_preprocessed_280824.Rdata"))

rm(df_agg, df)
# create data set only with the matched data ids
# ids structure in matched_data differ, make it as in pic variable in df_trial
matched_data <- matched_data %>%
  mutate(ID_Studie = sprintf("%03d", as.numeric(ID_Studie)))

d<-df_trial[df_trial$pic %in% matched_data$ID_Studie,]

# check unique timepoint names
table(d$timepoint)
#FU2  FU3    K K_FU   T2   T4   T6 
#704  409 2689  554 3280 2816 2656 
# K - Control Group T2
# K_FU - Control Group FU
# FU2 - ASD Group FU

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


summary_table <- vs %>%
  group_by(group, timepoint) %>%
  summarize(count = n(), .groups = 'drop')

print(summary_table)
# # A tibble: 4 × 3
#group timepoint count
#<fct> <fct>     <int>
# 1 ASD   T2         3280
#2 ASD   FU          704
#3 TD    T2         2689
#4 TD    FU          554


# for Consort Flow Chart
unique_ids_vs <- vs %>%
  group_by(group, timepoint) %>%  
  summarize(unique_ids_vs= n_distinct(pic), .groups = "drop")
print(unique_ids_vs)
# 1 ASD   T2                      54
# 2 ASD   FU                      20
# 3 TD    T2                      48
# 4 TD    FU                      15

# IDs Timepoint T2, Both Groups
ids_list_T2<- vs %>%
  filter( timepoint == "T2") %>%  
  distinct(pic) 

# IDs Timepoint FU, Both Groups
ids_list_FU<- vs %>%
  filter( timepoint == "FU") %>%  
  distinct(pic) 

## Data Sets for Sample Description----
# T2
filtered_df_T2<- matched_data %>% 
  filter(ID_Studie %in% ids_list_T2$pic)
# FU
filtered_df_FU <- df_FU %>%
  filter(ID_Studie %in% ids_list_FU$pic)

## T2 Tests ----

# Age
#Check Assumptions
shapiro.test(filtered_df_T2$Age[filtered_df_T2$group == "1"])
#---> not significant
shapiro.test(filtered_df_T2$Age[filtered_df_T2$group == "0"]) 
#---> not significant 

# Levene Test
car::leveneTest(Age~group, filtered_df_T2)
#--> significant

# Welch's t-test
t.test(Age~group,filtered_df_T2, var.equal = F)
##--> significant

# Test age
# check assumptions
shapiro.test(filtered_df_T2$test_age[filtered_df_T2$group == "1"])
#---< significant
shapiro.test(filtered_df_T2$test_age[filtered_df_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(test_age ~ group, data = filtered_df_T2)
# ---> not significant

# RBS-R
# check assumptions
shapiro.test(filtered_df_T2$IQ[filtered_df_T2$group == "1"])
#---< significant
shapiro.test(filtered_df_T2$IQ[filtered_df_T2$group == "0"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(IQ ~ group, data = filtered_df_T2)
# ---> significant

# CBCL
# check assumptions
shapiro.test(filtered_df_T2$CBCL[filtered_df_T2$group == "1"])
#---< not significant
shapiro.test(filtered_df_T2$CBCL[filtered_df_T2$group == "0"]) 
# ---> not significant

# Levene Test
car::leveneTest(CBCL~group, filtered_df_T2)
#--> not significant

# t-test
t.test(CBCL~group,filtered_df_T2)
##--> significant

# SRS
# check assumptions
shapiro.test(filtered_df_T2$SRS_Score[filtered_df_T2$group == "1"])
#---< not significant
shapiro.test(filtered_df_T2$SRS_Score[filtered_df_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(SRS_Score ~ group, data = filtered_df_T2)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(filtered_df_T2$RBSR_Score[filtered_df_T2$group == "1"])
#---< significant
shapiro.test(filtered_df_T2$RBSR_Score[filtered_df_T2$group == "0"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(RBSR_Score ~ group, data = filtered_df_T2)
# ---> significant

# chi-square test 
sex_group_T2<- table(filtered_df_T2$gender, filtered_df_T2$group)
chi_square_T2 <- chisq.test(sex_group_T2)
chi_square_T2
#---> not significant
print(chi_square_T2$expected) #--> over 5


##  Sample summary statistics T2 ----
stats_T2<- stby(
  data = filtered_df_T2,
  INDICES = filtered_df_T2$group, 
  FUN = descr, 
  stats = c("common")
)

print(stats_T2)

##  FU Tests ----
# t-tests

# Age
#Check Assumptions
shapiro.test(filtered_df_FU$Age[filtered_df_FU$group == "ASD"])
#---> not significant
shapiro.test(filtered_df_FU$Age[filtered_df_FU$group == "TD"]) 
#---> not significant 

# Levene Test
leveneTest(Age~factor(group), filtered_df_FU)
#--> not significant

# t-test
t.test(Age~group,filtered_df_FU)
##--> significant

# Test age
# check assumptions
shapiro.test(filtered_df_FU$test_age[filtered_df_FU$group == "ASD"])
#---< significant
shapiro.test(filtered_df_FU$test_age[filtered_df_FU$group == "TD"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(test_age ~ group, data = filtered_df_FU, exact = F)
# ---> not significant

# IQ
# check assumptions
shapiro.test(filtered_df_FU$IQ[filtered_df_FU$group == "ASD"])
#---< significant
shapiro.test(filtered_df_FU$IQ[filtered_df_FU$group == "TD"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(IQ ~ group, data = filtered_df_FU, exact = F)
# ---> significant

# CBCL
# check assumptions
shapiro.test(filtered_df_FU$CBCL[filtered_df_FU$group == "ASD"])
#---< significant
shapiro.test(filtered_df_FU$CBCL[filtered_df_FU$group == "TD"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test( CBCL ~ group, data = filtered_df_FU, exact = F)
# ---> significant

# SRS
# check assumptions
shapiro.test(filtered_df_FU$SRS_Score[filtered_df_FU$group == "ASD"])
#---< not significant
shapiro.test(filtered_df_FU$SRS_Score[filtered_df_FU$group == "TD"]) 
# ---> not significant

# Wilcoxon Rank-Sum test 
wilcox.test(SRS_Score ~ group, data = filtered_df_FU, exact = F)
# ---> significant

# RBS-R
# check assumptions
shapiro.test(filtered_df_FU$RBSR_Score[filtered_df_FU$group == "ASD"])
#---< not significant
shapiro.test(filtered_df_FU$RBSR_Score[filtered_df_FU$group == "TD"]) 
# ---> significant

# Wilcoxon Rank-Sum test 
wilcox.test(RBSR_Score ~ group, data = filtered_df_FU, exact = F)
# ---> significant


# chi-square test 
sex_group_T2<- table(filtered_df_FU$gender, filtered_df_FU$group)
chi_square_T2 <- chisq.test(sex_group_T2)
chi_square_T2
#---> not significant

#sex
# chi-square test 
sex_group_FU<- table(filtered_df_FU$gender, filtered_df_FU$group)
chi_square_FU <- chisq.test(sex_group_FU)
chi_square_FU
#--> not significant
print(chi_square_FU$expected) #--> over 5


##  Sample summary statistics FU-----
stats_FU<- stby(
  data = filtered_df_FU,
  INDICES = filtered_df_FU$group, 
  FUN = descr, 
  stats = c("common")
)

print(stats_FU)



