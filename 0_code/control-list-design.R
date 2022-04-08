# This script run analysis and create figures and tables for;
#Agerberg, Mattias, and Marcus Tannenberg. 
#"Dealing With Measurement Error in List Experiments:Choosing the Right Control List Design" 
#Research & Politics (2021) 


#load required libraries
library(tidyverse)
library(estimatr)
library(stargazer)
library(texreg)
library(rio)

set.seed(111)

###
#Our survey
###

#import dataset from /data folder
df <- rio::import("data/survey-ds.csv")

# create treatment, control and control_placebo indicator and outcome variable for the list
df <- df %>%
  mutate(treat_b = ifelse(!is.na(b_treatment), 1, 0), # get 1 if presented treatment list
         y_b = ifelse(treat_b==1, b_treatment,
                      ifelse(treat_b==0, b_control, NA)), # list outcome for treatment and (conventional) control
         y_b_placebo =  ifelse(treat_b==1, b_treatment,
                               ifelse(treat_b==0, b_placebo, NA)),# list outcome for treatment and placebo control
         sampler_dummy = sample(0:1, length(df[,1]), replace= TRUE),
         y_b_mixed = ifelse(sampler_dummy==1, b_control,            #uses half conventional and half placebo control
                            ifelse(sampler_dummy==0, b_placebo, NA)),
         y_b_mixed = ifelse(treat_b==1, b_treatment, y_b_mixed) # list outcome for treatment and mixed control
         )
         
  
###
# Estimate the "sensitive" item with the DiM estimator
###

#using conventional control group (j=4)
m <- lm_robust(y_b ~ treat_b, data = df) %>%
  tidy() %>%
  slice(2) %>%
  mutate(model = "Conventional")

n <- lm_robust(y_b ~ treat_b, data = df) %>% nobs()

conventional <- cbind(m,n) %>%
  mutate(model = paste(model, "\nN = ", as.character(n)))


#using placebo control group (j=5)
m <- lm_robust(y_b_placebo ~ treat_b, data = df) %>%
  tidy() %>%
  slice(2) %>%
  mutate(model = "Placebo")

n <- lm_robust(y_b_placebo ~ treat_b, data = df) %>% nobs()

placebo <- cbind(m,n) %>%
  mutate(model = paste(model, "\nN = ", as.character(n)))

#using mixed approach by aprox half sampled from conventional and other from placebo control group
m <- lm_robust(y_b_mixed ~ treat_b, data = df) %>%
  tidy() %>%
  slice(2) %>%
  mutate(model = "Mixed")

n <- lm_robust(y_b_mixed ~ treat_b, data = df) %>% nobs()

mixed <- cbind(m,n) %>%
  mutate(model = paste(model, "\nN = ", as.character(n)))


# lets compare them visually
estimates <- rbind(mixed, placebo, conventional)

# order as we want  
estimates$model <- factor(estimates$model, levels = c("Mixed \nN =  3317", "Placebo \nN =  3310", "Conventional \nN =  3284")) 

#Create Figure 2 
estimates %>%
  ggplot() +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high,
                      y = model, shape = rev(model), linetype = rev(model))) +
  geom_vline(aes(xintercept = (1/6)), linetype= 2, alpha=.5) +
  labs(y="", x="Estimated prevalence", shape="", linetype = "",
     title = "") +
  theme_classic() +
  theme(legend.position = "none") +
  NULL

ggsave("output/fig_2.pdf", width = 5, height = 2.5)


# how large is the bias for each? (info for ms text)
estimates %>% mutate(bias = (1/6) - estimate)


#create sum stats table

df %>% 
  select(age, income, education, female, urban = urban_hukou, treatment = b_treatment,
         conventional = b_control, placebo = b_placebo, mixed = y_b_mixed) %>% 
  mutate(mixed = ifelse(!is.na(treatment), NA, mixed)) %>% 
  stargazer(., digits = 1, omit.summary.stat = c("p25", "p75"),
            out = "output/sum_stats.tex")

# Manually add a column for means from a genpop sample (Asian barometer Wave 4)
# That data must be applied for which can be done here: http://www.asianbarometer.org/data/data-release

#recodings to harmonize with levels from our survey

 #rio::import("..data/asb_china.dta") %>% 
 # transmute(age = ifelse(se3_2 %in% 0:23, 1, 
 #              ifelse(se3_2 %in% 24:30, 2,
 #              ifelse(se3_2 %in% 31:40, 3, 
 #              ifelse(se3_2 %in% 41:55, 4, 5)))),
 #        education = ifelse(se5==1, 1,
 #                   ifelse(se5 %in% 2:6, se5 - 1, 
 #                   ifelse(se5 %in% 7:8, 5,
 #                   ifelse(se5 %in% 9:10, 7,NA)))),
 #        urban = ifelse(level==2, 1, 0), 
 #        female = ifelse(se2==2, 1, 0) 
 # ) %>% 
 # summarise_all(mean, na.rm=TRUE)


# Create table of  freqency distributions of responses for each group
treatment <- df %>%
  group_by(b_treatment) %>%
  summarise(N=n()) %>%
  slice(1:6) %>%
  mutate(prop = round(N/sum(N), digits = 3),
         treatment = paste(N, " (", prop, ")", sep="")) %>%
  select(treatment)

control <- df %>%
  group_by(b_control) %>%
  summarise(N=n()) %>%
  slice(1:6) %>%
  mutate(prop = round(N/sum(1648), digits = 3),
         control = paste(N, " (", prop, ")", sep="")) %>%
  select(control)

placebo <- df %>%
  group_by(b_placebo) %>%
  summarise(N=n()) %>%
  slice(1:6) %>%
  mutate(prop = round(N/sum(N), digits = 3),
         placebo = paste(N, " (", prop, ")", sep="")) %>%
  select(placebo)

cbind(treatment, control, placebo) %>%
  mutate(items = 0:5) %>%
  select(items, treatment, control, placebo) %>%
  stargazer(summary = FALSE, out = "output/freq_dist.tex")


# Create results table
m1 <- lm_robust(y_b ~ treat_b, data = df)
m2 <- lm_robust(y_b_placebo ~ treat_b, data = df)
m3 <- lm_robust(y_b_mixed ~ treat_b, data = df)

texreg(list(m1, m2, m3), file = "output/reg_table.tex",  include.ci = FALSE,
       custom.model.names = c("Conventional", "Placebo", "Mixed"))


###
#Impact illustration using data from Blair, Coppock and Moore's (APSR 2020) metaanalysis 
### 

#their data and replication files can be accessed here: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/YUXHZT

# import and recode hypothesised bias direction and assign optimal control list design
# based on table 1 in the paper
df <- rio::import("data/apsr_2020.RDS") %>%
  mutate(hypo = ifelse(polarity=="Prediction: Underreporting", "Under-\nreporting", "Over-\nreporting"),
         optimal = ifelse(hypo=="Under-\nreporting" & list_est <= .25, "Placebo", 
                   ifelse(hypo=="Under-\nreporting" & list_est > .25 & list_est <= .5, "Mixed",
                   ifelse(hypo=="Over-\nreporting" & list_est <= .25, "Mixed", "Conventional"))))

#order
df$optimal <- factor(df$optimal, levels = c("Conventional", "Placebo", "Mixed"))

# share of optimal control list design in past list experiments. (info for ms text)
df %>% 
  group_by(optimal) %>% 
  summarise(N=n(),
            share= n()/264)

# create figure 1, vizualising recommended control list design by hypothesized sensitivity
df %>%  
  ggplot() +
  geom_jitter(aes(y= list_est, x= hypo, shape = optimal), height = 0, width = .3, alpha=.8) +
  ylab("Estimated prevalence") +
  geom_hline(aes(yintercept = .25), linetype= 2, alpha=.5) +
  geom_hline(aes(yintercept = .5), linetype= 2, alpha=.5) +
  xlab("") +
  labs(shape="Recommended control list:") +
  coord_flip() + 
  theme_classic() +
  theme(legend.position="top", legend.text=element_text(size=10), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  guides(shape = guide_legend(override.aes = list(size=2.5)))
  NULL

ggsave("output/fig_1.pdf", width = 5.5, height = 3)

#for which subfield are our recommendations most relevant? (info for ms text)
df %>% group_by(final_category, optimal) %>% summarise(n())
