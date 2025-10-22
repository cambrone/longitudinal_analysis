#########################################################################################################################################
# Title: Exploratory Analysis of Framingham Heart Study 
# Author: Andres Cambronero Sanchez
# Date: 09/22/2025
# Purpose:  
#   - Working through Frees book on longitudinal analysis chapters on fixed and random effect models
#   - Using Framingham Heart Study Data to better understand methods
#   - Exercise: What is relationship between systolic blood pressure and taking blood pressure meds?
#########################################################################################################################################

set.seed(123)

#load libraries
library(dplyr)
library(ggplot2)
library(lattice)
library(corrplot)
library(gridExtra)

#load dataset
frmg_df<-read.csv("~/Desktop/stats_projects/longitudinal_analysis/data/frmgham2.csv")
colnames(frmg_df) <- toupper(colnames(frmg_df))

#####################################################################################################################
# Data cleaning
#####################################################################################################################
#original 1=male, 2=female
frmg_df$SEX2 = ifelse(frmg_df$SEX==1, 0, 1)

#####################################################################################################################
# EDA - Basic Checks
#####################################################################################################################
length(unique(frmg_df$RANDID))

# number of records per subject
frmg_df %>%
  group_by(RANDID) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(count = n()/length(unique(frmg_df$RANDID))) 

# percent of BPMEDS in the sample? 
mean(frmg_df$BPMEDS, na.rm=TRUE)

# do people's smoking status change over time? Yes
bpmeds_status_ot<- frmg_df %>%
  group_by(RANDID) %>%
  summarise(bpmeds_status = mean(BPMEDS, na.rm=TRUE))

nrow(bpmeds_status_ot[(bpmeds_status_ot$bpmeds_status>0) &
                  (bpmeds_status_ot$bpmeds_status<1),])/nrow(bpmeds_status_ot)

#####################################################################################################################
# EDA - univariate checks 
#####################################################################################################################
vars <- c("SYSBP", "BPMEDS", "CURSMOKE", "CIGPDAY", "EDUC", "SEX2", "AGE", "BMI", "TOTCHOL", "GLUCOSE")
for(i in vars){
  cat('*****************', i ,'**********************\n')
  # missing counts and unique values by BPMEDS status
  print(paste("mising cnt non-BPMEDS: ", sum(is.na(frmg_df[frmg_df$BPMEDS==0,c(i)]))))
  print(paste("Unique values non-BPMEDS:", paste(sort(unique(frmg_df[frmg_df$BPMEDS==0,c(i)])), collapse = ", ")))
  
  print(paste("mising cnt BPMEDS: ", sum(is.na(frmg_df[frmg_df$BPMEDS==1,c(i)]))))
  print(paste("Unique values BPMEDS:",
              paste(sort(unique(frmg_df[frmg_df$BPMEDS==1,c(i)])), collapse = ", ")))
  
  # distributions and distribution by smoking status
  print(paste('min: ', min(frmg_df[, c(i)], na.rm=TRUE)))
  print(paste('mean: ',  mean(frmg_df[, c(i)], na.rm=TRUE)))
  print(paste('max: ', max(frmg_df[, c(i)], na.rm=TRUE)))
  print(paste('sd: ', sd(frmg_df[, c(i)], na.rm=TRUE)))
  hist(frmg_df[, c(i)],
       main = paste("Distribution of ", i),
       xlab = "Value",
       ylab = "Frequency",
       breaks = 15, # Specify 10 bins
       col = "lightblue",
       border = "darkblue")
  
  print(paste('min nonBPMEDS: ', min(frmg_df[frmg_df$BPMEDS==0,c(i)], na.rm=TRUE)))
  print(paste('mean nonBPMEDS: ',  mean(frmg_df[frmg_df$BPMEDS==0,c(i)], na.rm=TRUE)))
  print(paste('max nonBPMEDS: ', max(frmg_df[frmg_df$BPMEDS==0,c(i)], na.rm=TRUE)))
  print(paste('sd nonBPMEDS: ', sd(frmg_df[frmg_df$BPMEDS==0,c(i)], na.rm=TRUE)))
  
  print('')
  
  print(paste('min BPMEDS: ', min(frmg_df[frmg_df$BPMEDS==1, c(i)], na.rm=TRUE)))
  print(paste('mean BPMEDS: ',  mean(frmg_df[frmg_df$BPMEDS==1, c(i)], na.rm=TRUE)))
  print(paste('max BPMEDS: ', max(frmg_df[frmg_df$BPMEDS==1, c(i)], na.rm=TRUE)))
  print(paste('sd BPMEDS: ', sd(frmg_df[frmg_df$BPMEDS==1, c(i)], na.rm=TRUE)))
  
  minx <- min(frmg_df[frmg_df$BPMEDS==0,c(i)], frmg_df[frmg_df$BPMEDS==1,c(i)], na.rm = TRUE)
  maxx <- max(frmg_df[frmg_df$BPMEDS==0,c(i)], frmg_df[frmg_df$BPMEDS==1,c(i)], na.rm = TRUE)
  
  # basic histogram
  hist(frmg_df[frmg_df$BPMEDS==0,c(i)],
       main=paste("Distribution of ", i),
       xlab="", 
       ylab="", 
       col = rgb(0, 0, 1, alpha = 0.5),
       xlim = c(minx, maxx))
  
  hist(frmg_df[frmg_df$BPMEDS==1,c(i)],
       xlab="", 
       ylab="",
       col = rgb(0, 1, 0, alpha = 0.5),
       add=TRUE)
  
  legend("topright", legend=c("non-BPMEDS", "BPMEDS"), fill=c(rgb(0, 0, 1, alpha = 0.5), 
                                                              rgb(0, 1, 0, alpha = 0.5)))
  
}

#####################################################################################################################
# EDA bivariate checks
#####################################################################################################################
vars <- c( "BPMEDS", "CIGPDAY", "EDUC", "SEX2", "AGE", "BMI", "TOTCHOL", "GLUCOSE")
complete_frmg_df <- frmg_df[complete.cases(frmg_df[vars]),]

for(i in vars){
  cat('*****************', i ,'**********************\n')
  # relationship with outcome
  # Calculate the quantile breaks
  breaks <- unique(quantile(frmg_df[,c(i)], probs = seq(0, 1, by = 1/10), na.rm = TRUE))

  # Cut the data into quantile-based bins
  frmg_df$binned_data <- cut(frmg_df[,c(i)], breaks = c(breaks, Inf), include.lowest = TRUE, labels = FALSE)
  frmg_df$binned_data <- factor(frmg_df$binned_data)
  
  boxplot(frmg_df$SYSB ~ frmg_df$binned_data, 
          data = frmg_df[!is.na(i)],
          main = "Boxplots by Groups of 10",
          xlab = paste('bkt', i),
          ylab = "SYSBP")
  
  # Partial regression plots
  var_interest <- vars[vars == i]
  vars_reduced <- vars[vars != i]
  
  form_var_intest <- as.formula(paste(var_interest, " ~ ", paste(vars_reduced, collapse= "+")))
  print(form_var_intest)
  resid_var_interest <-resid(lm(data=complete_frmg_df, form_var_intest))
  
  form_target <-as.formula(paste("SYSBP ~ ", paste(vars_reduced, collapse= "+")))
  print(form_target)
  resid_target <-resid(lm(data=complete_frmg_df, form_target))
  
  gg2<- ggplot(data=NULL) + 
    geom_point(aes(x=resid_var_interest, y=resid_target),size = 1, alpha = 0.1)+
    labs(title = paste("Partial Reg. Plot SYSBP v. ", i, " corr:", cor(resid_target, resid_var_interest)),
         x = paste("residual ",i),
         y = "residual SYSBP")
  
  print(gg2)
}

corr_mat<- cor(frmg_df[,vars],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", order="hclust", tl.col="black", tl.srt=45)


#####################################################################################################################
# EDA - longitudinal structure
#####################################################################################################################
#Individual SYSBOP trajectories
n_sample <- 100
sample_id <- sample(frmg_df[frmg_df$BPMEDS==0, "RANDID"], n_sample)
non_BPMEDS <- frmg_df[(frmg_df$RANDID %in% sample_id),]
non_BPMEDS_plot<- ggplot(non_BPMEDS, aes(x = non_BPMEDS$PERIOD, 
                        y = non_BPMEDS$SYSB, 
                        group = non_BPMEDS$RANDID, 
                        color = factor(non_BPMEDS$BPMEDS))) +
                    geom_line() +
                    geom_point() +
                    ylim(c(0, 210))+
                    labs(title = paste("SYSBP for",n_sample,"nonBPMEDS"),
                         x = "Period",
                         y = "SYSBP",
                         color = "Subject") +
                    theme_minimal()  + theme(legend.position="none")

sample_id <- sample(frmg_df[frmg_df$BPMEDS==1, "RANDID"], n_sample)
BPMEDS_df <- frmg_df[(frmg_df$RANDID %in% sample_id),]
BPMEDS_plot<- ggplot(BPMEDS_df, aes(x = BPMEDS_df$PERIOD, 
                                 y = BPMEDS_df$SYSB, 
                                 group = BPMEDS_df$RANDID, 
                                 color = factor(BPMEDS_df$BPMEDS))) +
              geom_line() +
              geom_point() +
              ylim(c(0, 210))+
              labs(title =  paste("SYSBP for ",n_sample,"BPMEDS"),
                   x = "Period",
                   y = "SYSBP",
                   color = "Subject") +
              theme_minimal() + theme(legend.position="none")

grid.arrange(non_BPMEDS_plot, BPMEDS_plot, ncol=2)
