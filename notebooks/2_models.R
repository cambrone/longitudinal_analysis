#########################################################################################################################################
# Title: Model Comparison
# Author: Andres Cambronero Sanchez
# Date: 09/30/2025
# Purpose:  
#   - Working through Frees book on longitudinal analysis chapters on fixed and random effect models
#   - Using Framingham Heart Study Data to better understand methods
#   - Exercise: What is relationship between systolic blood pressure and BPMEDS?
#########################################################################################################################################

set.seed(123)

#load libraries
library(car)

library(corrplot)
library(clubSandwich)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lattice)
library(lme4)
library(lmtest)
library(plm)
library(sandwich)
library(nlme)
library(sjPlot)
library(tibble)
library(tidyr)


#####################################################################################################################
# Introduction 
#####################################################################################################################
#  Using framington dataset, describe it.
# interesting things about this data:
#    - unbalanced. What are implications for FE and RE models?
#    - small number of obs per subjects. What implications for FE and RE?
#    - not all predictors are time varying or time invariante for each subject
#    - Focus on changes to coefficients, size of SE, different/ similar conclusions (interpretation) and model performance in sample and out
#   - Is AIC a good way to compare OLS, FE and RE? Yes https://uncdependlab.github.io/MLM_Tutorial/06_ModelSelection/model_comparison.html#comparing-relative-evidence-of-candidate-mlms
#   - Explain why this focus is importnat


# The main question is "what is the relationship between BPMEDS and SYSBP?" is deeper than initially looks.
# There isnt 1 single relationship in the data: There is the overall relationship between subjects and What there is within subjects where BPMEDS are used
# Show the EDA of these relatiosnhips.
# The research question is unclear. Is the interested in within subject effect or between subjects?
# it is important to separate these relationships because otherwise this coefficinet will end up being a weighted average of both effects which is uninterpretatble

#####################################################################################################################
#  Data Preparation 
#####################################################################################################################
#load dataset
df_frmg<-read.csv("~/Desktop/stats_projects/longitudinal_analysis/data/frmgham2.csv")
colnames(df_frmg) <- toupper(colnames(df_frmg))
df_frmg$RANDID = as.factor(df_frmg$RANDID)

#for this simple analysis, keep complete cases only
#variables chosen partly xbased on https://www.ahajournals.org/doi/10.1161/01.hyp.37.2.187
nrow(df_frmg)
df_frmg <- na.omit(df_frmg[,c("RANDID", 
                              "PERIOD",
                              "SYSBP",
                              "CIGPDAY",
                              "AGE",
                              "BMI",
                              "TOTCHOL",
                              "SEX",
                              'GLUCOSE',
                              "BPMEDS",
                              "EDUC")])
nrow(df_frmg)

length(unique(df_frmg$RANDID))

# number of records per subject
df_frmg %>%
  group_by(RANDID) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(count = n()/length(unique(df_frmg$RANDID))) 



#original 1=male, 2=female
df_frmg$SEX2 = ifelse(df_frmg$SEX==1, 0, 1)

# count number of observations per subject. 
# some subjects appear in period 1 and 3 when they only have 2 observations. The ordering doenst impact model estimatoion
n_periods <- df_frmg %>%
  group_by(RANDID) %>%
  summarise(period_count = n())

print(paste("% of subjects with 1 obs:", nrow(n_periods[n_periods$period_count==1,])/nrow(n_periods)))

nrow(df_frmg)
df_frmg <- left_join(df_frmg, n_periods, by="RANDID") 
nrow(df_frmg)

# Count if someone is always BPMEDS=0, =1 or changes throughout
temp <- df_frmg %>%
  group_by(RANDID) %>%
  summarise(mean_BPMED = mean(BPMEDS),
            max_BPMED = max(BPMEDS)) 

no_meds <- temp[temp$mean_BPMED==0,]
yes_meds <- temp[temp$mean_BPMED==1,]
some_meds <- temp[(temp$mean_BPMED!=0) & 
                  (temp$mean_BPMED!=1),]

print(paste('% of subjects with no BP meds: ', length(unique(no_meds$RANDID))/length(unique(temp$RANDID)))) 
print(paste('% of subjects with yes BP meds: ', length(unique(yes_meds$RANDID))/length(unique(temp$RANDID)))) 
print(paste('% of subjects with some BP meds: ', length(unique(some_meds$RANDID))/length(unique(temp$RANDID)))) 

#Merge with df_frmg
no_meds$no_meds <- 1
no_meds <- no_meds[,c("RANDID","no_meds")]

yes_meds$yes_meds <- 1
yes_meds <- yes_meds[,c("RANDID","yes_meds")]

some_meds$some_meds <- 1
some_meds <- some_meds[,c("RANDID","some_meds")]

nrow(df_frmg)
df_frmg<- left_join(df_frmg, no_meds, by="RANDID")
df_frmg<- left_join(df_frmg, yes_meds, by="RANDID")
df_frmg<- left_join(df_frmg, some_meds, by="RANDID")
nrow(df_frmg)

df_frmg$no_meds <- ifelse(is.na(df_frmg$no_meds), 0, df_frmg$no_meds)
df_frmg$yes_meds <- ifelse(is.na(df_frmg$yes_meds), 0, df_frmg$yes_meds)
df_frmg$some_meds <- ifelse(is.na(df_frmg$some_meds), 0, df_frmg$some_meds)

# is the direction always no meds -> yes meds? 87% yes
# for each subject choose their latest period
max_rows_by_group <- df_frmg %>%
  group_by(RANDID) %>%
  slice(which.max(PERIOD))

#for subjects with some meds, is their BPMEDS in last period ==1 
max_rows_by_group <- max_rows_by_group[max_rows_by_group$some_meds==1,]
mean(max_rows_by_group$BPMEDS)

# Create column with direction of BPMEDS
no_yes_meds<- max_rows_by_group[max_rows_by_group$BPMEDS==1,] #these subjects start with 0 then end in 1
no_yes_meds$no_yes_meds<- 1

no_yes_meds <- no_yes_meds[,c("RANDID","no_yes_meds")]

nrow(df_frmg)
df_frmg<- left_join(df_frmg, no_yes_meds, by="RANDID")
nrow(df_frmg)

df_frmg$no_yes_meds <- ifelse(is.na(df_frmg$no_yes_meds), 0, df_frmg$no_yes_meds)

df_frmg$yes_no_meds <- ifelse((df_frmg$some_meds==1) &
                                (df_frmg$no_yes_meds==0), 1, 0)

print(paste('no meds', mean(df_frmg$no_meds)))
print(paste('yes meds', mean(df_frmg$yes_meds)))
print(paste('no_yes meds', mean(df_frmg$no_yes_meds)))
print(paste('yes_no meds', mean(df_frmg$yes_no_meds)))

df_frmg$status_meds <-NA 
df_frmg$status_meds <- ifelse(df_frmg$no_meds==1,'no_change_meds', df_frmg$status_meds )
df_frmg$status_meds <- ifelse(df_frmg$yes_meds==1,'no_change_meds', df_frmg$status_meds)
df_frmg$status_meds <- ifelse(df_frmg$yes_no_meds==1,'yes_no_meds', df_frmg$status_meds)
df_frmg$status_meds <- ifelse(df_frmg$no_yes_meds==1,'no_yes_meds', df_frmg$status_meds)
df_frmg$status_meds <- factor(df_frmg$status_meds )


#####################################################################################################################
# Create Test sample
#####################################################################################################################
sample_test_id <-sample(df_frmg$RANDID, length(unique(df_frmg$RANDID))*0.1, replace=FALSE)

df_train <- df_frmg[!(df_frmg$RANDID %in% sample_test_id),]
df_test <- df_frmg[(df_frmg$RANDID %in% sample_test_id),]


#####################################################################################################################
# Within Training create chart with BPMED per group. Univariate view
#####################################################################################################################
# BPMEDS a) no meds entire time, b) yes meds entire time, c) no_yes_meds, d) yes_no_meds
summary_df <- df_train %>%
  group_by(status_meds, BPMEDS) %>%
  summarise(mean_sysbp = mean(SYSBP),
            count_id = n_distinct(RANDID),
            .groups = "drop")

# Calculate overall means for BPMEDS = 0 and 1
overall_df <- df_train %>%
  group_by(BPMEDS) %>%
  summarise(mean_sysbp = mean(SYSBP),
            count_id = n_distinct(RANDID),
            .groups = "drop") %>%
  mutate(status_meds = "overall")


summary_df <- bind_rows(summary_df, overall_df)


ggplot(summary_df, aes(x = factor(BPMEDS),
                       y = mean_sysbp, 
                       group = status_meds,
                       color = status_meds, 
                       size = count_id)) +
  geom_point() +
  geom_line(data = subset(summary_df, status_meds == "no_yes_meds"), linewidth=1, color="#1b9e77", linetype='dashed') +
  geom_line(data = subset(summary_df, status_meds == "yes_no_meds"), linewidth=1, color='#d95f02', linetype='dashed') +
  geom_line(data = subset(summary_df, status_meds == "no_change_meds"), linewidth=1, color="#7570b3", linetype='dashed') +
  geom_line(data = subset(summary_df, status_meds == "overall"), linewidth=2, color='black', linetype='solid') +
  scale_x_discrete(labels = c("0" = "BPMEDS=0", "1" = "BPMEDS=1"))+
  scale_color_manual(values = c(
    "no_yes_meds" = "#1b9e77",
    "yes_no_meds" = "#d95f02",
    "no_change_meds" = "#7570b3",
    "overall" = "black" )
    )
      
#####################################################################################################################
# create lists to score performance metrics
#####################################################################################################################
aic_list<-data.frame()
rmse_list<-data.frame()

#####################################################################################################################
# Model 1) OLS - Ignore the Repeated Measures (not recommended)
#####################################################################################################################
#   - If you could guarantee that there is no omitted variable bias (all factors of importance are present in model),
#      ols estimate for variable of interest would be unbiased. However, it coul still lead to incorrect inference because model assumes
#     indenpedent observations. As a result, SE would be too small. https://web.vu.lt/mif/a.buteikis/wp-content/uploads/2019/11/MultivariableRegression_4.pdf
#    Intuitively, using OLS, we build a model where we have less data than the model expects so standard errors are smaller than they should be. 
#     Model 1) OLS: SYSBP ~ cigs + time varying controls.
#              Interpretation 
#              Intime performance 
#              out of time performance
#              obtained standard errors follow OLS assumptions
#              diagnostics
#              standard errors can be adjusted with HAC standard errors https://www.fsb.miamioh.edu/lij14/672_s15.pdf

# set up model 
target <- "SYSBP"
predictors <-c('CIGPDAY', 'AGE', 'BMI', 'TOTCHOL', 'GLUCOSE', 'BPMEDS') # SEX and EDUC left out for now as they are time invariant
model_formula<- as.formula(paste(target, " ~ ", paste(predictors, collapse= "+")))

# a) run model
model_ols <- lm(model_formula, data = df_train)
summary(model_ols) # doesn't account for unobserved heterogeneity of subjects or repeated measures
                   # coefficients represent the average change in BPMEDS of records with BPMEDS=0 to BPMEDS=1
                   # cig has a small positive effect

# BPMEDS is positive in OLS model on whole training sample. Is relationship the same if we investigate those who switch
model_ols2 <- lm(model_formula, data = df_train[df_train$some_meds==1,]) # if we focus on those who switch, the relatioship is negative
summary(model_ols2) 

# b-i) Diagnostics: normally distributed errors
qqnorm(model_ols$residuals, ylab = 'Residuals')
qqline(model_ols$residuals)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
df_train$pred_model_ols <- predict(model_ols, newdata = df_train)
ggplot(df_train, aes(x = df_train$pred_model_ols, 
                     y = model_ols$residuals)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm") 

## b-iii) Diagnostics: serial correlation - present
df_train$resid_model_ols <- df_train$SYSBP - df_train$pred_model_ols
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_ols")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_ols.1", "resid_model_ols.2", "resid_model_ols.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)
acf(resid(model_ols))

# results of diagnostics suggest the need for SE correction
coeftest(model_ols, vcovCR(model_ols, type = "CR1",cluster = df_train$RANDID))

# c) evaluate model performance - metrics
model_aic<-AIC(model_ols)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_ols", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

# d) evaluate model performance - in-sample v. out-sample
df_train$pred_model_ols<-predict(model_ols, newdata = df_train)
df_test$pred_model_ols<-predict(model_ols, newdata = df_test)

insample_rmse<- sqrt(mean((df_train$SYSBP - df_train$pred_model_ols)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_ols)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "ols", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- rbind(rmse_list, new_row)




####################################################
# Model 2) Fixed effects
####################################################
# What does FE try to solve? 
#   - removes effects of unobserved time invariant heterogeneity from estimate of variable of interest. decreases omitted variable bias
#   - it does not remove unobsered bias from time varying heterogeniety. Limited to what is included in model
#   - there are two ways to fit a fixed effect model. 
#        a) OLS w dummy for clusters (SE account for degrees of freedom)
#        b) OLS w demeaning (Do not SE account for degrees of freedom)
#   - plm uses demeaning and adjusts SE for degress of freedom
#     Model 2) Fixed: SYSBP ~ BPMEDS + time varying controls
#             Differences in estimate of BPMEDS and other variables due to bias present in model 1 due to ommitting unobserved time invariant variables
#             Since FX estimates the association within subject we really only use those subjects with BPMEDS variation to estimate coefficient
#             The standard errors can be adjusted with HAC standard errors. This would be a way to keep model coefficeints the same if model is informed by theory
# general interpretation of fixed effect model coefficient, "your coefficient estimate is the marginal change expected in Y after changing one unit of X, controlling for all time-invariant heterogeneity in your groups. " within subject
# It is the average effect of a predictor within subjects
# "https://jmummolo.scholar.princeton.edu/sites/g/files/toruqf3341/files/mummolo_peterson_2018.pdf unit fixed effects are often used since they eliminate all between-unit variation, producing an estimate of a variable’s average effect within units over time"
# a) run model
model_fx <- plm(model_formula,
                data = df_train, 
                index = c("RANDID"),
                model = "within") # estimate accounting for subject heterogeneity
                                  # standard errors adjusted for implicit dummy variables used in fixed effect
                                  # standard errors assume no heteroskedasticity, no serial correlation
                                  # single observation subjects are not used for estimation of coefficients or SE.
summary(model_fx)

# show that same model is obtained if singl obs subjects are excluded
model_fx2 <- plm(model_formula,
                 data = df_train[df_train$period_count!=1,], 
                 index = c("RANDID"),
                 model = "within")
summary(model_fx2)      #shows same results as when we exclude single observation subjects
                        # single observation subjects would not get a dummy variable
                      
df_train$pred_model_fx <- predict(model_fx)  #plm with normal dataframe uses observed as prediction for single observation subjects
newdata.p <- pdata.frame(df_train) # make pdata.framepdata.frame simply leaves these single observation subjects. no prediction
df_train$pred_model_fx2<-predict(model_fx2, newdata = newdata.p)

# Does model use observations where BPMEDS is time invariant within subject? No 
model_fx3 <- plm(model_formula,
                  data = df_train[df_train$some_meds==0,], 
                  index = c("RANDID"),
                  model = "within")

# Would a model with only BPMED change if we include or exclude observations where BPMEDS is time in-variant?
# If only BPMEDS is used, the fixed effect only takes into estimation the observations where BPMEDS changes overtime.
# Same results if we use the entire dataset or use only those where some_meds==1 
model_fx4 <- plm(SYSBP ~ BPMEDS ,
                 data = df_train, 
                 index = c("RANDID"),
                 model = "within")
summary(model_fx4)

model_fx5 <- plm(SYSBP ~ BPMEDS,
                 data = df_train[df_train$some_meds==1,], 
                 index = c("RANDID"),
                 model = "within")
summary(model_fx5)

# If other time varying predictors included in formula, you get a different estimate for BPMEDS if whole sample is used or only those with some_meds==1
model_fx6 <- plm(SYSBP ~ BPMEDS + AGE,
                 data = df_train, 
                 index = c("RANDID"),
                 model = "within")
summary(model_fx6)

model_fx7 <- plm(SYSBP ~ BPMEDS + AGE,
                 data = df_train[df_train$some_meds==1,], 
                 index = c("RANDID"),
                 model = "within")
summary(model_fx7)


#show that ols without intercept but using dummies for a small subset is the same as fixed effect
sample_id <-sample(df_train[df_train$some_meds==1,'RANDID'], 30, replace=FALSE)
df_fx_subset <- df_train[(df_train$RANDID %in% sample_id),]

model_fx_dummy <- lm(SYSBP ~ -1 + CIGPDAY + AGE + BMI + TOTCHOL + GLUCOSE + BPMEDS + RANDID, df_fx_subset)
summary(model_fx_dummy)

model_fx9 <-  plm(model_formula,
                  data = df_fx_subset, 
                  index = c("RANDID"),
                  model = "within")
summary(model_fx9)


# b-i) Diagnostics: normally distributed errors
qqnorm(model_fx$residuals, ylab = 'Residuals')
qqline(model_fx$residuals)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
ggplot(df_train, aes(x = df_train$pred_model_fx, 
                     y = model_fx$residuals)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm") 

# b-iii) Diagnostics: serial correlation - present, but with only 1 to 3 observations this might not be a big issue or at least hard to conclude decisively
df_train$resid_model_fx <- df_train$SYSBP - df_train$pred_model_fx
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_fx")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_fx.1", "resid_model_fx.2", "resid_model_fx.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)

# Adjust SE given clustered nature of data, heteroskedasticity and auto correlation
# from my understanding vcovCL (cluster robust standard errors will account for autocorrelation within cluster)
coeftest(model_fx, vcovCR(model_fx, type = "CR1",cluster = 'RANDID'))


# c) evaluate model performance - metrics
# AIC adjusted for PLM takes into account implicit dummis and only time varying predictors used in model
# taken from https://stackoverflow.com/questions/46186527/how-to-calculate-bic-and-aic-for-a-gmm-model-in-r-using-plm
aicbic_plm <- function(object, criterion) {
  sp = summary(object)
  
  if((class(object)[1]=="plm") | (class(object)[1]=="pggls")){
    u.hat <- residuals(sp) # extract residuals
    df <- cbind(as.vector(u.hat), attr(u.hat, "index"))
    names(df) <- c("resid", "Country", "Time")
    c = length(levels(df$Country)) # extract country dimension 
    t = length(levels(df$Time)) # extract time dimension 
    np = length(sp$coefficients[,1]) # number of parameters
    n.N = nrow(sp$model) # number of data
    s.sq  <- log( (sum(u.hat^2)/(n.N))) # log sum of squares
    
    # effect = c("individual", "time", "twoways", "nested"),
    # model = c("within", "random", "ht", "between", "pooling", "fd")
    
    # I am making example only with some of the versions:
    
    if (sp$args$model == "within" & sp$args$effect == "individual"){
      n = c
      np = np+n+1 # update number of parameters
    }
    
    if (sp$args$model == "within" & sp$args$effect == "time"){
      T = t
      np = np+T+1 # update number of parameters
    }
    
    if (sp$args$model == "within" & sp$args$effect == "twoways"){
      n = c
      T = t
      np = np+n+T # update number of parameters
    }
    aic <- round(       2*np  +  n.N * (  log(2*pi) + s.sq  + 1 ),1)
    bic <- round(log(n.N)*np  +  n.N * (  log(2*pi) + s.sq  + 1 ),1)
    
    if(criterion=="AIC"){
      names(aic) = "AIC"
      return(aic)
    }
    if(criterion=="BIC"){
      names(bic) = "BIC"
      return(bic)
    }
  }
}

model_aic<-aicbic_plm(model_fx, "AIC")
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_fx", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

# d) evaluate model performance - in-sample v. out-sample
# Since fixed effect models are used primarily for adjusting heterogeneity between subjects and usually do not want to make claim about subject sout of sample
# we wont make out of sample predictions
df_train$pred_model_fx <- predict(model_fx)

insample_rmse<- sqrt(mean((df_train$SYSBP - df_train$pred_model_fx)^2))
print(paste('insample rmse: ', insample_rmse))

new_row <- data.frame(method = "fx", 
                      insample_rmse = insample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)


#####################################################################################################################
# Model 3) Feasible Generalized Least Squares
#####################################################################################################################
# OLS and FX both assume a no autocorrelation and no heteroskedasticity.
# when errors dont follow assiumption, many people siply change the SE to make them robust, without changes to estimates
# However to me and the author of this blog, these residuals imply a mispecified model and a more flexible approach might be needed
#    -https://fxdiebold.blogspot.com/2015/10/the-hac-emperor-has-no-clothes.html

# These flexible approaches include FGLS. FGLS is an estimation method that tries to find an appropriate variance covariance matrix suitable for your model 
# The approach multiplies this prescpedific matrix by the residuals to (hopefully) get residuals without any signale left in them https://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch5.pdf
# With GLS we have to prespecify a variance covariance matrix
# Based on our understanding of the data, we could attempt an AR(1)
# from documentation, General FGLS is based on a two-step estimation process: first a model is estimated by, fixed effects (model = "within") or first differences (model = "fd"), then its residuals are used to estimate an error covariance matrix for use in a feasible-GLS analysis.
# https://search.r-project.org/CRAN/refmans/plm/html/pggls.html

# a) run model
model_fgls <- pggls(model_formula,
                    data = df_train, 
                    index = c("RANDID"),
                    model = "within")


# b-i) Diagnostics: normally distributed errors
df_train$pred_model_fgls <- df_train$SYSBP - model_fgls$residuals
qqnorm(model_fgls$residuals, ylab = 'Residuals')
qqline(model_fgls$residuals)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
ggplot(df_train, aes(x = df_train$pred_model_fgls, 
                     y = model_fgls$residuals)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm") 

# b-iii) Diagnostics: serial correlation - present, again small sample per cluster. maybe had to determine
df_train$resid_model_fgls <- df_train$SYSBP - df_train$pred_model_fgls
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_fgls")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_fgls.1", "resid_model_fgls.2", "resid_model_fgls.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)


# c) performance metrics


# d) evaluate model performance - in-sample v. out-sample
# Since fixed effect models are used primarily for adjusting heterogeneity between subjects and usually do not want to make claim about subject sout of sample
# we wont make out of sample predictions
insample_rmse<- sqrt(mean((df_train$SYSBP - df_train$pred_model_fgls)^2))
print(paste('insample rmse: ', insample_rmse))

new_row <- data.frame(method = "fgls", 
                      insample_rmse = insample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)


#####################################################################################################################
# Model 4) Random Effects Model (random intercepts) with time varying variables only
#####################################################################################################################
# What are random effect models?
# One of the main differences between FX and RF is that RF assumes that the subject specific intercepts are uncorrelated with the predictors in the model.
# If that assumption does not hold and RF is used, we will end up with biased estimates of the regressors
# To test for this assumption we use the hausman test. It basically compares the coefficients of a FX model with time varying variables and a FX with time varying variables
# If H0: Random effect is consistent, Ha: FX model is consistent. IF p<0.05 we reject null and go with FX model
# Basically we are testing whether the coeffieints are uncorrelated with errors. IF they change between model they are correlated and should use FX
#  - Main limitations of FX:
#         - cannot make conclusions about observations out of sample because we know nothing about those groups effects
#           In contrast in random effects, we can make an make new predictions and conclusions by assuming the new group is similar to the other groups since they are realizations from the same distribution. 
#           Since we have the estimated mean and estimated variance, we can describe this distribution
#           Additionally in FX, we are not estimating an overall mean (we only have a mean for each subject in the sample) and the average within subject variation
#         Other ideas for benefits of partial pool that come with RF https://khakieconomics.github.io/2016/11/09/Random-effects,-partial-pooling,-and-exchangability.html
#        Shrinkage relates to the random intercept and slopes only. It is the result of assumping that these observations comes from a shared distribution with a shared mean and variance parameters
#         - do not extrapolate model predictions to new observations
#         - cannot estimate coefficients for time invariant variables 
#         - we are estimating 1 effect for each subject regardless of size. In random effects, with partial pooling we get subject averages that deviate from overall average depending on how much data they have. 
#           This might be another reason to use random effects for this data
# With random effects you have fewer parameters to estimate becuase we dont estimate the effect, but the parameters that describe the distribution of random effect https://stats.stackexchange.com/questions/658926/easy-way-to-understand-the-difference-between-a-cluster-variable-and-a-random-va#:~:text=A%20fixed%20effect%20is%20an,levels%20within%20the%20same%20school.
# how are random effects estimated? https://econ.pages.code.wm.edu/407/notes/docs/panel_re.html

# a) run model
# this fits model with unstructured variance covariance matrix
# this model has only time varying covariates to compare results against fixed effects fit before
# https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
model_rf_intercept_timevar <-lmer(SYSBP ~ CIGPDAY + AGE + BMI + TOTCHOL + GLUCOSE + BPMEDS + (1| RANDID), 
                             data=df_train)
summary(model_rf_intercept_timevar) # in RF effect model, the fixed effect coefficient represents the average effect of variable across subjects 
confint(model_rf_intercept_timevar) # represents the average effect of X over Y when X changes across time and between countries by one unit https://libguides.princeton.edu/R-Panel
                                      # RF use partial pooling. This means that the estimates for RF are closer to the overall average when subjects have few observations (unlikely/extreme observations with little weight)
                                      # results will be something between completely pooled models (ols) or no pooled models (models built on each subject or FX)





# What are the random interepts for subjects with only 1 observation (about 18% of subjects)?
# these subjects still receive a subject specific estimate. It seems like we might not be able to distinguish subject specific intercerpt from the error
# https://stats.stackexchange.com/questions/242821/how-will-random-effects-with-only-1-observation-affect-a-generalized-linear-mixe
# Since its only 1 obs, these estimates are likely unstable. This is another good reason why we are not interested in interpreting the random effcts themselves
# keep these observations in the model still allow use to make more precise estimates of fixed effect coefficients and ariance parameters of the individual- and group-level regressions https://stats.stackexchange.com/questions/51620/can-i-run-a-glmm-model-when-i-have-one-observation-for-most-subjects
ids_one_obs = df_train[df_train$period_count==1, "RANDID"]
random_eff_df <- ranef(model_rf_intercept_timevar)$RANDID 
random_eff_df <- cbind(RANDID = rownames(random_eff_df), random_eff_df)

random_eff_df_one_obs <- random_eff_df[random_eff_df$RANDID %in% ids_one_obs,] 

#does model change a lot if we remove subjects with 1 obseration only? It does change, but the residual variance doesnt change much
model_rf_intercept_timevar2 <-lmer(SYSBP ~ CIGPDAY + AGE + BMI + TOTCHOL + GLUCOSE + BPMEDS + (1| RANDID), 
                             data=df_train[df_train$period_count!=1,])
summary(model_rf_intercept_timevar)


ranef(model_rf_intercept_timevar)$RANDID %>% head(5)
coef(model_rf_intercept_timevar)$RANDID %>% head(5)


#intraclass correlation
# https://www.theanalysisfactor.com/the-intraclass-correlation-coefficient-in-mixed-models/
icc_df <- as.data.frame(VarCorr(model_rf_intercept_timevar))
(icc_df[icc_df$grp=='RANDID','vcov'])/(sum(icc_df[icc_df$grp=='RANDID','vcov'], icc_df[icc_df$grp=='Residual','vcov']))

df_train$pred_model_rf_intercept_timevar <-predict(model_rf_intercept_timevar, newdata = df_train)
df_train$pred_model_rf_intercept_timevar_norandeff<-predict(model_rf_intercept_timevar, newdata = df_train, re.form = ~0)

df_train$resid_model_rf_intercept_timevar <- residuals(model_rf_intercept_timevar)

# b-i) Diagnostics: normally distributed errors
qqnorm(df_train$resid_model_rf_intercept_timevar, ylab = 'Residuals')
qqline(df_train$resid_model_rf_intercept_timevar)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
ggplot(df_train, aes(x = df_train$pred_model_rf_intercept_timevar, 
                     y = df_train$resid_model_rf_intercept_timevar)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm",, formula = y ~ poly(x, 3))

# b-iii) Diagnostics: serial correlation - present
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_rf_intercept_timevar")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_rf_intercept_timevar.1", "resid_model_rf_intercept_timevar.2", "resid_model_rf_intercept_timevar.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)

#these resuls again might warrant a robust standard error
se_model_rf_intercept_timevar<-vcovCR(model_rf_intercept_timevar, type = "CR1",cluster = df_train$RANDID)
sqrt(diag(se_model_rf_intercept_timevar))

# c) evaluate model performance - metrics
model_aic<-AIC(model_rf_intercept_timevar)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_rf_intercept_timevar", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

# d) evaluate model performance - in-sample v. out-sample
#  - to make predictions for new observations we extract the fixed effect part only
#  - this implies that the new prediction will have the 0 random effect (mean of assumed distribution)
# https://pmarchand1.github.io/ECL7102/notes_cours/12E-Mixed_models_Part2.html
df_test$pred_model_rf_intercept_timevar<-predict(model_rf_intercept_timevar, newdata = df_test, re.form = ~0)

insample_rmse<-sqrt(mean((df_train$SYSBP - df_train$pred_model_rf_intercept_timevar)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_rf_intercept_timevar)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "rf_intercept_timevar", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)

#####################################################################################################################
# Model 5) Random slopes and intercepts.
#####################################################################################################################
# The fixed effect model and the random intercept model provide different direction of coeffcient for BPMEDS.
# The fixed effect model can only estimate the coefficient for BPMEDS from those subjects whose value is time varying. Therefore it represents the average within subject effect from no meds to yes meds
# The random intercept model can estimate the coefficient for BPMEDS from all subjects, both those whose value changes and those who doesnt. 
# as a result the BPMEDS coefficient is positive, similar to ols. 
# I want to see if a random slope model can capture that there are two types of subjects in the data: those whose value of BPMEDS doesnt change and those who it does (no meds to yes meds)
# a) run model 
model_rf_interceptslope_timevar <-lmer(SYSBP ~ CIGPDAY + AGE + BMI + TOTCHOL + GLUCOSE + BPMEDS + (1 + BPMEDS| RANDID), data=df_train)
summary(model_rf_interceptslope_timevar)

# What is the random slope for subjects with only 1 observation? We get an subject specific intercept and slope)
ids_one_obs = df_train[df_train$period_count==1, "RANDID"]
random_eff_df <- ranef(model_rf_interceptslope_timevar)$RANDID 
random_eff_df <- cbind(RANDID = rownames(random_eff_df), random_eff_df)

random_eff_df_one_obs <- random_eff_df[random_eff_df$RANDID %in% ids_one_obs,]
mean(random_eff_df_one_obs$BPMEDS) # the mean slope for these subjects is negative

# What is the average random slope for subjects that are placed on meds and those who dont get placed on meds?
ids_some_meds = df_train[df_train$some_meds==1, "RANDID"]
ids_no_meds = df_train[df_train$some_meds!=1, "RANDID"]

random_eff_df <- ranef(model_rf_interceptslope_timevar)$RANDID 
random_eff_df <- cbind(RANDID = rownames(random_eff_df), random_eff_df)

random_eff_df_some_meds<- random_eff_df[random_eff_df$RANDID %in% ids_some_meds,]
random_eff_df_no_meds<- random_eff_df[random_eff_df$RANDID %in% ids_no_meds,]

mean(random_eff_df_some_meds$BPMEDS) # negative slope if you are placed on meds (shallower increase than overall effect)
mean(random_eff_df_no_meds$BPMEDS) # small positive slope if you are not on meds, maybe increase overtime as people age?
                                   # for these subjects the "random slope" isnt really a slope since their value doesnt change over time
                                   # the random slope actually a further specification that the effect of BPMEDS is different for each subject. I am interpreting it as allowing 
 

# b) evaluate model performance - metrics
model_aic<- AIC(model_rf_interceptslope_timevar)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_rf_interceptslope_timevar", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

ranef(model_rf_interceptslope_timevar)$RANDID %>% head(5)
coef(model_rf_interceptslope_timevar)$RANDID %>% head(5)

df_train$pred_model_rf_interceptslope_timevar <-predict(model_rf_interceptslope_timevar, newdata = df_train)
df_train$resid_model_rf_interceptslope_timevar <- residuals(model_rf_interceptslope_timevar)

# b-i) Diagnostics: normally distributed errors
qqnorm(df_train$resid_model_rf_interceptslope_timevar, ylab = 'Residuals')
qqline(df_train$resid_model_rf_interceptslope_timevar)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
ggplot(df_train, aes(x = df_train$pred_model_rf_interceptslope_timevar, 
                     y = df_train$resid_model_rf_interceptslope_timevar)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm",, formula = y ~ poly(x, 3))

# b-iii) Diagnostics: serial correlation - present
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_rf_interceptslope_timevar")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_rf_interceptslope_timevar.1", "resid_model_rf_interceptslope_timevar.2", "resid_model_rf_interceptslope_timevar.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)

#these resuls again might warrant a robust standard error
se_model_rf_interceptslope_timevar<-vcovCR(model_rf_interceptslope_timevar, type = "CR1",cluster = df_train$RANDID)
sqrt(diag(se_model_rf_interceptslope_timevar))


# d) evaluate model performance - in-sample v. out-sample
#  - to make predictions for new observations we extract the fixed effect part only
#  - this implies that the new prediction will have the 0 random effect (mean of assumed distribution)
# https://pmarchand1.github.io/ECL7102/notes_cours/12E-Mixed_models_Part2.html
df_test$pred_model_rf_interceptslope_timevar<-predict(model_rf_interceptslope_timevar, newdata = df_test, re.form = ~0)

insample_rmse<-sqrt(mean((df_train$SYSBP - df_train$pred_model_rf_interceptslope_timevar)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_rf_interceptslope_timevar)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "rf_interceptslope_timevar", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)


############################################################################
# # Model 6) Within and Between effect model (no random slope)
############################################################################
create_cols_rewb <- function(df) {
  means_df <- df %>%
    group_by(RANDID) %>%
    summarise(mean_CIGPDAY = mean(CIGPDAY),
              mean_AGE = mean(AGE),
              mean_BMI = mean(BMI),
              mean_TOTCHOL = mean(TOTCHOL),
              mean_GLUCOSE = mean(GLUCOSE),
              mean_BPMEDS = mean(BPMEDS))
  
  nrow(df)  
  df<-left_join(df, means_df, by= 'RANDID')
  nrow(df)  
  
  df$diff_CIGPDAY <- df$CIGPDAY - df$mean_CIGPDAY
  df$diff_AGE <- df$AGE - df$mean_AGE
  df$diff_BMI <- df$BMI - df$mean_BMI 
  df$diff_TOTCHOL <- df$TOTCHOL - df$mean_TOTCHOL
  df$diff_GLUCOSE <- df$GLUCOSE - df$mean_GLUCOSE
  df$diff_BPMEDS <- df$BPMEDS - df$mean_BPMEDS
  
  return (df)
}

# prep data for rewb model
df_train<- create_cols_rewb(df_train)
df_test<- create_cols_rewb(df_test)

# a) run model 
model_rewb_intercept_timevar <-lmer(SYSBP ~ diff_CIGPDAY + mean_CIGPDAY +
                                      diff_AGE + mean_AGE + 
                                      diff_BMI + mean_BMI + 
                                      diff_TOTCHOL + mean_TOTCHOL+
                                      diff_GLUCOSE + mean_GLUCOSE +
                                      diff_BPMEDS + mean_BPMEDS + 
                                      (1 | RANDID), data=df_train) # this allows for getting the within effect, sme as fixed effect 

# following articles recommednation we can check if the within effect is different from the between effect
# the model with random intercepts only, assumes that the within effect=the between effect
# this model doesnt make that assumption. 
# Anova shows the rewb fits better so the effects are different and we should use this model instead
anova(model_rewb_intercept_timevar, model_rf_intercept_timevar)


# b) evaluate model performance - metrics
model_aic<- AIC(model_rewb_intercept_timevar)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_rewb_intercept_timevar", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

# d) evaluate model performance - in-sample v. out-sample
#  - to make predictions for new observations we extract the fixed effect part only
#  - this implies that the new prediction will have the 0 random effect (mean of assumed distribution)
# https://pmarchand1.github.io/ECL7102/notes_cours/12E-Mixed_models_Part2.html
df_train$pred_model_rewb_intercept_timevar <-predict(model_rewb_intercept_timevar, newdata = df_train)
df_test$pred_model_rewb_intercept_timevar <-predict(model_rewb_intercept_timevar, newdata = df_test, re.form = ~0)

insample_rmse<-sqrt(mean((df_train$SYSBP - df_train$pred_model_rewb_intercept_timevar)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_rewb_intercept_timevar)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "rewb_intercept_timevar", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)

#####################################################################################################################
# Model 7) REWB with time varying predictors (intercepts and slopes)
#####################################################################################################################
# a) run model
model_rewb_interceptslopes_timevar <- lmer(SYSBP ~ diff_CIGPDAY + mean_CIGPDAY +
                                             diff_AGE + mean_AGE + 
                                             diff_BMI + mean_BMI + 
                                             diff_TOTCHOL + mean_TOTCHOL+
                                             diff_GLUCOSE + mean_GLUCOSE +
                                             diff_BPMEDS + mean_BPMEDS +
                                             (1 + diff_BPMEDS| RANDID), data = df_train)

# b) evaluate model performance - metrics
model_aic<- AIC(model_rewb_interceptslopes_timevar)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_rewb_interceptslopes_timevar", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

# d) evaluate model performance - in-sample v. out-sample
#  - to make predictions for new observations we extract the fixed effect part only
#  - this implies that the new prediction will have the 0 random effect (mean of assumed distribution)
# https://pmarchand1.github.io/ECL7102/notes_cours/12E-Mixed_models_Part2.html
df_train$pred_model_rewb_interceptslopes_timevar <-predict(model_rewb_interceptslopes_timevar, newdata = df_train)
df_test$pred_model_rewb_interceptslopes_timevar <-predict(model_rewb_interceptslopes_timevar, newdata = df_test, re.form = ~0)

insample_rmse<-sqrt(mean((df_train$SYSBP - df_train$pred_model_rewb_interceptslopes_timevar)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_rewb_interceptslopes_timevar)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "rewb_interceptslopes_timevar", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)

#####################################################################################################################
# Model 8) REWB with time variant and invariany predictors (intercepts and slopes only)
#####################################################################################################################
model_rewb_interceptslopes_all <- lmer(SYSBP ~ diff_CIGPDAY + mean_CIGPDAY +
                                   diff_AGE + mean_AGE + 
                                   diff_BMI + mean_BMI + 
                                   diff_TOTCHOL + mean_TOTCHOL+
                                   diff_GLUCOSE + mean_GLUCOSE +
                                   diff_BPMEDS + mean_BPMEDS + 
                                   SEX + EDUC + (1 + diff_BPMEDS| RANDID), data = df_train)


summary(model_rewb_interceptslopes_all) # we interpret the fixed effects for a mixed model in the same way as a regression


ids_yes_meds  = df_train[df_train$yes_meds==1, "RANDID"]
ids_no_meds = df_train[df_train$no_meds==1, "RANDID"]
ids_no_yes_meds = df_train[df_train$no_yes_meds==1, "RANDID"]
ids_yes_no_meds = df_train[df_train$yes_no_meds==1, "RANDID"]

random_eff_df <- ranef(model_rewb_interceptslopes_all)$RANDID 
random_eff_df <- cbind(RANDID = rownames(random_eff_df), random_eff_df)

random_eff_df_yes_meds<- random_eff_df[random_eff_df$RANDID %in% ids_yes_meds,]
random_eff_df_no_meds<- random_eff_df[random_eff_df$RANDID %in% ids_no_meds,]
random_eff_df_no_yes_meds<- random_eff_df[random_eff_df$RANDID %in% ids_no_yes_meds,]
random_eff_df_yes_no_meds<- random_eff_df[random_eff_df$RANDID %in% ids_yes_no_meds,]

mean(random_eff_df_yes_meds$diff_BPMEDS)# negative slope if you are placed on meds (shallower increase than overall association)
var(random_eff_df_yes_meds$diff_BPMEDS)

mean(random_eff_df_no_meds$diff_BPMEDS) 
var(random_eff_df_no_meds$diff_BPMEDS) 

mean(random_eff_df_no_yes_meds$diff_BPMEDS) 
var(random_eff_df_no_yes_meds$diff_BPMEDS) 

mean(random_eff_df_yes_no_meds$diff_BPMEDS) 
var(random_eff_df_yes_no_meds$diff_BPMEDS) 

#b) evaluate model performance - metrics
model_aic<-AIC(model_rewb_interceptslopes_all)
print(paste('AIC: ', model_aic))

new_row <- data.frame(method = "model_rewb_interceptslopes_all", aic = model_aic)
aic_list<- rbind(aic_list, new_row)

df_train$pred_model_rewb_interceptslopes_all <-predict(model_rewb_interceptslopes_all, newdata = df_train)
df_train$resid_model_rewb_interceptslopes_all <- residuals(model_rewb_interceptslopes_all)

# b-i) Diagnostics: normally distributed errors
qqnorm(df_train$resid_model_rewb_interceptslopes_all, ylab = 'Residuals')
qqline(df_train$resid_model_rewb_interceptslopes_all)

# b-ii) Diagnostics: Constant variance - violated, fanning out pattern
ggplot(df_train, aes(x = df_train$pred_model_rewb_interceptslopes_all, 
                     y = df_train$resid_model_rewb_interceptslopes_all)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm",, formula = y ~ poly(x, 3))

# b-iii) Diagnostics: serial correlation - present
df_train_wide <- df_train[,c("RANDID","PERIOD","resid_model_rewb_interceptslopes_all")]
df_train_wide <- reshape(df_train_wide,
                         idvar = "RANDID",
                         timevar = "PERIOD",
                         direction = "wide")

corr_mat<-cor(df_train_wide[,c("resid_model_rewb_interceptslopes_all.1", "resid_model_rewb_interceptslopes_all.2", "resid_model_rewb_interceptslopes_all.3")],  use = "pairwise.complete.obs")
corrplot(corr_mat, type="upper", tl.col="black", tl.srt=45)

#these resuls again might warrant a robust standard error
se_model_rewb_interceptslopes_all<-vcovCR(model_rewb_interceptslopes_all, type = "CR1",cluster = df_train$RANDID)
sqrt(diag(se_model_rewb_interceptslopes_all))


# d) evaluate model performance - in-sample v. out-sample
#  - to make predictions for new observations we extract the fixed effect part only
#  - this implies that the new prediction will have the 0 random effect (mean of assumed distribution)
# https://pmarchand1.github.io/ECL7102/notes_cours/12E-Mixed_models_Part2.html
df_test$pred_model_rewb_interceptslopes_all<-predict(model_rewb_interceptslopes_all, newdata = df_test, re.form = ~0)

insample_rmse<-sqrt(mean((df_train$SYSBP - df_train$pred_model_rewb_interceptslopes_all)^2))
print(paste('insample rmse: ', insample_rmse))

outsample_rmse<-sqrt(mean((df_test$SYSBP - df_test$pred_model_rewb_interceptslopes_all)^2))
print(paste('outsample rmse: ', outsample_rmse))

new_row <- data.frame(method = "rewb_interceptslopes_all", 
                      insample_rmse = insample_rmse,
                      outsample_rmse = outsample_rmse)
rmse_list <- bind_rows(rmse_list, new_row)




####################################################################################
# Summary: Comparison of coefficients
####################################################################################
# ols with regular SE
# summary_model_ols <- as.data.frame(summary(model_ols)$coefficients)
# summary_model_ols <- cbind(predictor = rownames(summary_model_ols), summary_model_ols)
# rownames(summary_model_ols) <- 1:nrow(summary_model_ols)
# summary_model_ols$method_name <- c(rep("ols", nrow(summary_model_ols)))

# ols with robust SE
summary_model_ols_robustse <- coeftest(model_ols, vcovCR(model_ols, type = "CR1", cluster = df_train$RANDID))
summary_model_ols_robustse <- summary_model_ols_robustse[,] %>% 
  as_tibble() %>%
  mutate(predictor = rownames(summary_model_ols_robustse))
summary_model_ols_robustse$method_name <- c(rep("ols_robustse", nrow(summary_model_ols_robustse)))

# fx with regular SE
summary_model_fx <- as.data.frame(summary(model_fx)$coefficients)
summary_model_fx <- cbind(predictor = rownames(summary_model_fx), summary_model_fx)
rownames(summary_model_fx) <- 1:nrow(summary_model_fx)
summary_model_fx$method_name <- c(rep("fx", nrow(summary_model_fx)))
summary_model_fx <- summary_model_fx %>%
  rename("t value" = "t-value")

# random effects intercept only time varying
summary_model_rf_intercept_timevar <- as.data.frame(summary(model_rf_intercept_timevar)$coefficients)
summary_model_rf_intercept_timevar <- cbind(predictor = rownames(summary_model_rf_intercept_timevar), summary_model_rf_intercept_timevar)
rownames(summary_model_rf_intercept_timevar) <- 1:nrow(summary_model_rf_intercept_timevar)
summary_model_rf_intercept_timevar$method_name <- c(rep("rf_intercept_timevar", nrow(summary_model_rf_intercept_timevar)))

# random effects intercept and slope time varying
summary_model_rf_interceptslope_timevar <- as.data.frame(summary(model_rf_interceptslope_timevar)$coefficients)
summary_model_rf_interceptslope_timevar <- cbind(predictor = rownames(summary_model_rf_interceptslope_timevar), summary_model_rf_interceptslope_timevar)
rownames(summary_model_rf_interceptslope_timevar) <- 1:nrow(summary_model_rf_interceptslope_timevar)
summary_model_rf_interceptslope_timevar$method_name <- c(rep("rf_interceptslope_timevar", nrow(summary_model_rf_interceptslope_timevar)))

# RFWB intercept time varying 
summary_model_rewb_intercept_timevar <- as.data.frame(summary(model_rewb_intercept_timevar)$coefficients)
summary_model_rewb_intercept_timevar <- cbind(predictor = rownames(summary_model_rewb_intercept_timevar), summary_model_rewb_intercept_timevar)
rownames(summary_model_rewb_intercept_timevar) <- 1:nrow(summary_model_rewb_intercept_timevar)
summary_model_rewb_intercept_timevar$method_name <- c(rep("rewb_intercept_timevar", nrow(summary_model_rewb_intercept_timevar)))

# RFWB intercept and slope time varying 
summary_model_rewb_interceptslopes_timevar <- as.data.frame(summary(model_rewb_interceptslopes_timevar)$coefficients)
summary_model_rewb_interceptslopes_timevar <- cbind(predictor = rownames(summary_model_rewb_interceptslopes_timevar), summary_model_rewb_interceptslopes_timevar)
rownames(summary_model_rewb_interceptslopes_timevar) <- 1:nrow(summary_model_rewb_interceptslopes_timevar)
summary_model_rewb_interceptslopes_timevar$method_name <- c(rep("rewb_interceptslope_timevar", nrow(summary_model_rewb_interceptslopes_timevar)))


# RFWB intercept and slope time varying and time invariant
summary_model_rewb_interceptslopes_all <- as.data.frame(summary(model_rewb_interceptslopes_all)$coefficients)
summary_model_rewb_interceptslopes_all <- cbind(predictor = rownames(summary_model_rewb_interceptslopes_all), summary_model_rewb_interceptslopes_all)
rownames(summary_model_rewb_interceptslopes_all) <- 1:nrow(summary_model_rewb_interceptslopes_all)
summary_model_rewb_interceptslopes_all$method_name <- c(rep("rewb_interceptslope_all", nrow(summary_model_rewb_interceptslopes_all)))


coeff_summary_df <- bind_rows(summary_model_ols_robustse,
                          summary_model_fx, 
                          summary_model_rf_intercept_timevar,
                          summary_model_rf_interceptslope_timevar,
                          summary_model_rewb_intercept_timevar,
                          summary_model_rewb_interceptslopes_timevar,
                          summary_model_rewb_interceptslopes_all)


coeff_summary_df <- coeff_summary_df %>%
  rename("estimate" = "Estimate",
         "se" = "Std. Error") %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

#keep only BPMEDS
coeff_summary_df <- coeff_summary_df[coeff_summary_df$predictor %in% c('BPMEDS', 'diff_BPMEDS', 'mean_BPMEDS'),]

# Plot 
predictor_order <- c("BPMEDS", "diff_BPMEDS", "mean_BPMEDS")
method_order<- c("ols_robustse", 
                 "fx",
                 "rf_intercept_timevar",
                 "rf_interceptslope_timevar",
                 "rewb_intercept_timevar",
                 "rewb_interceptslope_timevar",
                 "rewb_interceptslope_all")

coeff_summary_df <- coeff_summary_df %>%
  dplyr::mutate(
    predictor   = factor(predictor, levels = predictor_order),
    method_name = factor(method_name, levels = method_order)
  )

ggplot(coeff_summary_df, aes(x = estimate, y = predictor,
                             color = method_name)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(x = "Coefficient Estimate", y = "Predictor",
       color = "Method") +
  # Flip axis order if you want top-to-bottom
  scale_y_discrete(limits = rev(predictor_order)) +
  scale_color_discrete(limits = method_order) +
  theme_minimal()


####################################################################################
# Summary: Intime AIC
####################################################################################
aic_list$method <- substr(aic_list$method, 7, nchar(aic_list$method))

method_order <- c("ols", "fx", "rf_intercept_timevar",
                  "rf_interceptslope_timevar","rewb_intercept_timevar",
                  "rewb_interceptslopes_timevar","rewb_interceptslopes_all")

aic_list <- aic_list %>%
  dplyr::mutate(method = factor(method, levels = method_order))

ggplot(aic_list, aes(x=method,  y=aic))+ 
  geom_bar(stat = "identity", width=0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(label = round(aic,2) ), vjust=3, size=3, color="white") +
  coord_cartesian(ylim=c(30000,80000)) +
  ggtitle("AIC")
  

####################################################################################
# Summary: in sample RMSE v out of sample RMSE
####################################################################################
rmse_list <- rmse_list[rmse_list$method!='fgls',]
rmse_list_long <- rmse_list %>%
  pivot_longer(cols = c(insample_rmse, outsample_rmse), 
               names_to = "sample", 
               values_to = "Value")
rmse_list_long$sample <- substr(rmse_list_long$sample, 1, nchar(rmse_list_long$sample) - 5)

method_order <- c("ols", "fx", "rf_intercept_timevar",
                  "rf_interceptslope_timevar","rewb_intercept_timevar",
                  "rewb_interceptslopes_timevar","rewb_interceptslopes_all")
rmse_list_long <- rmse_list_long %>%
  dplyr::mutate(method = factor(method, levels = method_order))

ggplot(rmse_list_long, aes(x = method, y = Value, fill = sample)) +
  geom_bar(stat = "identity", 
           width=0.4,
           position = position_dodge()) +
  labs(x = "Method", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = round(Value,2), y = Value + 1),
            position = position_dodge(width = 0.9),
            vjust = -0.2,
            size = 3,
            angle =45) +
  coord_cartesian(ylim=c(0,25)) +
  ggtitle("RMSE")






