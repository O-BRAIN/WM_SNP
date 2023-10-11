#load packages
{library(lme4)       # for gmler()
library(car)         # for Anova()
library(sjPlot)      # for plot_model()
library(readxl)      # for read_excel()
library("lsmeans")   # for posthoc tests
library(performance) # for model evaluations
library(ggplot2)     # for plotting
library(tidyverse)   # for re-arranging data
library(broom.mixed) # for printing results tables
library(report)      # for reporting R environment
  
  #for saving Anova() output into APA-style table
  library(officer)
  library(jtools)
  library(huxtable)
}

#check R background/environment
sessionInfo()                                     # should be R version 4.2.0 (2022-04-22 ucrt)

setwd("Work/WM_SNP/Source_data/")                 # set your working directory to where the data is stored
data <- read_excel("~/Work/WM_SNP/Source_data/BEDOB_WC_GREADT_long.xlsx")   #read in data


#some checks bevor we start ----
length(unique(data$ID))                           # N = 324    
#check data 
{sublist = unique(data$ID)

check = matrix(ncol = 2, nrow = length(sublist))
library(tidyverse)
for (k in 1:length(sublist)){ 
  filter= filter(data,  data$ID == sublist[[k]])
  check[k,1] = sublist[k]
  check[k,2] = nrow(filter)
}

check = as.data.frame(check)
check2 = check[c(which(check$V2 < 128)), ] #check2 and check 3 should be 0 obs., as they indicate whether subjects had less or more than 128 trials, which should not be the case
check3 = check[c(which(check$V2 > 128)), ]
}

# format data 
{
data$ID = as.factor(data$ID)
data$COMT = as.factor(data$COMT_rs4680)
data$DRD2 = as.factor(data$DRD2)
data$Gender = as.factor(data$Gender)
data$Age = as.numeric(data$Age)
data$BMI = as.numeric(data$BMI)
data$zDFS = scale(as.numeric(data$DFS))
data$IQ = as.numeric(data$IQ)
data$BED = as.factor(data$BED)
data$project = as.factor(data$project)
data$condition = as.factor(data$condition)
data$correct = as.factor(data$correct)
data$WM_concentration = as.numeric(data$WM_concentration)
data$WM_tiredness = as.numeric(data$WM_tiredness)
data$DARPP = as.factor(data$DARPP)
data$zAge      <- as.numeric(scale(data$Age))
data$zBMI      <- as.numeric(scale(data$BMI))
data$zIQ       <- as.numeric(scale(data$IQ))
data$zWM_conc  <- as.numeric(scale(data$WM_concentration))
data$zWM_tired <- as.numeric(scale(data$WM_tiredness))
}

# check if there are participants who performed with less than 50% on each condition
{
  sublist = unique(data$ID)
  perform = as.data.frame(matrix(ncol = 3, nrow = length(sublist)*4))
  colnames(perform) = c("ID","condition","meanACC")
  
  #filter according to condition (ignore = 0, m1 = 1, update = 2, m2 = 3) and ID
  subcount = nrow(perform)
  fill = subcount - 1
  for (k in 1:length(sublist)){
    for (j in 0:3){
      filter= filter(data, data$condition == j, data$ID == sublist[[k]])
      meanACC = sum(filter$correct == 1,na.rm=TRUE)/32
      
      #fill in values for filter condition
      perform$condition[subcount-fill+j] = as.character(filter$condition[1])
      perform$meanACC[subcount-fill+j] = meanACC
      perform$ID[subcount-fill+j] = as.character(sublist[k])
    }
    fill = fill - 4
  }
  
  lowperform = perform$ID[which(perform$meanACC < 0.5)]  #203; 230; 235 all have < 50% in all 4 conditions
  lowperform
}  
#remove low performing subjects
{
  data = data[-c(which(data$ID == "203")),]
  data = data[-c(which(data$ID == "230")),]
  data = data[-c(which(data$ID == "235")),]
  data = data[-c(which(data$ID == "542")),] #reported wrong button use 
}

#set miss trials to incorrect and remove trials with 2ms < RT < 0.2ms
data$correct[data$correct == 666] <- 0
data = data[data$RT > 0.2, ] #false alarm
data = data[data$RT < 2, ]   #miss trials

#check if classes are right now
classes <- lapply(data, class)
classes

#set contrasts
options(contrasts = c("contr.sum","contr.poly")) 
#https://faculty.nps.edu/sebuttre/home/r/contrasts.html

# check contrasts - is it set to -1 -1 -1 (poly because there is no "real" reference level)?
contrasts(data$condition)  

#sample descriptives
{
  #overall 
  {
  keep = c("ID","Age","BMI","IQ", "DFS", "Gender", "project")
  data_wide = data[ , (names(data) %in% keep)]
  data_wide = unique(data_wide)
  
  #populate table
  
  mean(data_wide$BMI) # 26.38158
  sd(data_wide$BMI)   # 6.347586
  min(data_wide$BMI)  # 17.51
  max(data_wide$BMI)  # 45.54
  
  length(which(is.na(data_wide$Age)==TRUE)) #how many NA in Age
  mean(data_wide$Age) # 26.93036
  sd(data_wide$Age)   # 6.816784
  min(data_wide$Age)  # 12.16667
  max(data_wide$Age)  # 49.75
  
  length(which(is.na(data_wide$IQ)==TRUE)) #how many NA in IQ
  mean(as.numeric(data_wide$IQ))  # 105.4156
  sd(as.numeric(data_wide$IQ))    # 10.60545
  min(as.numeric(data_wide$IQ))   # 71
  max(as.numeric(data_wide$IQ))   # 122
  
  length(unique(data$ID[(which(data$Gender == 0))])) #how many male?
  
  length(which(is.na(data_wide$DFS)==TRUE)) #how many NA in DFS
  mean(na.omit(data_wide$DFS))
  sd(na.omit(data_wide$DFS))
  min(na.omit(data_wide$DFS))
  max(na.omit(data_wide$DFS))
  }
  
  #per project 
  {
    # for BEDOB
    {
    data_wide2 = data_wide[data_wide$project == "BEDOB",]
    
    mean(data_wide2$BMI)
    sd(data_wide2$BMI)
    min(data_wide2$BMI)
    max(data_wide2$BMI)
    
    length(which(is.na(data_wide2$Age)==TRUE)) #how many NA in Age
    mean(data_wide2$Age)
    sd(data_wide2$Age)
    min(data_wide2$Age)
    max(data_wide2$Age)
    
    length(which(is.na(data_wide2$IQ)==TRUE)) #how many NA in IQ
    mean(as.numeric(data_wide2$IQ))
    sd(as.numeric(data_wide2$IQ))
    min(as.numeric(data_wide2$IQ))
    max(as.numeric(data_wide2$IQ))
    
    length(which(data_wide2$Gender == 0)) #how many female?
    
    length(which(is.na(data_wide2$DFS)==TRUE)) #how many NA in DFS
    mean(na.omit(data_wide2$DFS))
    sd(na.omit(data_wide2$DFS))
    min(na.omit(data_wide2$DFS))
    max(na.omit(data_wide2$DFS))
    }
    # for GREADT
    {
      data_wide2 = data_wide[data_wide$project == "GREADT",]
      
      mean(data_wide2$BMI)
      sd(data_wide2$BMI)
      min(data_wide2$BMI)
      max(data_wide2$BMI)
      
      length(which(is.na(data_wide2$Age)==TRUE)) #how many NA in DFS
      mean(data_wide2$Age)
      sd(data_wide2$Age)
      min(data_wide2$Age)
      max(data_wide2$Age)
      
      length(which(is.na(data_wide2$IQ)==TRUE)) #how many NA in DFS
      mean(as.numeric(data_wide2$IQ))
      sd(as.numeric(data_wide2$IQ))
      min(as.numeric(data_wide2$IQ))
      max(as.numeric(data_wide2$IQ))
      
      length(which(data_wide2$Gender == 1)) #how many female?
      
      length(which(is.na(data_wide2$DFS)==TRUE)) #how many NA in DFS
      mean(na.omit(data_wide2$DFS))
      sd(na.omit(data_wide2$DFS))
      min(na.omit(data_wide2$DFS))
      max(na.omit(data_wide2$DFS))
    }
    # for WORMCRI
    {
      data_wide2 = data_wide[data_wide$project == "WORMCRI",]
      
      mean(data_wide2$BMI)
      sd(data_wide2$BMI)
      min(data_wide2$BMI)
      max(data_wide2$BMI)
      
      length(which(is.na(data_wide2$Age)==TRUE)) #how many NA in DFS
      mean(data_wide2$Age)
      sd(data_wide2$Age)
      min(data_wide2$Age)
      max(data_wide2$Age)
      
      length(which(is.na(data_wide2$IQ)==TRUE)) #how many NA in DFS
      mean(as.numeric(data_wide2$IQ))
      sd(as.numeric(data_wide2$IQ))
      min(as.numeric(data_wide2$IQ))
      max(as.numeric(data_wide2$IQ))
      
      length(which(data_wide2$Gender == 1)) #how many female?
      
      length(which(is.na(data_wide2$DFS)==TRUE)) #how many NA in DFS
      mean(na.omit(data_wide2$DFS))
      sd(na.omit(data_wide2$DFS))
      min(na.omit(data_wide2$DFS))
      max(na.omit(data_wide2$DFS))
    }
  }
  }

# clear workspace
rm(list= ls()[!(ls() %in% c('data'))])

# WM gating x BMI? ----
#select best model first
#select data, so that model comparison can work (because there is no AAratio in BEDOB, and then BED only has one level if)
{
keep = c("ID", "correct", "condition", "zBMI", "zAge","zIQ","BED","project","Gender","zWM_tired","zWM_conc","zDFS", "COMT", "DRD2", "DARPP","DRD2_rs6277") 
data1 = data[ ,(names(data) %in% keep)]

#model comparison

model_BMI_full <- glmer(correct ~ condition * zBMI + zIQ + zWM_tired + zWM_conc + Gender + BED + project + zDFS + zAge + (0+condition|ID), 
                   family = "binomial", 
                   nAGQ = 1,
                   data = na.omit(data1), 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_BMI_full, type = 3)
#BED, project, zDFS, zAge insignificant

model_BMI_reduced <- glmer(correct ~ condition * zBMI + zIQ + zWM_tired + zWM_conc  + Gender + (0+condition|ID), 
                   family = "binomial", 
                   nAGQ = 1,
                   data = na.omit(data1), 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

anova(model_BMI_full, model_BMI_reduced)
#AIC & BIC model reduced smaller --> go with reduced model
}
#main analysis on full sample with best model:
{
model_BMI <- glmer(correct ~ condition * zBMI + zIQ + zWM_tired + zWM_conc + Gender + (1|ID), 
                             family = "binomial", 
                             nAGQ = 1,
                             data = data, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_BMI, type = 3)
summary(model_BMI)
tab_model(model_BMI)

# plot
{
  keep = c("ID","correct","BMI","condition", "zIQ", "zWM_tired", "Gender", "zWM_conc")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_BMI)
  plotdata = unique(plotdata)
  
  ggplot(data = plotdata, aes(BMI, fitted, colour = condition)) +
    geom_point(alpha = .5, size = 3) +
    scale_color_manual(values = c("black","grey30","grey45","grey60"), 
                       labels = c("ignore", "m1", "update", "m2")) +
    geom_smooth(method = lm) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    labs(y = "probability correct") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial"))
    
  
}
}


# COMT x DRD2 x BMI? ----
#model comparisons on subset (excluding AAratio for na.omit to work)
{
  model_COMT_Taq1A_BMI_full <- glmer(correct ~ condition * COMT * DRD2 * zBMI + zIQ + Gender + zWM_tired + zWM_conc + zAge + project + BED + zDFS + (1|ID), 
                        family = "binomial", 
                        nAGQ = 1,
                        data = na.omit(data1), 
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_COMT_Taq1A_BMI_full, type = 3)

model_COMT_Taq1A_BMI_reduced <- glmer(correct ~ condition * COMT * DRD2 * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                                      family = "binomial", 
                                      nAGQ = 1,
                                      data = na.omit(data1), 
                                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

anova(model_COMT_Taq1A_BMI_full,model_COMT_Taq1A_BMI_reduced)
}
#main model on full dataset, 0+condition doesn't work... model is singular
{
  model_COMT_Taq1A_BMI <- glmer(correct ~ condition * COMT * DRD2 * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                                        family = "binomial", 
                                        nAGQ = 1,
                                        data = data, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))  
  
Anova(model_COMT_Taq1A_BMI, type = 3)
tab_model(model_COMT_Taq1A_BMI)
}
#save into APA-table
{
  anova_results <- Anova(model_COMT_Taq1A_BMI, Type = 3)
  anova_df <- data.frame(
    Term = rownames(anova_results),
    Chisq = anova_results$Chisq,
    Df = anova_results$Df,
    p_value = anova_results$Pr(>Chisq)
  )
  anova_df$Chisq <- round(anova_df$Chisq, 3)
  anova_df$p_value <- sapply(anova_df$p_value, function(p) {
    if (p < .001) {
      return("< .001")
    } else {
      return(sprintf("%.3f", p))
    }
  })
  
  hux <- as_hux(anova_df)
  hux <- set_bold(hux, 1, everywhere, TRUE)   # Make header row bold
  hux <- set_bottom_border(hux, 1, everywhere, 1)   # Border below header
  
  quick_docx(hux, file = "apa_table.docx")  
}

#check if model assumptions are met
{
#collinearity of the predictors
x = check_collinearity(model_COMT_Taq1A_BMI)  
plot(x) # all main effects below 5 

plot_model(model_COMT_Taq1A_BMI, type = "diag") #looks good
}

#how to make sense of 3-way interaction condition x BMI x DRD2
#post hoc test
{posthoc1 <- emtrends(model_COMT_Taq1A_BMI, ~ DRD2 | condition, var = "zBMI")
pairs(posthoc1) #only update sig (p = 0.0017)
#extract trends 
posthoc1
}

# plot
{
keep = c("ID","correct","BMI","condition", "COMT","DRD2", "zIQ", "zWM_tired", "Gender", "zWM_conc", "DSBack")
plotdata = data[ , (names(data) %in% keep)]
plotdata = na.omit(plotdata)
plotdata$fitted = fitted(model_COMT_DRD2_BMI)
plotdata = unique(plotdata)

ggplot(data = plotdata, aes(BMI, fitted, colour = DRD2, facets = condition)) +
        geom_point(alpha = .3) +
        scale_color_manual(values = c("deepskyblue4", "darkgoldenrod2")) +
        geom_smooth(method = lm) +
        theme_classic() +
        facet_grid(.~ condition,
                   labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))

ggplot(data = plotdata, aes(BMI, fitted, colour = condition)) +
  geom_point(alpha = .5, size = 3) +
  scale_color_manual(values = c("black","grey30","grey45","grey60")) +
  geom_smooth(method = lm) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
  theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
  labs(y = "probability correct") +
  theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
  theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) 
  

#each condition on its own
for (var in unique(plotdata$condition)){
  ggplot(data = plotdata[plotdata$condition==var,], aes(BMI, fitted, colour = DRD2)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = .5, size = 3) +
    scale_color_manual(values = c("deepskyblue4","darkgoldenrod2")) +
    scale_y_continuous(limits = c(0.5, 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    labs(y = "probability correct") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial"))
  +
    facet_wrap(.~ condition,
               labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
}

}
        
# DARPP x BMI ----
#model comparison on data1
{
  model_DARPP_BMI_full <- glmer(correct ~ DARPP * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + zAge + project + BED + zDFS + (1|ID), 
                                family = "binomial", 
                                nAGQ = 1,
                                data = na.omit(data1), 
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_DARPP_BMI_full, type = 3)
  model_DARPP_BMI_reduced <- glmer(correct ~ DARPP * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + (1|ID), 
                                family = "binomial", 
                                nAGQ = 1,
                                data = na.omit(data1), 
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  anova(model_DARPP_BMI_full,model_DARPP_BMI_reduced)
  
  }
# main model on full data
{
model_DARPP_BMI <- glmer(correct ~ DARPP * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + (1|ID), 
                           family = "binomial", 
                           nAGQ = 1,
                           data = data, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
Anova(model_DARPP_BMI, type = 3)
tab_model(model_DARPP_BMI)
}
#posthoc tests
{posthoc2 <- emtrends(model_DARPP_BMI, ~ DARPP | condition, var = "zBMI")
pairs(posthoc2) #only update shows sig. difference p = 0.0115
posthoc2
}

# plot
{
  keep = c("ID","correct","BMI","condition", "DARPP")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_DARPP_BMI)
  plotdata = unique(plotdata)
  
  ggplot(data = plotdata, aes(BMI, fitted, colour = DARPP, facets = condition)) +
    geom_point(alpha = .3) +
    scale_color_manual(values = c("deepskyblue4", "darkgoldenrod2")) +
    geom_smooth(method = lm) +
    theme_classic() +
    facet_grid(.~ condition,
               labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
  
  for (var in unique(plotdata$condition)){
    ggplot(data = plotdata[plotdata$condition==var,], aes(BMI, fitted, colour = DARPP)) +
      geom_smooth(method = "lm") +
      geom_point(alpha = .5, size = 3.5) +
      scale_color_manual(values = c("deepskyblue4","darkgoldenrod2")) +
      theme_classic() +
      theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
      theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
      labs(y = "probability correct") +
      theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
      theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
      facet_wrap(.~ condition,
                 labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
  }
  
  }


# C957T  x BMI ----
#model comparison on data1
{
  model_C957T_BMI_full <- glmer(correct ~ DRD2_rs6277 * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + zAge + project + BED + zDFS + (1|ID), 
                                family = "binomial", 
                                nAGQ = 1,
                                data = na.omit(data1), 
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_C957T_BMI_full, type = 3)
  model_C957T_BMI_reduced <- glmer(correct ~ DRD2_rs6277 * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + (1|ID), 
                                   family = "binomial", 
                                   nAGQ = 1,
                                   data = na.omit(data1), 
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  anova(model_C957T_BMI_full,model_C957T_BMI_reduced)
  
} 

#analysis with main model
{model_C957T_BMI <- glmer(correct ~ DRD2_rs6277 * zBMI * condition + zIQ + zWM_conc + zWM_tired + Gender + (1|ID), 
                        family = "binomial", 
                        nAGQ = 1,
                        data = data, 
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
Anova(model_C957T_BMI, type = 3)
tab_model(model_C957T_BMI)
}




# AAratio x BMI  ----
# no data available for BEDOB
# BMI range for AA sample
{AAdata = data[data$project != "BEDOB",]
mean(AAdata$BMI)
sd(AAdata$BMI)
min(AAdata$BMI)
max(AAdata$BMI)}

# model comparisons with na.omit, because there are NaN in DFS
# exclude factor BED, because BEDOB project doesn't have AAratio
# scale AAratio because of model convergence problems
{
model_AA_BMI_full <- glmer(correct ~ scale(AAratio) * zBMI * condition + zIQ + zWM_conc + Gender + zWM_tired + zDFS + project + zAge + (1|ID), 
                           family = "binomial", 
                           nAGQ = 1,
                           data = na.omit(data), 
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
Anova(model_AA_BMI_full, type = 3)
#zWM_tired, zDFS, project, zAge not sig. --> remove from model, 
#Gender trend sig. --> leave in model, AIC better if I leave it in

model_AA_BMI_reduced <- glmer(correct ~ scale(AAratio) * zBMI * condition + zIQ + zWM_conc + Gender + (1|ID), 
                              family = "binomial", 
                              nAGQ = 1,
                              data = na.omit(data), 
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

anova(model_AA_BMI_reduced, model_AA_BMI_full) #AIC & BIC model_reduced smaller, go with reduced model
}
#model on full dataset
{model_AA_BMI <- glmer(correct ~ AAratio * zBMI * condition + zIQ + zWM_conc + Gender + (1|ID), 
                              family = "binomial", 
                              nAGQ = 1,
                              data = data, 
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_AA_BMI, type = 3)
tab_model(model_AA_BMI)

plot_model(model_AA_BMI, type = "pred", terms = c("zBMI","AAratio","condition"), show.data = F)
plot_model(model_AA_BMI_reduced, type = "pred", terms = c("zIQ"))
plot_model(model_AA_BMI_reduced, type = "pred", terms = c("zWM_conc"))
}

#post hoc test
#for interpretation purposes, subset the data
{
  subset = filter(data,data$condition == 0 | data$condition == 2)
  model_AA_BMI_subset1 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset1, type = 3) #0.000263 ***
  summary(model_AA_BMI_subset1)
  tab_model(model_AA_BMI_subset1)
  #check VIF on model without interactions
  {
  model = glmer(correct ~ zBMI + condition + AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                 family = "binomial", 
                 nAGQ = 1,
                 data = subset, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    
  check_collinearity(model)
  }
  
  
  subset = filter(data,data$condition == 0 | data$condition == 1)
  model_AA_BMI_subset2 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset2, type = 3) #0.246173
  #summary(model_AA_BMI_subset2)
  #tab_model(model_AA_BMI_subset)
  
  subset = filter(data,data$condition == 2 | data$condition == 3)
  model_AA_BMI_subset3 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset3, type = 3) #p = 0.0917754 .
  
  subset = filter(data,data$condition == 1 | data$condition == 3)
  model_AA_BMI_subset4 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset4, type = 3) #p = 0.987743
  
  subset = filter(data,data$condition == 1 | data$condition == 2)
  model_AA_BMI_subset5 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset5, type = 3) #p = 0.0576240 .
  
  subset = filter(data,data$condition == 0 | data$condition == 3)
  model_AA_BMI_subset6 <- glmer(correct ~ zBMI * condition * AAratio + zIQ + Gender + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_subset6, type = 3) #p = 0.2524662
  
  #so there is a difference in update vs. ignore of BMI x AAratio effect
  #but not so in ignore vs. m1 or update vs. m2
  
}

#is the sig. 3-way interaction due to BMI outlier?
{
  #subset data
  library(tidyverse)
  subset = filter(data, data$project != "BEDOB")
  boxplot(subset$BMI)
  unique(subset$ID[subset$BMI > 35])
  subset = filter(subset, subset$BMI < 35)
  
  #run model again on subset
  model_AA_BMI_outl <- glmer(correct ~ AAratio * zBMI * condition + zIQ + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = subset, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  
  Anova(model_AA_BMI_outl, type = 3)
  #tab_model(model_AA_BMI_outl)
  #interaction becomes trend significant, because effect in ignore is not as big anymore
  plot_model(model_AA_BMI_subset, type = "pred", terms = c("zBMI","AAratio","condition"))
}

# plot AA_BMI_condition
{
  keep = c("ID","correct","BMI","condition", "AAratio" ,"AA_group")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_AA_BMI)
  plotdata = unique(plotdata)
  #plotdata = plotdata[-c(which(plotdata$correct == 0)),]
  
  ggplot(data = plotdata, aes(BMI, fitted, colour = AA_group)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = .4) +
    scale_color_manual(values = c( "darkblue","deepskyblue4","darkgoldenrod2")) +
    theme_classic() +
    facet_grid(.~ condition,
               labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
  
  #each condition on its own
  for (var in 0:4){
    ggplot(data = plotdata[plotdata$condition==var,], aes(BMI, fitted, colour = AA_group)) +
      geom_smooth(method = "lm") +
      geom_point(alpha = .5, size = 3.5) +
      scale_color_manual(values = c( "darkblue","deepskyblue4","darkgoldenrod2")) +
      theme_classic() +
      theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
      theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
      labs(y = "probability correct") +
      theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
      theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
      facet_wrap(.~ condition,
                 labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
    }
  
  
#plot w0 BMI outlier
{plotdata = plotdata[-c(which(plotdata$BMI > 35)),]
  ggplot(data = plotdata, aes(BMI, fitted, colour = AA_group)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = .4) +
    scale_color_manual(values = c( "darkblue","deepskyblue4","darkgoldenrod2")) +
    theme_classic() +
    facet_grid(.~ condition,
               labeller = labeller(condition = c(`0` = "ignore", `1` = "m1", '2' = "update", '3' = "m2")))
}
}



# effect of updating dependent on probetype? ----

data_update = data[data$condition == 2, ]
data_update$probetype <- as.factor(data_update$probetype)
#probetype: 
# 0 = novel
# 1 = distractor
# 2 = target

#Taq1A
{
model_Taq1A_BMI_probe <- glmer(correct ~ probetype * DRD2 *zBMI + zIQ + zWM_tired + zWM_conc + (1|ID), 
                             family = "binomial", 
                             nAGQ = 1,
                             data = data_update, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
#summary(model_DRD2_BMI_probe)
Anova(model_Taq1A_BMI_probe, type = 'III')

plot_model(model_Taq1A_BMI_probe, type = "pred", terms = c("zBMI", "DRD2","probetype"), show.data = F) + theme_classic()
plot_model(model_Taq1A_BMI_probe, type = "pred", terms = c("probetype"))

#plot nicer
{
  keep = c("ID","correct","BMI","probetype", "DRD2")
  plotdata = data_update[ , (names(data_update) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_Taq1A_BMI_probe)
  plotdata = unique(plotdata)
  
  gen = plotdata$DRD2 #DRD2
  
  ggplot(data = plotdata, aes(BMI, fitted, colour = gen, facets = probetype)) +
    geom_point(alpha = .3) +
    scale_color_manual(values = c("deepskyblue4", "darkgoldenrod2")) +
    geom_smooth(method = lm) +
    theme_classic() +
    theme(text = element_text(size=20)) +
    facet_grid(.~ probetype,
               labeller = labeller(probetype = c(`0` = "novel", `1` = "distractor", '2' = "target")))
}

}

# DARPP
{
model_DARPP_BMI_probe <- glmer(correct ~ probetype * DARPP *zBMI + zIQ + zWM_tired + zWM_conc + (1|ID), 
                              family = "binomial", 
                              nAGQ = 1,
                              data = data_update, 
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
#summary(model_DARPP_BMI_probe)
Anova(model_DARPP_BMI_probe, type = 'III')

plot_model(model_DARPP_BMI_probe, type = "pred", terms = c("probetype","DARPP")) + theme_classic()

#post hoc test
{
  posthoc <- emtrends(model_DARPP_BMI_probe, ~ DARPP | probetype, var = "zBMI")
  pairs(posthoc)
}

#plot nicer
{
  keep = c("ID","correct","BMI","probetype", "DARPP")
  plotdata = data_update[ , (names(data_update) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_DARPP_BMI_probe)
  plotdata = unique(plotdata)
  
  gen = plotdata$DARPP
  
  ggplot(data = plotdata, aes(BMI, fitted, colour = gen, facets = probetype)) +
    geom_point(alpha = .3) +
    scale_color_manual(values = c("deepskyblue4", "darkgoldenrod2")) +
    geom_smooth(method = lm) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    theme(axis.text.y = element_text(face = "bold", size = 22, color = "black", family = "arial")) +
    labs(y = "probability correct") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 22, face = "bold", family = "Arial")) +
    facet_grid(.~ probetype,
               labeller = labeller(probetype = c(`0` = "novel", `1` = "distractor", '2' = "target")))
}

}


# supplements: SNP interactions ----
model_Taq1A_C957T_BMI <- glmer(correct ~ condition * DRD2_rs6277 * DRD2 * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                             family = "binomial", 
                             nAGQ = 1,
                             data = data, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_Taq1A_C957T_BMI, type = 3)


model_Taq1A_AA_BMI <- glmer(correct ~ condition * AAratio * DRD2 * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                              family = "binomial", 
                              nAGQ = 1,
                              data = data, 
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_Taq1A_AA_BMI, type = 3)


model_DARPP_AA_BMI <- glmer(correct ~ condition * AAratio * DARPP * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                            family = "binomial", 
                            nAGQ = 1,
                            data = data, 
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_DARPP_AA_BMI, type = 3)
plot_model(model_DARPP_AA_BMI, type = "pred", terms = c("zBMI", "condition"), show.data = F) + theme_classic()


model_COMT_C957T_BMI <- glmer(correct ~ condition * DRD2_rs6277 * COMT * zBMI + zIQ + Gender + zWM_tired + zWM_conc + (1|ID), 
                               family = "binomial", 
                               nAGQ = 1,
                               data = data, 
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

Anova(model_COMT_C957T_BMI, type = 3)

