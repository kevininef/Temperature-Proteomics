##################################################################################
 #### predict biological aging ####
library(tidyverse)
library(lightgbm)
library(Boruta)
library(caret)
library(lme4)

data_age <- protein_temp_long %>% 
            dplyr:: select(idnr,phase,sex,age,ends_with('.raw'))

data_age <- data_age[,!colnames(data_age)%in%protein_lod] # excluding the proteins below LOD
data_age$sex <- as.factor(data_age$sex)


data_age <- data_age[complete.cases(data_age),]

# perfrom linear mixed model to remove the repeated measurment 

residuals_df <- data.frame(idnr= data_age$idnr)

proteins_columns <- grep('.raw',names(data_age),value = T)

residuals_df <- data.frame(idnr = data_age$idnr)



# Loop through each protein column and calculate residuals
for (i in 1:282) {
  
  p <- proteins_columns[i]
  # Fit the linear mixed model
  lmm <- lmer(as.formula(paste0(p, "~ 1 + (1 | idnr)")), data = data_age)
  
  # Extract residuals and store them in the residuals_df
  residuals_df <- cbind(residuals_df,resid(lmm))
  colnames(residuals_df)[i+1] <- p
  print(i)
}

colnames(residuals_df)[2:283] <- gsub(pattern = '.raw',
                               replacement = '.res',
                               x=colnames(residuals_df)[2:283])

data_age_res <- cbind(residuals_df, data_age[,c('sex','age','phase')])
data_age_res$id <- paste0(data_age_res$idnr,'_',data_age_res$phase)
#### fit c

set.seed(0816)


trainIndex <- createDataPartition(data_age_res$age,p=0.7,list = FALSE)

trainData <- data_age_res[trainIndex,]
testData <- data_age_res[-trainIndex,]
rownames(trainData) <- trainData$id
rownames(testData) <- testData$id
trainData$id <- NULL
testData$id <- NULL

trainData$idnr <- NULL
trainData$sex <- NULL
trainData$phase <- NULL

testData$idnr <- NULL
testData$sex <- NULL
testData$phase <- NULL



train_matrix <- as.matrix(trainData %>% dplyr::select(ends_with('.res'))) # only include proteins 
train_label <- trainData$age

test_matrix <- as.matrix(testData %>% dplyr::select(ends_with('.res'))) # only include proteins
test_label <- testData$age              

# convert data to lightGBM dataset 
dtrain <- lgb.Dataset(data = train_matrix, label = train_label)


# Define a custom feature importance function using LightGBM
lightgbm_importance <- function(X, y) {
  dtrain <- lgb.Dataset(data = as.matrix(X), label = y)
  
  # Train a LightGBM model
  model <- lgb.train(
    params = list(objective = "regression", boosting = "gbdt"),
    data = dtrain,
    nrounds = 100,
    verbose = 0
  )
  
  # Get feature importance
  importance <- lgb.importance(model)
  # Match importance to feature names and return as a named vector
  imp_vec <- setNames(importance$Gain, importance$Feature)
  return(imp_vec[match(colnames(X), names(imp_vec))])
}


# Use Boruta with the custom LightGBM importance function
boruta_result <- Boruta(
  x = trainData[,-283],
  y = train_label, 
  getImp = lightgbm_importance, 
  maxRuns = 500,
  doTrace = 2
  
)

# Print the results
print(boruta_result)

# Get the final decision on features
final_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
print(final_features)

# Plot the results
plot(boruta_result)

# use the proteins selected by boruta to reconstruct the lightGBM model 

trainData <- trainData[,colnames(trainData)%in%final_features]
train_matrix <- as.matrix(trainData %>% dplyr::select(ends_with('.res'))) # remove outcome age

testData <- testData[,colnames(testData)%in%final_features]
test_matrix <- as.matrix(testData %>% dplyr::select(ends_with('.res')))

# rebuild the lgb datast 
dtrain <- lgb.Dataset(data = train_matrix, label = train_label)


# Set parameters for LightGBM
params <- list(
  objective = "regression",
  metric = "rmse",        # using customed function to evaluate 
  boosting_type = "gbdt", # Gradient Boosting Decision Tree
  num_leaves = 31,        # Maximum number of leaves in one tree
  learning_rate = 0.05,   # Learning rate
  feature_fraction = 0.9  # Proportion of features to be used in each iteration
)


cv_results <- lgb.cv(
  params = params,
  data = dtrain,
  nfold = 5,                   # 5-fold cross-validation
  nrounds = 500,               # Number of boosting iterations
  early_stopping_rounds = 10,  # Stop early if validation error doesn't improve
  verbose = 1
)

best_nrounds <- cv_results$best_iter

model_age <- lgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 1
)

preds <- predict(model_age, test_matrix)
rmse <- sqrt(mean((preds - test_label)^2))
cat("RMSE on test data:", rmse, "\n")

# importance

importance <- lgb.importance(model_age)
View(importance)

write.csv(importance,file = 'LightGBM_featureselection_age_importance.csv')

rownames(data_age_res) <- data_age_res$id
data_age_res$id <- NULL

full_matrix <- as.matrix(data_age_res[,colnames(data_age_res)%in%final_features])

full_preds <- predict(model_age,full_matrix)

View(data_age_res)
protein_temp_long2 <- protein_temp_long[protein_temp_long$id%in%rownames(data_age_res),]

data_age_res$proteomic_age <- full_preds
data_age_res$id <- rownames(data_age_res)
cor.test(data_age_res$age,data_age_res$proteomic_age)

protein_temp_long2 <- merge(protein_temp_long2,
                            data_age_res[,c('id','proteomic_age')],
                            by='id',
                            all.x = T)

library(hrbrthemes)

plot_agingclock_a <- ggplot(data=data_age_res,
       aes(x=age, y=proteomic_age)) +
  geom_point() +
 # geom_abline(slope = 0.899, intercept = 5,linewidth=2,color='red') +
  geom_smooth(method = 'lm',linewidth=1.5,color='red',se=T) +
  theme_minimal() +
  xlab('Chronological age')+
  annotate("text", x = 25.5, y = 27.5, label = "Correlation= 0.93, R-square = 0.85", size = 4, color = "black")+
  theme(plot.margin = margin(t=1,
                             r=0,
                             b=0,
                             l=0))

plot_agingclock_b <- importance %>% 
             mutate(Feature= gsub(pattern='.res',
                                  replacement='',
                                  Feature),
                    Feature= factor(Feature,levels=Feature)) %>% 
             ggplot(aes(x=Feature,y=Gain))+
             geom_bar(stat = 'identity',fill='#3585b4')+
             theme_minimal()+
             theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
             xlab('Proteins')+ylab('Gain')+
  theme(plot.margin = margin(t=1,
                             r=0,
                             b=0,
                             l=0))

plot_agingclock_c <- protein_temp_long2 %>% 
                 mutate(phase= gsub(pattern='A',
                                    replacement='',
                                    x= phase)) %>% 
               ggplot(aes(x=phase,y=age_acce,fill=phase))+
                              geom_violin(width=0.8)+
               geom_boxplot(width=0.4, color="grey", alpha=0.2)+
               theme_minimal()+
               xlab('Phase')+ylab('Proteomic Age Accerelation')+
               theme(legend.position = 'none')+
               scale_color_nejm()+
  theme(plot.margin = margin(t=0,
                             r=0,
                             b=1,
                             l=0))

      

ss_total <- sum((protein_temp_long2$proteomic_age - mean(protein_temp_long2$age))^2)         # Total sum of squares
ss_residual <- sum((protein_temp_long2$age - protein_temp_long2$proteomic_age)^2)  # Residual sum of squares
r2 <- 1 - (ss_residual / ss_total)                   # R² calculation
r2                                            


cor.test(protein_temp_long2$proteomic_age,protein_temp_long2$age)
summary(lm(age~proteomic_age,data = protein_temp_long2))

# ggplot(data=protein_temp_long,
#        aes(x=factor(phase),y=age,color=phase))+
#         geom_violin()+
#         theme_ipsum()

protein_temp_long2$age_acce <- protein_temp_long2$proteomic_age-protein_temp_long2$age

# calcuate the proteomic age accerelation 
mod_age <- lmer(age~proteomic_age+(1|idnr),data = protein_temp_long2)

protein_temp_long2$age_residual <- -residuals(mod_age)


ggplot(data=protein_temp_long2,
       aes(x=factor(phase),y=age_residual,fill=phase))+
  geom_boxplot()+
  facet_wrap(vars(sex))+
  theme_ipsum()


# association between temperature exposure and proteomic age accerelation
library(dlnm)
library(splines)
cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long2[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),
                           lag =4,argvar = list(fun='ns',knots=c(25,50,75)),arglag = list(fun='poly'))





test <- lmer(age_residual~cbtemp.lag04+ns(time,df=6)+sex+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr),data = protein_temp_long2)
predict <- crosspred(cbtemp.lag04,test,cen=68,at=c(1:99),cumul = T)

predict$allfit
predict$allhigh
predict$alllow
plot(predict,cex=1.5,col='blue',lag=4)
plot(predict)

# get the MMT from the model 


MMT0 <- predict$predvar[which.min(predict$cumfit[,1])]
MMT1 <- predict$predvar[which.min(predict$cumfit[,2])]  
MMT3 <- predict$predvar[which.min(predict$matfit[,4])]
MMT4 <- predict$predvar[which.min(predict$matfit[,5])]

cold_ageacce_lag04 <- cbind(predict$predvar,
                      as.numeric(predict$matfit[,5]),
                      as.numeric(predict$matse[,5]),
                      as.numeric(predict$mathigh[,5]),
                      as.numeric(predict$matlow[,5]))
cold_ageacce_lag04 <- as.data.frame(cold_ageacce_lag04)
colnames(cold_ageacce_lag04) <- c('temperature','allfit','allse','allhigh','alllow')


heat_ageacce_lag01 <- cbind(predict$predvar,
                            as.numeric(predict$cumfit[,1]),
                            as.numeric(predict$cumse[,1]),
                            as.numeric(predict$cumhigh[,1]),
                            as.numeric(predict$cumlow[,1])
                            )
heat_ageacce_lag01 <- as.data.frame(heat_ageacce_lag01)
colnames(heat_ageacce_lag01) <- c('temperature','allfit','allse','allhigh','alllow')



cold_ageacce_lag04$alllow <- cold_ageacce_lag04$allfit-1.96*cold_ageacce_lag04$allse
cold_ageacce_lag04$allhigh <- cold_ageacce_lag04$allfit+1.96*cold_ageacce_lag04$allse

heat_ageacce_lag01$alllow <- heat_ageacce_lag01$allfit-1.96*heat_ageacce_lag01$allse
heat_ageacce_lag01$allhigh <- heat_ageacce_lag01$allfit+1.96*heat_ageacce_lag01$allse



# scale on MMT 

mmt <- cold_ageacce_lag04$temperature[which.min(cold_ageacce_lag04$allfit)] 

cold_ageacce_lag04$allfit <- cold_ageacce_lag04$allfit- cold_ageacce_lag04$allfit[cold_ageacce_lag04$temperature==mmt]
cold_ageacce_lag04$alllow <- cold_ageacce_lag04$alllow- cold_ageacce_lag04$alllow[cold_ageacce_lag04$temperature==mmt]
cold_ageacce_lag04$allhigh <- cold_ageacce_lag04$allhigh- cold_ageacce_lag04$allhigh[cold_ageacce_lag04$temperature==mmt]




plot_cold04_aging <- ggplot(data = cold_ageacce_lag04,aes(x=temperature,y=allfit))+
  geom_line(aes(color=range),linewidth=1.3)+
  geom_ribbon(aes(ymin=alllow,ymax=allhigh),alpha=0.3,fill='grey70')+
  xlab('Temperature percentile')+
  ylab('Proteomic Age Accerelation')+
  theme_minimal()+
  scale_color_manual(values = c('#0072B5','#BC3C29'))+
  geom_abline(intercept = 0,slope = 0,lty=2)+
  theme(legend.position = 'none')
plot_cold04_aging


plot_heat01_aging <- ggplot(data = heat_ageacce_lag01,aes(x=temperature,y=allfit))+
  geom_line(aes(color=range),linewidth=1.3)+
  geom_ribbon(aes(ymin=alllow,ymax=allhigh),alpha=0.3,fill='grey70')+
  xlab('Temperature percentile')+
  ylab('Proteomic Age Accerelation')+
  theme_minimal()+
  scale_color_manual(values = c('#0072B5','#BC3C29'))+
  geom_abline(intercept = 0,slope = 0,lty=2)+
  theme(legend.position = 'none')

plot_heat01_aging

plot_agingclock_d <-  plot_heat01_aging %>% 
          mutate(allse= case_when(
            range=='heat'~allse,
            range=='cold'~allse
          )) %>% 
          mutate(alllow=allfit-1.96*allse,
                 allhigh=allfit+1.96*allse) %>% 
          ggplot(aes(x=temperature,y=allfit))+
          geom_line(aes(color=range),linewidth=1.3)+
          geom_ribbon(aes(ymin=alllow,ymax=allhigh),alpha=0.3,fill='grey70')+
          xlab('Temperature percentile')+
          ylab('Proteomic Age Accerelation')+
  theme_minimal()+
  scale_color_manual(values = c('#0072B5','#BC3C29'))+
  geom_abline(intercept = 0,slope = 0,lty=2)+
  theme(legend.position = 'none')+
  theme(plot.margin = margin(t=0,
                             r=0,
                             b=1,
                             l=0))

plot_agingclock_d



library(cowplot)
plot_temp_aging <- plot_grid(plot_agingclock_a,plot_agingclock_b,plot_agingclock_c,plot_agingclock_d, nrow = 2,labels = c('A','B','C','D'))
ggsave(filename = 'temp_aging.tiff',plot = plot_temp_aging,width = 14,height = 10,units = 'in',dpi = 300,compression='lzw')

protein_temp_long %>% 
      distinct(idnr,.keep_all = T) %>% 
      summarise(number=n())
