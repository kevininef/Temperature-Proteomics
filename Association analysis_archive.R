# Analysis corrected 20240827 ##### 

rm(list = ls())

# load the library
library(tidyverse)
library(lubridate)
library(haven)
library(splines)
library(lme4)
library(writexl)
library(haven)
library(readxl)
library(dlnm)


# inverse transformation for all the proteins ####
proteins <- protein_temp_long[,6:370]


# check the normality for NPX values
results.test <- data.frame()
for(i in 1:length(proteins)){	
  norm.test <- ks.test(proteins[,i], 'pnorm')
  results.normKS.stat <- as.data.frame(norm.test$statistic)
  results.normKS.p <- as.data.frame(norm.test$p.value)
  results.norm <- cbind(results.normKS.stat,results.normKS.p)
  names(results.norm) <- c("statistic","pvalue")
  results.norm$test <- "Kolmogorov-Smirnov"
  results.norm
  # results.norm$Protein <- protname$Protein_biomarker[i]
  results.test <- rbind(results.test,results.norm)
}

results.test$protein <- colnames(proteins)

results.test$normality <- ifelse(results.test$pvalue<=0.05,'non-normal','normal')
table(results.test$normality)
# non-normal 365 

BiocManager::install("FRGEpistasis")
library(FRGEpistasis)

# use inverse normal transformation 

protein_temp_long$id <- paste(protein_temp_long$idnr,protein_temp_long$phase,sep = '_') # 0 duplicates
protein_temp_long <- protein_temp_long[!duplicated(protein_temp_long$id),]

names(protein_temp_long)[6:370] <- paste(names(protein_temp_long)[6:370],'raw',sep = '.')

rownames(protein_temp_long) <- protein_temp_long$id

prot.int <- list()
for (i in colnames(protein_temp_long)[7:371]){
  x <- protein_temp_long[,i]
  print(head(x))
  prot <- x
  c <- 3/8
  prot.int[[i]] <- rankTransPheno(prot,c)
}

int.all <- as.data.frame(do.call(rbind,prot.int))
int.all <- t(int.all)
int.all <- as.data.frame(int.all)
rownames(int.all) <- rownames(protein_temp_long)

names(int.all) <- gsub(x = names(int.all), pattern = "raw", replacement = "int")

int.all <- int.all %>%  # add the id column
  add_column(id = rownames(int.all),
             .before = "IL2.int")

# add the transformed to the dataset

protein_temp_long <- merge(protein_temp_long,
                           int.all,by = 'id',all.x = T)

# check the missing vaules in the covariates 
summary(protein_temp_long$smoke_24)
summary(protein_temp_long$vacc_1st)
summary(protein_temp_long$vacc_2nd)
summary(protein_temp_long$vacc_3rd)
table(protein_temp_long$season,useNA = 'always')
table(protein_temp_long$covid,useNA = 'always')

protein_temp_long <- merge(protein_temp_long,
                           BAMSE[,c('idnr','q24ifdat')],
                           by='idnr',all.x = T)

protein_temp_long$time <- difftime(protein_temp_long$date_blood,protein_temp_long$q24ifdat)
protein_temp_long$time <- as.numeric(protein_temp_long$time)







# step 1 : comparing the model with or without temperature, to select the signficant proteins and the relevant P-value ####
protein_temp_long$temp_lag0 <- protein_temp_long$temp
protein_temp_long$temp_lag04 <- apply(protein_temp_long[,c('temp_lag0','temp_lag1','temp_lag2','temp_lag3','temp_lag4')],1,mean)
protein_temp_long$temp_lag04_pct <- ntile(protein_temp_long$temp_lag04,100)

result_sig <- matrix(NA,ncol = 2,nrow = ncol(proteins))
colnames(result_sig) <- c('protein','p_value')

covar <- 'ns(time,df=6)+weekday+factor(phase)+sex+age+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'
covar_temp <- 'ns(temp_lag04_pct,df=2)+factor(phase)+sex+age+ns(time,df=6)+weekday+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'


proteins <- protein_temp_long %>% dplyr::select(ends_with('.int'))

for (k in 1:ncol(proteins)){
  x <- colnames(proteins)[k]
  fmla_notemp <- as.formula(paste(x,'~',covar,sep = ''))
  fmla_withtemp <- as.formula(paste(x,'~',covar_temp,sep = ''))
  model_notemp <- lmer(fmla_notemp,data = protein_temp_long)
  model_withtemp <- lmer(fmla_withtemp,data = protein_temp_long)
  p_results <- lmtest::lrtest(model_notemp,model_withtemp)
  result_sig[k,1] <- x
  result_sig[k,2] <- p_results$`Pr(>Chisq)`[2]
  print(k)
  k <- k+1
}

result_sig <- as.data.frame(result_sig)
result_sig$p <- gsub(pattern = '.int',replacement = '',x =  result_sig$protein)


result_sig$p.adj <- p.adjust(result_sig$p_value,method = 'BH')
sum(result_sig$p.adj<0.05,na.rm = T)

result_significant <- result_sig[result_sig$p.adj<0.05|result_sig$p=='SCGB1A1',]

write.csv(result_sig,file = 'step1_protein_tempvsnotemp.csv')
write.csv(result_significant,file = 'step1_Temp_proteins.csv')

#### step 2 make the ER curve for significant proteins ####

covar_nl <- 'ns(temp_lag04_pct,df=5)+ns(time,df=6)+phase+sex+age+weekday+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

protein_predicted <- matrix(NA,ncol = 7,nrow = 1)
colnames(protein_predicted) <- c('x','predicted','std.error','conf.low','conf.high','group','protein')
for (k in 1:length(result_significant$protein)){
  x <- result_significant$protein[k]
  fmla_nonlinear <- as.formula(paste(x,'~',covar_nl,sep = ''))
  model_nonlinear <- lmer(fmla_nonlinear,data = protein_temp_long)
  df <- predict_response(model_nonlinear,terms = 'temp_lag04_pct[all]',margin = 'empirical')
  df$protein <- x
  protein_predicted <- rbind(protein_predicted,df)
  #  ggsave(plot=p_plot,filename = paste0(x,'_ER.tiff'),width = 6,height = 4,units = 'in',dpi = 300,compression='lzw')
  print(k)
  k <- k+1
}

# plot by facet_grid 

View(protein_predicted)
protein_predicted <- protein_predicted[2:length(protein_predicted$x),]
protein_predicted$protein <- gsub(pattern = '.int',replacement = '',x= protein_predicted$protein)

protein_ercurve <- ggplot(data = protein_predicted,aes(x=x,y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill='grey70',alpha=0.7)+
  facet_wrap(vars(protein),ncol = 7,scales = 'free_y')+
  xlab('Temperature percentile')+
  ylab('Predicted protein levels')

ggsave(plot = protein_ercurve,filename = 'protein_ercurve_df5.tiff',units = 'in',width = 12,height = 15,dpi = 300,compression='lzw')


##### step 3 quantify the cold and heat effect using the DLNM function #######

library(dlnm)

# using the most frequent temperature levels as cutoff # 68 percentile

# change to percentile 
protein_temp_long <- protein_temp_long %>% 
            mutate(temp_lag0_pct= ntile(temp_lag0,100),
                   temp_lag1_pct= ntile(temp_lag1,100),
                   temp_lag2_pct= ntile(temp_lag2,100),
                   temp_lag3_pct= ntile(temp_lag3,100),
                   temp_lag4_pct= ntile(temp_lag4,100))


cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

hist(protein_temp_long$temp_lag04_pct)

table(protein_temp_long2$temp_lag04_pct)

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 65,at=95)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 65,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
}


cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[protein_temp_long$phase==c('A2','A3'),c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)),arglag = list(fun='ns'))

model <- lmer(CXCL1.int~cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr),data=protein_temp_long[protein_temp_long$phase==c('A2','A3'),])
fit <- crosspred(cbtemp.lag04,model,cen=75,at=c(1:99))
plot(fit,lag=3)


# add the p_value and adjusted p_value 

res_temp_proteins <- merge(res_temp_proteins,
                           result_sig,
                           by.x = 'p',by.y='protein',
                           all.x = T)
res_temp_mft <- res_temp_proteins



# check the proteins 
res_temp_mft$cold_p_adj <- p.adjust(res_temp_mft$cold_pvalue,method = 'fdr')
sum(res_temp_mft$cold_p_adj<0.05&res_temp_mft$p.adj<0.05,na.rm = T)
  cold_proteins <- res_temp_mft$protein[res_temp_mft$cold_p_adj<0.05&res_temp_mft$p.adj<0.05]

res_temp_mft$heat_p_adj <- p.adjust(res_temp_mft$heat_pvalue,method = 'fdr')
sum(res_temp_mft$heat_p_adj<0.05&res_temp_mft$p.adj<0.05,na.rm=T)
heat_proteins <- res_temp_mft$protein[res_temp_mft$heat_p_adj<0.05&res_temp_mft$p.adj<0.05]

p_temp <-  intersect(cold_proteins,heat_proteins) # overlap among the cold- and heat- proteins 

p_cold <- res_temp_mft$p[res_temp_mft$cold_p_adj<0.05&res_temp_mft$heat_p_adj>0.05&res_temp_mft$p.adj<0.05] # cold only
p_heat <- res_temp_mft$p[res_temp_mft$cold_p_adj>0.05&res_temp_mft$heat_p_adj<0.05&res_temp_mft$p.adj<0.05] # heat only 

library(lubridate)
 temp_daily %>% 
       mutate(month=as.numeric(month(Datum)),
              season=ifelse(month%in%c(4,5,6,7,8,9),'warm','cold')) %>% 
       group_by(season) %>% 
       summarise(mean_temp = mean(temp, na.rm = TRUE))

### step 4 sensitivity analysis ####
 
protein_temp_long$temp_lag01_pct <- ntile(protein_temp_long$temp_lag01,100) 
protein_temp_long$temp_lag07_pct <- ntile(protein_temp_long$temp_lag07,100) 
protein_temp_long$temp_lag014_pct <- ntile(protein_temp_long$temp_lag014,100)
protein_temp_long$temp_lag021_pct <- ntile(protein_temp_long$temp_lag021,100)
# add the additional temperature exposure 
protein_temp_long <- protein_temp_long %>% 
  mutate(temp_lag8_pct=ntile(temp_lag8,100),
         temp_lag9_pct=ntile(temp_lag9,100),
         temp_lag10_pct=ntile(temp_lag10,100),
         temp_lag11_pct=ntile(temp_lag11,100),
         temp_lag12_pct=ntile(temp_lag12,100),
         temp_lag13_pct=ntile(temp_lag13,100),
         temp_lag14_pct=ntile(temp_lag14,100),
         temp_lag15_pct=ntile(temp_lag15,100),
         temp_lag16_pct=ntile(temp_lag16,100),
         temp_lag17_pct=ntile(temp_lag17,100),
         temp_lag18_pct=ntile(temp_lag18,100),
         temp_lag19_pct=ntile(temp_lag19,100),
         temp_lag20_pct=ntile(temp_lag20,100),
         temp_lag21_pct=ntile(temp_lag21,100))

protein_temp_long <- protein_temp_long %>% 
   mutate(temp_lag5_pct= ntile(temp_lag5,100),
          temp_lag6_pct= ntile(temp_lag6,100),
          temp_lag7_pct= ntile(temp_lag7,100))

# use lag01 exposure 
 
cbtemp.lag01 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct')]),
                            lag = 1,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag01+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag01,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag01,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_lag01 <- res_temp_proteins
res_sensi_lag01$type <- 'Lag01'

# use lag 0-7 day exposure 

cbtemp.lag07 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct',
                                                           'temp_lag2_pct','temp_lag3_pct',
                                                           'temp_lag4_pct','temp_lag5_pct',
                                                           'temp_lag6_pct','temp_lag7_pct')]),
                            lag = 7,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag07+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag07,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag07,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_lag07 <- res_temp_proteins
res_sensi_lag07$type <- 'Lag07'

# use lag 0-14 day exposure 

cbtemp.lag014 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct',
                                                           'temp_lag2_pct','temp_lag3_pct',
                                                           'temp_lag4_pct','temp_lag5_pct',
                                                           'temp_lag6_pct','temp_lag7_pct',
                                                           'temp_lag8_pct','temp_lag9_pct',
                                                           'temp_lag10_pct','temp_lag11_pct',
                                                           'temp_lag12_pct','temp_lag13_pct',
                                                           'temp_lag14_pct')]),
                            lag = 14,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag014+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag014,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag014,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_lag014 <- res_temp_proteins
res_sensi_lag014$type <- 'Lag014'

# use lag 0-21 day exposure 

cbtemp.lag021 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct',
                                                            'temp_lag2_pct','temp_lag3_pct',
                                                            'temp_lag4_pct','temp_lag5_pct',
                                                            'temp_lag6_pct','temp_lag7_pct',
                                                            'temp_lag8_pct','temp_lag9_pct',
                                                            'temp_lag10_pct','temp_lag11_pct',
                                                            'temp_lag12_pct','temp_lag13_pct',
                                                            'temp_lag14_pct','temp_lag15_pct',
                                                            'temp_lag16_pct','temp_lag17_pct',
                                                            'temp_lag18_pct','temp_lag19_pct',
                                                            'temp_lag20_pct','temp_lag21_pct')]),
                             lag = 21,
                             argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag021+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag021,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag021,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_lag021 <- res_temp_proteins
res_sensi_lag021$type <- 'Lag021'

## adjusting for PM2.5 


cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+pm25_lag01+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_adjpm25 <- res_temp_proteins
res_sensi_adjpm25$type <- 'Adjusting for PM2.5'

# adjusting for NO2 

cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+no2_lag01+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_adjno2 <- res_temp_proteins
res_sensi_adjno2$type <- 'Adjusting for NO2'

# adjusting for Ozone 
cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+Ozone_lag01+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_adjozone <- res_temp_proteins
res_sensi_adjozone$type <- 'Adjusting for Ozone'

# adjust for RH 
cbtemp.lag04 <- crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+rh_lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_adjrh <- res_temp_proteins
res_sensi_adjrh$type <- 'Adjusting for RH'

# restrict to phase 2&3 only 


protein_p23 <-  subset(protein_temp_long,protein_temp_long$phase%in%c('A2','A3'))

cbtemp.lag04 <- crossbasis(as.matrix(protein_p23[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_p23)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_phase23 <- res_temp_proteins
res_sensi_phase23$type <- 'Phase2&3 only'

# restrict to 1&2 only 

protein_p12 <-  subset(protein_temp_long,protein_temp_long$phase%in%c('A1','A2'))

cbtemp.lag04 <- crossbasis(as.matrix(protein_p12[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_p12)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_phase12 <- res_temp_proteins
res_sensi_phase12$type <- 'Phase1&2 only'

# excluding participants with asthma medication intake 

protein_asthmamed <-  subset(protein_temp_long,protein_temp_long$asthma_medication==0)

cbtemp.lag04 <- crossbasis(as.matrix(protein_asthmamed[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),lag = 4,argvar = list(fun='ns',knots=c(25,50,75)))
covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_asthmamed)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_asthmamed <- res_temp_proteins
res_sensi_asthmamed$type <- 'Excluding people taking asthma medications'

## using cummulative exposure lag window ##### 

cbtemp.lag04 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag04_pct')]),
                            lag = 0,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 68,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_cumlag04 <- res_temp_proteins
res_sensi_cumlag04$type <- 'Lag04 without DLNM'

## change the referent temperature to medium 

cbtemp.lag04 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),
                            lag = 4,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 50,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 50,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_medium <- res_temp_proteins
res_sensi_medium$type <- '50th as reference'

# change the reference temperature to 25th for cold, 75th for heat 
cbtemp.lag04 <-  crossbasis(as.matrix(protein_temp_long[,c('temp_lag0_pct','temp_lag1_pct','temp_lag2_pct','temp_lag3_pct','temp_lag4_pct')]),
                            lag = 4,
                            argvar = list(fun='ns',knots=c(25,50,75)))

covar <- 'cbtemp.lag04+ns(time,df=6)+phase+age+weekday+sex+vacc_1st+vacc_2nd+vacc_3rd+covid_combined+season+(1|idnr)'

res_temp_proteins <-as.data.frame(matrix(NA,ncol=7,nrow = ncol(proteins)))
colnames(res_temp_proteins) <- c('p','heat_est','heat_se','heat_pvalue','cold_est','cold_se','cold_pvalue')

for (k in 1:ncol(proteins)) {
  x <- colnames(proteins)[k]
  fmla <- as.formula(paste(x,'~',covar,sep = ''))
  mymodel <- lmer(fmla,data = protein_temp_long)
  heatpredict <- crosspred(cbtemp.lag04,mymodel,cen = 75,at=99)
  coldpredict <- crosspred(cbtemp.lag04,mymodel,cen = 25,at=1)
  a1 <- heatpredict$allfit
  b1 <- heatpredict$allse
  df1 <- df.residual(mymodel)
  p_heat <-  2*pt(-abs(a1/b1),df1)
  A3 <- coldpredict$allfit
  b2 <- coldpredict$allse
  df2 <- df.residual(mymodel)
  p_cold <-  2*pt(-abs(A3/b2),df2)
  res_temp_proteins[k,1] <- x
  res_temp_proteins[k,2] <- a1
  res_temp_proteins[k,3] <- b1
  res_temp_proteins[k,4] <- p_heat
  res_temp_proteins[k,5] <- A3
  res_temp_proteins[k,6] <- b2
  res_temp_proteins[k,7] <- p_cold
  k <- k+1
}
res_sensi_quartile <- res_temp_proteins
res_sensi_quartile$type <- '25th as reference'






## combine the sensitivity analysis results #### 
library(tidyverse)
sensitivty_results <- res_temp_mft %>% 
  dplyr::select(p,heat_est,heat_se,heat_pvalue,cold_est,cold_se,cold_pvalue) %>% 
  mutate(type='Main analysis') %>% 
  bind_rows(res_sensi_lag01) %>% 
  bind_rows(res_sensi_lag07) %>% 
  bind_rows(res_sensi_lag014) %>% 
  bind_rows(res_sensi_lag021) %>%
  bind_rows(res_sensi_cumlag04) %>% 
  bind_rows(res_sensi_medium) %>% 
  bind_rows(res_sensi_quartile) %>% 
  bind_rows(res_sensi_adjpm25) %>% 
  bind_rows(res_sensi_adjno2) %>% 
  bind_rows(res_sensi_adjozone) %>% 
  bind_rows(res_sensi_adjrh) %>% 
  bind_rows(res_sensi_phase12) %>% 
  bind_rows(res_sensi_phase23) %>% 
  bind_rows(res_sensi_asthmamed) 


write.csv(sensitivty_results,file = 'sensitvity_fulllist_202400905.csv')






