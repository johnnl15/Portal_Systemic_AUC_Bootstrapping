#We thank Jinsu Park for help with bootstrapping and statistical methods. Method based on 
# Bootstrapping, Using Efron, B. and Tibshirani, R. (1993). An Introduction to the Bootstrap. 

# The below code is adapted from this publication which contains Johnny Le's github: https://pubmed.ncbi.nlm.nih.gov/37337122/

library(stringi)
library(stringr)
library(purrr)
library(ggpubr)
library(reshape2)
library(cowplot)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr,warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(gtools)
library(RColorBrewer)
library(tidyr)

mycurrentdirectory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mycurrentdirectory)

list.files("../Datafiles")
excel_sheets("../Datafiles/BDF_C1_C2_water_vs_HFCS_PvsS_glucose.xlsx")

Dataframe<-read_excel("../Datafiles/BDF_C1_C2_water_vs_HFCS_PvsS_glucose.xlsx",sheet = "Merged")

colnames(Dataframe)[1]<-"time"

Dataframemeltnas<-melt(Dataframe,id.vars = "time") %>% 
  mutate(genetics = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$1"),
         treatment = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$2"),
         vessel = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$3"),
         ID = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$4")) %>%
  group_by(time,genetics,treatment,ID) %>%
  filter(any(value == "na")) #%>%
 # mutate(value = as.numeric(value))

write.csv(Dataframemeltnas,file = "../Datafiles/Glucose_NA_values.csv",row.names = FALSE)

Dataframemelt<-melt(Dataframe,id.vars = "time") %>% 
  mutate(genetics = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$1"),
         treatment = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$2"),
         vessel = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$3"),
         ID = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)_(.*)","$4")) %>%
  group_by(time,genetics,treatment,ID) %>%
  filter(!any(value == "na")) %>%
  mutate(value = as.numeric(value))
  #summarise(sum = n()) #just to check whether we've eliminated all with na's and their pair %>%
  
Systemicmelt<-Dataframemelt %>% filter(vessel == "systemic") %>% select(!c(vessel,variable)) %>% rename("serumvalue" = "value")
Portalmelt<-Dataframemelt %>% filter(vessel == "portal") %>% select(!c(vessel,variable)) %>% rename("batvalue" = "value")

b<-merge(Systemicmelt,Portalmelt)

b$Compound_temp<-paste("glucose",b$genetics,b$treatment, sep = "_")

compound_temp_type<-unique(b$Compound_temp)

b2<-b

bootmet<-data.frame(compound_temp=(rep(NA,length(compound_temp_type))),
                    bataverage=(rep(NA,length(compound_temp_type))),
                    serumaverage=(rep(NA,length(compound_temp_type))),
                    portalminussystemic=(rep(NA,length(compound_temp_type))),
                    tstat=(rep(NA,length(compound_temp_type))),
                    CILower97.5=(rep(NA,length(compound_temp_type))),
                    CIHigher2.5=(rep(NA,length(compound_temp_type))),
                    pvalue_boot=(rep(NA,length(compound_temp_type))),
                    bootstrappingtrialsnumber=(rep(NA,length(compound_temp_type))))

for (i in 1:length(compound_temp_type)) { 
  b<-b2%>% filter(Compound_temp==compound_temp_type[i])

alpha<-0.05 # significance level

time_set<-unique(b$time)[order(unique(b$time),decreasing = FALSE)]

n<-length(time_set)

bbat_mean<-b %>% group_by(time) %>% summarise(average=mean(batvalue)) 
colnames(bbat_mean)[colnames(bbat_mean)=="average"]<-"bataverage"
bserum_mean<-b %>% group_by(time) %>% summarise(average=mean(serumvalue)) 
colnames(bserum_mean)[colnames(bserum_mean)=="average"]<-"serumaverage"

b_mean_sep<-merge(bbat_mean,bserum_mean)

bbat_sd<-b %>% group_by(time) %>% summarise(average=sd(batvalue)) 
bserum_sd<-b %>% group_by(time) %>% summarise(average=sd(serumvalue)) 

################################################################
# Setting the Time
################################################################
a_set<-rep(NA,n)

a_set[1]<-(time_set[2]-time_set[1])/2
for(j in 2:n) {
  a_set[j]<-(time_set[j+1]-time_set[j-1])/2
}
a_set[n]<-(time_set[n]-time_set[n-1])/2

timetable<-data.frame(time=time_set,
                      a_set=a_set)

################################################################

b_mean_sep_merged<-merge(b_mean_sep,timetable)

AUC_mean<-b_mean_sep %>% summarise(BAT_AUC=sum(a_set*bataverage),Serum_AUC=sum(a_set*serumaverage))

AUC_sd<-b %>% 
  group_by(time) %>% 
  summarise(BATsd=sd(batvalue)^2/n(),
            Serumsd=sd(serumvalue)^2/n())
AUC_sd<-merge(AUC_sd,timetable)
AUC_sd$a_set<-AUC_sd$a_set^2
AUC_sd <- AUC_sd %>% group_by(time) %>% mutate(BATsdaset=BATsd*a_set,
                                                   Serumsdaset=Serumsd*a_set)
AUC_sd <- AUC_sd %>% ungroup() %>%
  summarise(BAT_AUC_sd=sqrt(sum(BATsdaset)),Serum_AUC_sd=sqrt(sum(Serumsdaset)))# %>%

BatSerum_mean<-b %>% group_by(time) %>% summarise(BatSerumAve=mean(c(batvalue,serumvalue)))

Boot_B<-merge(b,BatSerum_mean,by="time")
Boot_B<-merge(Boot_B,b_mean_sep)
Boot_BSmooth<-Boot_B %>% mutate(BATSmoothValue=batvalue-bataverage+BatSerumAve)
Boot_BSmooth<-Boot_BSmooth %>% mutate(SerumSmoothValue=serumvalue-serumaverage+BatSerumAve)
 
Boot_BSmooth<-merge(Boot_BSmooth,timetable,by="time")
Boot_BSmooth$BATSmoothValueAset<-Boot_BSmooth$BATSmoothValue*Boot_BSmooth$a_set #not using the smoothened value
Boot_BSmooth$SerumSmoothValueAset<-Boot_BSmooth$SerumSmoothValue*Boot_BSmooth$a_set #not using the smoothened value

library(foreach)
library(doParallel)
library(boot)
numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores

trials <- 10000

system.time({
Bootstrapping_stat <- foreach(icount(trials), .combine=c) %dopar% {
    samples<-Boot_BSmooth %>% 
      group_by(time) %>%
      do(sample_n(.,n(),replace = TRUE))
    samples2<-samples %>% group_by(time) %>% 
      summarise(BATaverage=mean(BATSmoothValueAset),Serumaverage=mean(SerumSmoothValueAset)) %>%
      summarise(BATAUC=sum(BATaverage),SerumAUC=sum(Serumaverage)) %>%
      summarise(Bootstrapping_difference=BATAUC-SerumAUC)
    samples3<-samples %>% 
      group_by(time) %>% summarise(BATsd=sd(BATSmoothValueAset)^2/n(),Serumsd=sd(SerumSmoothValueAset)^2/n())
    samples3<-merge(samples3,timetable)
    samples3$a_set<-samples3$a_set^2
    samples3 <- samples3 %>% group_by(time) %>% 
      mutate(BATsdaset=BATsd*a_set,Serumsdaset=Serumsd*a_set) %>% 
      ungroup() %>%
      summarise(BAT_sd=sqrt(sum(BATsdaset)),Serum_sd=sqrt(sum(Serumsdaset))) %>%
      summarise(Bootstrapping_sd=sqrt(BAT_sd^2+Serum_sd^2))
    samples2[[1]]/samples3[[1]] 
  }
})

t_stat<-(AUC_mean$BAT_AUC-AUC_mean$Serum_AUC)/sqrt(AUC_sd$BAT_AUC_sd^2+AUC_sd$Serum_AUC_sd^2)
CI_Boot<-AUC_mean$BAT_AUC-AUC_mean$Serum_AUC-sqrt(AUC_sd$BAT_AUC_sd^2+AUC_sd$Serum_AUC_sd^2)*quantile(Bootstrapping_stat,probs=c(1-alpha/2,alpha/2),na.rm = TRUE)
pvalue_Boot<-2*min(mean(t_stat>=Bootstrapping_stat),mean(t_stat<=Bootstrapping_stat))

bootmet[i,1]<-compound_temp_type[i]
bootmet[i,2]<-AUC_mean$BAT_AUC
bootmet[i,3]<-AUC_mean$Serum_AUC
bootmet[i,4]<-AUC_mean$BAT_AUC-AUC_mean$Serum_AUC
bootmet[i,5]<-t_stat
bootmet[i,6]<-CI_Boot[1]
bootmet[i,7]<-CI_Boot[2]
bootmet[i,8]<-pvalue_Boot
bootmet[i,9]<-trials
}

#stopped here and will continue. 

bootmet <- bootmet %>% rename("portalAUCaverage"="bataverage" ,
                              "systemiAUCcaverage"="serumaverage",
                              "AUCportalminusAUCsystemicaverage"="portalminussystemic")

write.csv(bootmet,file = "../Datafiles/Bootstrapped_Result_Glucose.csv",row.names = FALSE)


