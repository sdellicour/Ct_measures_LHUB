###############################################################################
# Ct: math modelling study by age group with comparison with estimated prevalence
# Nicolas Franco - UHasselt/UNamur
# Version: June 25 2021
###############################################################################

# cleaning if needed
rm(list = ls()); gc() 

#path to working directory
path_to_wd <- "."  
#path to csv file with confirmed ct
CT_file <- "./ctages.csv"
#path to txt from comparmental model (with percentiles)
model_file <- "./p50.txt"
model_file_low <- "./p5.txt"
model_file_up <- "./p95.txt"
#number of bootstrap
BS_num = 999
#lenght window
Window = 14

setwd(path_to_wd)

#libraries
library(tidyverse)
library(readr)
library(zoo)
library(scales)

#bootstrap function
bootstrapsCI = function(vS, type="median", conf.level=0.95, R=999){
  differences <- rep(NA, R)
  vS <- vS[which(!is.na(vS))]
  if(length(vS)==0){
    return(c(NA,NA,NA))
  }
  if(length(vS)==1){
    return(c(NA,NA,vS))
  }
  for (i in 1:R)
  {
    if (type == "median") differences[i] <- median(sample(vS,length(vS),replace=T))-median(vS)
    if (type == "mean") differences[i] <- mean(sample(vS,length(vS),replace=T))-median(vS)
  }
  return(median(vS)+c(quantile(differences, probs=c(0.5-(conf.level/2),0.5+(conf.level/2))),0))
}

bootstrapsCI_tribble = function(vS){
  output <- bootstrapsCI(vS,R=BS_num)
  tibble(medianCt = output[3],median_lower = output[1], median_upper = output[2])
}


#read data
Modelbrut <- as.data.frame(read.delim(model_file, header = FALSE,sep = " ")) 
Model <- data.frame(date=as.Date(character()),prevalence=numeric(),agegroup=character())
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(3:6,8)])/(11500000),agegroup=" All ages")
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(11:14,16)])/3250000,agegroup="0-24")
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(19:22,24)])/3000000,agegroup="25-44")
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(27:30,32)])/3080000,agegroup="45-64")
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(35:38,40)])/1150000,agegroup="65-74")
Model <- Model %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence=rowSums(Modelbrut[,c(43:46,48,51:54,56)])/(1020000),agegroup="75+")

Modelbrut <- as.data.frame(read.delim(model_file_low, header = FALSE,sep = " ")) 
Model_low <- data.frame(date=as.Date(character()),prevalence_low=numeric(),agegroup=character())
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(3:6,8)])/(11500000),agegroup=" All ages")
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(11:14,16)])/3250000,agegroup="0-24")
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(19:22,24)])/3000000,agegroup="25-44")
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(27:30,32)])/3080000,agegroup="45-64")
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(35:38,40)])/1150000,agegroup="65-74")
Model_low <- Model_low %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_low=rowSums(Modelbrut[,c(43:46,48,51:54,56)])/(1020000),agegroup="75+")

Modelbrut <- as.data.frame(read.delim(model_file_up, header = FALSE,sep = " ")) 
Model_up <- data.frame(date=as.Date(character()),prevalence_up=numeric(),agegroup=character())
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(3:6,8)])/(11500000),agegroup=" All ages")
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(11:14,16)])/3250000,agegroup="0-24")
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(19:22,24)])/3000000,agegroup="25-44")
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(27:30,32)])/3080000,agegroup="45-64")
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(35:38,40)])/1150000,agegroup="65-74")
Model_up <- Model_up %>% add_row(date=as.Date("2020-02-29")+Modelbrut[,1],prevalence_up=rowSums(Modelbrut[,c(43:46,48,51:54,56)])/(1020000),agegroup="75+")

#clean data
Ctbrut <- read.csv(CT_file,colClasses = c("Age" = "integer","Ct" = "real"),na.strings = "0")
Ctbrut <- filter(Ctbrut,Ct>0,!is.na(Age),grepl("M 2000 Real Time",Analyseur))
Ctwindow <- Ctbrut
Ctwindow$date <- as.Date(Ctbrut$Date.prelv,"%d/%m/%Y")
for (i in 1:(Window-1)) {
  Cttemp <- Ctbrut
  Cttemp$date <- as.Date(Ctbrut$Date.prelv,"%d/%m/%Y")+i
  Ctwindow <- full_join(Ctwindow,Cttemp)
}

#agregate data 
Ctagregages <- as.data.frame(Ctwindow) %>% filter(date <= as.Date("2021-05-17"),date >= as.Date("2020-05-01")) %>% mutate(agegroup = case_when(Age >= 75 ~ "75+",Age >= 65 ~ "65-74",Age >= 45 ~ "45-64",Age >= 25 ~ "25-44",Age >= 0 ~ "0-24"))
Ctagregall <- as.data.frame(Ctwindow) %>% filter(date <= as.Date("2021-05-17"),date >= as.Date("2020-05-01")) %>% mutate(agegroup = case_when(Age >= 0 ~ " All ages"))
Ctagreg <- Ctagregall
Ctagreg <-  full_join(Ctagreg,Ctagregages)
Ctagreg <- Ctagreg %>% group_by(agegroup,date)  %>% summarise(bootstrapsCI_tribble(Ct))
Ctagreg <- Ctagreg %>% filter(date <= as.Date("2021-05-17"),date >= as.Date("2020-05-14")) 

Commondataframe <- inner_join(Ctagreg,Model)
Commondataframe <- inner_join(Commondataframe,Model_low)
Commondataframe <- inner_join(Commondataframe,Model_up)

Commondataframe$agegroup <- factor(Commondataframe$agegroup)

#plot data
Coef <- 0.05
Trans <- 30
plot1 <- ggplot(Commondataframe) +
  geom_line(aes(x = date, y=Coef*(-medianCt+Trans)), colour="red",linetype = "solid",size=0.5)   +
  geom_ribbon(aes(x = date,ymin=Coef*(-median_lower+Trans), ymax=Coef*(-median_upper+Trans)) ,fill="red", alpha=0.3) +
  geom_line(aes(x = date, y=100*prevalence),colour="gray30",linetype = "solid",size=0.5)   +
 geom_ribbon(aes(x = date,ymin=100*prevalence_low, ymax=100*prevalence_up,group = 1) ,fill="gray30", alpha=0.2) +
  facet_wrap(~agegroup, strip.position="top", ncol=2) +
  scale_y_continuous(name = "Daily estimated prevalence",labels = function(x) paste0(x, "%"),sec.axis = sec_axis(~./-Coef+Trans, name="Ct median over 14 days")) +
  scale_x_date(date_minor_breaks = "1 month",name="",labels = date_format("%d/%m/%Y")) +
  coord_cartesian(ylim = c(0, 1.2),xlim=as.Date(c("2020-05-15","2021-05-17"))) +
  theme_minimal() +
  theme(axis.text.y  = element_text(color = 'gray30'),axis.title.y = element_text(color='gray30'),axis.text.y.right =  element_text(color = 'red'),axis.title.y.right = element_text(color='red'))

plot1

#pdf output
pdf("Figure_ct_prevalence.pdf", width=10, height=8)
print(plot1);
dev.off() 

