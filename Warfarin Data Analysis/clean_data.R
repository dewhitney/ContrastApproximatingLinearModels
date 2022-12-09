################################################################################
###                                                                          ###
###                             clean_data.R                                 ###
###                                                                          ###
################################################################################
###                                                                          ###
###  File modified: 18 Jul 2022                                              ###
###  File creator: David Whitney (from Wallace, Moodie, Stephens 2018)       ###
###  Purpose: to create an analysis data for comparison against IWPC dosing  ###
###           algorithm as well as other optimal treatment rules             ###
###                                                                          ###
################################################################################

# Begin Wallace et al code, with two changes to match (n=1732)
# 1. DW commented out a block that imputed Enzyme/Amiodarone..Cordarone
# 2. DW commented out filter on Warfarin Dose < 6

# The above changes yield a match on descriptive statistics reported in e.g.
# Schulz and Moodie (2021), Wallace et al (2018), Chen et al (2016)

# read in data

setwd("~/Documents/Rcode/WarfarinDataAnalysis")

# Note: download the 'IWPC dataset' from https://www.pharmgkb.org/downloads
# Save the required sheet from the Excel data download
dose <- read.delim('IWPC.csv',sep=",",header=TRUE,as.is=TRUE)

index1 <- which(dose$Subject.Reached.Stable.Dose.of.Warfarin==1)
index2 <- which(is.na(dose$INR.on.Reported.Therapeutic.Dose.of.Warfarin)==FALSE)
dose.tmp <- dose[intersect(index1,index2),]
index3 <- which(dose.tmp$Race..OMB.=='Unknown')
dose.tmp1 <-  dose.tmp[-index3,]
index4 <- which(is.na(dose.tmp1$Age))
index5 <- which(is.na(dose.tmp1$Weight..kg.))
index6 <- which(is.na(dose.tmp1$Height..cm.))
index7 <- which(is.na(dose.tmp1$Gender))
dose.tmp2 <-  dose.tmp1[-union(union(union(index4,index5),index6),index7),]
#### genetics vs clinical ####


##### imputation based on page 13 of supplement of NEJM 2009 #########
index8 = which(is.na(dose.tmp2$CYP2C9.consensus))
dose.tmp3 <- dose.tmp2[-index8,]
index9 = which(is.na(dose.tmp3$VKORC1..1639.consensus))
#### IMPUTE rs9923231
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1173.consensus=="C/C"),index9)]="G/G"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1173.consensus=="T/T"),index9)]="A/A"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1173.consensus=="C/T"),index9)]="A/G"

index9 = intersect(which(is.na(dose.tmp3$VKORC1..1639.consensus)),which(dose.tmp3$Race..OMB.=='White'))
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1542.consensus=="G/G"),index9)]="G/G"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1542.consensus=="C/C"),index9)]="A/A"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.1542.consensus=="C/G"),index9)]="A/G"

dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.2255.consensus=="C/C"),index9)]="G/G"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.2255.consensus=="T/T"),index9)]="A/A"
dose.tmp3$VKORC1..1639.consensus[intersect(which(dose.tmp3$VKORC1.2255.consensus=="C/T"),index9)]="A/G"

index9 = which(is.na(dose.tmp3$VKORC1..1639.consensus))

dose.tmp4 = dose.tmp3[-index9,]

## one type is taken then enzyme is one, none of taken is zero, otherwise missing.
dose.tmp4$Enzyme = NA
dose.tmp4$Enzyme[which(dose.tmp4$Rifampin.or.Rifampicin==1)] = 1
dose.tmp4$Enzyme[which(dose.tmp4$Carbamazepine..Tegretol.==1)] = 1
dose.tmp4$Enzyme[which(dose.tmp4$Phenytoin..Dilantin.==1)] = 1
dose.tmp4$Enzyme[(dose.tmp4$Rifampin.or.Rifampicin==0)+(dose.tmp4$Phenytoin..Dilantin.==0)+(dose.tmp4$Carbamazepine..Tegretol.==0) > 2] = 0
#dose.tmp4$Rifampin.or.Rifampicin   dose.tmp4$Carbamazepine..Tegretol. dose.tmp4$Phenytoin..Dilantin.

### DW commented out -- it resulted in a different sample from comparison papers

## Updated version: We believe one would not use two type of drug at the same time.

# index10 = intersect(which(is.na(dose.tmp4$Enzyme)),which(dose.tmp4$Amiodarone..Cordarone. == 1))
# index11 = intersect(which(is.na(dose.tmp4$Amiodarone..Cordarone.)),which(dose.tmp4$Enzyme == 1))
# dose.tmp4$Enzyme[index10] = 0
# dose.tmp4$Amiodarone..Cordarone.[index11] = 0

### end DW comment out

### get rid of missing data. 
dose.tmp5 = dose.tmp4[-which(is.na(dose.tmp4$Amiodarone..Cordarone.)),]
dose.tmp6 = dose.tmp5[-which(is.na(dose.tmp5$Enzyme)),]


### Final dataset construction
traindata = data.frame(Age = dose.tmp6$Age,Weight = dose.tmp6$Weight..kg.,Height = dose.tmp6$Height..cm.,
                       Enzyme = dose.tmp6$Enzyme, Amiodarone= dose.tmp6$Amiodarone..Cordarone.)

### 0 for 10-19
traindata$Age =  as.numeric(as.factor(dose.tmp6$Age)) - 1
### 0 for female
traindata$Gender =  as.numeric(as.factor(dose.tmp6$Gender)) - 1

## 0 for White
traindata$Black = as.numeric(dose.tmp6$Race..OMB. ==  'Black or African American')
traindata$Asian = as.numeric(dose.tmp6$Race..OMB. ==  'Asian')
traindata$VKORC1.AG = as.numeric(dose.tmp6$VKORC1..1639.consensus =='A/G')
traindata$VKORC1.AA = as.numeric(dose.tmp6$VKORC1..1639.consensus =='A/A')
## 0 for CYP2C9.consensus =='*1/*1'
traindata$CYP2C9.12 = as.numeric(dose.tmp6$CYP2C9.consensus =='*1/*2')
traindata$CYP2C9.13 = as.numeric(dose.tmp6$CYP2C9.consensus =='*1/*3')
traindata$CYP2C9.other = 1 - as.numeric(dose.tmp6$CYP2C9.consensus =='*1/*1') - traindata$CYP2C9.12 - traindata$CYP2C9.13
traindata$Dose = dose.tmp6$Therapeutic.Dose.of.Warfarin
traindata$INR = dose.tmp6$INR.on.Reported.Therapeutic.Dose.of.Warfarin

### DW commented out Dose < 6 as 5.81 is reported in ODTR comparison papers
# indexoutlier = which(traindata$Dose<6| traindata$Dose>95)
indexoutlier = which(traindata$Dose>95)
traindata = traindata[-indexoutlier,]
X = data.matrix(traindata[,-15])

X = data.frame(X)

# End Wallace et al code

rm(list = setdiff(ls(), "traindata"))
####





