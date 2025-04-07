##This analysis tests for differences in slope between two groups (wild vs insectary) and subsequently the presence of elevational or major-axis shifts if slopes are not significantly different

#load libraries
library(smatr)
library('openxlsx')
library(ggplot2)
#read in data, data here should be with the volume of individual samples
data<-read.csv("Table S10- wild vs insectary.csv",header=T)

#structuring data for analysis
data$wild.insectary<-as.factor(data$wild.insectary)
data$sp.abb<-as.factor(data$sp.abb)
data$sex<-as.factor(data$sex)
data$clade<-as.factor(data$clade)
#note: sp.app= species name abbrevation, sim.clade=simplified clade name 
#omit NAs
data<-na.omit(data)
#add logs to AL ALH, GL and rCBR
data$AL.log <-log10(data$AL)
data$ALH.log<-log10(data$ALH)
data$GL.log<-log10(data$GL)
data$rCBR.log<-log10(data$rCBR)

#check data structure 
str(data)
#subset data to their species group
melp<-subset(data,clade=="MELPOMENE")
cydno<-subset(data,clade=="CYDNO")
hecale<-subset(data,sp.abb=="Hhel")
ism<-subset(data,sp.abb=="Hism")
erato<-subset(data,clade=="ERATO")
wild<-subset(data,wild.insectary=="wild")
insectary<-subset(data,wild.insectary=="insectary")

##test for differences in the AL of H.erato,melpomene,cydno,ismenius and hecale
#erato 
#test for differences in slope
#run sma 
eALslope<-sma(AL.log~rCBR.log*wild.insectary, robust=F, data=erato)
summary(eALslope)
#write to file
write.xlsx(eALslope$coef, "erato coefficients AL slope.xlsx")
#slope not significantly different, test for elevational shifts and major-axis shifts 
eALele<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=erato)
summary(eALele)
#write to file
write.xlsx(eALele$coef, "erato coefficients ele shift.xlsx")
#elevation not significant, test for major-axis shift 
eALma<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=erato)
summary(eALma)
#write to file
write.xlsx(eALma$coef, "erato coefficients major axis shift.xlsx")
#no significant major-axis shift 


#melpomene
#run sma
#test for differences in slope
mALslope<- sma(AL.log~rCBR.log*wild.insectary, robust=F, data=melp)
summary(mALslope)
#write to file
write.xlsx(mALslope$coef, "melp coefficients AL slope.xlsx")
#slope is significantly different 

#cydno
#run sma
#test for differences in slope
cALslope<- sma(AL.log~rCBR.log*wild.insectary, robust=F, data=cydno)
summary(cALslope)
#write to file
write.xlsx(cALslope$coef, "cydno coefficients AL slope.xlsx")
#slope not significant, test for elevation and major-axis shift 
cALele<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=cydno)
summary(cALele)
#write to file 
write.xlsx(cALele$coef, "cydno coefficients ele shift.xlsx")
#elevation not significant, test for major-axis shift
cALma<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=cydno)
summary(cALma)
#write to file 
write.xlsx(cALma$coef, "cydno coefficients major axis shift.xlsx")
#no significant major axis shift as well 

#ismenius
#run sma 
#test for differences in slope
iALslope<- sma(AL.log~rCBR.log*wild.insectary, robust=F, data=ism)
summary(iALslope)
#write to file
write.xlsx(iALslope$coef, "ismenius coefficients AL slope.xlsx")
#slope not significant, test for elevation and major-axis shift 
iALele<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=ism)
summary(iALele)
#write to file
write.xlsx(iALele$coef, "ismenius coefficients ele shift.xlsx")
#elevation not significant, test for major-axis shift
iALma<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=ism)
summary(iALma)
#write to file
write.xlsx(iALma$coef, "ismenius coefficients major axis shift.xlsx")
#no significant major axis shift as well 

#hecale
#run sma
#test for differences in slope
hALslope<- sma(AL.log~rCBR.log*wild.insectary, robust=F, data=hecale)
summary(hALslope)
#write to file
write.xlsx(hALslope$coef, "hecale coefficients AL slope.xlsx")
#slopes not significant, test for elevation and major-axis shift 
hALele<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=hecale)
summary(hALele)
#write to file
write.xlsx(hALele$coef, "hecale coefficients ele shift.xlsx")
#elevational shifts significant, test for major-axis shifts
hALma<- sma(AL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=hecale)
summary(hALma)
#write to file
write.xlsx(hALma$coef, "hecale coefficients major axis shift.xlsx")
#major-axis shifts significant 


##test for differences in the ALH of H.erato,melpomene,cydno,ismenius and hecale
#erato 
#run sma 
#test for differences in slope
eALHslope<-sma(ALH.log~rCBR.log*wild.insectary, robust=F, data=erato)
summary(eALHslope)
#write to file
write.xlsx(eALHslope$coef, "erato coefficients ALH slope.xlsx")
#slopes not significant, test for elevation and major-axis shift 
eALHele<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=erato)
summary(eALHele)
#write to file
write.xlsx(eALHele$coef, "erato coefficients ALH ele shift.xlsx")
#elevational shifts not significant, test for major-axis shifts
eALHma<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="shift", data=erato)
summary(eALHma)
#write to file
write.xlsx(eALHma$coef, "erato ALH coefficients major axis shift.xlsx")
#sisgnificant major-axis shifts

#melpomene
#run sma 
#test for differences in slope
mALHslope<- sma(ALH.log~rCBR.log*wild.insectary, robust=F, data=melp)
summary(mALHslope)
#write to file
write.xlsx(mALHslope$coef, "melp coefficients ALH slope.xlsx")
#slopes not significant, test for elevation and major-axis shift 
mALHele<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=melp)
summary(mALHele)
#write to file
write.xlsx(mALHele$coef, "melp coefficients ALH ele shift.xlsx")
#elevational shifts not significant, test for major-axis shifts
mALHma<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="shift", data=melp)
summary(mALHma)
#write to file
write.xlsx(mALHma$coef, "melp ALH coefficients major axis shift.xlsx")
#significant major-axis shifts

#cydno
#run sma 
#test for differences in slope
cALHslope<-sma(ALH.log~rCBR.log*wild.insectary, robust=F, data=cydno)
summary(cALHslope)
#write to file
write.xlsx(cALHslope$coef, "cydno coefficients ALH slope.xlsx")
#slopes not significant, test for elevation and major-axis shift 
cALHele<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=cydno)
summary(cALHele)
#write to file
write.xlsx(cALHele$coef, "cydno coefficients ALH ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
cALHma<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="shift", data=cydno)
summary(cALHma)
write.xlsx(cALHma$coef, "cydno ALH coefficients major axis shift.xlsx")
#major-axis shifts not significant

#ismenius
#run sma 
#test for differences in slope
iALHslope<-sma(ALH.log~rCBR.log*wild.insectary, robust=F, data=ism)
summary(iALHslope)
#write to file
write.xlsx(iALHslope$coef, "ismenius coefficients ALH slope.xlsx")
#slopes not significant, test for elevation and major-axis shift 
iALHele<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=ism)
summary(iALHele)
#write to file
write.xlsx(iALHele$coef, "ismenius coefficients ALH ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
iALHma<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="shift", data=ism)
summary(iALHma)
#write to file
write.xlsx(iALHma$coef, "ismenius ALH coefficients major axis shift.xlsx")
#significant major-axis shifts

#hecale
#run sma 
#test for differences in slope
hALHslope<-sma(ALH.log~rCBR.log*wild.insectary, robust=F, data=hecale)
summary(hALHslope)
#write to file
write.xlsx(hALHslope$coef, "hecale coefficients ALH slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
hALHele<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=hecale)
summary(hALHele)
#write to file
write.xlsx(hALHele$coef, "hecale coefficients ALH ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
hALHma<- sma(ALH.log~rCBR.log+wild.insectary, robust=F,type="shift", data=hecale)
summary(hALHma)
#write to file
write.xlsx(hALHma$coef, "hecale ALH coefficients major axis shift.xlsx")
#major-axis shifts not significant


##test for differences in the GL of H.erato,melpomene,cydno,ismenius and hecale
#erato 
#run sma 
#test for differences in slope
eGLslope<-sma(GL.log~rCBR.log*wild.insectary, robust=F, data=erato)
summary(eGLslope)
#write to file
write.xlsx(eGLslope$coef, "erato coefficients GL slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
eGLele<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=erato)
summary(eGLele)
#write to file
write.xlsx(eGLele$coef, "erato coefficients GL ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
eGLma<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=erato)
summary(eGLma)
#write to file
write.xlsx(eGLma$coef, "erato GL coefficients major axis shift.xlsx")
#major-axis shifts not significant

#melpomene
#test for differences in slope
mGLslope<-sma(GL.log~rCBR.log*wild.insectary, robust=F, data=melp)
summary(mGLslope)
#write to file
write.xlsx(mGLslope$coef, "melpomene coefficients GL slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
mGLele<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=melp)
summary(mGLele)
#write to file
write.xlsx(mGLele$coef, "melpomene coefficients GL ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
mGLma<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=melp)
summary(mGLma)
#write to file
write.xlsx(mGLma$coef, "melpomene GL coefficients major axis shift.xlsx")
#major-axis shifts not significant

#cydno
#test for differences in slope
cGLslope<-sma(GL.log~rCBR.log*wild.insectary, robust=F, data=cydno)
summary(cGLslope)
#write to file
write.xlsx(cGLslope$coef, "cydno coefficients GL slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
cGLele<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="elevation", data=cydno)
summary(cGLele)
#write to file
write.xlsx(cGLele$coef, "cydno coefficients GL ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
cGLma<- sma(GL.log~rCBR.log+wild.insectary, robust=F,type="shift", data=cydno)
summary(cGLma)
#write to file
write.xlsx(cGLma$coef, "cydno GL coefficients major axis shift.xlsx")
#major-axis shifts not significant

#ismenius
#test for differences in slope
iGLslope<-sma(GL.log~rCBR.log*wild.insectary, robust=F, data=ism)
summary(iGLslope)
#write to file
write.xlsx(iGLslope$coef, "ismenius coefficients GL slope.xlsx")
#slopes significantly different

#hecale
#test for differences in slope
hGLslope<-sma(GL.log~rCBR.log*wild.insectary, robust=F, data=hecale)
summary(hGLslope)
#write to file
write.xlsx(hGLslope$coef, "hecale coefficients GL slope.xlsx")
#slopes significantly different


##test for differences in the GL~ALH of H.erato,melpomene,cydno,ismenius and hecale
#erato 
#run sma 
#test for differences in slope
ebothslope<-sma(GL.log~ALH.log*wild.insectary, robust=F, data=erato)
summary(ebothslope)
#write to file 
write.xlsx(ebothslope$coef, "erato coefficients both slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
ebothele<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="elevation", data=erato)
summary(ebothele)
#write to file 
write.xlsx(ebothele$coef, "erato coefficients both ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
ebothma<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="shift", data=erato)
summary(ebothma)
#write to file 
write.xlsx(ebothma$coef, "erato both coefficients major axis shift.xlsx")
#major-axis shifts not significant

#melpomene
#run sma 
#test for differences in slope
mbothslope<-sma(GL.log~ALH.log*wild.insectary, robust=F, data=melp)
summary(mbothslope)
#write to file 
write.xlsx(mbothslope$coef, "melpomene coefficients both slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
mbothele<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="elevation", data=melp)
summary(mbothele)
#write to file 
write.xlsx(mbothele$coef, "melpomene coefficients both ele shift.xlsx")
#elevational shifts not significant, test for major-axis shifts
mbothma<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="shift", data=melp)
summary(mbothma)
#write to file 
write.xlsx(mbothma$coef, "melpomene both coefficients major axis shift.xlsx")
#major-axis shifts not significant

#cydno
#run sma 
#test for differences in slope
cbothslope<-sma(GL.log~ALH.log*wild.insectary, robust=F, data=cydno)
summary(cbothslope)
#write to file 
write.xlsx(cbothslope$coef, "cydno coefficients both slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
cbothele<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="elevation", data=cydno)
summary(cbothele)
#write to file 
write.xlsx(cbothele$coef, "cydno coefficients both ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
cbothma<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="shift", data=cydno)
summary(cbothma)
#write to file 
write.xlsx(cbothma$coef, "cydno both coefficients major axis shift.xlsx")
#major-axis shifts not significant

#ismenius 
#run sma 
#test for differences in slope
ibothslope<-sma(GL.log~ALH.log*wild.insectary, robust=F, data=ism)
summary(ibothslope)
#write to file 
write.xlsx(ibothslope$coef, "ismenius coefficients both slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
ibothele<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="elevation", data=ism)
summary(ibothele)
#write to file 
write.xlsx(ibothele$coef, "ismenius coefficients both ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
ibothma<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="shift", data=ism)
summary(ibothma)
#write to file 
write.xlsx(ibothma$coef, "ismenius both coefficients major axis shift.xlsx")
#major-axis shifts not significant

#hecale
#run sma 
#test for differences in slope
hbothslope<-sma(GL.log~ALH.log*wild.insectary, robust=F, data=hecale)
summary(hbothslope)
#write to file 
write.xlsx(hbothslope$coef, "hecale coefficients both slope.xlsx")
#slopes not significant, test for elevation and major-axis shift
hbothele<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="elevation", data=hecale)
summary(hbothele)
#write to file 
write.xlsx(hbothele$coef, "hecale coefficients both ele shift.xlsx")
#significant elevational shifts, test for major-axis shifts
hbothma<- sma(GL.log~ALH.log+wild.insectary, robust=F,type="shift", data=hecale)
summary(hbothma)
#write to file
write.xlsx(hbothma$coef, "hecale both coefficients major axis shift.xlsx")
#major-axis shifts not significant