


## -----------------------------------------------------
library(TraMineR)
library(tidyverse)
library(knitr)
library(ggsci)

df<-read_csv("data_EE.csv")
colors=pal_npg("nrc", alpha = .9)(9)[c(1:4,7)]
labels=c("No Sex-No Union","Sex- No Union","Union-No Sex","No Children","Children")
scode=c("NoSex-NoUn","Sex-NoUn","Un-NoSex","NoChild", "Child")

seqdata=seqdef(df[,19:37], #<<
               cnames=12:30,
               cpal=colors,
               missing=NA,
               states=scode,
               labels=labels,
               missing.color="white", 
               weights=df$aw)




## -----------------------------------------------------
example <- seqdata[1:5,1:10]
kable(example)



## -----------------------------------------------------
seqformat(example, from = "STS", to = "SPS",compress = TRUE) #<<


## -----------------------------------------------------
seqdss(example)



## -----------------------------------------------------
seqIplot(example)


## -----------------------------------------------------
seqdplot(seqdata)


## -----------------------------------------------------
seqmtplot(seqdata)


## ----echo=FALSE---------------------------------------
duration<-seqistatd(seqdata)
df$age2firstevent<-duration[,1]+9
df$age2firstevent[df$S30==0]<-NA
df$age2sex<-duration[,1]+duration[,3]+9
df$age2sex[df$S30==0 | df$S30==3]<-NA

df$age2union<-duration[,1]+duration[,2]+9
df$age2union[df$S30<2]<-NA

df$age2child<-duration[,1]+duration[,2]+duration[,3]+duration[,4]+9
df$age2child[df$S30<4]<-NA

library(data.table)
library(Hmisc)
TableAge2Event<-df %>% 
  group_by(countryname, birth_5y) %>% 
  summarise(AgeFirstEvent=mean(age2firstevent,na.rm=T),
    AgeFirstSex=mean(age2sex, na.rm=T),
    AgeFirstUnion=mean(age2union, na.rm=T),
    AgeFirstChild=mean(age2child, na.rm=T),
    SDAgeFirstSex=sd(age2sex, na.rm=T),
    SDAgeFirstUnion=sd(age2union, na.rm=T),
    SDAgeFirstChild=sd(age2child, na.rm=T),
    SDAgeFirstEvent=sd(age2firstevent, na.rm=T),
    NAgeFirstEvent=sum((is.na(age2firstevent)==0) ),
    NAgeFirstSex=sum((is.na(age2sex)==0) ),
    NAgeFirstUnion=sum((is.na(age2union)==0) ),
    NAgeFirstChild=sum((is.na(age2child)==0) )
  )  %>% 
  arrange(countryname, birth_5y) %>% 
  ungroup()

library(RColorBrewer)
mypal <-  c(brewer.pal(9,"OrRd")[c(3,5,7,9)], brewer.pal(9,"Greens")[c(7,9)], brewer.pal(9,"Blues")[c(5,7,9)], brewer.pal(9,"Purples")[c(5,9)])

child=ggplot(data=TableAge2Event, aes(x=birth_5y, y=AgeFirstChild, group=countryname, col=countryname, fill=countryname)) +
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=AgeFirstChild-2*SDAgeFirstChild/sqrt(NAgeFirstChild), ymax=AgeFirstChild+2*SDAgeFirstChild/sqrt(NAgeFirstChild),fill=countryname), linetype=0, alpha=0.1)+
  theme_light()+ xlab("Birth Cohort")+ylab("Age at first child")+ theme(legend.position=c(.85,.85))+   scale_fill_manual(values = mypal)   +   scale_color_manual( values = mypal)

child




## ----echo=FALSE---------------------------------------

df$sex2union<-duration[,2]
df$sex2union[df$S30<2]<-NA
df$sex2child<-duration[,2]+duration[,4]
df$sex2child[df$S30<4]<-NA
df$union2child<-duration[,3]+duration[,4]
df$union2child[df$S30<4]<-NA

TableYears2Event<- df %>% 
  group_by(countryname, birth_5y) %>% 
  summarise(Time_sex_union=mean(sex2union,na.rm=T),
            Time_sex_child=mean(sex2child, na.rm=T),
            Time_union_child=mean(union2child, na.rm=T),
            SDTime_sex_union=sd(sex2union,na.rm=T),
            SDTime_sex_child=sd(sex2child, na.rm=T),
            SDTime_union_child=sd(union2child, na.rm=T),
            NTime_sex_union=sum((is.na(sex2union)==0) ),
            NTime_sex_child=sum((is.na(sex2child)==0) ),
            NTime_union_child=sum((is.na(union2child)==0) )
  )
union2child=ggplot(data=TableYears2Event, aes(x=birth_5y, y=Time_union_child, col=countryname, fill=countryname)) +
  geom_point()+
  geom_line()+
  theme_light()+
  xlab("Birth Cohort")+
  geom_ribbon(aes(ymin=Time_union_child-2*SDTime_union_child/sqrt(NTime_union_child), ymax=Time_union_child+2*SDTime_union_child/sqrt(NTime_union_child),
                  fill=countryname), linetype=0, alpha=0.1)+
  ylab("Time union to child")+ theme(legend.position=c(.85,.85))+   scale_fill_manual(values = mypal)   +   scale_color_manual(values = mypal)
union2child




## -----------------------------------------------------
submat<-seqsubm(seqdata, method="TRATE", with.missing=TRUE)

kable(submat)


## -----------------------------------------------------
dist.om1=seqdist(seqdata, #<<
                 method="OM",
                 indel=1,
                 sm=submat,
                 with.missing=TRUE)
kable(dist.om1[1:8,1:8], digits=2)


## ----eval=FALSE, include=TRUE-------------------------
 library(cluster)
 seq.clusterward = agnes(dist.om1, diss = T, method = "ward")


## ----eval=FALSE, include=TRUE-------------------------
 cl4 = cutree(seq.clusterward, k = 4)


## ----eval=FALSE, include=TRUE-------------------------
 
 wardRange <- as.clustrange(seq.clusterward, diss = dist.om1, ncluster = 10)
 summary(wardRange, max.rank = 2)
 plot(wardRange, stat = c("ASW", "HG", "PBC", "HC"))

