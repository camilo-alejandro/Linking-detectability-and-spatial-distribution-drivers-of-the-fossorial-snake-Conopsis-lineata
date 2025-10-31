#occupancy modelling C. lineata - Camilo Cruz-Arroyave - 2025

#### detection history ####
library(dplyr)
setwd("D:/OneDrive - Universidad de Antioquia/Postgrado/Articulos/Herpetological_Journal/data/")
dethist<-read.table("dethist_final.csv",sep = ",",header = T)%>% arrange(id)
#rownames(dethist)<-dethist$id 

#### detection covariates  ####

n.cov.obj<-read.csv("Detect_#NatCoverObjects.csv",sep=";",header = F)
colnames(n.cov.obj)<-c("transect","v1","v2","v3","v4","v5","total.F")
rownames(n.cov.obj)<-n.cov.obj$transect
n.cov.obj<-n.cov.obj[-c(1,2,33),-c(1,7)]
n.cov.obj$v1<-as.integer(n.cov.obj$v1)
str(n.cov.obj)

Tmin<-read.csv("Detect_Tmin.csv",sep=";",header = F)
colnames(Tmin)<-c("transect","v1","v2","v3","v4","v5","total.F")
rownames(Tmin)<-Tmin$transect
Tmin<-Tmin[-c(1,2,33),-c(1,7)]
Tmin$v1<-as.numeric(Tmin$v1)
str(Tmin)

Tmax<-read.csv("Detect_Tmax.csv",sep=";",header = F)
colnames(Tmax)<-c("transect","v1","v2","v3","v4","v5","total.F")
rownames(Tmax)<-Tmax$transect
Tmax<-Tmax[-c(1,2,33),-c(1,7)]
Tmax$v1<-as.numeric(Tmax$v1)
str(Tmax)

Hour<-read.csv("Detect_Hour.csv",sep=";",header = F)
colnames(Hour)<-c("transect","v1","v2","v3","v4","v5","total.F")
rownames(Hour)<-Hour$transect
Hour<-Hour[-c(1,2,33),-c(1,7)]
Hour[which(Hour == "",arr.ind=T)]<-NA
Hour[] <- lapply(Hour, function(col) as.POSIXct(col, format = "%I:%M:%S %p"))
str(Hour)
# Extract fraction of the day
hour_fraction <- lapply(Hour, function(col) {
  as.numeric(format(col, "%H")) / 24 + 
    as.numeric(format(col, "%M")) / (24 * 60) + 
    as.numeric(format(col, "%S")) / (24 * 3600)
})
hour_fraction <- as.data.frame(hour_fraction)

obsCovs<-list(n.cov.obj=n.cov.obj,
              Tmin=Tmin,
              Tmax=Tmax,
              hour=hour_fraction)

obsCovs.str<-list(n.cov.obj=as.data.frame(scale(n.cov.obj)),
              Tmin=as.data.frame(scale(Tmin)),
              Tmax=as.data.frame(scale(Tmax)),
              hour=as.data.frame(scale(hour_fraction)))
str(obsCovs)
str(obsCovs.str)


#### site covariates ####

#Import sitecovs
sitecovs<-read.table("siteCovs_julio2025.csv",sep=";",header = T)%>% arrange(id)

rownames(sitecovs)<-sitecovs$id

# Creates dummy variables (0/1) for all levels of veg type
vegtype_dummi <- model.matrix(~vegtype - 1, sitecovs)
siteCovs<-cbind(sitecovs,vegtype_dummi)

str(siteCovs)
names(siteCovs)

#Ensuring correct reading of dataframe
siteCovs$vegtype<-as.factor(siteCovs$vegtype)
siteCovs$agro<-as.factor(siteCovs$agro)
siteCovs$vegtypecrops<-as.factor(siteCovs$vegtypecrops)
siteCovs$vegtypeforest<-as.factor(siteCovs$vegtypeforest)
siteCovs$vegtypegrass<-as.factor(siteCovs$vegtypegrass)
siteCovs$vegtypemixed<-as.factor(siteCovs$vegtypemixed)
siteCovs$vegtypescrub<-as.factor(siteCovs$vegtypescrub)

siteCovs<-select(siteCovs,id,vegh,veghmax,
                 vegtype,vegtypecrops,vegtypeforest,vegtypegrass,vegtypemixed,vegtypescrub,
                 slopemax,agro,
                 convergence,cti,slope,hfp,
                 evapotranspiration,productivity,
                 bio_01,bio_05,bio_06,bio_07,bio_08,bio_12,bio_15,bio_16,bio_17)
str(siteCovs)

#standardizing numeric variables
siteCovs.str <- siteCovs %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) %>% 
  arrange(id)

str(siteCovs.str)

#### detection covariates correlation analyses ####

n.obj<-c(n.cov.obj[,1],n.cov.obj[,2],n.cov.obj[,3],n.cov.obj[,4],n.cov.obj[,5])
T_min<-c(Tmin[,1],Tmin[,2],Tmin[,3],Tmin[,4],Tmin[,5])
T_max<-c(Tmax[,1],Tmax[,2],Tmax[,3],Tmax[,4],Tmax[,5])
Hour<-c(hour_fraction[,1],hour_fraction[,2],hour_fraction[,3],hour_fraction[,4],hour_fraction[,5])
covdet<-cbind(n.obj,T_min,T_max,Hour)
covdet<-as.data.frame(covdet)
str(covdet)

library(corrplot)
cormat.det<-cor(covdet,use = "pairwise.complete.obs",method = "spearman")
corrplot(cormat.det,method = c("number"),type="upper",
         tl.col = "black", tl.cex = 1.2,tl.srt = 30,number.cex = 1,
         mar=c(1, 1, 1, 2)) 
#there is not high correlated covariates of detection.

#### site covariates correlation analyses ####

#exploring database

plot(siteCovs$vegh~siteCovs$vegtype) 
#note: vegetation type is congruent with vegetation height, 
plot(siteCovs$vegh~siteCovs$agro)
#note: it makes sense

library(corrplot)
sitecovs_num <- sitecovs %>% select(where(is.numeric))
str(sitecovs_num)

cormat.site<-cor(sitecovs_num,use = "pairwise.complete.obs",method = "spearman")
corrplot(cormat.site,method = c("number"),type="upper",
         tl.col = "black", tl.cex = 0.6,tl.srt = 30,number.cex = 0.5,cl.cex = 0.5,
         mar=c(0,0,0,0)) 
#note: Considering that bio 5 was reported as important for C. lineata, bio1,4,6,7,12 will be omitted because showed high correlation index at diferent thresholds.  
#note: There are not representative changes in final site variables if we do not choose bio5 apriori

#to create corrplot for chooses variables 
cormat.site<-cor(siteCovs[,-c(3,6,16,17,18,19,20)],use = "pairwise.complete.obs",method = "spearman")
corrplot(cormat.site,method = c("number"),type="upper",
         tl.col = "black", tl.cex = 1,tl.srt = 30,number.cex = 0.8,cl.cex = 0.8,
         mar=c(0,0,0,0)) 


#### occupancy modelling ####
library(unmarked)

conopsis <- unmarkedFrameOccu(y=dethist[,2:6],obsCovs = obsCovs.str,siteCovs =siteCovs.str)
summary(conopsis)
plot(conopsis)

(m1 <- occu(~1 ~1, conopsis))
(m2<- occu(~n.cov.obj ~1, conopsis))
(m4<-occu(~Tmin ~1, conopsis))
(m5<-occu(~Tmax ~1, conopsis))
(m6<-occu(~hour ~1, conopsis))
(m6.1<-occu(~hour^2 ~1, conopsis))
(m8<-occu(~n.cov.obj+Tmin ~1, conopsis))
(m9<-occu(~n.cov.obj+Tmax ~1, conopsis))
(m10<-occu(~n.cov.obj+hour ~1, conopsis))
(m14<-occu(~Tmin+Tmax ~1, conopsis))
(m15<-occu(~Tmin+hour ~1, conopsis))
(m16<-occu(~Tmax+hour ~1, conopsis))

library(AICcmodavg)

models<-list("p(.)ψ(.)"=m1,
             "p(n.cov.obj)ψ(.)"=m2,
             "p(Tmin)ψ(.)"=m4,
             "p(Tmax)ψ(.)"=m5,
             "p(hour)ψ(.)"=m6.1,
             "p(n.cov.obj+Tmin)ψ(.)"=m8,
             "p(n.cov.obj+Tmax)ψ(.)"=m9,
             "p(n.cov.obj+hour)ψ(.)"=m10,
            "p(Tmin+Tmax)ψ(.)"=m14,
            "p(Tmin+hour)ψ(.)"=m15,
            "p(Tmax+hour)ψ(.)"=m16
             )
             

mods_det<-aictab(models)
mods_det
#write.table(mods_det,"graphs/mods_det_julio_extended.csv",sep=";",
#  row.names=F,col.names = T)

models2<-list("p(.)ψ(.)"=m1,
             "p(n.cov.obj)ψ(.)"=m2,
             "p(Tmin)ψ(.)"=m4,
             "p(Tmax)ψ(.)"=m5,
             "p(hour)ψ(.)"=m6
)
mods_det2<-aictab(models2)
mods_det2
#write.table(mods_det2,"graphs/mods_det_sep.csv",sep=";",
#            row.names=F,col.names = T)

#to analise relation between n.cov.obj and Tmax

n.cov.obj<-c(n.cov.obj$v1,n.cov.obj$v2,n.cov.obj$v3, n.cov.obj$v4,n.cov.obj$v5)
Tmax<-c(Tmax$v1,Tmax$v2,Tmax$v3,Tmax$v4, Tmax$v5)
plot(Tmax~n.cov.obj)
cor(n.cov.obj,Tmax,method = "spearman",use="pairwise.complete.obs")
#note: n.cov.obj and Tmax correlation = 0.011

## Modelling occupancy
#(m20<-occu(~n.cov.obj ~VegHmin, conopsis)) #SE > 4
(m21<-occu(~n.cov.obj ~veghmax, conopsis))
(m22<-occu(~n.cov.obj ~vegtype, conopsis)) #SE very high
(m23.1<-occu(~n.cov.obj ~vegtypegrass, conopsis))
(m23.2<-occu(~n.cov.obj ~vegtypescrub, conopsis)) #SE very high
(m23.3<-occu(~n.cov.obj ~vegtypemixed, conopsis)) 
(m23.4<-occu(~n.cov.obj ~vegtypeforest, conopsis))
(m23.5<-occu(~n.cov.obj ~vegtypecrops, conopsis)) 
(m24<-occu(~n.cov.obj ~slopemax, conopsis))#SE=2.1

(m25<-occu(~n.cov.obj ~hfp, conopsis))
(m26<-occu(~n.cov.obj ~agro, conopsis))
(m27<-occu(~n.cov.obj ~vegh, conopsis))#SE >3

(m28<-occu(~n.cov.obj ~convergence, conopsis))
(m29<-occu(~n.cov.obj ~cti, conopsis)) 
(m30<-occu(~n.cov.obj ~slope, conopsis))
#(m31<-occu(~n.cov.obj ~tri, conopsis))
 
(m33<-occu(~n.cov.obj ~evapotranspiration, conopsis))
(m34<-occu(~n.cov.obj ~productivity, conopsis))
(m35<-occu(~n.cov.obj ~bio_01, conopsis))#SE>5
(m36<-occu(~n.cov.obj ~bio_05, conopsis))
(m37<-occu(~n.cov.obj ~bio_06, conopsis))
(m38<-occu(~n.cov.obj ~bio_07, conopsis))
(m39<-occu(~n.cov.obj ~bio_08, conopsis))#SE>2
(m40<-occu(~n.cov.obj ~bio_12, conopsis))#SE>6
(m41<-occu(~n.cov.obj ~bio_15, conopsis))
(m42<-occu(~n.cov.obj ~bio_16, conopsis))#SE>5
(m43<-occu(~n.cov.obj ~bio_17, conopsis))

mods.occu<-list("p(.)ψ(.)"=m1,
                "p(n.cov.obj)ψ(.)"=m2,
                  # "p(n.cov.obj)ψ(VegHmin)"=m20,
                  "p(n.cov.obj)ψ(VegHmax)"=m21,
                 # "p(n.cov.obj)ψ(VegType)"=m22,
                  "p(n.cov.obj)ψ(Grass)"=m23.1,
                  #"p(n.cov.obj)ψ(Scrub)"=m23.2,
                  "p(n.cov.obj)ψ(Mixed)"=m23.3,
                  "p(n.cov.obj)ψ(Forest)"=m23.4,
                  "p(n.cov.obj)ψ(Crops)"=m23.5,
                  "p(n.cov.obj)ψ(slopeMax)_SE=2"=m24,
                  "p(n.cov.obj)ψ(hfp)"=m25,
                  "p(n.cov.obj)ψ(Agro)"=m26,
                  "p(n.cov.obj)ψ(vegh)-SE>3"=m27,
                  "p(n.cov.obj)ψ(convergence)"=m28,
                  "p(n.cov.obj)ψ(cti)"=m29,
                  "p(n.cov.obj)ψ(slope)"=m30,
                  #"p(n.cov.obj)ψ(tri)"=m31,
                  "p(n.cov.obj)ψ(evapotrans)"=m33,
                  "p(n.cov.obj)ψ(productivity)"=m34,
                  "p(n.cov.obj)ψ(bio1)-SE>5"=m35,
                "p(n.cov.obj)ψ(bio5)"=m36,
                "p(n.cov.obj)ψ(bio6)"=m37,
                "p(n.cov.obj)ψ(bio7)"=m38,
                "p(n.cov.obj)ψ(bio8)-SE>2"=m39,
                "p(n.cov.obj)ψ(bio12)-SE>6"=m40,
                "p(n.cov.obj)ψ(bio15)"=m41,
                "p(n.cov.obj)ψ(bio16)-SE>5"=m42,
                "p(n.cov.obj)ψ(bio17)"=m43
)


(mods_ocu<-aictab(mods.occu))
#write.table(mods_ocu,"graphs/mods_occu_julio.csv",sep=";",
#            row.names=F,col.names = T)
#note: SlopeMax,vegtype show good performance (delta aic<2), bio5 too but do not converge.

##Testing model set based on tmax

(m21t<-occu(~Tmax ~veghmax, conopsis))
(m22t<-occu(~Tmax ~vegtype, conopsis))
(m23.1t<-occu(~Tmax ~vegtypegrass, conopsis))
(m23.2t<-occu(~Tmax ~vegtypescrub, conopsis)) 
(m23.3t<-occu(~Tmax ~vegtypemixed, conopsis)) 
(m23.4t<-occu(~Tmax ~vegtypeforest, conopsis)) 
(m23.5t<-occu(~Tmax ~vegtypecrops, conopsis))
(m24t<-occu(~Tmax ~slopemax, conopsis))
(m25t<-occu(~Tmax ~hfp, conopsis))
(m26t<-occu(~Tmax ~agro, conopsis))
(m27t<-occu(~Tmax ~vegh, conopsis))

(m28t<-occu(~Tmax ~convergence, conopsis))
(m29t<-occu(~Tmax ~cti, conopsis)) 
(m30t<-occu(~Tmax ~slope, conopsis))
#(m31t<-occu(~Tmax ~tri, conopsis))

(m33t<-occu(~Tmax ~evapotranspiration, conopsis))
(m34t<-occu(~Tmax ~productivity, conopsis))
(m35t<-occu(~Tmax ~bio_01, conopsis))
(m36t<-occu(~Tmax ~bio_05, conopsis))
(m37t<-occu(~Tmax ~bio_06, conopsis))
(m38t<-occu(~Tmax ~bio_07, conopsis))
(m39t<-occu(~Tmax ~bio_08, conopsis))
(m40t<-occu(~Tmax ~bio_12, conopsis))
(m41t<-occu(~Tmax ~bio_15, conopsis))
(m42t<-occu(~Tmax ~bio_16, conopsis))
(m43t<-occu(~Tmax ~bio_17, conopsis))

mods.occu.t<-list("p(.)ψ(.)"=m1,
                  "p(n.cov.obj)ψ(.)"=m2,
                  "p(Tmax)ψ(.)"=m5,
                  # "p(n.cov.obj)ψ(VegHmin)"=m20,
                  "p(n.cov.obj)ψ(VegHmax)"=m21,
                  # "p(n.cov.obj)ψ(VegType)"=m22,
                  "p(n.cov.obj)ψ(Grass)"=m23.1,
                  #"p(n.cov.obj)ψ(Scrub)"=m23.2,
                  "p(n.cov.obj)ψ(Mixed)"=m23.3,
                  "p(n.cov.obj)ψ(Forest)"=m23.4,
                  "p(n.cov.obj)ψ(Crops)"=m23.5,
                  "p(n.cov.obj)ψ(slopeMax)_SE=2"=m24,
                  "p(n.cov.obj)ψ(hfp)"=m25,
                  "p(n.cov.obj)ψ(Agro)"=m26,
                  #"p(n.cov.obj)ψ(vegh)-SE>3"=m27,
                  "p(n.cov.obj)ψ(convergence)"=m28,
                  "p(n.cov.obj)ψ(cti)"=m29,
                  "p(n.cov.obj)ψ(slope)"=m30,
                  #"p(n.cov.obj)ψ(tri)"=m31,
                  "p(n.cov.obj)ψ(evapotrans)"=m33,
                  "p(n.cov.obj)ψ(productivity)"=m34,
                  #"p(n.cov.obj)ψ(bio1)-SE>5"=m35,
                  "p(n.cov.obj)ψ(bio5)"=m36,
                  "p(n.cov.obj)ψ(bio6)"=m37,
                  "p(n.cov.obj)ψ(bio7)"=m38,
                  #"p(n.cov.obj)ψ(bio8)-SE>2"=m39,
                  #"p(n.cov.obj)ψ(bio12)-SE>6"=m40,
                  "p(n.cov.obj)ψ(bio15)"=m41,
                  #"p(n.cov.obj)ψ(bio16)-SE>5"=m42,
                  "p(n.cov.obj)ψ(bio17)"=m43,
                
                  "p(Tmax)ψ(VegHmax)"=m21t,
                  "p(Tmax)ψ(VegType)"=m22t,
                  "p(Tmax)ψ(Grass)"=m23.1t,
                  "p(Tmax)ψ(Scrub)"=m23.2,
                  "p(Tmax)ψ(Mixed)"=m23.3t,
                  "p(Tmax)ψ(Forest)"=m23.4t,
                  "p(Tmax)ψ(Crops)"=m23.5t,
                  "p(Tmax)ψ(slopeMax)"=m24t,
                  "p(Tmax)ψ(hfp)"=m25t,
                  "p(Tmax)ψ(Agro)"=m26t,
                  "p(Tmax)ψ(vegh)"=m27t,
                  "p(Tmax)ψ(convergence)"=m28t,
                  "p(Tmax)ψ(cti)"=m29t,
                  "p(Tmax)ψ(slope)"=m30t,
                 # "p(Tmax)ψ(tri)"=m31t,
                  "p(Tmax)ψ(evapotrans)"=m33t,
                  "p(Tmax)ψ(productivity)"=m34t,
                  "p(Tmax)ψ(bio1)"=m35t,
                  "p(Tmax)ψ(bio5)"=m36t,
                  "p(Tmax)ψ(bio6)"=m37t,
                  "p(Tmax)ψ(bio7)"=m38t,
                  "p(Tmax)ψ(bio8)"=m39t,
                  "p(Tmax)ψ(bio12)"=m40t,
                  "p(Tmax)ψ(bio15)"=m41t,
                  "p(Tmax)ψ(bio16)"=m42t,
                  "p(Tmax)ψ(bio17)"=m43t
)
#Tmax improved many models performance

(mods_ocu.t<-aictab(mods.occu.t))
#write.table(mods_ocu.t,"graphs/mods_occu_total_septiembre.csv",sep=";",
#            row.names=F,col.names = T,fileEncoding = "UTF-8")

#m24 and m24t are the best models and consider SlopeMax as occuvar

#Parametric bootstrap

(pb1 <- parboot(m24, nsim=100, report=1))
pb1
plot(pb1, main="")
(pb2 <- parboot(m24t, nsim=100, report=1))
pb2
plot(pb2, main="")
#note: fails to reject Ho, meaning the data fit the model well

#### Creating response curves ####
conopsis2 <- unmarkedFrameOccu(y=dethist[,2:6],obsCovs = obsCovs,siteCovs =siteCovs)
summary(conopsis2)
plot(conopsis2)

(m24_units<-occu(~n.cov.obj ~slopemax, conopsis2))
(m24t_units<-occu(~Tmax ~slopemax, conopsis2))
(m25_units<-occu(~Tmax^2 ~hfp, conopsis2)) #NANs produced


newdata1<-data.frame(slopemax=seq(min(siteCovs$slopemax),
                                    max(siteCovs$slopemax),
                                    length=30),
                     hfp=seq(min(siteCovs$hfp),
                             max(siteCovs$hfp),
                             length=30))

statepredict1<-predict(m24t_units,type="state",newdata=newdata1,appendData=T)
str(statepredict1)

library(ggplot2)

ggplot(statepredict1, aes(x = slopemax, y = Predicted)) +
  ylim(0,1) +
  labs(x = "Slope Max (º)", y = "Occupancy") +
  geom_ribbon(aes(ymin = pmax(0, Predicted - SE), 
                  ymax = pmin(1, Predicted + SE)),
              alpha = 0.5,
              fill = "#cfe0e3") +
  geom_line(colour = "#2aa3bb", size = 1) +
  theme_classic()

newdata2<-data.frame(n.cov.obj=seq(min(n.cov.obj,na.rm=T),
                                  max(n.cov.obj,na.rm=T),
                                  length=30),
                     Tmax=seq(min(Tmax,na.rm=T),
                              max(Tmax,na.rm=T),
                              length=30))

statepredict2<-predict(m24t_units,type="det",newdata=newdata2,appendData=T)
str(statepredict2)

ggplot(statepredict2, aes(x = Tmax, y = Predicted)) +
  ylim(0,1) +
  labs(x = "T max (ºC)", y = "Detectability") +
  geom_ribbon(aes(ymin = pmax(0, Predicted - SE), 
                  ymax = pmin(1, Predicted + SE)),
              alpha = 0.5,
              fill = "#cfe0e3") +
  geom_line(colour = "#2aa3bb", size = 1) +
  theme_classic()

statepredict3<-predict(m24_units,type="det",newdata=newdata2,appendData=T)
str(statepredict3)

ggplot(statepredict3, aes(x = n.cov.obj, y = Predicted)) +
  ylim(0,1) +
  labs(x = "Number cover objects", y = "Detectability") +
  geom_ribbon(aes(ymin = pmax(0, Predicted - SE), 
                  ymax = pmin(1, Predicted + SE)),
              alpha = 0.5,
              fill = "#cfe0e3") +
  geom_line(colour = "#2aa3bb", size = 1) +
  theme_classic()

#### snake individual microhabitat analises ####

ind<-read.table("individual_microhabitat.csv",sep=";",header = T)
ind<-ind[-32,]
str(ind)

ind$hour <-as.POSIXct(ind$hour, format = "%I:%M:%S %p")

#analising cover object's weight 

summary(ind$weight.obj)
sd(ind$weight.obj)

#hypothesis test
summary(ind$Tmicro) #50% of individuals was founded betwwen 17 and 24°C
sd(ind$Tmicro)
shapiro.test(ind$Tmicro) #is not normal

#binomial test to probe if Tmicro within preferred temperatures 26±4°C is more than expected by chance
(inside <- sum(ind$Tmicro >= 22 & ind$Tmicro <= 30,na.rm = T))
(outside <- sum(ind$Tmicro < 22 | ind$Tmicro > 30,na.rm = T))

binom.test(inside, inside + outside, p = 0.5,alternative = "greater")

#binomial test to probe if Tmicro within preferred temperatures 24.9±4.1°C is more than expected by chance
(inside <- sum(ind$Tmicro >= 20.8 & ind$Tmicro <= 29,na.rm = T))
(outside <- sum(ind$Tmicro < 20.9 | ind$Tmicro > 29,na.rm = T))

binom.test(inside, inside + outside, p = 0.5,alternative = "greater")

#### graphs ####

#exporting corrplots

jpeg(filename = "graphs/corrplot_sitecovs_julio2025.jpeg", 
     width = 3000,        # Ancho en píxeles
     height = 3000,       # Alto en píxeles
     units = "px",        # Unidades en píxeles
     res = 300,           # Resolución en DPI (300 para alta calidad)
     quality = 90)        # Calidad JPEG (95 = muy alta)

corrplot(cormat.site,method = c("number"),type="upper",
         tl.col = "black", tl.cex = 0.6,tl.srt = 30,number.cex = 0.5,cl.cex = 0.5,
         mar=c(0,0,0,0)) 
dev.off()

#Tmax~Vegtype
T.vegtype<-data.frame(vegtype=rep(siteCovs$VegType,5),
           Tmax2=c(Tmax[,1],Tmax[,2],Tmax[,3],Tmax[,4],Tmax[,5]))
str(T.vegtype)

boxplot(T.vegtype$Tmax2~T.vegtype$vegtype)

anova_model <- aov(Tmax2 ~ vegtype, data = T.vegtype)
summary(anova_model)
TukeyHSD(anova_model)
plot(TukeyHSD(anova_model))

barplot(siteCovs$VegType)

#Tmicro  vs Tmax and min

# Creating row numbers for x-axis if dates are irregular or duplicated
ind$observation <- 1:nrow(ind)
library(ggplot2)

# Create the plot
ggplot(ind, aes(x=ind)) +
  geom_ribbon(aes(ymin=Tmin, ymax=Tmax), fill="lightblue", alpha=0.5) +
  geom_ribbon(aes(ymin=22, ymax=30), fill="#DAA520", alpha=0.5) +
  geom_line(aes(y=Tmicro), color="darkblue", size=1) +
  geom_point(aes(y=Tmicro), color="darkblue", size=2) +
  labs(x="Individual observation",
       y="Temperature (°C)") +
  annotate("text", x=max(ind$observation)*0.8, y=max(ind$Tmax)*0.9, 
           label="Tmicro = Blue line\nTmax/Tmin = Blue area", 
           hjust=1) +
  theme_minimal() +
  theme(
    axis.title = element_text(face="bold"),
    panel.grid.minor = element_blank()
  )

#is Tmicro and Tmin are related to hour?
plot(ind$Tmicro~ind$hour,xlab="Hour (24h)",ylab="Microhabitat temperature (°C)",type="p")
abline(reg<-lm(ind$Tmicro~ind$hour))
text(x=mean(ind$hour),y=25, labels="Tmicro = -7.18+ 4.11(hour)",cex=0.7)
text(x=mean(ind$hour),y=24, labels="Adjusted R-squared:  0.7458",cex=0.7)
text(x=mean(ind$hour),y=23, labels="p-value: 1.471e-11***",cex=0.7)

plot(ind$Tmin~ind$hour)
abline(reg2<-lm(ind$Tmin~ind$hour))
text(x=mean(ind$hour)+1,y=25, labels="Tmin = -3+ 1.77(hour)")
text(x=mean(ind$hour)+1,y=24, labels="Adjusted R-squared:  0.7458")
text(x=mean(ind$hour)+1,y=23, labels="p-value: 1.471e-11***") #no clear

#histograms
ggplot(siteCovs, aes(x = slopemax)) +
  geom_histogram(fill = "darkred",bins = 10) +
  labs(x = "Slope max", y = "Frecuency")+
  theme_minimal() 

ggplot(siteCovs, aes(x = cti)) +
  geom_histogram(fill = "#FF7F00",bins=10) +
  labs(x = "Compound Topographic Index", y = "Frecuency")+
  theme_minimal() 

ggplot(siteCovs, aes(x = vegtype)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Vegetation type", y = "Frecuency")+
  theme_minimal() 

ggsave(filename="hist_vegtype.tiff",
       path = "graphs/",
       plot = last_plot(),
       width =132, 
       height = 66, 
       dpi = 300,
       units = "mm")
#____________________________________________

ind$hour<-as.POSIXct(ind$hour, format = "%I:%M:%S %p")

hist(sitecovs$slope,xlab = "Slope (%)",main = "",col = "brown")

hist(ind$weight.obj,xlab = "Rock weight (kg)", main = "",col="darkgrey")
ggplot(ind, aes(x = weight.obj)) +
  geom_histogram(fill = "darkgrey",bins=7) +
  labs(x = "Rock weight (kg)", y = "Frecuency")+
  theme_minimal()

hist(ind$hour,breaks = 8,xlab = "Hour",main ="",col="lightblue")

hist(ind$Tmicro,xlab="Microhabitat temperature",main = "")
polygon(x=c(22.3, 29.5, 29.5, 22.3), y=c(0, 0, 10, 10)
        , col = rgb(red = 218/255, green = 165/255, blue = 32/255, alpha = 0.5)
        , border = NA, lwd = 2) #Figueroa et al.,2024
polygon(x=c(20.8, 29, 29, 20.8), y=c(0, 0, 10, 10)
        , col = rgb(red = 0, green = 0.2, blue = 0.8, alpha = 0.5)
        , border = NA, lwd = 2) #Castañeda-Gonzalez et al., 2011

plot(siteCovs.str$slopemax~siteCovs.str$slope)
cor.test(siteCovs.str$slopemax,siteCovs.str$slope) #there are correlation, both measures are congruent

cor.test(ind$Tmicro,ind$weight.obj)
plot(ind$Tmicro~ind$weight.obj)

par(las=2)
heatmap(as.matrix(dethist[,2:6]),scale="none",Rowv = NA, Colv = NA,
        col=c("lightgrey","black"),cexCol = 1.2)



ggsave(filename="predicted_det-ncovobj_julio.tiff",
       path = "graphs/",
       plot = last_plot(),
       width =66, 
       height = 66, 
       dpi = 300,
       units = "mm")

