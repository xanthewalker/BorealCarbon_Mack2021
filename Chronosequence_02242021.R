library(nlme)
library(AICcmodavg)
library(dplyr)
library(vegan)
library(MuMIn)
library(MASS)
library(ggplot2)
library(ggeffects)
library(gtools)
library(effects)
library(scales)
library(emmeans)
rm(list=ls())

### 
number_ticks <- function(n) {function(limits) pretty(limits, n)}

setwd("C:\\Users\\xjw5\\Documents\\NAU\\Lab Materials\\Michelle_Misc\\Archived Data")
#####################
#
data=read.csv(file="Chronosequence_02242021.csv")

data=mutate(data, above.below=aboveC/belowC)


names(data)

data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))
data$Age<-factor(data$Age.Class, levels=c("Seedlings", "> 20 years"))

############## new data
g=ggplot()+geom_point(data=data, size=3,  aes(x=rel.decid.dens*100, y=rel.decid.bio*100,fill=traj.20), pch=21)
g=g + theme_bw() + labs(x="Relative Deciduous Density (%)", y ="Proportion Deciduous Biomass (%)")
g=g+stat_smooth(data=data, aes(x=rel.decid.dens*100, y=rel.decid.bio*100),color="black", se=T, alpha=0.2, method=lm, formula = y ~ x + I(x^2), size = 1) 
#g = g + geom_ribbon(data=v, aes(x=rel.decid.dens,ymin=lower,ymax=upper),  alpha=0.25, size=1)
#g=g+geom_line(data=v, aes(y=fit, x=rel.decid.dens),color="black", size=1.3)
g = g +scale_fill_manual(values=c("#443983FF" ,"#77D153FF","gold1")) +scale_color_viridis_d()
g = g+ scale_y_continuous(breaks=c(0,10,20 , 30, 40, 50, 60, 70, 80, 90, 100))
g= g+ scale_x_continuous(breaks=c(0,10,20 , 30, 40, 50, 60, 70, 80, 90, 100))
g = g +facet_wrap(~Age) + labs(fill = "Trajectory") +theme(legend.position = "top")
g

data1=mutate(data, outlier=ifelse(treesnagsolC>14000, "outlier", "not"))

names(data1)

data2=filter(data1, outlier=="outlier")


############
# total C
mod20=lme(data=data1, treesnagsolC ~ age * traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod20)
plot(mod20)


mod21=lme(data=data1, treesnagsolC ~ age + traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod21)
plot(mod21)

anova(mod20, mod21)

mod20=lme(data=data1, treesnagsolC ~ age * traj.20,  random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod20)
plot(mod20)
r.squaredGLMM(mod20)

g=ggplot(data=data1, aes(y=treesnagsolC, x=traj.20)) +geom_boxplot()
g

g=ggplot(data=data1, aes(y=aboveC, x=traj.20)) +geom_boxplot()
g

g=ggplot(data=data1, aes(y=belowC, x=traj.20)) +geom_boxplot()
g

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod20, add.smooth = FALSE, which = 1)
E <- resid(mod20)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

par(op)
plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")




emm_options(opt.digits = FALSE)
head(data)
emt.intercept1 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=0))
summary(emt.intercept1,infer=c(TRUE,TRUE))
pairs(emt.intercept1)
inter1=as.data.frame(emt.intercept1)
head(inter1)
emt.intercept100 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=100))
summary(emt.intercept100,infer=c(TRUE,TRUE))
pairs(emt.intercept100)
inter100=as.data.frame(emt.intercept100)
head(inter100)


emt.slope1 <- emtrends(mod20, "traj.20", var = "age")
summary(emt.slope1,infer=c(TRUE,TRUE))
pairs(emt.slope1)
slope1=as.data.frame(emt.slope1)
head(slope1)
emmip(mod20, traj.20 ~ age, cov.reduce = range, CIs=T)
summary(mod20)

interact1 <- Effect(c("age", "traj.20"),mod20, xlevels=100) 
names(interact1)
interact1<-as.data.frame(interact1)
head(interact1)
interact1$traj.20<-factor(interact1$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))
data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))


########
#######above C
mod20=lme(data=data1, aboveC ~ age * traj.20, weights=varComb(varExp(form=~age), varIdent(form=~1|traj.20)), random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod20)
plot(mod20)


mod21=lme(data=data1, aboveC ~ age + traj.20, weights=varComb(varExp(form=~age), varIdent(form=~1|traj.20)),random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod21)
plot(mod21)

anova(mod20, mod21)

mod20=lme(data=data1, aboveC ~ age * traj.20, weights=varComb(varExp(form=~age), varIdent(form=~1|traj.20)),random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod20)
plot(mod20)
r.squaredGLMM(mod20)


op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod20, add.smooth = FALSE, which = 1)
E <- resid(mod20)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")
par(op)


library(emmeans)
emm_options(opt.digits = FALSE)
head(data)
emt.intercept1 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=0), type="response")
summary(emt.intercept1,infer=c(TRUE,TRUE))
pairs(emt.intercept1)
inter1=as.data.frame(emt.intercept1)
head(inter1)
emt.intercept100 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=100), type="response")
summary(emt.intercept100,infer=c(TRUE,TRUE))
pairs(emt.intercept100)
inter100=as.data.frame(emt.intercept100)
head(inter100)


library(emmeans)
emt.slope1 <- emtrends(mod20, "traj.20", var = "age", type="response")
summary(emt.slope1,infer=c(TRUE,TRUE))
pairs(emt.slope1)
slope1=as.data.frame(emt.slope1)
head(slope1)
emmip(mod20, traj.20 ~ age, cov.reduce = range, CIs=T, type="response")
summary(mod20)
library(ggeffects)
library(gtools)
library(effects)
interact1a <- Effect(c("age", "traj.20"),mod20, xlevels=100, back.transform = TRUE) 
names(interact1a)
interact1a<-as.data.frame(interact1a)
head(interact1a)
interact1a$traj.20<-factor(interact1a$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))
data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))

#

##########
## below C
mod20=lme(data=data1, belowC ~ age * traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod20)
plot(mod20)


mod21=lme(data=data1, belowC ~ age + traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod21)
plot(mod21)

anova(mod20, mod21)

mod20=lme(data=data1, belowC ~ age + traj.20, random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod20)
plot(mod20)
r.squaredGLMM(mod20)


op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod20, add.smooth = FALSE, which = 1)
E <- resid(mod20)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")
par(op)


library(emmeans)
emm_options(opt.digits = FALSE)
head(data)
emt.intercept1 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=0))
summary(emt.intercept1,infer=c(TRUE,TRUE))
pairs(emt.intercept1)
inter1=as.data.frame(emt.intercept1)
head(inter1)
emt.intercept100 <- emmeans(mod20, ~ traj.20, var = "age",adjust="tukey", at=list(age=100))
summary(emt.intercept100,infer=c(TRUE,TRUE))
pairs(emt.intercept100)
inter100=as.data.frame(emt.intercept100)
head(inter100)


library(emmeans)
emt.slope1 <- emtrends(mod20, "traj.20", var = "age")
summary(emt.slope1,infer=c(TRUE,TRUE))
pairs(emt.slope1)
slope1=as.data.frame(emt.slope1)
head(slope1)
emmip(mod20, traj.20 ~ age, cov.reduce = range, CIs=T)
summary(mod20)
library(ggeffects)
library(gtools)
library(effects)
interact1b <- Effect(c("age", "traj.20"),mod20, xlevels=100) 
names(interact1b)
interact1b<-as.data.frame(interact1b)
head(interact1b)
interact1b$traj.20<-factor(interact1b$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))
data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))

####### 
#nitrogen
names(data)
levels(data$type2)

mod1=lme(data=data, treesnagsolN ~ age * traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod1)
plot(mod1)
r.squaredGLMM(mod1)

mod10=lme(data=data, treesnagsolN ~ age + traj.20, random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod10)

anova(mod1, mod10)

###
mod2=lme(data=data, treesnagsolN ~ age * traj.20, random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod2)
plot(mod2)
r.squaredGLMM(mod2)



op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod2, add.smooth = FALSE, which = 1)
E <- resid(mod2)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")
par(op)




library(emmeans)
head(data)
emt.intercept2 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=0))
summary(emt.intercept2,infer=c(TRUE,TRUE))
pairs(emt.intercept2)
inter2=as.data.frame(emt.intercept2)
head(inter2)

emt.intercept200 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=100))
summary(emt.intercept200,infer=c(TRUE,TRUE))
pairs(emt.intercept200)
inter200=as.data.frame(emt.intercept200)
head(inter200)


emt.slope2= emtrends(mod2, "traj.20", var = "age")
summary(emt.slope2,infer=c(TRUE,TRUE))
pairs(emt.slope2)
slope2=as.data.frame(emt.slope2)
emmip(mod2, traj.20 ~ age, cov.reduce = range, CIs=T)
summary(mod2)


interact2 <- Effect(c("age", "traj.20"),mod2, xlevels=100) 
names(interact2)
interact2<-as.data.frame(interact2)
head(interact2)
interact2$traj.20<-factor(interact2$traj.20, levels=c("Spruce", "Mixed", "Deciduous"))


#############
####### C:N

names(data)
mod2=lme(data=data,treesnagsolCN ~ (age) * traj.20,  random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod2)
plot(mod2)

mod3=lme(data=data,treesnagsolCN ~ (age) + traj.20,  random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod3)
plot(mod3)

anova(mod2, mod3)

mod3=lme(data=data,treesnagsolCN ~ (age) * traj.20, weights=varIdent(form=~1|traj.20), random=~1|project/fire, na.action=na.omit, method="ML")
summary(mod3)
AIC(mod3)
plot(mod3)

mod2=lme(data=data,treesnagsolCN ~ (age) +traj.20, weights=varIdent(form=~1|traj.20),
         random=~1|project/fire, na.action=na.omit, method="REML")

plot(mod2)
anova(mod2, mod3)


summary(mod2)
plot(mod2)
AIC(mod2)
r.squaredGLMM(mod2)


anova(mod1, mod2)


op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod2, add.smooth = FALSE, which = 1)
E <- resid(mod2)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")
par(op)


head(data)
emt.intercept3 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=0))
summary(emt.intercept3,infer=c(TRUE,TRUE))
pairs(emt.intercept3)
inter3=as.data.frame(emt.intercept3)
head(inter3)

head(data)
emt.intercept300 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=100))
summary(emt.intercept300,infer=c(TRUE,TRUE))
pairs(emt.intercept300)
inter300=as.data.frame(emt.intercept300)
head(inter300)


emt.slope3= emtrends(mod2, "traj.20", var = "age")
summary(emt.slope3,infer=c(TRUE,TRUE))
pairs(emt.slope3)
slope3=as.data.frame(emt.slope3)
emmip(mod2, traj.20 ~ age, cov.reduce = range, CIs=T)
summary(mod2)


interact3<- Effect(c("age", "traj.20"),mod2, xlevels=100) 
names(interact3)
interact3<-as.data.frame(interact3)
head(interact3)
interact3$traj.20<-factor(interact3$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))
levels(interact3$type2)

##########
##above:below
names(data)

mod1=lme(data=data, above.below ~ age * traj.20, weights=varIdent(form=~1|traj.20),random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod1)
plot(mod1)
r.squaredGLMM(mod1)

mod3=lme(data=data, above.below ~ age + traj.20, weights=varIdent(form=~1|traj.20),random=~1|project/fire, na.action=na.omit, method="REML")
summary(mod3)
plot(mod3)
r.squaredGLMM(mod3)


mod2=lme(data=data,above.below ~ (age) * traj.20, weights=varIdent(form=~1|traj.20),
         random=~1|project/fire, na.action=na.omit, method="REML")


summary(mod2)
AIC(mod2)
summary(mod2)
plot(mod2)
r.squaredGLMM(mod2)


op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod2, add.smooth = FALSE, which = 1)
E <- resid(mod2)
hist(E, xlab = "Residuals", main = "")
plot(data$age, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire, E, 
     ylab = "Residuals")
plot(data$project, E, 
     ylab = "Residuals")
par(op)




head(data)
emt.intercept4 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=0))
summary(emt.intercept4,infer=c(TRUE,TRUE))
pairs(emt.intercept4)
inter4=as.data.frame(emt.intercept4)
head(inter4)

head(data)
emt.intercept400 <- emmeans(mod2, ~ traj.20, var = "age",adjust="tukey", at=list(age=100))
summary(emt.intercept400,infer=c(TRUE,TRUE))
pairs(emt.intercept400)
inter400=as.data.frame(emt.intercept400)
head(inter400)


emt.slope4 <- emtrends(mod2, "traj.20", var = "age")
summary(emt.slope4,infer=c(TRUE,TRUE))
pairs(emt.slope4)
slope4=as.data.frame(emt.slope4)
emmip(mod1, traj.20 ~ age, cov.reduce = range, CIs=T)
summary(mod2)


interact4 <- Effect(c("age", "traj.20"),mod2, xlevels=100) 
names(interact4)
interact4<-as.data.frame(interact4)
head(interact4)
interact4$traj.20<-factor(interact4$traj.20, levels=c( "Deciduous", "Mixed", "Spruce"))


######
#figure
names(data)
dat2<-dplyr::select(data, site, traj.20, treesnagsolC, treesnagsolN, treesnagsolCN, above.below, aboveC, belowC, age)

library(reshape)
mdata <- melt(dat2, id=c("site","traj.20", "age"))
head(mdata)
names(mdata)

library(ggplot2)
names(new7)
str(new7)


#devtools::install_github("zeehio/facetscales")
library(g)
library(facetscales)

my_breaks=mutate(mdata, ifelse)
scales_y <- list(
  '\"Carbon (g C m\"^-2 * \")\"' = scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 2500)),
  '\"Above Carbon (g C m\"^-2 * \")\"' = scale_y_continuous(limits = c(0, 13000), breaks = seq(0, 14000, 2000)),
  '\"Below Carbon (g C m\"^-2 * \")\"' = scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 12000, 2000)),
  '\"Nitrogen (g N m\"^-2 * \")\"' = scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)),
  'C:N' = scale_y_continuous(limits = c(20, 140), breaks = seq(20, 140,  20)), 
  'Aboveground:Belowground'=scale_y_continuous(limits=c(-1, 6), breaks=seq(0, 6, 1)))
head(scales_y)

names(mdata)
names(interact1)
head(interact1)
interact1=mutate(interact1, variable="Carbon")
head(interact1a)
interact1a=mutate(interact1a, variable="Aboveground Carbon")
head(interact1b)
interact1b=mutate(interact1b, variable="Belowground Carbon")
head(interact2)
interact2=mutate(interact2, variable="Nitrogen")
interact3=mutate(interact3, variable="C:N")
interact4=mutate(interact4, variable="Above:Below")
head(interact4)

interact.all=rbind(interact1, interact1a)
head(interact.all)
interact.all=rbind(interact.all, interact1b)
interact.all=rbind(interact.all, interact2)
interact.all=rbind(interact.all, interact3)
interact.all=rbind(interact.all, interact4)

str(interact.all)
interact.all$variable=as.factor(interact.all$variable)

levels(interact.all$variable)[levels(interact.all$variable)=="Carbon"] <- "\"Carbon (g C m\"^-2 * \")\""
levels(interact.all$variable)[levels(interact.all$variable)=="Aboveground Carbon"] <- "\"Above Carbon (g C m\"^-2 * \")\""
levels(interact.all$variable)[levels(interact.all$variable)=="Belowground Carbon"] <- "\"Below Carbon (g C m\"^-2 * \")\""
levels(interact.all$variable)[levels(interact.all$variable)=="Nitrogen"] <- "\"Nitrogen (g N m\"^-2 * \")\""
levels(interact.all$variable)[levels(interact.all$variable)=="C:N"]<- "C:N"
levels(interact.all$variable)[levels(interact.all$variable)=="Above:Below"]= "\"Above C : Below C\""


levels(mdata$variable)[levels(mdata$variable)=="treesnagsolC"] <- "\"Carbon (g C m\"^-2 * \")\""
levels(mdata$variable)[levels(mdata$variable)=="aboveC"] <- "\"Above Carbon (g C m\"^-2 * \")\""
levels(mdata$variable)[levels(mdata$variable)=="belowC"] <- "\"Below Carbon (g C m\"^-2 * \")\""
levels(mdata$variable)[levels(mdata$variable)=="treesnagsolN"] <- "\"Nitrogen (g N m\"^-2 * \")\""
levels(mdata$variable)[levels(mdata$variable)=="treesnagsolCN"]<- "C:N"
levels(mdata$variable)[levels(mdata$variable)=="above.below"]= "\"Above C : Below C\""

levels(interact.all$variable)
levels(mdata$variable)


dat_text <- data.frame(
  label = c("A", "B", "C", "E", "F", "D"),
  variable   = c("\"Carbon (g C m\"^-2 * \")\"" , "\"Above Carbon (g C m\"^-2 * \")\"" ,"\"Below Carbon (g C m\"^-2 * \")\"" , "\"Nitrogen (g N m\"^-2 * \")\"", "C:N", "\"Above C : Below C\"")
)

mdata=mdata%>%
  group_by(variable)%>%
  mutate(y_max=ifelse(variable=="C:N", 140, 
                      ifelse(variable=="\"Above C : Below C\"", 6, 
                             ifelse(variable=="\"Carbon (g C m\"^-2 * \")\"", 20000, 
                                    ifelse(variable=="\"Above Carbon (g C m\"^-2 * \")\"", 14000, 
                                           ifelse(variable=="\"Below Carbon (g C m\"^-2 * \")\"", 12000, 
                                                  500))))), 
         y_min=ifelse(variable=="C:N", 0, 
                      ifelse(variable=="\"Above C : Below C\"", 0, 
                             ifelse(variable=="\"Carbon (g C m\"^-2 * \")\"", 0, 0))))
head(mdata)
head(interact.all)
levels(interact.all$variable)
levels(mdata$traj.20)
library(scales)
interact.all$traj.20<-factor(interact.all$traj.20, levels=c("Spruce",  "Mixed", "Deciduous"))
mdata$traj.20<-factor(mdata$traj.20, levels=c("Spruce",  "Mixed", "Deciduous"))
interact.all$variable<-factor(interact.all$variable, levels=c("\"Carbon (g C m\"^-2 * \")\"", "\"Above Carbon (g C m\"^-2 * \")\"", "\"Below Carbon (g C m\"^-2 * \")\"",  "\"Nitrogen (g N m\"^-2 * \")\"", 
                                                              "C:N", "\"Above C : Below C\""))
levels(interact.all$variable)

interact.all$variable<-factor(interact.all$variable, 
                              levels=c("\"Carbon (g C m\"^-2 * \")\"", 
                                       "\"Above Carbon (g C m\"^-2 * \")\"", 
                                       "\"Below Carbon (g C m\"^-2 * \")\"",
                                       "\"Above C : Below C\"",
                                       "\"Nitrogen (g N m\"^-2 * \")\"", 
                                       "C:N"))

dat_text$variable<-factor(dat_text$variable, 
                          levels=c("\"Carbon (g C m\"^-2 * \")\"", 
                                   "\"Above Carbon (g C m\"^-2 * \")\"", 
                                   "\"Below Carbon (g C m\"^-2 * \")\"",
                                   "\"Above C : Below C\"",
                                   "\"Nitrogen (g N m\"^-2 * \")\"", 
                                   "C:N"))

mdata$variable<-factor(mdata$variable, 
                       levels=c("\"Carbon (g C m\"^-2 * \")\"", 
                                "\"Above Carbon (g C m\"^-2 * \")\"", 
                                "\"Below Carbon (g C m\"^-2 * \")\"",
                                "\"Above C : Below C\"",
                                "\"Nitrogen (g N m\"^-2 * \")\"", 
                                "C:N"))

p=ggplot()
p=p+ geom_point(data=mdata,aes(x=age,y=value, fill=traj.20, shape=traj.20), size=1.5, color="black") 
p=p+facet_wrap(~variable, scales="free_y", strip.position="left", labeller = label_parsed) + 
  theme(strip.background = element_blank(),  strip.placement = "outside") 
p


p=ggplot()

p=p+ geom_line(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Spruce"),
               aes(x=age,y=fit), color="#443983FF",size=0.5)
p=p+ geom_line(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Mixed"&age<125),
               aes(x=age,y=fit) ,  color="#77D153FF",size=0.5) 
p=p+ geom_line(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Deciduous"&age<125),
               aes(x=age,y=fit) ,  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Spruce"), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#443983FF",alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Nitrogen (g N m\"^-2 * \")\""&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="gold1", alpha=0.3)

p=p+ geom_line(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"),
               aes(x=age,y=fit), color="#443983FF",size=0.5)
p=p+ geom_line(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125),
               aes(x=age,y=fit) ,  color="#77D153FF",size=0.5) 
p=p+ geom_line(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125),
               aes(x=age,y=fit) ,  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#443983FF",alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="gold1", alpha=0.3)

p=p+ geom_line(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"),
               aes(x=age,y=fit), color="#443983FF",size=0.5)
p=p+ geom_line(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125),
               aes(x=age,y=fit) ,  color="#77D153FF",size=0.5) 
p=p+ geom_line(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125),
               aes(x=age,y=fit) ,  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#443983FF",alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="gold1", alpha=0.3)

p=p+ geom_line(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"),
               aes(x=age,y=fit), color="#443983FF",size=0.5)
p=p+ geom_line(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125),
               aes(x=age,y=fit) ,  color="#77D153FF",size=0.5) 
p=p+ geom_line(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125),
               aes(x=age,y=fit) ,  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Spruce"), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#443983FF",alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Below Carbon (g C m\"^-2 * \")\""&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="gold1", alpha=0.3)


p=p+ geom_line(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Deciduous"&age<125), 
               aes(x=age,y=fit),  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="gold1", alpha=0.3)
p=p+ geom_line(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Mixed"&age<125), 
               aes(x=age,y=fit),  color="#77D153FF",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+ geom_line(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Spruce"), 
               aes(x=age,y=fit),linetype="dotted", color="#443983FF",size=0.5)
p=p+geom_ribbon(data=subset(interact.all,variable=="\"Above C : Below C\""&traj.20=="Spruce"), 
                aes(x=age, ymin=lower ,ymax=upper),fill="#443983FF", alpha=0.3)

p=p+ geom_line(data=subset(interact.all,variable=="C:N"&traj.20=="Deciduous"&age<125), 
               aes(x=age,y=fit) , linetype="dotted",  color="gold1",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="C:N"&traj.20=="Deciduous"&age<125), 
                aes(x=age, ymin=lower,ymax=upper),fill="gold1", alpha=0.3)
p=p+ geom_line(data=subset(interact.all,variable=="C:N"&traj.20=="Mixed"&age<125), 
               aes(x=age,y=fit) , linetype="dotted",  color="#77D153FF",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="C:N"&traj.20=="Mixed"&age<125), 
                aes(x=age, ymin=lower,ymax=upper),fill="#77D153FF", alpha=0.3)
p=p+ geom_line(data=subset(interact.all,variable=="C:N"&traj.20=="Spruce"), 
               aes(x=age,y=fit),linetype="dotted",  color="#443983FF",size=0.5) 
p=p+geom_ribbon(data=subset(interact.all,variable=="C:N"&traj.20=="Spruce"), 
                aes(x=age, ymin=lower,ymax=upper),fill="#443983FF", alpha=0.3)

p=p+ geom_point(data=mdata,aes(x=age,y=value, fill=traj.20, shape=traj.20), size=1.8, stroke=0, color="black") 


p=p+scale_x_continuous(breaks=pretty_breaks(6)) +scale_y_continuous(breaks=pretty_breaks(7))+theme_bw()
p=p+scale_shape_manual(name="Regeneration Trajectory", values=c(21, 22, 24))
p=p+scale_fill_manual(name="Regeneration Trajectory", values=c("#443983FF" ,"#77D153FF", "gold1")) #+geom_blank(data=mdata, aes(y = y_max)) +geom_blank(data=mdata, aes(y=y_min))
p=p+scale_color_manual(name="Regeneration Trajectory", values=c("#443983FF" ,"#77D153FF", "gold1"), guide=F)
p=p+xlab("Time After Fire (years)") + theme(strip.text.y = element_text(size = 8)) +theme(legend.title=element_text(size=8))
p=p + geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf, label = label),  hjust   = -0.3,  vjust   = 1.2, size=4, fontface=2)
p=p+facet_wrap(~variable, scales="free_y", strip.position="left", labeller = label_parsed, ncol=2, dir="h") + 
  theme(strip.background = element_blank(),  strip.placement = "outside") 
p=p+theme(legend.position="top") + theme(legend.text=element_text(size=8)) + 
  theme(legend.background = element_rect(fill="transparent"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.25))+
  theme(legend.margin=margin(0.15,0,-0.15,0), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.0, 'cm'))+
  theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.18))+
  theme(strip.text.y = element_text(margin = margin(-1,-0.05,-1,-1), "cm")) + theme(plot.margin = unit(c(0.15,0.15,0.05,0), "cm"))+
  theme(legend.box.margin=margin(0,-1,-5,-1))+
  theme(axis.text.x  = element_text(size=7)) +theme(axis.text.y  = element_text(size=7))+
  theme(axis.title.x  = element_text(size=8, vjust=2))+theme(axis.title.y  = element_blank())+
  geom_blank(data=mdata, aes(y = y_max)) +geom_blank(data=mdata, aes(y=y_min))
p

ggsave("Figure_3.pdf", plot = p, width =12, height = 15, units = "cm",
       dpi = 600 , device=cairo_pdf)




