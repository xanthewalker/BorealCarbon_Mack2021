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
library(extrafont)
library(reshape)
library(plotrix)
library(multcompView)
loadfonts(device = "win")

rm(list=ls())
### 
number_ticks <- function(n) {function(limits) pretty(limits, n)}

setwd("C:\\Users\\xjw5\\Documents\\NAU\\Lab Materials\\Michelle_Misc\\Archived Data")

#################### 
#####

data<-read.csv("JFSP_02242021.csv")
names(data)


data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed", "Deciduous"))


mod2=lme(data=data,postsol.c ~ presol.c * traj.20, random=~1|fire.x, na.action=na.omit, method="ML")
AIC(mod2)
plot(mod2)
summary(mod2)
r.squaredGLMM(mod2)

mod1=lme(data=data,postsol.c ~ presol.c + traj.20, random=~1|fire.x, na.action=na.omit, method="ML")
anova(mod1, mod2)


mod2=lme(data=data,postsol.c ~ presol.c * traj.20 ,random=~1|fire.x, na.action=na.omit, method="REML")
summary(mod2)
r.squaredGLMM(mod2)
plot(mod2)



mod78=lme(data=subset(data, !site=="BF78"), postsol.c ~ presol.c * traj.20 ,random=~1|fire.x, na.action=na.omit, method="REML")
summary(mod78)
r.squaredGLMM(mod78)
plot(mod78)

emt.slope<- emtrends(mod78, "traj.20", var = "presol.c")
summary(emt.slope,infer=c(TRUE,TRUE))
pairs(emt.slope)



op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(mod2, add.smooth = FALSE, which = 1)
E <- resid(mod2)
hist(E, xlab = "Residuals", main = "")
plot(data$presol.c, E, 
     ylab = "Residuals")
g=ggplot(data, aes(x=traj.20, y=E)) + geom_boxplot()
g

plot(data$fire.x, E, 
     ylab = "Residuals")

par(op)
################
###########
library(emmeans)
head(data)
emt.intercept <- emmeans(mod2, ~ traj.20, var = "presol.c",adjust="tukey", at=list(presol.c=0))
summary(emt.intercept,infer=c(TRUE,TRUE))
pairs(emt.intercept)
inter=as.data.frame(emt.intercept)
head(inter)
emt.slope<- emtrends(mod2, "traj.20", var = "presol.c")
summary(emt.slope,infer=c(TRUE,TRUE))
pairs(emt.slope)
slope=as.data.frame(emt.slope)
emmip(mod2, traj.20 ~ presol.c, cov.reduce = range, CIs=T)
summary(mod2)


interact <- Effect(c("presol.c", "traj.20"),mod2, xlevels=1000) 
names(interact)
interact<-as.data.frame(interact)
head(interact)
interact$traj.20<-factor(interact$traj.20, levels=c( "Spruce", "Mixed" ,"Deciduous"))

spruce=lme(data=subset(data, traj.20=="Spruce"),postsol.c ~ presol.c, random=~1|fire.x, na.action=na.omit, method="ML")
r.squaredGLMM(spruce)
mixed=lme(data=subset(data, traj.20=="Mixed"),postsol.c ~ presol.c, random=~1|fire.x, na.action=na.omit, method="ML")
r.squaredGLMM(mixed)
Deciduous=lme(data=subset(data, traj.20=="Deciduous"),postsol.c ~ presol.c, random=~1|fire.x, na.action=na.omit, method="ML")
r.squaredGLMM(Deciduous)

interact1=filter(interact, traj.20=="Spruce"&presol.c>2999)
interact2=filter(interact, traj.20=="Mixed"&presol.c>2999)
interact3=filter(interact, traj.20=="Deciduous")

interact4=rbind(interact1, interact2, interact3)
interact4$traj.20<-factor(interact4$traj.20, levels=c( "Spruce", "Mixed" ,"Deciduous"))
head(interact4)
head(interact)
interact5=subset(interact, fit>0)

int=mutate(interact, new.upper=ifelse(upper<0, 0, upper))
int=mutate(int, new.lower=ifelse(lower<0, 0, lower))
int=mutate(int, new.fit=ifelse(fit<0, NA, fit))


data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed" ,"Deciduous"))

#####

##########PRE C
names(data)

data2=dplyr::select(data, traj.20, preabove.c, prebelow.c, pretot.c.new, 
                    preabove.n, prebelow.n, pretot.n.new, 
                    preabove.cn, prebelow.cn, pretot.cn.new, 
                    postabove.c, postbelow.c, posttot.c.new, 
                    postabove.n, postbelow.n, posttot.n.new, 
                    postabove.cn, postbelow.cn, posttot.cn.new, 
                    lostabove.c, lostbelow.c, losttot.c.new, 
                    lostabove.n, lostbelow.n, losttot.n.new, 
                    lostabove.cn, lostbelow.cn, losttot.cn.new)


data3=melt(data2, id="traj.20")
head(data3)
dat.sum=data3%>%
  dplyr::group_by(variable, traj.20)%>%
  dplyr::summarise(mean=mean(value), 
            std.error=std.error(value))

data$traj.20<-factor(data$traj.20, levels=c("Spruce", "Mixed", "Deciduous"))


###preabove C
preabove.c <- lme(preabove.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(preabove.c)
plot(preabove.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(preabove.c, add.smooth = FALSE, which = 1)
E <- resid(preabove.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
           ylab = "Residuals")
par(op)


leastsquare = emmeans(preabove.c,
                      pairwise ~ traj.20,
                      adjust="tukey")
CLD = CLD(leastsquare, alpha=0.05, Letters=letters)
CLD


#### pre below C
prebelow.c <- lme(prebelow.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(prebelow.c)
plot(prebelow.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(prebelow.c, add.smooth = FALSE, which = 1)
E <- resid(prebelow.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


leastsquare = emmeans(prebelow.c,
                      pairwise ~ traj.20,
                      adjust="tukey")
CLD = CLD(leastsquare, alpha=0.05, Letters=letters)
CLD


#### pretot C
pretot.c.new <- lme(pretot.c.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(pretot.c.new)
plot(pretot.c.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(pretot.c.new, add.smooth = FALSE, which = 1)
E <- resid(pretot.c.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(pretot.c.new,
                      pairwise ~ traj.20,
                      adjust="tukey")

##################POST C

### posta above C
postabove.c <- lme(postabove.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postabove.c)
plot(postabove.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postabove.c, add.smooth = FALSE, which = 1)
E <- resid(postabove.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(postabove.c,pairwise ~ traj.20,
                      adjust="tukey")

#### pot below C
postbelow.c <- lme(postbelow.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postbelow.c)
plot(postbelow.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postbelow.c, add.smooth = FALSE, which = 1)
E <- resid(postbelow.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


leastsquare = emmeans(postbelow.c,
                      pairwise ~ traj.20,
                      adjust="tukey")
summary(leastsquare)


###post tot C
posttot.c.new <- lme(posttot.c.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(posttot.c.new)


plot(posttot.c.new)
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(posttot.c.new, add.smooth = FALSE, which = 1)
E <- resid(posttot.c.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(posttot.c.new,
                      pairwise ~ traj.20,
                      adjust="tukey")


################## LOST C

###lost above C
lostabove.c <- lme(lostabove.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostabove.c)
plot(lostabove.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostabove.c, add.smooth = FALSE, which = 1)
E <- resid(lostabove.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(lostabove.c,
                      pairwise ~ traj.20,
                      adjust="tukey")


####lost below C

lostbelow.c <- lme(lostbelow.c ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostbelow.c)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostbelow.c, add.smooth = FALSE, which = 1)
E <- resid(lostbelow.c)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(lostbelow.c, pairwise ~ traj.20,
                      adjust="tukey")


#### l0st total

losttot.c.new <- lme(losttot.c.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(losttot.c.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(losttot.c.new, add.smooth = FALSE, which = 1)
E <- resid(losttot.c.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


emmeans(losttot.c.new,
                      pairwise ~ traj.20,
                      adjust="tukey")



##########PRE N
### pre above N

names(data)
preabove.n <- lme(preabove.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(preabove.n)
plot(preabove.n)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(preabove.n, add.smooth = FALSE, which = 1)
E <- resid(preabove.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(preabove.n, pairwise ~ traj.20,
                      adjust="tukey")

### pre below N
prebelow.n <- lme(prebelow.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(prebelow.n)
plot(prebelow.n)


op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(prebelow.n, add.smooth = FALSE, which = 1)
E <- resid(prebelow.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)



emmeans(prebelow.n,pairwise ~ traj.20,
                      adjust="tukey")

#### pre tot N

pretot.n.new <- lme(pretot.n.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(pretot.n.new)
plot(pretot.n.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(pretot.n.new, add.smooth = FALSE, which = 1)
E <- resid(pretot.n.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(pretot.n.new, pairwise ~ traj.20,
                      adjust="tukey")

##################POST N

#post above N
postabove.n <- lme(postabove.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postabove.n)
plot(postabove.n)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postabove.n, add.smooth = FALSE, which = 1)
E <- resid(postabove.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(postabove.n, pairwise ~ traj.20,
                      adjust="tukey")

# post below N
postbelow.n <- lme(postbelow.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postbelow.n)
plot(postbelow.n)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postbelow.n, add.smooth = FALSE, which = 1)
E <- resid(postbelow.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(postbelow.n,pairwise ~ traj.20,
                      adjust="tukey")

# post total N

post.total.n.new <- lme(posttot.n.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(post.total.n.new)
plot(post.total.n.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(post.total.n.new, add.smooth = FALSE, which = 1)
E <- resid(post.total.n.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(post.total.n.new, pairwise ~ traj.20,
                      adjust="tukey")



################## LOST N

#lost above N
lostabove.n <- lme(lostabove.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostabove.n)
plot(lostabove.n)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostabove.n, add.smooth = FALSE, which = 1)
E <- resid(lostabove.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(lostabove.n,  pairwise ~ traj.20,
                      adjust="tukey")


names(data)

# lost below N
lostbelow.n <- lme(lostbelow.n ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostbelow.n)
plot(lostbelow.n)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostbelow.n, add.smooth = FALSE, which = 1)
E <- resid(lostbelow.n)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(lostbelow.n,  pairwise ~ traj.20,
                      adjust="tukey")


#lost total N
losttotal.n.new <- lme(losttot.n.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(losttotal.n.new)
plot(losttotal.n.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(losttotal.n.new, add.smooth = FALSE, which = 1)
E <- resid(losttotal.n.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(losttotal.n.new,
                      pairwise ~ traj.20,
                      adjust="tukey")



##########PRE CN
names(data)

# pre above CN
preabove.cn <- lme(preabove.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(preabove.cn)
plot(preabove.cn)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(preabove.cn, add.smooth = FALSE, which = 1)
E <- resid(preabove.cn)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)



emmeans(preabove.cn,pairwise ~ traj.20,
                      adjust="tukey")

#prebelow CN 

prebelow.cn <- lme(prebelow.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(prebelow.cn)

plot(prebelow.cn)
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(prebelow.cn, add.smooth = FALSE, which = 1)
E <- resid(prebelow.cn)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(prebelow.cn, pairwise ~ traj.20,
                      adjust="tukey")


#pre tot CN
pretot.cn.new <- lme(pretot.cn.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(pretot.cn.new)
plot(pretot.cn.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(pretot.cn.new, add.smooth = FALSE, which = 1)
E <- resid(pretot.cn.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(pretot.cn.new, pairwise ~ traj.20,
                      adjust="tukey")


##################POST CN

#post above CN
postabove.cn <- lme(postabove.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postabove.cn)
plot(postabove.cn)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postabove.cn, add.smooth = FALSE, which = 1)
E <- resid(pretot.c.new)
hist(E, xlab = "postabove.cn", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(postabove.cn,
                      pairwise ~ traj.20,
                      adjust="tukey")


## post below CN
postbelow.cn <- lme(postbelow.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(postbelow.cn)
plot(postbelow.cn)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(postbelow.cn, add.smooth = FALSE, which = 1)
E <- resid(postbelow.cn)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(postbelow.cn,pairwise ~ traj.20,
                      adjust="tukey")


## post tot CN

posttot.cn.new <- lme(posttot.cn.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(posttot.cn.new)
plot(posttot.cn.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(posttot.cn.new, add.smooth = FALSE, which = 1)
E <- resid(posttot.cn.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

 emmeans(posttot.cn.new,
                      pairwise ~ traj.20,
                      adjust="tukey")


################## LOST CN
 
 #lost above CN
lostabove.cn <- lme(lostabove.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostabove.cn)
plot(lostabove.cn)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostabove.cn, add.smooth = FALSE, which = 1)
E <- resid(lostabove.cn)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(lostabove.cn,
                      pairwise ~ traj.20,
                      adjust="tukey")


# lost below CN
lostbelow.cn <- lme(lostbelow.cn ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(lostbelow.cn)
plot(lostbelow.cn)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lostbelow.cn, add.smooth = FALSE, which = 1)
E <- resid(lostbelow.cn)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(lostbelow.cn, pairwise ~ traj.20,
                      adjust="tukey")


## lost total CN
losttot.cn.new <- lme(losttot.cn.new ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(losttot.cn.new)
plot(losttot.cn.new)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(losttot.cn.new, add.smooth = FALSE, which = 1)
E <- resid(losttot.cn.new)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)

emmeans(losttot.cn.new, pairwise ~ traj.20,
                      adjust="tukey")



##################### post SOL depth and % SOL combusted
#postsol.depth
names(data)
data=mutate(data, perc.SOL.C.lost= (lostsol.c/presol.c)*100)
data$perc.SOL.C.lost

datasum=data%>%
  group_by(traj.20)%>%
  summarize(meanburndepth=mean(lostsol.depth), 
            meanpercentlost=mean(perc.SOL.C.lost), 
            err.burn=std.error(lostsol.depth), 
            err.perc=std.error(perc.SOL.C.lost))


m1 <- lme(lostsol.depth ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(m1)
plot(m1)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(m1, add.smooth = FALSE, which = 1)
E <- resid(m1)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


library(lsmeans)
emmeans(m1,pairwise ~ traj.20,
                      adjust="tukey")




m2 <- lme(perc.SOL.C.lost ~traj.20 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(m2)
plot(m2)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(m2, add.smooth = FALSE, which = 1)
E <- resid(m2)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


library(lsmeans)
leastsquare = emmeans(m2,
                      pairwise ~ traj.20,
                      adjust="tukey")
CLD = CLD(leastsquare, alpha=0.05, Letters=letters)
CLD


m2 <- lme(postsol.depth ~traj.new17 , random=~1|fire.x,method="REML", data=data, na.action=na.omit)
summary(m2)
plot(m2)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(m2, add.smooth = FALSE, which = 1)
E <- resid(m2)
hist(E, xlab = "Residuals", main = "")
plot(data$fire.x, E, 
     ylab = "Residuals")
plot(data$traj.20, E, 
     ylab = "Residuals")
par(op)


library(lsmeans)
leastsquare = emmeans(m2,
                      pairwise ~ traj.new17,
                      adjust="tukey")
CLD = CLD(leastsquare, alpha=0.05, Letters=letters)
CLD


######
#figure 1
p2 <- ggplot() 
p2= p2+theme_bw()+ylab(expression('Post-fire SOL Carbon (g C m'^-2*')'))+ xlab(expression('Pre-fire SOL Carbon (g C m'^-2*')'))
#p2= p2+geom_point(data=data,aes(x=presol.c,y=postsol.c, fill=traj.20, shape=traj.20), size=1.5, color="black")
p2=p2+scale_shape_manual( name="Regeneration Trajectory", values=c(21, 22, 24), labels=c(expression('Spruce (R'^2*'= 0.67)'), expression('Mixed (R'^2*'= 0.84)'),  expression('Deciduous (R'^2*'= 0.03)')))
p2=p2+ scale_fill_manual( name="Regeneration Trajectory", values=c("#443983FF" ,"#77D153FF","gold1"), labels=c(expression('Spruce (R'^2*'= 0.67)'), expression('Mixed (R'^2*'= 0.84)'),  expression('Deciduous (R'^2*'= 0.03)')))
p2=p2+  theme(legend.position=c(0.02, 0.98), legend.spacing.y = unit(0.1, 'cm'), legend.spacing.x = unit(0.1, 'cm')) + theme(legend.text=element_text(size=6)) + theme(legend.background = element_rect(fill="transparent"))                                                                                                                                                                              
p2= p2 +  theme(legend.text.align = 0) +scale_x_continuous(breaks=number_ticks(5)) + scale_y_continuous(breaks=number_ticks(5))
p2= p2+ theme(legend.title=element_text(size=7)) + geom_abline(intercept = 0, slope = 1, size=0.18)
p2=p2 + geom_line(data=int, aes(color=int$traj.20, x=presol.c, y=new.fit), size=0.18)
p2 = p2 + scale_color_manual(values=c("#443983FF" ,"#77D153FF","gold1"), guide=FALSE)
p2 = p2+geom_ribbon(data=int, aes(fill=int$traj.20, x=presol.c, ymin=new.lower,ymax=new.upper), alpha=0.2) 
p2= p2+geom_point(data=data,aes(x=presol.c,y=postsol.c, fill=traj.20, shape=traj.20), stroke=0.1,size=2, color="black")
p2 = p2+theme(panel.grid.minor = element_line(size = 0.18), panel.grid.major = element_line(size = 0.18))
p2= p2+theme(axis.title.x = element_text( size=8, vjust=1.5), axis.text.x  = element_text(size=7)) 
p2= p2+theme(axis.title.y = element_text( size=8, vjust=1.5), axis.text.y  = element_text(size=7))
p2= p2+theme(axis.title.x = element_text( size=8, vjust=1.5), axis.text.x  = element_text(size=7)) 
p2= p2+theme(axis.title.y = element_text( size=8, vjust=1.5), axis.text.y  = element_text(size=7))
p2 = p2 + theme(legend.margin=margin(0,0.2,0.2,0.2), legend.justification = c("left", "top"),legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.1, 'cm'))
p2 = p2 + theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18))
p2 = p2 + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.18))
p2 = p2 + theme(plot.margin = unit(c(0.1,0.15,-0.1,0), "cm"))
p2= p2 + theme(axis.text.x = element_text(vjust = 1), axis.title.x = element_text(vjust = 2))  
p2= p2 + theme(axis.text.y = element_text(hjust = 1), axis.title.y = element_text(vjust = 0)) 
p2 = p2 +theme(legend.box.margin=margin(-1,-1,-5,-1))
p2
ggsave("Figure_2.pdf", plot = p2, device = NULL, path = NULL,
       scale = 1, width =5.5, height = 6, units = "cm",
       dpi = 600)

############
##### and a barplot figure 1
names(data)
dat2<-dplyr::select(data, site, traj.20, prebelow.c, preabove.c, postbelow.c, postabove.c, lostbelow.c, lostabove.c)

library(reshape)
mdata <- melt(dat2, id=c("site","traj.20"))
head(mdata)
names(mdata)


levels(mdata$variable)
mdata2<-mutate(mdata, class=ifelse(variable=="prebelow.c"|variable=="preabove.c", "Pre-fire C",
                                   ifelse(variable=="postbelow.c"|variable=="postabove.c", "Post-fire C",
                                          ifelse(variable=="lostbelow.c"|variable=="lostabove.c", "C Loss", NA))))

mdata2<-mutate(mdata2, above.below=ifelse(variable=="prebelow.c"|variable=="postbelow.c"|variable=="lostbelow.c", "Belowground",
                                          ifelse(variable=="preabove.c"| variable=="postabove.c"|variable=="lostabove.c", "Aboveground", NA)))




names(mdata2)
head(mdata2)
mdata2$class<-as.factor(mdata2$class)
mdata2$above.below<-as.factor(mdata2$above.below)


mdata2$class<-factor(mdata2$class, levels=c("Pre-fire C", "Post-fire C", "C Loss"))

mdata2$traj.20<-factor(mdata2$traj.20, levels=c("Spruce", "Mixed", "Deciduous"))

library(Rmisc)
tgc <- summarySE(mdata2, measurevar="value", groupvars=c("traj.20", "variable", "class","above.below"))
head(tgc)


head(tgc)
head(tgc)
tgc$class<-factor(tgc$class, levels=c("Pre-fire C", "Post-fire C", "C Loss"))

tgc$traj.20<-factor(tgc$traj.20, levels=c("Spruce", "Mixed", "Deciduous"))

data$traj.20<-factor(data$traj.20, levels=c( "Spruce", "Mixed" ,"Deciduous"))

mus3 <- within(tgc,var2 <- ave(value,traj.20,class, FUN=cumsum))
head(mus3)

#write.csv(mus3, "err1.csv")
########## write it - then run the models below to add in teh letters ##
### then read it in and add in the ltext

tj=read.csv("err1.csv")
head(tj)
mus3$class<-factor(mus3$class, levels=c("Pre-fire C", "Post-fire C", "C Loss"))

mus3$traj.20<-factor(mus3$traj.20, levels=c("Spruce", "Mixed", "Deciduous"))
names(mus3)
mus3$variable<-factor(mus3$variable, levels=c( "preabove.c", "prebelow.c", "lostabove.c", "lostbelow.c", "postabove.c", "postbelow.c"))
head(mus3)

tgc$variable<-factor(tgc$variable, levels=c( "preabove.c", "prebelow.c", "lostabove.c", "lostbelow.c", "postabove.c", "postbelow.c"))
head(tgc)
names(mdata2)
head(mus3)

#write.csv(tgc, file="traj.error.bars.csv")
head(tj)
tj=dplyr::select(tj, traj.20,variable,class, above.below, N, text)
head(tj)
mus4=merge(mus3, tj, by=c("traj.20", "variable", "class", "above.below", "N"))
head(mus4)

head(mus3)
g=ggplot(data=tgc, aes(x=traj.20, y=value, fill=traj.20)) + scale_alpha_discrete(range = c(0.3, 1))+
  geom_bar(stat="identity", size=0.18, aes(alpha=above.below), color="black") + facet_wrap(~class) +
  scale_fill_manual(values=c("#443983FF" ,"#77D153FF","gold1"), guide=F)+
  geom_errorbar(data=mus3,aes(ymin=var2-se, ymax=var2+se),width=0.1, size=0.18) + theme_bw() +
  scale_x_discrete(labels=c("spruce" = "Spruce", "mixed" = "Mixed", "decid" = "Deciduous"))+ 
  theme(legend.title=element_blank()) + theme(legend.position="top", legend.spacing.x = unit(0.1, 'cm'))+ 
  theme(axis.title.x = element_text( size=8, vjust=1.5), axis.text.x  = element_text(size=7)) +
  theme(axis.title.y = element_text( size=8, vjust=1.5), axis.text.y  = element_text(size=7))+
  ylab(expression('Carbon (g C m'^-2*')'))+ xlab("Regeneration Trajectory")+
  theme(legend.background=element_blank()) + 
  scale_y_continuous(breaks=number_ticks(8), limits=c(0,9100)) + theme(strip.text.x = element_text(size = 8)) +theme(strip.background = element_rect(fill="white", size=0.18))+
  theme(legend.text=element_text(size=7)) +
  geom_text(data=subset(mus4,above.below=="Belowground"), aes(x=traj.20, y=var2-2*se,label=text),size=2.5)+ 
  geom_text(data=subset(mus4,above.below=="Aboveground"&class=="Pre-fire C"), aes(x=traj.20, y=var2+4*se,label=text), size=2.5)+
  geom_text(data=subset(mus4,above.below=="Aboveground"&class=="C Loss"), aes(x=traj.20, y=var2+6*se,label=text), size=2.5)+
  geom_text(data=subset(mus4,above.below=="Aboveground"&class=="Post-fire C"), aes(x=traj.20, y=var2+4*se,label=text), size=2.5)+
  theme(axis.title.x = element_blank())+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.18))+
  theme(legend.margin=margin(0.15,0,-0.15,0), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.0, 'cm'))+
  theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.18))+
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm"))) + theme(plot.margin = unit(c(0.15,0.15,0.05,0), "cm"))+
  theme(axis.text.y = element_text(hjust = 1), axis.title.y = element_text(vjust = 0)) +
  theme(legend.box.margin=margin(0,-1,-5,-1))
g


ggsave("Figure_1.pdf", plot = g, device = NULL, path = NULL,
scale = 1, width =12, height = 6.5, units = "cm",
dpi = 600)


