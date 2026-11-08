# Example 1
Before=c(268,225,252,192,307,228,246,298,231,185)
After=c(106,186,223,110,203,101,211,176,194,203)
Subject=1:10
data1=data.frame(Before,After,Subject)
str(data1)
data1$Subject=as.factor(data1$Subject)
head(data1)

par(mfrow=c(2,3))
hist(Before);hist(After);hist(After-Before)
boxplot(Before);boxplot(After);boxplot(After-Before)

model1<-lm(After~Before,data1)
head(model1)
model.matrix(model1) #Compute model designed matrix X.
summary(model1)

par(mfrow=c(2,2))
plot(model1)

data2<-read.csv("/Users/tang/Documents/before_after.csv")
str(data2)
model1<-lm(Trt~Msr,data=data2)
model.matrix(model1)
summary(model1)
plot(model1)
data2$Subj<-as.factor(data2$Subj)
data2$Subj<-relevel(data2$Subj,"subj1")
head(data2$Subj);is.factor(data2$Subj)
model2<-lm(Trt ~ ., data=data2) #Model Trt with all other variables.
par(mfrow=c(1,1))
plot(model2,which=2)

# Example 1B
var1<-c(12,13,12,42,14,12,32,44,12,4,21,44)
length(var1)
var2<-c("yes","no","yes","yes","yes","no","yes","yes","yes","no","yes","yes")
length(var2)
set<-data.frame(var1,var2)
head(set)
set$var2<-as.factor(set$var2) #Making R to treat its values as categories
set$var2<-relevel(set$var2,"no") #chooses one factor category as the reference level.
model1<-lm(var1~var2,data=set)
summary(model1)
model.matrix(model1)
par(mfrow=c(2,2))
plot(model1)

# Example 2
library("GLMsData")
data("turbines")
head(turbines)
names(turbines)
is.numeric(turbines$Hours)
is.numeric(turbines$Turbines)
is.numeric(turbines$Fissures)
#All variables are quantitative.
summary(turbines$Hours)
summary(turbines$Turbines)
summary(turbines$Fissures)
prop<-turbines$Fissures/turbines$Hours
head(prop)
plot(turbines$Hours,prop,
xlab="Proportion of failure",ylab="Number of hours")
model1<- lm(turbines$Hours~prop)
summary(model1)

# Example 3
library("GLMsData")
data("humanfat")
head(humanfat)
names(humanfat)
humanfat$Gender<-as.factor(humanfat$Gender)
humanfat$Gender<-relevel(humanfat$Gender,"F")
head(humanfat$Gender)
summary(humanfat$Age)
summary(humanfat$Percent.Fat)
summary(humanfat$Gender)
summary(humanfat$BMI)

par(mfrow=c(2,3))
plot(humanfat$Age,humanfat$BMI,xlab="Age (Years)",ylab="BMI")
plot(humanfat$Age,humanfat$Percent.Fat,xlab="Age (Years)",ylab="Body fat (%)")
plot(humanfat$Percent.Fat,humanfat$BMI,xlab="Body fat (%)",ylab="BMI")
plot(humanfat$Gender,humanfat$BMI,xlab="Gender",ylab="BMI")
plot(humanfat$Gender,humanfat$Age,xlab="Gender",ylab="Age (Years)")
plot(humanfat$Gender,humanfat$Percent.Fat,xlab="Gender",ylab="Body fat (%)")

model.matrix(~Gender-1,data=humanfat)
model1<-lm(BMI~Gender,data=humanfat)
summary(model1)$adj.r.squared
model2<-lm(BMI~Age,data=humanfat)
summary(model2)$adj.r.squared
model3<-lm(BMI~Percent.Fat,data=humanfat)
summary(model3)$adj.r.squared
model1<-lm(BMI~Percent.Fat+Age,data=humanfat)
summary(model1)$adj.r.squared
model1<-lm(BMI~Percent.Fat+Gender,data=humanfat)
summary(model1)$adj.r.squared
model1<-lm(BMI~Age+Gender,data=humanfat)
summary(model1)$adj.r.squared
model1<-lm(BMI~Age+Gender+Percent.Fat,data=humanfat)
summary(model1)$adj.r.squared

# Example 4
y<-c(0,1)
mu0<-c(0.1,1)
z<-log(mu0)+(y-mu0)*(1/mu0)
w<-matrix(0,nrow=2,ncol=2)
w[1,1]<-mu0[1]
w[2,2]<-mu0[2]

irwis<-function(mu0,x,y,maxn,epsilon){
    beta<-0
    for (i in 1:maxn){
        z<-log(mu0)+(y-mu0)*(1/mu0)
        w<-matrix(0,nrow=length(x),ncol=length(x))
        for (i in 1:length(x)){
            w[i,i]<-mu0[i]
        }
        beta1<-solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%z
        if (abs(beta-beta1)<epsilon){
            beta<-beta1
            break
        }
        else {
           beta<-beta1
           mu0<-exp(x%*%beta1)
        }
    }
    return(beta)
}
x0<-matrix(c(1,2),nrow=2)
y0<-matrix(c(0,1),nrow=2)
m0<-matrix(c(0.5,0.75),nrow=2)
irwis(m0,x0,y0,100,0.001)

summary(glm(y0~x0-1,family=poisson))

x0<-matrix(c(1,2,3,4),nrow=4)
y0<-matrix(c(0,1,3,2),nrow=4)
m0<-matrix(c(0.5,0.75,0.5,0.75),nrow=4)
irwis(m0,x0,y0,100,0.001)

summary(glm(y0~x0-1,family=poisson))

# Example 4B
irwis<-function(mu0,x,y,n,epsilon){
    beta<-0
    for (i in 1:n){
        w<-matrix(0,nrow=length(x),ncol=length(x))
        z<-log(mu0/(1-mu0))+(y-mu0)/(mu0*(1-mu0))
        for (k in 1:length(x)){
            w[k,k]<-mu0[k]*(1-mu0[k])
        }
        beta1<-solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%z
        if (abs(beta1-beta)<epsilon){
            beta<-beta1
            break
        }
        else{
            beta<-beta1
            mu0<-exp(x%*%beta)/(1+exp(x%*%beta))
        }
    }
    return(beta)
}

x0<-matrix(c(1,2),nrow=2)
y0<-matrix(c(0,1),nrow=2)
m0<-matrix(c(0.5,0.75),nrow=2)
irwis(m0,x0,y0,100,0.001)

summary(glm(y0~x0-1,family=binomial(link="logit")))

x0<-matrix(c(1,2,3,4),nrow=4)
y0<-matrix(c(0,1,1,1),nrow=4)
m0<-matrix(c(0.5,0.75,0.5,0.75),nrow=4)
irwis(m0,x0,y0,100,0.001)

summary(glm(y0~x0-1,family=binomial(link="logit")))

# Example 5
data1<-read.csv("/Users/tang/Documents/fish_length(1).csv")
head(data1)
data1$River<-as.factor(data1$River)
data1$River<-relevel(data1$River,ref="Wear")

model1<-glm(Length~River,data=data1,family=gaussian)
model2<-glm(Length~River,data=data1,family=Gamma)
model3<-glm(Length~River,data=data1,family=poisson)

summary(model1)
summary(model2)
summary(model3)

# Example 6
generator<-function(n){
    y<-matrix(0,nrow=n,ncol=1)
    x<-runif(n,0,100)
    temp<-exp(-2+0.04*x)
    pi1<-temp/(1+temp)
    for (i in 1:n){
    y[i]<-rbinom(1,1,pi1)
    }
    model<-fitted(glm(y~x,family=binomial))
    return(cor(y-model,model))
}
set.seed(10)
generator(25);generator(100);generator(1000);generator(10000)

# Part b
pchisq(c(0.602,14.83,2.83),df=1,lower.tail=FALSE)

beta<-c(0.143,1.247,-0.706)
se<-c(0.19,0.45,0.45)
wald<-beta/se
pnorm(-abs(wald),0,1)*2

alpha<-0.025
zvalue<-qnorm(1-alpha,0,1)
lower<-beta-zvalue*se
upper<-beta+zvalue*se
cbind(lower,upper)

# Part c
temp<-rnorm(100,1,21)*1.1
temp2<-rbinom(100,1,0.67)
x0<-matrix(temp,nrow=100)
x1<-matrix(temp2,nrow=100)
y<-round(runif(100,0,1),0)
y0<-matrix(y,nrow=100)
model1<-glm(y0~x0+x1,family=binomial)
names(model1)
temp<-summary(model1)
temp<-temp$coefficients
temp
estimate<-temp[,1]
se<-temp[,2]
estimate
se
wald<-estimate/se
wald
lower<-estimate-qnorm(0.975,0,1)*se
upper<-estimate+qnorm(0.975,0,1)*se
cbind(lower,upper)

library(MASS)
confint(model1)

wald<-function(design,y){
    model1<-glm(y~design,family=binomial)
    temp<-summary(model1)$coefficients
    estimate<-temp[,1]
    sd<-temp[,2]
    wald<-estimate/sd
    z<-qnorm(0.975,0,1)
    lower<-estimate-z*sd
    upper<-estimate+z*sd
    print(cbind(lower,upper))
    return(confint(model1))
}
temp<-rnorm(100,1,21)*1.1
temp2<-rbinom(100,1,0.67)
y<-round(runif(100,0,1),0)
y0<-matrix(y,nrow=100)
x<-as.matrix(cbind(temp,temp2))
wald(x,y0)


# Part d
simulation<-function(n,r,beta0,beta1){
    set.seed(10)
    temp1<-matrix(0,nrow=n)
    temp15<-matrix(0,nrow=n)
    temp2<-matrix(0,nrow=n)
    temp3<-matrix(0,nrow=n)
    for (i in 1:n){
        x<-rbinom(r,1,0.37)
        eta<-beta0+beta1*x
        mu0=exp(eta)
        y<-rpois(r,lambda=mu0)
        model1<-summary(glm(y~x,family=poisson))
        temp1[i]<-model1$coefficients[1,1]
        temp15[i]<-model1$coefficients[1,2]
        temp2[i]<-model1$coefficients[2,1]
        temp3[i]<-model1$coefficients[2,2]
    }
    z<-(temp2-beta1)/temp3
    z1<-(temp1-beta0)/temp15
    temp4<-cbind(temp1,temp15,z1,temp2,temp3,z)
    colnames(temp4)<-c("estimated beta0","s.d. of beta0","z for beta0",
    "estimated beta1","s.d. of beta1","z for beta1")
    par(mfrow=c(2,2))
    hist(temp1,breaks=50,main="Histogram of beta0",xlab="beta0")
    hist(temp2,breaks=50,main="Histogram of beta1",xlab="beta1")
    qqnorm(temp1,main="Normal QQ plot for beta0")
    qqline(temp1)
    qqnorm(temp2,main="Normal QQ plot for beta1")
    qqline(temp2)
    return(head(temp4))
}
simulation(1000,100,0.1,0.4)

# Part dii
set.seed(1)
simulation<-function(r,lambda,beta0,beta1){
    temp1<-numeric(1000)
    temp2<-numeric(1000)
    temp3<-numeric(1000)
    for (i in 1:1000){
        x<-rpois(r,lambda)
        mu0<-1/(1+exp(-(beta0+beta1*x)))
        y<-rbinom(r,1,mu0)
        model1<-summary(glm(y~x,family=binomial))
        temp1[i]<-model1$coefficients[1,2]
        temp2[i]<-model1$coefficients[2,2]
    }
    par(mfrow=c(1,2))
    hist(temp1,breaks=30,xlab="beta1",main="Histogram for beta1")
    qqnorm(temp1,main="Normal QQ plot for beta1")
    qqline(temp1)
}
simulation(100,1,0,1)

# Part e
set.seed(1)
simulation<-function(x,beta0,beta1,shape,alpha,x0){
    n<-length(x)
    mu=exp(beta0+beta1*x)
    y<-rgamma(n,shape=shape,rate=shape/mu)
    temp<-glm(y~x,family=Gamma(link="log"))
    model<-summary(temp)
    constant<-model$coefficients[1,1]
    estimate<-model$coefficients[2,1]
    sd<-model$coefficients[2,2]
    test<-qnorm((1-alpha/2),0,1)
    lower<-estimate-sd*test
    upper<-estimate+sd*test
    ci<-cbind(lower,upper)
    colnames(ci)<-c("Lower bound","Upper bound")
    mean1<-predict(temp,newdata=data.frame(x=x0),type="response")
    print(paste("Mean response at a given covariate value:", mean1))
    return(ci)
}
x<-rnorm(1000,0,1)
simulation(x,1,2,1,0.05,2)

# Question 7
beta<-c(-1.2,0.45,0.8)
se<-c(0.25,0.12,0.35)
z<-beta/se
2*pnorm(-abs(z),0,1)
dfull<-112.4;dreduced<-118.9
testst<- dreduced-dfull
pchisq(testst,df=1,lower.tail=FALSE)

# Part b
da<-210.5;db<-175.2;dc<-162.8
dfa<-99;dfb<-98;dfc<-97
test1<-da-db
pchisq(test1,df=dfa-dfb,lower.tail=FALSE)
test1<-db-dc
pchisq(test1,df=dfb-dfc,lower.tail=FALSE)

# Part c
est<-0.3;se<-0.07
test<-qnorm(0.95,0,1)
est-se*test; est+se*test

beta<-c(0.23,0.04,0.06,0.01,0.09,0.05,0.3)
se<-c(0.13,0.04,0.04,0.03,0.06,0.02,0.07)
beta/se
2*(1-pnorm(abs(beta/se),0,1))

# Part d
dr<-245.3; df=230.1
dfr<-98;dff<-97
diff<-dr-df
pchisq(diff,df=dfr-dff,lower.tail=FALSE)
z<-1.9
2*(1-pnorm(abs(z),mean=0,sd=1))