##############################################################################
#GLMM with BNP prior on the random effects  
#############  GLMM -  JAGS -  240 patients - 17 hospitals        #################
##############################################################################
setwd("C:/LocalDiskAlessandra/DocuLavoro/Didattica/StatisticaBayesiana/mat2425/SLIDES_2425/R-Notebooks_2425")
rm(list=ls())

# Continuous covariates are centered
# Import the dateset
input = read.table("db_centrato_xDP.txt",header=T)
names(input) 

head(input)
summary(input)

n.data=dim(input)[1] # total num of patients
J=17 # num of hospitals
M=40 #  truncation level of the stick-breaking sum
## eta and logOB have been standardized

#######################################
library(rjags)      #   R to JAGS
library(plotrix)    #  
set.seed(1)         #  set seed to reproduce output


# list of all data needed in JAGS 
data <- list(Y=input$vivo,
		AGE=input$eta, KILLIP=input$killip,LOGOB = input$logOB,
		CENTRO=input$centro,
		n=n.data,J=J,M=M)

###########  
#  theta's are the location points of the DP 
# (of the truncated Sethuraman's series)

r=rep(0.5,M)
theta=rep(0,M)
S=rep(1,J)

inits = list(beta=rep(0,3),a = 1,  
            #lambda.bb = 1,
		r=r,
		theta= theta,
		S = S,
		Snew = 1,
		.RNG.seed = 2,
		.RNG.name = 'base::Wichmann-Hill'
)

#########################
### CREATE JAGS MODEL ###
#########################

modelGLMM_DP=jags.model("GLMM_DP_model.bug",data=data,inits=inits,n.adapt=1000,n.chains=1) 
# I run 19K more iterations 
update(modelGLMM_DP,19000)


variable.names=c("bb", "beta",  "tau.bb", "newcentro", "alpha","K")
# Parameter to save for posterior inference 
n.iter=50000 
thin=10  

outGLMM_DP=coda.samples(model=modelGLMM_DP,variable.names=variable.names,n.iter=n.iter,thin=thin)
#  save the MCMC of the parameters of interest


save(outGLMM_DP,file='GLMMDP_output_23.Rdata')

##############################
### OUTPUT SAVED in WeBeep ###
##############################
setwd("C:/LocalDiskAlessandra/DocuLavoro/Didattica/StatisticaBayesiana/mat2324/SLIDES_NOTES/Rnotebook")
rm(list=ls())

# Import data  
input = read.table("db_centrato_xDP.txt",header=T)
names(input) 

n.data=dim(input)[1] # total number of patients 
J=17 # number of hospitals
M=40 # upper value of truncated stick-braeking sum

#install.packages("plotrix")

library(coda)        # pacchetto per analizzare catene
library(plotrix)     # per fare plot CIs

load('GLMMDP_output_23.Rdata') 
  # upload output  from  coda.samples

#  summary(output)

data=as.matrix(outGLMM_DP) #  dataframe into a matrix 
data=data.frame(data)
attach(data)
n.chain=dim(data)[1]   # lunghezza catena (final sample size)

##################
### TRACEPLOTS ###
##################
x11()
par(mfrow=c(1,3))
plot(data[,'beta.1.'],main='età/age',type='l')
plot(data[,'beta.2.'],main='log(OB)',type='l')
plot(data[,'beta.3.'],main='killip',type='l')

x11()
par(mfrow=c(1,3))
plot(density(data[,'beta.1.']),main='age')
plot(density(data[,'beta.2.']),main='log(OB)')
plot(density(data[,'beta.3.']),main='killip')
# age and killip   significative;  log (OB) less significative  
mean(data[,'beta.1.']<0)
mean(data[,'beta.3.']<0)
mean(data[,'beta.2.']<0)

x11()
par(mfrow=c(3,3)) #  traceplot  of some random effects parameters
plot(data[,'bb.1.'],main='b1',type='l')
plot(data[,'bb.2.'],main='b2',type='l')
plot(data[,'bb.3.'],main='b3',type='l')
plot(data[,'bb.7.'],main='b7',type='l')
plot(data[,'bb.10.'],main='b10',type='l')
plot(data[,'bb.13.'],main='b13',type='l')
plot(data[,'bb.15.'],main='b15',type='l')
plot(data[,'bb.16.'],main='b16',type='l')
plot(data[,'bb.17.'],main='b17',type='l')

x11()
par(mfrow=c(1,2))
plot(data[,'alpha'],main='total mass',type='l') # total mass parameter of the DP
plot(data[,'K'],main='K_n',type='l')

x11()
par(mfrow=c(1,2))
plot(density(data[,'alpha']),main='total mass')
plot(table(data[,'K']),main='K_n')
mean(alpha);var(alpha)
mean(data[,'K']);var(data[,'K'])

x11()
par(mfrow=c(1,1))
plot(data[,'tau.bb'],main='tau',type='l')
#plot(1/sqrt(data[,'tau.bb']),main='lambda',type='l')
#plot(density(1/sqrt(data[,'tau.bb'])),main='lambda')


######################################
##  Plot IC for the random intecepts##  
######################################

### ordered by number of patients ###
# creo un vettore contenente il numero di pazienti per ogni ospedale
patients=rep(0,J)
for(i in 1:J){
	patients[i]=length(which(input$centro==i))
}
sort_patients=sort(patients,index.return=T)   # ordino


# ricalcolo i quantili per ogni ospedale in ordine di numero di pazienti
Q=matrix(nrow=J+1, ncol=3) 
for (j in 1:J){
	Q[j,]=quantile(data[, 2 + sort_patients$ix[j]  ],probs=c(0.025,0.5,0.975))
}
Q[J+1,]=quantile(data$newcentro,probs=c(0.025,0.5,0.975))
colnames(Q) <- c("2.5","median","97.5")


# symbol for point estimate: all points are round, 
#    the last point is x (new random hospital)
pch=c(rep(21,J),4)   
x11()
plotCI(x=seq(1,J+1),y=Q[,2],uiw=(Q[,3]-Q[,2]) ,liw=(Q[,2]-Q[,1]),pch=pch,
       scol=1 , xlab="hospitals (sort by increasing number of patients)", 
       ylab="b_j", main="CIs of the Random Intercept")  
abline(h=mean(Q[,2]))
##################################
# Cluster estimate di Lau&Green
################################
#load('M1.Rdata') # load coda output 
#data=as.matrix(M1) 
#data=data.frame(data)  # convert into a dataframe
J=17

label.mat = as.matrix(data[,3:19]) # extract cluster labels

m=J 
G=n.chain
pihat <- matrix(0,ncol=m,nrow=m)
for(i in 1:G){
  ss <- label.mat[i,]
  cij <- outer(ss,ss,'==')
  pihat <- pihat+cij
}

pihat <- pihat/G

#####Binder loss function
FF <- vector("numeric")
K <- 0.7 #prova con K=0.5, >=0.7
for(i in 1:G){
  ss <- label.mat[i,] 
  cij <- outer(ss,ss,'==')
  pluto <- (pihat-K)*as.matrix(cij)
  pluto <-  pluto[upper.tri(pluto)]
  FF[i] <- sum(pluto)
}

plot(FF)

ind.bind <- which.max(FF)[1]
label.mat[ind.bind,]#leti(ind.bind,L,G,m)
#plot(FF)
ll.bind <- label.mat[ind.bind,] #leti(ind.bind,L,G,m)
unici <- unique(ll.bind)
unici
l.uni <- length(unici)#  estimated number of clusters
l.uni

ncl=l.uni
for(i in 1:ncl){
  print(as.numeric(which(ll.bind==unici[i])))
}
##########################################################################
### PLot of the random intercepts with different color by estimated cluster
Sest=rep(0,J)
for(i in 1:ncl){
  Sest[as.numeric(which(ll.bind==unici[i]))] = i
}

Sest
table(Sest)

nclusLG=length(unique(Sest))
nclusLG

bins = as.numeric(names(table(Sest)))
freqs = as.vector(table(Sest))

J=17
label=rep(0,J)

mylist=list()

for(i in 1:nclusLG)
{
  in_gr=which(Sest==bins[i])
  label[in_gr]=i
  mylist[[i]]=in_gr
}

# "label" is the estimate of the vector of allocation variables

label
mylist
#This is the estimated random partition \rho

##################################################
# Figure of the posterior CI of the random effects
# each color denotes a cluster
# the grey CI's are singletons, i.e. cluster made 
#    by one observation only

gr1 <- mylist[[1]]
gr2 <- mylist[[2]]
gr3 <- mylist[[3]]
gr4 <- mylist[[4]]
gr5 <- mylist[[5]] 
gr6 <- mylist[[6]]
gr7 <- mylist[[7]]
gr8 <- mylist[[8]]
gr9 <- mylist[[9]] 
gr10 <- mylist[[10]] 
gr11 <- mylist[[11]] 




colore = rep('ciao',J)
colore[gr2]='grey' # grey: color of the singletons 
colore[gr5]='grey'
colore[gr6]='grey'
colore[gr8]='grey'
colore[gr10]='grey'
colore[gr11]='grey'
colore[gr1]='red'    
colore[gr3]='green'   
colore[gr4]='blue' 
colore[gr7]='yellow' 
colore[gr9]='magenta'

pch=rep(20,J)
# Nella matrice 'data' dalla colonna 3 in poi ci sono i b,
# gli effetti casuali (hospital-specific) al primo livello
# Per verificare il contenuto delle colonne di 'data': names(data) 
# ricalcolo i quantili per ogni ospedale in ordine di numero di pazienti
Q=matrix(nrow=J, ncol=3) #matrix(nrow=J, ncol=3) 
for (j in 1:J){
	Q[j,]=quantile(data[, 2 + j ],probs=c(0.025,0.5,0.975))
}

x11()
plotCI(x=seq(1,J),y=Q[,2],uiw=(Q[,3]-Q[,2]) ,liw=(Q[,2]-Q[,1]),pch=pch,,xaxt='n',
       col=colore,scol=colore,xlab="hospital",ylab=' ', main=" ", lwd=2)  
axis(side=1,at=seq(1,J),labels=seq(1,J),cex.axis=0.8)

#mean of all posterior medians
media <- mean(Q[,2])
abline(h= media)




