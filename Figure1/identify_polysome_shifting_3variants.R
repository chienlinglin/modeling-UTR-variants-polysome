
library(MASS)
library(qvalue)

utr3<-read.csv("20190502_3UTR.csv", row.names=1)
gtable<-read.csv("20181220_grouptable.csv")

myplot<-function(wtc, mtc){
	wtc<-log(wtc+1)
	mtc<-log(mtc+1)
	plot(wtc[1:3], x=1:3, pch=15, ylim=c(min(wtc, mtc), max(wtc, mtc)), main=as.character(gtable[candidate[k],1]) )
	points(wtc[4:6], x=1:3, pch=16)
	points(wtc[7:9], x=1:3, pch=17)
	points(wtc[1:3], x=1:3, type="l")
	points(wtc[4:6], x=1:3, type="l")
	points(wtc[7:9], x=1:3, type="l")
	
	points(mtc[1:3], x=1:3, pch=15, col=2)
	points(mtc[4:6], x=1:3, pch=16, col=2)
	points(mtc[7:9], x=1:3, pch=17, col=2)
	points(mtc[1:3], x=1:3, type="l", col=2)
	points(mtc[4:6], x=1:3, type="l", col=2)
	points(mtc[7:9], x=1:3, type="l", col=2)
}

################################################################################

frac<-rep(1:3, 6)
expt<-rep(rep(1:3, each=3), 2)
mut<-rep(0:1, each=9)
pvals<-NULL
chisqs1<-NULL; chisqs2<-NULL; chisqs3<-NULL
for (i in 1:3168){
	wtname<-as.character(gtable[i,2])
	mtname<-as.character(gtable[i,3])
	utr3i<-utr3[c(wtname, mtname), c(2:4)]
	for (j in 2:3) utr3i<-cbind(utr3i, utr3[c(wtname, mtname), c((j-1)*4+2:4)])
	count<-as.numeric(as.matrix(t(utr3i)))

	if (sum(apply(matrix(count, nr=3), 2, sum)<20)>0){
		chisqs1<-c(chisqs1, NA)
		chisqs2<-c(chisqs2, NA)
		chisqs3<-c(chisqs3, NA)
	} else{ 
	chisq1<-NA; chisq2<-NA; chisq3<-NA
	try({
	fit<-glm.nb(count~as.factor(frac)+as.factor(expt)+mut*as.factor(expt)+mut+mut*as.factor(frac))
	chisq1<-t(fit$coef[9:10])%*%solve(summary(fit)$cov.scaled[9:10, 9:10])%*%fit$coef[9:10]
	chisq2<-t(fit$coef[4:5])%*%solve(summary(fit)$cov.scaled[4:5, 4:5])%*%fit$coef[4:5]
	chisq3<-t(fit$coef[7:8])%*%solve(summary(fit)$cov.scaled[7:8, 7:8])%*%fit$coef[7:8]
	}, silent=TRUE)
	chisqs1<-c(chisqs1, chisq1)
	chisqs2<-c(chisqs2, chisq2)
	chisqs3<-c(chisqs3, chisq3)
	}
	print(i); flush.console()
}

temp1<-cbind(chisqs1, chisqs2, chisqs3)
chisqs1.filter<-chisqs1[chisqs2<5.99 & chisqs3<5.99] 
index.filter<-(1:3168)[chisqs2<5.99 & chisqs3<5.99]
index.filter.temp3<-(1:3168)[as.numeric(row.names(temp3))]
table(pchisq(chisqs1.filter, df=2, lower.tail=FALSE)<0.05/sum(!is.na(chisqs1.filter)))

hist(pchisq(chisqs1.filter, df=2, low=FALSE))

candidate<-index.filter[order(chisqs1.filter, decreasing=TRUE)[1:27]]

gtable[candidate,]
filter3<-na.omit(gtable[index.filter,])

par(mfrow=c(5, 6))
for (k in c(1:27)){

wtname.top<-as.character(gtable[candidate[k],2])
mtname.top<-as.character(gtable[candidate[k],3])

wtc<-as.numeric(as.matrix(utr3[wtname.top, c(2:4, 6:8, 10:12)]))
mtc<-as.numeric(as.matrix(utr3[mtname.top, c(2:4, 6:8, 10:12)]))

myplot(wtc, mtc)
}
