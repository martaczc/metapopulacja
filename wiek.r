library(fields)

data<-read.csv('wiek.txt',sep=' ')

time<-data$t
as<-unique(data$as)
ar<-unique(data$ar)
N<-data[[6]]
varN<-data[[7]]
age<-data[[8]]

NSpring<-matrix(N[time==0],nrow=101)
NAtumn<-matrix(N[time==3],nrow=101)
varNSpring<-matrix(varN[time==0],nrow=101)
varNAtumn<-matrix(varN[time==3],nrow=101)
ageSpring<-matrix(age[time==0],nrow=101)
ageAtumn<-matrix(age[time==3],nrow=101)


r<-rev(1:101)
ar<-ar[r]
NSpring<-NSpring[r,]
NAtumn<-NAtumn[r,]
ageSpring<-ageSpring[r,]
ageAtumn<-ageAtumn[r,]

#scale<-round((log(max(varN)-min(varN)))/log(10))-1
#start<-ceiling(min(varN)/10^scale)
#stop<-floor(max(varN)/10^scale)
#brk<-c(1:9*0.1,1:9,seq(start,stop,by=1)*10^scale)+0.001
#print(min(varN))
#print(max(varN))
#print(brk)
labels<-c(1,2,5)
brk<-c(labels*0.001,labels*0.01,labels*0.1,labels,labels*10)+0.001

#pdf(file="wiek.pdf")
par(mfrow=c(2,3),pty="s")
image.plot(ar,as,NSpring,main='średnie zagęszczenie populacji wiosną (t mod 5 = 0)',xlab='ar',ylab='as',zlim=c(min(N),max(N)))
image.plot(ar,as,log(varNSpring+0.001),main='wariancja zagęszczenia populacji wiosną (t mod 5 = 0)',xlab='ar',ylab='as',zlim=c(min(log(varN+0.001)),max(log(varN+0.001))),axis.args=list(at=log(brk), labels=brk-0.001))
image.plot(ar,as,ageSpring,main='średni wiek [miesiące] osobników w populacji wiosną (t mod 5 = 0)',xlab='ar',ylab='as',zlim=c(min(age),max(age)))

image.plot(ar,as,NAtumn,main='średnie zagęszczenie populacji jesienią (t mod 5 = 3)',xlab='ar',ylab='as',zlim=c(min(N),max(N)))
image.plot(ar,as,log(varNAtumn+0.001),main='wariancja zagęszczenia populacji jesienią (t mod 5 = 3)',xlab='ar',ylab='as',zlim=c(min(log(varN+0.001)),max(log(varN+0.001))),axis.args=list( at=log(brk), labels=brk-0.001))
image.plot(ar,as,ageAtumn,main='średni wiek [miesiące] osobników w populacji jesienią (t mod 5 = 3)',xlab='ar',ylab='as',zlim=c(min(age),max(age)))
#dev.off()

