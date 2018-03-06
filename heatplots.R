# read and plot output with different capacity terms
allread=F
obs=list()
k=1
while(!allread){
  filedr1=paste('out/heat_temperature-',k,'.dat',sep="")
  filedr2=paste('drutes-dev/out/heat_temperature-',k,'.dat',sep="")
  if(file.exists(filedr1)){
    obs[[k]]=read.table(paste('out/heat_temperature-',k,'.dat',sep=""))
  }else if(file.exists(filedr2)){
      obs[[k]]=read.table(paste('out/heat_temperature-',k,'.dat',sep=""))
  }
  else{
      allread=TRUE
  }
  k=k+1
}

dev.new()
plot(obs1$V2,obs1$V3,type='l',ylab='Temperature [°C]',xlab='Wall [cm]',col=grey(.9),ylim=c(0,30))
par(new=T)
plot(obs2$V2,obs2$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(.8),ylim=c(0,30))
par(new=T)
plot(obs3$V2,obs3$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(0.7),ylim=c(0,30))
par(new=T)
plot(obs4$V2,obs4$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(0.6),ylim=c(0,30))
par(new=T)
plot(obs5$V2,obs5$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(0.5),ylim=c(0,30))
par(new=T)
plot(obs6$V2,obs6$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(0.4),ylim=c(0,30))
par(new=T)
plot(obs7$V2,obs7$V3,type='l',ylab='',xlab='',xaxt='n',yaxt='n',col=grey(0.3),ylim=c(0,30))