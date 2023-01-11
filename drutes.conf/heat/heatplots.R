#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# read and plot output with different capacity terms
plotname=args[1]
if(is.na(plotname)){
  plotname=""
}

mycol=colorRampPalette(c("darkred", "darkorange","goldenrod1","deepskyblue","royalblue2","darkblue"))

read_data_plot=function(filedr1,filedr2,ext,col1,col2,xlabs,ylabs,whatisplotted,skip=0,k=0,legpos="topleft",idfix=0){
  allread=F
  obs=list()
  n=1
  while(!allread){
    filedr1x=paste(filedr1,k,ext,sep="")
    filedr2x=paste(filedr2,k,ext,sep="")
    if(file.exists(filedr1x)){
      obs[[n]]=read.table(filedr1x,skip=skip)
    }else if(file.exists(filedr2x)){
      obs[[n]]=read.table(filedr2x,skip=skip)
    }else{
      allread=TRUE
    }
    if(!allread){
      if(n==1){
        mins=min(obs[[n]][,col2])
        maxs=max(obs[[n]][,col2])
      }else{
        if(min(obs[[n]][,col2]) < mins){
          mins=min(obs[[n]][,col2])
        }
        if(max(obs[[n]][,col2]) > maxs){
          maxs=max(obs[[n]][,col2])
        }
      }
    }
    
    k=k+1
    n=n+1
  }
  ln_obs=length(obs)
  if(ln_obs>0){
    mycolors=mycol(ln_obs)
    png(paste("obs_",whatisplotted,"_temp_",plotname,".png",sep=""),width = 800, height = 600, units = "px")
    par(cex=1.7,mar=c(5,4.5,4,2))
    leg.txt=c()
    for(i in 1:ln_obs){
      if(i == 1){
        plot(obs[[i]][,col1],obs[[i]][,col2],type="l",pch=i,col=mycolors[i],ylab=ylabs,xlab=xlabs,ylim=c(mins,maxs),lwd=1.5)
      }else{
        par(new=T)
        plot(obs[[i]][,col1],obs[[i]][,col2],type="l",pch=i,col=mycolors[i],ylab="",xlab="",ylim=c(mins,maxs),axes=F,lwd=1.5)
      } 
      leg.txt[i]=paste("obs.", whatisplotted,i-idfix)
    }
    ncols=1
    if(ln_obs>8){
      ncols=2
    }
    if(ln_obs>12){
      ncols=3
    }
    legend(legpos,leg.txt,ncol=ncols,col=mycolors,seg.len=0.5,lty=1,lwd=2,bty="n")
    invisible(dev.off())
  }
  return(ln_obs)
  }
  
ln_obs=read_data_plot('out/heat_temperature-','drutes-dev/out/heat_temperature-'
                 ,'.dat',col1=2,col2=3,ylabs=expression(paste("temperature [",degree*C,"]",sep=""))
                 ,xlabs="wall depth [m]",'time',idfix=1)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: Temperature vs. wall depth "))
}
ln_obs=read_data_plot('out/obspt_heat-','drutes-dev/out/obspt_heat-','.out',col1=1,col2=4,ylabs=expression(paste("heat flux [W",~m^2,"]",sep=""))
                 ,xlabs="time [h]",'point',skip=10,k=1,"topright")
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation points created: heat flux vs. time"))
}
ln_obs=read_data_plot('out/obspt_heat-','drutes-dev/out/obspt_heat-','.out',col1=1,col2=5,ylabs=expression(paste("cumulative heat flux [W",~m^2,"]",sep=""))
                 ,xlabs="time [h]",'point_cum',skip=10,k=1)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation points created: cumulative heat flux vs. time"))
}
