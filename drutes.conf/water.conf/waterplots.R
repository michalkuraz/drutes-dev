#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# read and plot output with different capacity terms
plotname=args[1]
if(is.na(plotname)){
  plotname=""
}

mycol=colorRampPalette(c("darkred", "darkorange","goldenrod1","deepskyblue","royalblue2","darkblue"))

read_data_plot=function(filedr1,filedr2,ext,col1,col2,xlabs,ylabs,whatisplotted,skip=0,k=0,legpos="topleft",idfix=0,lims=1:1000){
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
    png(paste("obs_",whatisplotted,"_", plotname,".png",sep=""),width = 800, height = 600, units = "px")
    par(cex=1.7,mar=c(5,4.5,4,2))
    leg.txt=c()
    for(i in 1:ln_obs){
      if(i == 1){
        plot(obs[[i]][lims,col1],obs[[i]][lims,col2],type="l",pch=i,col=mycolors[i],ylab=ylabs,xlab=xlabs,ylim=c(mins,maxs),lwd=1.5)
      }else{
        par(new=T)
        plot(obs[[i]][lims,col1],obs[[i]][lims,col2],type="l",pch=i,col=mycolors[i],ylab="",xlab="",ylim=c(mins,maxs),axes=F,lwd=1.5)
      } 
      leg.txt[i]=paste("obs. time ",i-idfix)
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
  
ln_obs=read_data_plot('out/RE_matrix_theta-','drutes-dev/out/RE_matrix_theta-'
                 ,'.dat',col1=2,col2=3,ylabs=expression(paste("vol. water content [-]",sep=""))
                 ,xlabs="depth [cm]",'water',idfix=1,lims=300:400)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: water content vs. depth"))
}

ln_obs=read_data_plot('out/RE_matrix_press_head-','drutes-dev/out/RE_matrix_press_head-','.dat',col1=2,col2=3,ylabs="pressure head [cm]"
                 ,xlabs="depth [cm]",'press_head',idfix=0,lims=300:400)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: pressure head vs. depth"))
}
