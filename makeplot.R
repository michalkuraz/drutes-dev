#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# read and plot output with SHP

ln_args=length(args)
isNAmin=T
isNAmax=T
isNAplot=T
isNAleg=T
legopts=c("bottomright", "bottom", "bottomleft", 
          "left", "topleft", "top", "topright", "right", "center","none")
assign_args=T
if(ln_args>0){
  options(warn=-1)
  i=1
  while(assign_args){
    switch(args[i],
           "-name"={plotname=args[i+1];isNAplot=F},
           "-legpos"={legpos=args[i+1];isNAleg=F;if(!any(legopts==legpos)){stop(paste("-legpos only takes",legopts," as keyword. You entered:",args[i+1]))}},
           "-min"={if(is.na(as.numeric(args[i+1])))
           {stop(paste("-min only takes numeric values. You entered:",args[i+1]))}
             else
             {lims_min=as.numeric(args[i+1]);isNAmin=F};
             if(lims_min<0){stop(paste("-min only takes positive values. You entered:",args[i+1]))}
             if(lims_min%%1!=0){stop(paste("-min only takes integer values. You entered:",args[i+1]))}
           },
           "-max"={if(is.na(as.numeric(args[i+1])))
           {stop(paste("-max only takes numeric values. You entered:",args[i+1]))}
             else
             {lims_max=as.numeric(args[i+1]);isNAmax=F};
             if(lims_max<0){stop(paste("-max only takes positive values. You entered:",args[i+1]))}
             if(lims_max%%1!=0){stop(paste("-max only integer values. You entered:",args[i+1]))}
           },
           "-man"={readline("Manual:      -name - defines plotname 
                            -min - first index of data to be plotted
                            -max - last index of data to be plotted
                            -legpos - position of legend in the plot"); i=i-1
           },
           stop(paste("no valid argument, options are -name, -min, -legpos or -man. You entered:",args[i])))
    i=i+2
    if(i>ln_args){
      assign_args=F
    }
  }
  
}else{
  plotname=""
  lims_min=1
  lims_max=100
  isNAmin=T
  isNAmax=T
  isNAplot=T
  legpos="topleft"
}

if(isNAplot){
  plotname=""
}
if(isNAleg){
  legpos="topleft"
}

if(isNAmin){
  lims_min=1
}

if(isNAmax){
  lims_max=100
}

lims=lims_min:lims_max

mycol=colorRampPalette(c("darkred", "darkorange","goldenrod1","deepskyblue","royalblue2","darkblue"))

read_data_plot=function(filedr1,filedr2,ext,col1,col2,xlabs,ylabs,whatisplotted,skip=0,k=0,legpos=legpos,idfix=0,lims=lims){
  allread=F
  obs=list()
  n=1
  if(any(unlist(strsplit(unlist(strsplit(filedr1,split = "/"))[2],"_"))=="obspt")){
    legend_entry="obs. point"
  }else{
    legend_entry=" obs. time"
  }
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
    if(isNAmax){
      if(lims_min>1){
        if(lims[length(lims)]!=length(obs[[1]][,col1])){
          lims=lims_min:length(obs[[1]][,col1])
        }
      }else{
        if(lims[length(lims)]!=length(obs[[1]][,col1])){
          lims=1:length(obs[[1]][,col1])
          print("Entire domain is plotted")
        }
      }
      
      
    }
    
    mycolors=mycol(ln_obs)
    pname=paste("obs_",whatisplotted,"_", plotname,".png",sep="")
    idname=2
    while(file.exists(pname)){
      pname2=paste("obs_",whatisplotted,"_", plotname,idname,".png",sep="")
      #print(paste(pname, " already exists"))
      pname=pname2
      idname=idname+1
    }
    png(pname,width = 800, height = 600, units = "px")
    print(paste("plotting", pname))
    par(cex=1.9,mar=c(5,4.5,4,2))
    leg.txt=c()
    for(i in 1:ln_obs){
      if(i == 1){
        plot(obs[[i]][lims,col1],obs[[i]][lims,col2],type="l",pch=i,col=mycolors[i],ylab=ylabs,xlab=xlabs,ylim=c(mins,maxs),lwd=1.5,main=plotname)
      }else{
        par(new=T)
        plot(obs[[i]][lims,col1],obs[[i]][lims,col2],type="l",pch=i,col=mycolors[i],ylab="",xlab="",ylim=c(mins,maxs),axes=F,lwd=1.5)
      } 
      leg.txt[i]=paste(legend_entry,i-idfix)
    }
    ncols=1
    if(ln_obs>8){
      ncols=2
    }
    if(ln_obs>12){
      ncols=3
    }
    if(legpos!="none"){
      legend(legpos,leg.txt,ncol=ncols,col=mycolors,seg.len=0.5,lty=1,lwd=2,bty="n",cex=0.8)
    }
    invisible(dev.off())
  }
  return(ln_obs)
}



ln_obs=read_data_plot('out/heat_temperature-','drutes-dev/out/heat_temperature-'
                      ,'.dat',col1=2,col2=3,ylabs=expression(paste("temperature [",degree*C,"]",sep=""))
                      ,xlabs="wall depth [m]",'time',idfix=1,skip=0,k=0,legpos=legpos,lims=lims)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: Temperature vs. depth "))
}

ln_obs=read_data_plot('out/obspt_heat-','drutes-dev/out/obspt_heat-'
                      ,'.out',col1=1,col2=4,ylabs=expression(paste("heat flux [W",~m^-2,"]",sep=""))
                      ,xlabs="time [h]",'point',skip=10,k=1,legpos="topright",lims=lims)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation points created: heat flux vs. time"))
}

# ln_obs=read_data_plot('out/obspt_heat-','drutes-dev/out/obspt_heat-'
#                       ,'.out',col1=1,col2=2,ylabs=expression(paste("temperature [°C]",sep=""))
#                       ,xlabs="time [h]",'temp_point',skip=10,k=1,legpos="topright",lims=lims)
# if(ln_obs>0){
#   print(paste("plot of",ln_obs,"observation points created: temperature vs. time"))
# }
ln_obs=read_data_plot('out/obspt_heat-','drutes-dev/out/obspt_heat-'
                      ,'.out',col1=1,col2=5,ylabs=expression(paste("cumulative heat flux [W",~m^-2,"]",sep=""))
                      ,xlabs="time [h]",'point_cum',skip=10,k=1,legpos="topright",lims=lims)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation points created: cumulative heat flux vs. time"))
}



ln_obs=read_data_plot('out/obspt_RE_matrix-','drutes-dev/out/obspt_RE_matrix-'
                      ,'.out',col1=1,col2=3,ylabs=expression(paste("vol. water content [-]",sep=""))
                      ,xlabs="time [T]",'water_point',idfix=0,lims=lims,legpos = legpos ,skip=5,k=1)
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation points created: water content vs. time"))
}


# ln_obs=read_data_plot('out/obspt_ADER_in_liquid-','drutes-dev/out/obspt_ADER_in_liquid-'
#                       ,'.out',col1=1,col2=3,ylabs=expression(paste("concentration in liquid [M ",L^-3,"]",sep=""))
#                       ,xlabs="time [T]",'conc_point',idfix=0,lims=lims,legpos = legpos ,skip=5,k=1)
# if(ln_obs>0){
#   print(paste("plot of",ln_obs,"observation points created: concentration in liquid vs. time"))
# }
# 
# ln_obs=read_data_plot('out/RE_matrix_theta-','drutes-dev/out/RE_matrix_theta-'
#                        ,'.dat',col1=2,col2=3,ylabs=expression(paste("vol. water content [-]",sep=""))
#                        ,xlabs="depth [L]",'water',idfix=1,lims=lims,legpos = legpos )
# if(ln_obs>0){
#   print(paste("plot of",ln_obs,"observation times created: water content vs. depth"))
# }
# 
# ln_obs=read_data_plot('out/ADER_in_liquid_solute_concentration-','drutes-dev/out/ADER_in_liquid_solute_concentration-'
#                       ,'.dat',col1=2,col2=3,ylabs=expression(paste("concentration in liquid [M ",L^-3,"]",sep=""))
#                       ,xlabs="depth [L]",'conta',idfix=1,lims=lims,legpos = legpos )
# if(ln_obs>0){   
#   print(paste("plot of",ln_obs,"observation times created: concentrationvs. depth"))
# }
# 
# ln_obs=read_data_plot('out/obspt_RE_matrix-','drutes-dev/out/obspt_RE_matrix-'
#                       ,'.out',col1=1,col2=2,ylabs=expression(paste("pressure head [cm]",sep=""))
#                       ,xlabs="time [days]",'press_point',idfix=0,lims=lims,legpos = legpos ,skip=5,k=1)
# if(ln_obs>0){
#   print(paste("plot of",ln_obs,"observation points created: pressure head vs. time"))
# }

ln_obs=read_data_plot('out/RE_matrix_theta-','drutes-dev/out/RE_matrix_theta-'
                      ,'.dat',col1=2,col2=3,ylabs=expression(paste("vol. water content [-]",sep=""))
                      ,xlabs="depth [cm]",'water',idfix=1,lims=lims,legpos = legpos )
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: water content vs. depth"))
}

ln_obs=read_data_plot('out/RE_matrix_press_head-','drutes-dev/out/RE_matrix_press_head-','.dat',col1=2,col2=3,ylabs="pressure head [cm]"
                      ,xlabs="depth [cm]",'press_head',idfix=0,lims=lims,legpos = legpos )
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: pressure head vs. depth"))
}

ln_obs=read_data_plot('out/RE_matrix_flux-','drutes-dev/out/RE_matrix_flux-','.dat',col1=2,col2=3,ylabs=expression(paste("flux [cm ",d^-1,"]")),
                      xlabs="depth [cm]",'flux',idfix=0,lims=lims,legpos = legpos )
if(ln_obs>0){
  print(paste("plot of",ln_obs,"observation times created: flux vs. depth"))
}