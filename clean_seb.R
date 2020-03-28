
# replace filename as necessary
filename <- "~/Desktop/CS_07/surface_energy.out"


clean_seb <- function(filename){
  dta <- read.table(filename, skip = 6)
  dta_new <- as.data.frame(matrix(ncol = ncol(dta), nrow = length(unique(dta$V1))))
  k <- 1
  for(rows in 1:(nrow(dta)-1)){
    t_old <- dta[rows,1]
    t_new <- dta[rows+1,1]
    
    if(t_new == t_old){
      
    }else{
      dta_new[k,] <- dta[rows,]
      k <- k+1
    }
    if(rows == (nrow(dta)-1)){
      if(t_new == t_old){
        dta_new[k,] <- dta[rows+1,]
      }else{
        dta_new[k,] <- dta[rows+1,]
      }
    }
  }
  return(dta_new)
}

data_clean<- clean_seb(filename)

names(data_clean) <- c("time","Rad", "Hs", "LE","E", "G")



plot(data_clean$time/86400, data_clean$E*86400*100,
     xlab = "time [d]", ylab = "Evaproation rate [cm/d]",
     pch = 4)


ylims <- c(-100, 300)
plot(data_clean$time/86400, data_clean$LE,
     xlab = "time [d]", ylab = "Energy flux[W/m2]",
     type = "l", ylim = ylims, col = "royalblue2", lwd = 2)
par(new = T)
plot(data_clean$time/86400, data_clean$G,
     xlab = "", ylab = "", col = "darkgreen", lwd = 2,
     type = "l", ylim = ylims, xaxt = "n", yaxt = "n")
par(new = T)
plot(data_clean$time/86400, data_clean$Hs,
     xlab = "", ylab = "", col = "goldenrod1", lwd = 2,
     type = "l", ylim = ylims, xaxt = "n", yaxt = "n")
par(new = T)
plot(data_clean$time/86400, data_clean$Rad,
     xlab = "", ylab = "", col = "firebrick1", lwd = 2,
     type = "l", ylim = ylims, xaxt = "n", yaxt = "n")

leg.txt <- c("LE", "G", "Hs", "Rn")
legend("topright", leg.txt, col = c("royalblue2","darkgreen","goldenrod1","firebrick1"), lty = 1, bty ="n", lwd = 2)
