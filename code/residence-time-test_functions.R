################################################################################
#                                                                              #
# residence-time-test Functions Source Code                                    #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Emmi Mueller                                                     #
#                                                                              #
#                                                                              #
################################################################################
#                                                                              #
# Notes:                                                                       #
#                                                                              #
# Issues:                                                                      #
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#                                                                              #
#                                                                              #
################################################################################

BP_fxn <- function(CPMs, Tau, set){
  ##extract whole experiment info from top of csv
  #date experiment was run
  date <- as.Date(as.character(CPMs[1,2]), "%m/%d/%Y")
  #date the standard was produced
  date_std <- as.Date(as.character(CPMs[2,2]), "%m/%d/%Y")
  #DPM of the standard at date of production
  DPM_std <- as.double(as.character(CPMs[3,2]))
  #DPM of the standard based on scintillation counter on experiment date
  DPM_curr <- as.double(as.character(CPMs[4,2]))
  #half life of tritium - 12.346 years
  half_life <- as.double(as.character(CPMs[5,2]))
  #Mols of leucine in each tube based on hot leucine stock concentration
  M_Leu <- as.double(as.character(CPMs[6,2]))
  #CPMs of the voucher on experiment date
  Voucher <- as.double(as.character(CPMs[7,2]))
  
  ##remove whole experiment info from top of dataframe
  CPMs <- CPMs[-c(1:9),]
  colnames(CPMs) <- c("Sample", "CPM", "Kill")
  
  ##calculate time from the experiment date to the standard production date
  t <- as.numeric(date - date_std)/365
  
  ##calculate the expected DPMs of the standard based on t N(t) = N(0)*exp((ln(2)/t0.5)*t))
  DPM_exp <- (DPM_std)*exp((-0.693/half_life)*t)
  
  ##calculate scintillation efficiency as DPM ratio
  efficiency <- DPM_curr/DPM_exp
  #divide CPMs into kill and sample dataframes
  Kills <- subset(CPMs, Kill == "T")
  CPMs <- subset(CPMs, Kill == "F")
  #convert CPMs to DPMs, DPMs = CPMs/efficiency
  CPMs$CPM <- as.numeric(as.character(CPMs$CPM))
  CPMs$DPM <- CPMs$CPM / efficiency
  Kills$CPM <- as.numeric(as.character(Kills$CPM))
  Kills$DPM <- Kills$CPM / efficiency
  #average DPMs for each sample and add to Tau
  for(x in unique(CPMs$Sample)){
    Tau[Tau$Tau == x, "DPM"] <- as.numeric(mean(CPMs[CPMs$Sample == x, "DPM"]))
  }

  #for each sample, subtract the corresponding kill DPM
  for (x in unique(Tau$Tau)){
    if(Tau[Tau$Tau == x, "Set"] == set){
      Tau[Tau$Tau == x, "DPMKills"] <- Tau[Tau$Tau ==x, "DPM"] - (as.numeric(Kills[Kills$Sample == x, "CPM"])/efficiency)
    }
  }

  #Determine Mols Leucine based on MLeu_sample = MLeu * DPM/voucher
  Tau$molLeu <- Tau$DPMKills *(1/2.22e12)*(1/161)*(1/1000)
  Tau$molLeuLhr <- Tau$molLeu * 1 * (1/0.0015)
  #Convert molLeu/L/hr to gC/L/hr
  Tau$gCLhr <- Tau$molLeuLhr * 131.2 * (1/0.073)*0.86*2
  #Conversion of g to mol to umol to get umolC/hr
  Tau$uMChr <- Tau$gCLhr * 0.083333 * 100000
  Tau$log_uMChr <- log(Tau$uMChr, 10)
  return(Tau)
}

N_fxn <- function(N, Tau){
  for(x in unique(N$Sample)){
    Tau$N[Tau$Tau == x] <- mean(N$N[N$Sample == x])
  }
  Tau$log_N <- log(Tau$N, 10)
  
  return(Tau)
}

OT_fxn <- function(OT, Tau){
  for(x in unique(OT$Tau)){
    Tau$OT[Tau$Tau == x] <- mean(OT$OT[OT$Tau == x])-mean(OT$OT[OT$Tau == 0])
  }
  
  return(Tau)
}

EP_fxn <- function(EP, Tau){
  EP$Tau <- as.numeric(EP$Tau)
  for(x in unique(EP$Tau)){
    Tau$NumRes[Tau$Tau == x] <- EP$NumRes[EP$Tau == x]
    Tau$AvgRes[Tau$Tau == x] <- EP$Avg[EP$Tau ==x]
  }
  
  return(Tau)
}

S.cal <- function(x = ""){
  rowSums(x > 0) * 1
}

S_fxn <- function(S, Tau){
  for(x in unique(S$V1)){
    Tau$Day_0.S[Tau$Day_0_Seq == x] <- S$V2[S$V1 == x]
    Tau$Day_20.S[Tau$Day_20_Seq == x] <- S$V2[S$V1 == x]
  }
  return(Tau)
}

SimpE.cal <- function(x = ""){
  S <- S.cal(x)
  x = as.data.frame(x)
  D <- diversity(x, index = "inv")
  E <- (D)/S
  return(E)
}

SimpE_fxn <- function(SimpE, Tau){
  for(x in unique(SimpE$V1)){
    Tau$Day_0.SimpE[Tau$Day_0_Seq == x] <- SimpE$V2[SimpE$V1 == x]
    Tau$Day_20.SimpE[Tau$Day_20_Seq == x] <- SimpE$V2[SimpE$V1 == x]
  }
  return(Tau)
}

Evar.cal <-function(x){
  x <- as.vector(x[x > 0])
  return(1 - (2/pi)*atan(var(log(x))))
}

Evar_fxn <- function(OTU_r, Tau){
  for(x in unique(na.omit(Tau$Day_0_Seq))){
    Tau$Day_0.Evar[Tau$Day_0_Seq == x] <- Evar.cal(OTU_r[x,])
  }
  for(x in unique(na.omit(Tau$Day_20_Seq))){
    Tau$Day_20.Evar[Tau$Day_20_Seq == x] <- Evar.cal(OTU_r[x,])
  }
  return(Tau)
}

S.ace <- function(x = "", thresh = 10){
  x <- x[x>0]
  S.abund <- length(which(x > thresh))
  S.rare <- length(which(x <= thresh))
  singlt <- length(which(x ==1))
  N.rare <- sum(x[which(x <= thresh)])
  C.ace <- 1 - (singlt/N.rare)
  i <- c(1:thresh)
  count <- function(i, y){
    length(y[y == i])
  }
  a.1 <- sapply(i, count, x)
  f.1 <- (i * (i -1)) * a.1
  G.ace <- (S.rare/C.ace) * (sum(f.1)/(N.rare*(N.rare-1)))
  S.ace <- S.abund + (S.rare/C.ace) + (singlt/C.ace) * max(G.ace,0)
  return(S.ace)
}