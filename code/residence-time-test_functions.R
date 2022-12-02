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
  
  ##calculate the expected DPMs of the standard based on t
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
  Tau$MLeu <- (M_Leu * Tau$DPMKills)/Voucher
  #Convert MLeu to ug C/L/hr
  Tau$ugCLhr <- Tau$MLeu * 131.2 * (1/0.073)*0.86*2*1000000
  Tau$uMChr <- Tau$ugCLhr *0.083333
  Tau$log_uMChr <- log(Tau$uMChr, 10)
  return(Tau)
}

#BP[[i]]$Leucine <- ((BP[[i]]$DPM / 2.2e12) / 153) / 1000
# DPM*1Ci/2.2e12DPM *1mmolLeu/153Ci * 1molLeu/1000mmol
#BP[[i]]$Leucine.per <- (BP[[i]]$Leucine * (1/1) * (1/0.0015))
# Leu incorporated * 1/time (hrs) * 1/vol (L)
#BP[[i]]$Protein <- ( BP[[i]]$Leucine.per / 0.073 ) * 131.2
# mol leu * 1 mol protein/0.073 mol leu * 131.2 g protein/1mol protein
#BP[[i]]$Carbon <- BP[[i]]$Protein * (1/0.63) * (0.54/1) * (10^6)
# g protein * 1 g DW/0.63 g Pro * 0.54 g C/1g DW = gC/l/hr  * 10^6 = ugC/L/Hr

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