### load the necessary packages

library(tidyverse)
library(readxl)
library(PBSmodelling)
library(foreach)
library(doParallel)

## Read data base

EcoData <- read_csv("EcoData2019.csv")


##Make a list of the data input spearated by SEQ_NO

Data2 <- EcoData %>% split(.$TIME)


### make functions for SQS

############################################
##########Function for each simulation######
############################################

RichnessWhile <- function(Data, quorum = 0.95){
  #Reorder the dataset randomly
  Order <- sample(1:nrow(Data), nrow(Data))
  Data <- Data[Order,]
  # Add species richness, sample and goods U as new columns
    Data$Sample <- NA
  Data$GoodsU <- NA
    # Start Richness, sample number and Good's U as 0
  
  GoodsU <- 0
  i <- 0
  #loop until Goods'u is equal or bigger than the thresholded quorum
  while(GoodsU < quorum){
    i <- i + 1
    
    Data$Sample[i] <- i
    # Calculate Good's U
    GoodsU <- 1 - (nrow(Data[1:i,] %>% group_by(GENUS) %>% summarise(N = n()) %>% filter(N == 1))/nrow(Data[1:i,]))
    #Asign it to the proper row
    Data$GoodsU[i] <- GoodsU
    #### Give out a message every 100 samples
    if(i %% 100 == 0){
      message(paste("Sample", i, "of", nrow(Data)))
      message(paste("Good's U is:", round(Data$GoodsU[i], 3)))
    }
  }
  
  ## Calculate species richness
  Richness <- nrow(Data[1:i,] %>% group_by(GENUS) %>% summarise(NSpp = n()))
  Trop <- nrow(Data[1:i,] %>% filter(CLIM == "tropical") %>% group_by(GENUS) %>% summarise(NSpp = n()))
  Temp <- nrow(Data[1:i,] %>% filter(CLIM == "temperate")%>% group_by(GENUS) %>% summarise(NSpp = n()))
  ## Calculate richness of CODEs
  Functional <- Data[1:i,] %>% group_by(CODE) %>% summarise(NEco = n()) %>%  spread(key = CODE, value = NEco, fill = 0)
  Ecos <- nrow(Data[1:i,] %>% group_by(CODE) %>% summarise(NEco = n()))
  EcosTrop <- nrow(Data[1:i,]  %>% filter(CLIM == "tropical") %>% group_by(CODE) %>% summarise(NEco = n()))
  EcosTemp <- nrow(Data[1:i,]  %>% filter(CLIM == "temperate") %>% group_by(CODE) %>% summarise(NEco = n()))
  
  ## Asign them to the proper columns
  ### Make a plot of samples against Good's U
  g2 <- ggplot(Data[1:i,], aes(x = Sample, y = GoodsU)) + ylim(c(0,1))+ geom_line() + geom_hline(yintercept = quorum, lty = 2, color = "red") + theme_classic()
  #### Return the graph, the Ensemble of the samples, and a dataframe with the Good's U, the Number of samples, the species Richness and the number of CODES
  return(list(Graph = g2, Ensemble =  Data[1:i,], Data = data.frame(GoodsU = GoodsU, Sample = i, Richness = Richness, RichnessTemp = Temp, RichnessTrop = Trop, Ecos = Ecos, EcosTrop = EcosTrop, EcosTemp = EcosTemp), Functional = Functional))
}


#############################################################
#################Function to make many simulations of sqs####
#############################################################

ResampleWhile <- function(DF, nsim = 10, Quorum, method = "Geometric", ncores = NULL){
  
  ##Loop through various simulations of richnesswhile
  if(is.null(ncores)){
    ### Create list for samples and ensembles
    Resample <- list()
    EcoEnt <- list()
    for(i in 1:nsim){
      Temp <- RichnessWhile(DF, quorum = Quorum)
      Resample[[i]] <- Temp$Data
      Resample[[i]]$Simulation <- i
      EcoEnt[[i]] <- Temp$Functional
      message(paste("Finished simulation", i, "of", nsim))
    }
    ### Join all the data.frames with the number of data
    Resample <- bind_rows(Resample)
    EcoEnt <- bind_rows(EcoEnt)
    EcoEnt[is.na(EcoEnt)] <- 0
  }
  
  if(!is.null(ncores)){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    ### Create list for samples and ensembles
    
    res <- foreach(i = 1:nsim, .packages = c("tidyverse"), .export=c("RichnessWhile")) %dopar% {
      
      Temp <- RichnessWhile(DF, quorum = Quorum)
      Resample <- Temp$Data
      Resample$Simulation <- i
      EcoEnt <- Temp$Functional
      return(list(Resample, EcoEnt))
      
    }
    
    stopCluster(cl)
    
    Resample <- bind_rows(lapply(res,function(x){x[[1]]}))
    EcoEnt <- bind_rows(lapply(res,function(x){x[[2]]}))
    EcoEnt[is.na(EcoEnt)] <- 0
  }
  ### Calculate the median richness to select the best simulations
  if(method == "Median"){
    richness <- median(Resample$Richness)
    RichnessTemp <- median(Resample$RichnessTemp)
    RichnessTrop <- median(Resample$RichnessTrop)
    Ecos <- median(Resample$Ecos)
    EcosTrop <- median(Resample$EcosTrop)
    EcosTemp <- median(Resample$EcosTemp)
    
  }
  if(method == "Geometric"){
    richness <- calcGM(Resample$Richness, offset = 0.0001, exzero = FALSE)
    RichnessTemp <- calcGM(Resample$RichnessTemp, offset = 0.0001, exzero = FALSE)
    RichnessTrop <- calcGM(Resample$RichnessTrop, offset = 0.0001, exzero = FALSE)
    Ecos <- calcGM(Resample$Ecos, offset = 0.0001, exzero = FALSE)
    EcosTrop <- calcGM(Resample$EcosTrop, offset = 0.0001, exzero = FALSE)
    EcosTemp <- calcGM(Resample$EcosTemp, offset = 0.0001, exzero = FALSE)
    
  }
  

  gc()
  ### Plot an histogram of the richness frequency of simulations and add a doted red line on the median richness
  #hist(Resample$Richness, main = paste("Resampled richness after", nsim, "simulations"), xlab = "Richness")
  #abline(v = richness, lty = 2, col = "red")
  ## Calculate the delta of the selected richness to establish the colsest ensemble
  Resample <- Resample %>% mutate(Delta =abs(Richness - richness)) %>% arrange(Delta, desc(Sample))
  EcoEnt <- EcoEnt %>% summarise_all(funs(calcGM(., offset = 0.0001, exzero = FALSE)))
  
  
  
  ### Return the simulation reults (Resample), the calculated species richness, and the selected ensamble
  return(list(Resample = Resample, richness = richness, RichnessTemp = RichnessTemp, RichnessTrop = RichnessTrop, Ecos = Ecos, EcosTrop = EcosTrop, EcosTemp = EcosTemp, EcoEnt = EcoEnt))
}


Bootstraps = 500

A <- Sys.time()
Results <- list()
Results2 <- list()

set.seed(2019)


for(x in 1:length(Data2)){
  s<- data.frame(Richness = rep(NA, Bootstraps), RichnessTemp = rep(NA, Bootstraps), RichnessTrop = rep(NA, Bootstraps), Ecos = rep(NA, Bootstraps), EcosTrop = rep(NA, Bootstraps), EcosTemp = rep(NA, Bootstraps))
  b <- list()
  for(i in 1:nrow(s)){
    bootstrapped<-Data2[[x]][sample(1:nrow(Data2[[x]]),replace=TRUE),]
    a <- ResampleWhile(DF = bootstrapped ,nsim = 500, Quorum = 0.4, method = "Geometric", ncores = 40)
    s$Richness[i] <- a$richness
    s$RichnessTemp[i] <- a$RichnessTemp
    s$RichnessTrop[i] <- a$RichnessTrop
    s$Ecos[i] <- a$Ecos
    s$EcosTrop[i] <- a$EcosTrop
    s$EcosTemp[i] <- a$EcosTemp
    b[[i]] <- a$EcoEnt
    gc()
    if(i %% 50 == 0){
      print(paste("bootstrap", i , "of", nrow(s), "ready!!!", Sys.time()))
    }
    
  }
  
  Results[[x]] <- data.frame(Meassurement = c("Richness", "RichnessTemp", "RichnessTrop", "Ecos", "EcosTrop", "EcosTemp"), Mean = rep(NA,6), Lower = rep(NA,6), Upper = rep(NA,6), Time = unique(Data2[[x]]$TIME))
  Results[[x]][1,2:4] <- quantile(s$Richness,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]][2,2:4] <- quantile(s$RichnessTemp,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]][3,2:4] <- quantile(s$RichnessTrop,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]][4,2:4] <- quantile(s$Ecos,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]][5,2:4] <- quantile(s$EcosTrop,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]][6,2:4] <- quantile(s$EcosTemp,probs=c(0.5,0.025,0.975), na.rm = T)
  Results[[x]] <- bind_rows(Results[[x]])
  b <- bind_rows(b)
  b[is.na(b)] <- 0
  b <- b %>% gather(key = "Eco", value = "Abundance") %>% group_by(Eco) %>% summarise(Mean = quantile(Abundance, 0.5), Lower = quantile(Abundance, 0.025), Upper = quantile(Abundance, 0.975)) %>% mutate(Time = unique(Data2[[x]]$TIME))
  Results2[[x]] <- b
  message(paste("bin", x, "ready", Sys.time()))
  saveRDS(Results, "Results.rds")
  saveRDS(Results2, "Results2.rds")
}

B <- Sys.time()

Results <- bind_rows(Results)
saveRDS(Results, "Results.rds")

Results2 <- bind_rows(Results2)
saveRDS(Results2, "Results2.rds")


Results3 <- Results %>% mutate(Type = case_when(Meassurement %in% c("Ecos", "EcosTemp", "EcosTrop") ~ "Ecological entities",Meassurement %in% c("Richness", "RichnessTemp", "RichnessTrop") ~ "Richness")) %>% mutate(Meassurement = case_when(Meassurement %in% c("Ecos", "Richness") ~ "Total", Meassurement %in% c("EcosTrop", "RichnessTrop") ~ "Tropical", Meassurement %in% c("EcosTemp", "RichnessTemp") ~ "Temperate")) %>% mutate(Type = relevel(factor(Type), "Richness")) %>% mutate(Meassurement = relevel(factor(Meassurement), "Total"))


Results3 <- bind_rows(Results3)
saveRDS(Results3, "Results3.rds")


ggplot(Results3, aes(y = Mean, x = Time, group = Meassurement)) + geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = Meassurement), alpha = 0.3) + geom_path(aes(color = Meassurement)) + geom_point(aes(color = Meassurement)) + theme_bw() + facet_grid(Type ~ ., scales = "free_y") + theme(legend.position = "bottom") + scale_x_reverse() + xlab("Time [Ma]") + ylab("") + geom_vline(xintercept = 250, lty = 2, color = "red") + geom_vline(xintercept = 200, lty = 2, color = "red")

B-A
