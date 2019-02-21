### load the necessary packages

library(tidyverse)
library(readxl)
library(EnvStats)
library(foreach)
library(doParallel)

## Read data base
Data_Imput <- read_excel("Data.Imput.xlsx") %>% mutate(SEQ_NO = as.numeric(SEQ_NO)) %>% arrange(SEQ_NO)

##Make a list of the data input spearated by SEQ_NO

Data2 <- Data_Imput %>% split(.$SEQ_NO)


### make functions for SQS

############################################
##########Function for each simulation######
############################################

RichnessWhile <- function(Data, quorum = 0.95){
  #Reorder the dataset randomly
  Order <- sample(1:nrow(Data), nrow(Data))
  Data <- Data[Order,]
  # Add species richness, sample and goods U as new columns
  Data$Richness <- NA
  Data$Sample <- NA
  Data$GoodsU <- NA
  Data$Ecos <- NA
  # Start Richness, sample number and Good's U as 0
  Richness <- 0
  Sample <- 0
  GoodsU <- 0
  Ecos <- 0
  i <- 0
  #loop until Goods'u is equal or bigger than the thresholded quorum
  while(GoodsU < quorum){
    i <- i + 1
    ## Calculate species richness
    Richness <- nrow(Data[1:i,] %>% group_by(GENUS) %>% summarise(NSpp = n()))
    ## Calculate richness of ecocodes
    Ecos <- nrow(Data[1:i,] %>% group_by(ECOCODE) %>% summarise(NEco = n()))
    ## Asign them to the proper columns
    Data$Richness[i] <- Richness
    Data$Ecos[i] <- Ecos
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
  ### Make a plot of samples against Good's U
  g2 <- ggplot(Data, aes(x = Sample, y = GoodsU)) + ylim(c(0,1))+ geom_line() + geom_hline(yintercept = quorum, lty = 2, color = "red") + theme_classic()
  #### Return the graph, the Ensemble of the samples, and a dataframe with the Good's U, the Number of samples, the species Richness and the number of ECOCODES
  return(list(Graph = g2, Ensemble =  Data[1:i,], Data = data.frame(GoodsU = GoodsU, Sample = i, Richness = Richness, Ecos = Ecos)))
}

#############################################################
#################Function to make many simulations of sqs####
#############################################################

ResampleWhile <- function(DF, nsim = 10, Quorum, method = "Geometric", ncores = NULL){
  
  ##Loop through various simulations of richnesswhile
  if(is.null(ncores)){
    ### Create list for samples and ensembles
    Resample <- list()
    Ensembles <- list()
    for(i in 1:nsim){
      Temp <- RichnessWhile(DF, quorum = Quorum)
      Resample[[i]] <- Temp$Data
      Resample[[i]]$Simulation <- i
      Ensembles[[i]] <- Temp$Ensemble
      Ensembles[[i]]$Simulation <- i
      message(paste("Finished simulation", i, "of", nsim))
    }
    ### Join all the data.frames with the number of data
    Resample <- bind_rows(Resample)
    ### Join all the data.frames with the Ensembles
    Ensembles <- bind_rows(Ensembles)
  }
  
  if(!is.null(ncores)){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    ### Create list for samples and ensembles
    
    res <- foreach(i = 1:nsim, .packages = c("tidyverse"), .export=c("RichnessWhile")) %dopar% {
      
      Temp <- RichnessWhile(DF, quorum = Quorum)
      Resample <- Temp$Data
      Resample$Simulation <- i
      Ensembles <- Temp$Ensemble
      Ensembles$Simulation <- i
      return(list(Resample, Ensembles))
      
    }
    
    stopCluster(cl)
    
    Resample <- bind_rows(lapply(res,function(x){x[[1]]}))
    Ensembles <- bind_rows(lapply(res,function(x){x[[2]]}))
  }
  ### Calculate the median richness to select the best simulations
  if(method == "Median"){
    richness <- median(Resample$Richness)  
  }
  if(method == "Geometric"){
    richness <- geoMean(Resample$Richness)
  }
  gc()
  ### Plot an histogram of the richness frequency of simulations and add a doted red line on the median richness
  #hist(Resample$Richness, main = paste("Resampled richness after", nsim, "simulations"), xlab = "Richness")
  #abline(v = richness, lty = 2, col = "red")
  ## Calculate the delta of the selected richness to establish the colsest ensemble
  Resample <- Resample %>% mutate(Delta =abs(Richness - richness)) %>% arrange(Delta, desc(Sample))
  #### Select the closest ensemble
  Selected <- Ensembles %>% filter(Simulation == Resample$Simulation[1])
  ### Show the number of ecocodes in that simulation
  Ecos <- length(unique(Selected$Ecos))
  ### Return the simulation reults (Resample), the calculated species richness, and the selected ensamble
  return(list(Resample = Resample, richness = richness, Ecos = Ecos, Selected = Selected))
}


A <- Sys.time()
Results <- list()
set.seed(2019)
for(x in 1:10){
  s<- data.frame(Richness = rep(NA, 100), Eco = rep(NA, 100))
  for(i in 1:nrow(s)){
    bootstrapped<-Data2[[x]][sample(1:nrow(Data2[[x]]),replace=TRUE),]
    a <- ResampleWhile(DF = bootstrapped ,nsim = 500, Quorum = 0.4, method = "Geometric", ncores = 10)
    s$Richness[i] <- a$richness
    s$Eco[i] <- a$Ecos
    gc()
    print(paste("bootstrap", i , "of", nrow(s), "ready!!!"))
  }
  
  # To remove any NAs (bootstrapping can produce a lot of pseudo-singletons)
  # Returns the median and 95% confidence intervals
  
  Results[[x]] <- data.frame(Meassurement = c("Richness", "Eco"), Mean = rep(NA,2), Lower = rep(NA,2), Upper = rep(NA,2), Bin = unique(Data2[[x]]$SEQ_NAME), NumericBin = unique(Data2[[x]]$SEQ_NO))
  Results[[x]][1,2:4] <- quantile(s$Richness,probs=c(0.5,0.025,0.975))
  Results[[x]][2,2:4] <- quantile(s$Eco,probs=c(0.5,0.025,0.975))
  message(paste("bin", x, "ready" ))
  saveRDS(Results, "Results.rds")
}

B <- Sys.time()

Results <- bind_rows(Results)
saveRDS(Results, "Results.rds")
Rich <- dplyr::filter(Results, Meassurement == "Richness")

ggplot(Rich, aes(x = NumericBin, y = Mean)) + geom_ribbon(aes(ymax = Upper, ymin = Lower), alpha = 0.5, fill = "red")+ geom_line() + theme_classic() + ylab("Richness") + xlab("Bin") + scale_x_continuous(breaks = unique(Rich$NumericBin))
ggplot(dplyr::filter(Results, Meassurement == "Eco"), aes(x = NumericBin, y = Mean)) + geom_ribbon(aes(ymax = Upper, ymin = Lower), alpha = 0.5, fill = "blue")+ geom_line() + theme_classic() + ylab("Functional Diversity") + xlab("Bin") + scale_x_continuous(breaks = unique(Rich$NumericBin))

B-A