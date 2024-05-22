Data <- function(dataname){
  result <- list()
  
  #### Mixte datasets
  if(dataname=='abalone'){
    abalone <- read_csv("Data/mixte/abalone.data", col_names = FALSE)
    result$x <- abalone[, 2:ncol(abalone)]
    result$y <- as.numeric(as.factor(abalone$X1))
    result$name <- 'abalone'
    rm(abalone)
  }
  
  if(dataname=='ecoli'){
    ecoli <- read_table("Data/mixte/ecoli.data", col_names = FALSE)
    result$x <- ecoli[,2:(ncol(ecoli)-1)]
    result$y <- as.numeric(as.factor(ecoli$X9))
    result$name <- 'ecoli'
    rm(ecoli)
  }
  
  if(dataname=='wine'){
    wine <- read_csv("Data/mixte/wine.data", col_names = FALSE)
    result$x <- wine[,2:ncol(wine)]
    result$y <- as.numeric(as.factor(wine$X1))
    result$name <- 'wine'
    rm(wine)
  }
  
  if(dataname=='thyroid'){
    data("thyroid")
    result$x <- thyroid[,2:ncol(thyroid)]
    result$y <- as.numeric(as.factor(thyroid$Diagnosis))
    result$name <- 'thyroid'
  }
  
  #### multi-view
  if(dataname=='prokaryoti'){
    prokaryoti <- R.matlab::readMat('Data/multi-view/prokaryotic.mat')
    result$x <- list(prokaryoti$text,prokaryoti$proteome.comp,
                         prokaryoti$gene.repert) 
    result$y <- as.numeric(as.factor(prokaryoti$truth))
    result$name <- 'prokaryoti'
    rm(prokaryoti)
  }
  
  if(dataname=='mfeat'){
    mfeat_fac <- read_table("Data/multi-view/mfeat/mfeat-fac", col_names = FALSE)
    mfeat_fou <- read_table("Data/multi-view/mfeat/mfeat-fou", col_names = FALSE)
    mfeat_kar <- read_table("Data/multi-view/mfeat/mfeat-kar", col_names = FALSE)
    mfeat_y <- rep(1:10, each=200)
    
    result$x <- list(mfeat_fac, mfeat_fou, mfeat_kar)
    result$y <- mfeat_y
    result$name <- 'mfeat'
    rm(mfeat_fac, mfeat_fou, mfeat_kar, mfeat_y)
  }
  
  #### categorical datasets
  if(dataname=='car'){
    car <- read.csv("Data/categorical/car.data", header=FALSE)
    result$x <- car[,1:(ncol(car)-1)]
    result$y <- as.numeric(as.factor(car$V7))
    result$name <- 'car'
    rm(car)
  }
  
  if(dataname=='cancer'){
    cancer <- read.csv("Data/categorical/breast-cancer.data", header=FALSE)
    result$x <- cancer[,2:ncol(cancer)]
    result$y <- as.numeric(as.factor(cancer$V1))
    result$name <- 'cancer'
    rm(cancer)
  }
  
  if(dataname=='zoo'){
    zoo <- read.csv("Data/categorical/zoo.data", header=FALSE)
    result$x <- zoo[,setdiff(2:(ncol(zoo)-1), 14)]
    result$y <- as.numeric(as.factor(zoo$V18))
    result$name <- 'zoo'
    rm(zoo)
  }
  
  
  ### Time serie
  if(dataname=='biofam'){
    load("Data/time-series/biofam3c.rda")
    children_long <- transform_table(biofam3c$children, "children")
    married_long <- transform_table(biofam3c$married, "married")
    left_long <- transform_table(biofam3c$left, "left")
    biofam <- children_long %>%
      full_join(married_long, by = c("id", "key")) %>%
      full_join(left_long, by = c("id", "key"))
    biofam$time <- as.numeric(as.factor(biofam$key))
    
    result$x <- biofam
    result$y <- as.numeric(as.factor(biofam3c$covariates$sex))
    result$name <- 'biofam'
    rm(children_long,married_long,left_long,biofam3c)
  }
  
  if(dataname=='UWaveGestureLibrary'){
    UWaveGestureLibrary <- read_excel("Data/time-series/UWaveGestureLibrary.xlsx")
    result$x <- UWaveGestureLibrary
    df_unique <- UWaveGestureLibrary %>% distinct(id, label, .keep_all = FALSE)
    result$y <- as.numeric(as.factor(df_unique$label))
    result$name <- 'UWaveGestureLibrary'
    rm(UWaveGestureLibrary)
  }
  
  if(dataname=='StandWalkJump'){
    StandWalkJump <- read_excel("Data/time-series/StandWalkJump.xlsx")
    result$x <- StandWalkJump
    df_unique <- StandWalkJump %>% distinct(id, label, .keep_all = FALSE)
    result$y <- as.numeric(as.factor(df_unique$label))
    result$name <- 'StandWalkJump'
    rm(StandWalkJump)
  }
  
  return(result)
}

