## ## ## ## ## ## ## ## ##
# # # this script performs the initial analysis using the kluster_sim function....

source("kluster_dev.R")

#function to 
simit <- function (n,smpl, smpl_iter, iter_simu){
  #generate and save the data
  sepvalue = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7) 
  output = list()
  for (j in 1:length(sepvalue)){
    sepValu <- sepvalue[j]
    dat <- genRandomClust(numClust=n, 
                          sepVal=sepValu, 
                          numNonNoisy=2, 
                          numNoisy=0, 
                          numOutlier=0, 
                          numReplicate=1, 
                          fileName="test",  
                          clustszind=1, # equal cluster size
                          clustSizeEq=200, # total data points is nClust^10
                          rangeN=c(50,200), 
                          clustSizes=NULL, 
                          covMethod=c("eigen"), 
                          rangeVar=c(1, 10), 
                          lambdaLow=1, 
                          ratioLambda=10,  
                          alphad=1,
                          eta=1,
                          rotateind=TRUE, 
                          iniProjDirMethod=c("SL"), 
                          projDirMethod=c("newton"), 
                          alpha=0.05, 
                          ITMAX=20, 
                          eps=1.0e-10, 
                          quiet=TRUE, 
                          outputDatFlag=FALSE, 
                          outputLogFlag=FALSE, 
                          outputEmpirical=FALSE, 
                          outputInfo=FALSE)$datList[1]
    
    da <- data.frame(dat)
    rm(dat)
    colnames(da) <- c("x", "y")
    run <- stringi::stri_rand_strings(1,7)
    da2 <- da
    da2$run <- run
    # write.csv(da2, paste0("/data/size_",nrow(da),"_n_",n,"_sepval_",sepValu,".csv"))
    rm(da2)
    #store kluster, but need a modified version that also stores each run's results
    out <- kluster_sim(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = smpl)$sim
    out$sepVal = sepValu
    out$run = run
    out$sample = smpl
    out$sampling_iteration = smpl_iter
    out$simulation_iteration = iter_simu    
    # jpeg(paste0("/plots/size_",nrow(da),"_n_",n,"_sepval_",sepValu,".jpg"))
    # plot(da)
    # dev.off()
    output[[j]] <- out
    rm(da,run,out)
  }
  result <- do.call(rbind, lapply(output, data.frame, stringsAsFactors=FALSE))
  write.csv(result, paste0("results1/",n,smpl, smpl_iter, iter_simu,".csv"))
}

for (t in 3:15){
  simit(n=t, smpl = 100,smpl_iter = 100,iter_simu = 25)
}

