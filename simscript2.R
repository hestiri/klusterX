## ## ## ## ## ## ## ## ##
# # # this script performs the large data analysis using the kluster function....
source("kluster_dev.R")


# now only kluster on large data. Cool stuff! 
simit2 <- function (n, smpl_iter, iter_simu){
  sepvalue = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7) 
  output2 = list()
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
                          clustSizeEq=n*10000, # total data points is nClust^10
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
    run <- stringi::stri_rand_strings(1,5)
    da2 <- da
    da2$run <- run
    # write.csv(da2, paste0("/data/size_",nrow(da),"_n_",n,"_sepval_",sepValu,".csv"))
    rm(da2)
    
    out100 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 100)$sim
    out100$sample = 100
    out200 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 200)$sim
    out200$sample = 200
    out300 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 300)$sim
    out300$sample = 300
    out400 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 400)$sim
    out400$sample = 400
    out500 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 500)$sim
    out500$sample = 500
    out1000 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 1000)$sim
    out1000$sample = 1000
    out2 = rbind(out100,out200,out300,out400,out500,out1000)
    out2$sepVal = sepValu
    out2$run = run
    out2$sampling_iteration = smpl_iter
    out2$simulation_iteration = iter_simu    
    # jpeg(paste0("plots/size_",nrow(da),"_n_",n,"_sepval_",sepValu,".jpg"))
    plot(da)
    dev.off()
    output2[[j]] <- out2
    rm(da,run,out2)
  }
  result2 <- do.call(rbind, lapply(output2, data.frame, stringsAsFactors=FALSE))
  write.csv(result2, paste0("results2/",n, smpl_iter, iter_simu,".csv"))
}


for (t in 3:15){
  simit2(n=t,smpl_iter = 100,iter_simu = 1)
}


