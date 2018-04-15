###Experimental analyses of the results from the second set of datasets
source("load.R")
options(scipen = 999)
#function to normalize between 0 and 1
# Adding ... will allow you to pass through na.rm = T if you want to omit missing values from the calculation (they will still be present in the results):
norm1.0 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
##### binding all result tables.
path = "results2"
names <- list.files(path)
# names <- names[-1] #dropping the first name, which is plots folder
l <- length(names)

output <- list()
for (n in 1:l) {
  output[[n]] = data.frame(read.csv(paste0(path,"/",names[n],sep="")))
}
allresults.new.2 <- do.call(rbind, lapply(output, data.frame, stringsAsFactors=FALSE))
setDT(allresults.new.2)

##processing time summary
speed.test.wide <- dcast(subset(allresults.new.2,allresults.new.2$n == 2250000), data ~ method+sample , value.var="ptime", fun.aggregate = mean)
summary(speed.test.wide)


# #plot data
# ggplot(allresults.new.2) +
#   geom_point(aes(x=k_orig,y=sepVal, size=n),shape = 21, colour = "black", fill = "white")+
#   # theme_bw() +
#   theme(panel.grid.major.y = element_line(colour = "gray"),
#         panel.grid.minor.y = element_line(colour = "gray"),
#         panel.grid.major.x = element_line(colour = "gray"),
#         panel.grid.minor.x = element_line(colour = "gray"),
#         # axis.line = element_line(size=1, colour = "black"),
#         # panel.border = element_blank(), panel.background = element_blank(),
#         # plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
#         # text=element_text(family="Tahoma", face = "bold"),
#         axis.text.x=element_text(colour="black", size = 7, face="plain"),
#         axis.text.y=element_text(colour="black", size = 7, face="plain"),
#         legend.position="none") +
#   xlab("number of clusters") +
#   ylab("separation value")
# 
# 

##creating a percent error to compare overall performance
allresults.new.2$err.freq <- round((allresults.new.2$e_freq/allresults.new.2$k_orig),3)
allresults.new.2$err.mean <- round((allresults.new.2$e_mean/allresults.new.2$k_orig),3)
allresults.new.2$method <- as.character(allresults.new.2$method)
allresults.new.2$method <- ifelse(allresults.new.2$method == "BIC_kluster", "bic_kluster",allresults.new.2$method)



allresults.new.2$X <- seq.int(nrow(allresults.new.2))
allresults.new.2$data <- paste0(allresults.new.2$k_orig,"_",allresults.new.2$sepVal)
# allresults.new.2$err.absolute <- abs(allresults.new.2$err)
allresults.new.2$err.absolute.norm.mean <- norm1.0(abs(allresults.new.2$err.mean))
allresults.new.2$err.absolute.norm.freq <- norm1.0(abs(allresults.new.2$err.freq))

# for ranking, I have to create a new value because it ranks stuff by larger numbers
allresults.new.2$err.absolute.norm.rev.mean <- 1-norm1.0(abs(allresults.new.2$err.mean))
allresults.new.2$err.absolute.norm.rev.freq <- 1-norm1.0(abs(allresults.new.2$err.freq))


# write.csv(allresults.new.2, "longresults1.csv")

##making the wide format as well
library(data.table)
setDT(allresults.new.2)
allresults.new.2.wide.mean <- dcast(allresults.new.2, data ~ method, value.var="err.absolute.norm.rev.mean", fun.aggregate = mean)
allresults.new.2.wide.freq <- dcast(allresults.new.2, data ~ method, value.var="err.absolute.norm.rev.freq", fun.aggregate = mean)
colnames(allresults.new.2.wide.mean) <- c("data", "ap_kluster_mean","bic_kluster_mean","cal_kluster_mean","pam_kluster_mean")
colnames(allresults.new.2.wide.freq) <- c("data", "ap_kluster_frq","bic_kluster_frq","cal_kluster_frq","pam_kluster_frq")
allresults.new.2.wide <- merge(allresults.new.2.wide.mean,allresults.new.2.wide.freq, by="data")

##now doing the Nemenyi test to rank the kluster procedures on big data
# the Nemenyi test, which is not a recommended choice in practice
nm <- nemenyiTest(as.matrix(allresults.new.2.wide[,-1]))
nm
nm$diff.matrix
png(filename=paste0(getwd(),"NemenyiBIG.PNG"),
    units="in",
    width=8,
    height=5,
    pointsize=10,
    res=300)
nemenyi <- plotCD(results.matrix = as.matrix(allresults.new.2.wide[,-1]), alpha = 0.05)
dev.off()


## we then proceed with evaluating kluster procedures for the Friedman post-hoc test with Bergmann and Hommel’s correction
test.res <- postHocTest(data = as.matrix(allresults.new.2.wide[,-1]), test = 'friedman', correct = 'bergmann')
average.ranking <- colMeans(rankMatrix(as.matrix(allresults.new.2.wide[,-1])))
# drawAlgorithmGraph(pvalue.matrix = test.res$corrected.pval,mean.value = average.ranking)
# write.csv(test.res$corrected.pval, paste0(getwd(),"friedman.csv"))


##which one is faster?
# focusing on PAM and BIC
speed.test <- subset(allresults.new.2, allresults.new.2$method %in% c("bic_kluster","pam_kluster"))
# ggplot(speed.test)+
#   geom_boxplot(aes(x=method, y=ptime/60))



wilcox.speed.test <- pairwise.wilcox.test(speed.test$ptime,speed.test$method,paired = T)

setDT(speed.test)
speed.test.wide <- dcast(speed.test, data ~ method, value.var="ptime", fun.aggregate = mean)
# wst <- wilcoxonSignedTest(speed.test.wide$bic_kluster,speed.test.wide$pam_kluster)


speed.test.bic.faster <- subset(speed.test.wide,speed.test.wide$bic_kluster < speed.test.wide$pam_kluster)
wilcox.2 <- wilcox.test(speed.test.bic.faster$bic_kluster,speed.test.bic.faster$pam_kluster,alternative = "two.sided")
library(reshape2)
long1 <- melt(speed.test.bic.faster, id.vars = c("data"))
long1$faster <- "BIC"
wilcox.speed.test.bic <- pairwise.wilcox.test(long1$value,long1$variable,paired = T)

speed.test.pam.faster <- subset(speed.test.wide.mean,speed.test.wide.mean$bic_kluster > speed.test.wide.mean$pam_kluster)
wilcox.3 <- wilcox.test(speed.test.pam.faster$bic_kluster,speed.test.pam.faster$pam_kluster,alternative = "two.sided")
long2 <- melt(speed.test.pam.faster, id.vars = c("data"))
long2$faster <- "PAM"
wilcox.speed.test.pam <- pairwise.wilcox.test(long2$value,long2$variable,paired = T)

speed.test2 <- rbind(long1,long2)

# ### processing time over size of data
# p7 = ggplot(speed.test)+
#   # geom_boxplot(aes(x=variable, y=as.numeric(value)))+
#   # geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
#   # geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
#   # geom_line(aes(color=as.factor(method)),alpha = 0.7) +
#   # geom_point(size = 1, alpha = 0.4 ) +
#   # geom_boxplot(alpha = 0.7)+
#   geom_jitter(aes(x=n, y=scale(ptime/60),color=method,size=sepVal), fill=alpha = 0.4 ) + 
#   # geom_jitter(aes(x=n, y=scale(ptime/60)))+
#   # geom_boxplot(aes(x=factor(n), y=scale(ptime/60),color=method),  alpha = 0.4 ) +
#   # geom_smooth(aes(x=n, y=scale(ptime/60),color=method))+
#   # scale_y_continuous(limits=c(-2,2))+
#   # scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Procedure", nrow= 1)) + 
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(colour = "gray"),
#         panel.grid.minor.y = element_blank(),
#         axis.line = element_line(size=1, colour = "black"),
#         panel.border = element_blank(), panel.background = element_blank(),
#         plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
#         text=element_text(family="Tahoma", face = "bold"),
#         axis.text.x=element_text(colour="black", size = 7, face="plain"),
#         axis.text.y=element_text(colour="black", size = 7, face="plain"),
#         legend.position="bottom") +
#   xlab("Size of database (number of rows)") +
#   ylab("processing time (minutes)") +
#   facet_wrap(~method, ncol = 2)
# 
# ggsave(filename=paste0(getwd(),"/rev2-bigdataresearch/figures/fig7.PNG"), 
#        plot=p7, dpi = 300, width = 12, height = 7)





#now testing the sampling strategies on both BIC most frequent kluster products
# allresults.new.2.wide.sample.mean <- dcast(subset(allresults.new.2,allresults.new.2$method == "bic_kluster" ), data ~ sample, value.var="err.absolute.norm.rev.mean", fun.aggregate = mean)
# allresults.new.2.wide.sample.freq <- dcast(subset(allresults.new.2,allresults.new.2$method == "bic_kluster" ), data ~ sample, value.var="err.absolute.norm.rev.freq", fun.aggregate = mean)
# # colnames(allresults.new.2.wide.sample.mean) <- c("data", "ap_kluster_mean","bic_kluster_mean","cal_kluster_mean","pam_kluster_mean")
# # colnames(allresults.new.2.wide.sample.freq) <- c("data", "ap_kluster_frq","bic_kluster_frq","cal_kluster_frq","pam_kluster_frq")
# allresults.new.2.wide.sample <- rbind(allresults.new.2.wide.sample.mean,allresults.new.2.wide.sample.freq)

allresults.new.2.wide.sample.freq.BIC <- dcast(subset(allresults.new.2,allresults.new.2$method == "bic_kluster" ), data ~ sample, value.var="err.absolute.norm.rev.freq", fun.aggregate = mean)
##now doing the Nemenyi test to rank the sampling strategies on big data
nm2 <- nemenyiTest(as.matrix(allresults.new.2.wide.sample.freq.BIC[,-1]))
nm2
nm2$diff.matrix
png(filename=paste0(getwd(),"NemenyiBIGSample.PNG"),
    units="in",
    width=7,
    height=3,
    pointsize=10,
    res=300)
nemenyi <- plotCD(results.matrix = as.matrix(allresults.new.2.wide.sample.freq.BIC[,-1]), alpha = 0.05)
dev.off()

## we then proceed with evaluating kluster on BIC sample sizes for the Friedman post-hoc test with Bergmann and Hommel’s correction
test.res <- postHocTest(data = as.matrix(allresults.new.2.wide.sample.freq.BIC[,-1]), test = 'friedman', correct = 'bergmann')
average.ranking <- colMeans(rankMatrix(as.matrix(allresults.new.2.wide.sample.freq.BIC[,-1])))
test.res$corrected.pval
