###Experimental analyses of the results from the first set of datasets
source("load.R")
options(scipen = 999)
#function to normalize between 0 and 1
# Adding ... will allow you to pass through na.rm = T if you want to omit missing values from the calculation (they will still be present in the results):
norm1.0 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
##### binding all result tables.
path = "results1"
names <- list.files(path)
# names <- names[-1] #dropping the first name, which is plots folder
n <- length(names)

output <- list()
N <- length(names)
for (n in 1:N) {
  output[[n]] = data.frame(read.csv(paste0(path,"/",names[n],sep="")))
}
allresults.new <- do.call(rbind, lapply(output, data.frame, stringsAsFactors=FALSE))

##creating a percent error to compare overall performance
allresults.new$err <- round((allresults.new$e/allresults.new$k_orig),3)
allresults.new$method <- as.character(allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "BIC.best", "BIC",allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "pamk.best", "PAM",allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "calinski.best", "CAL",allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "apclus.best", "AP",allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "BIC_kluster_frq", "bic_kluster_frq",allresults.new$method)
allresults.new$method <- ifelse(allresults.new$method == "BIC_kluster_mean", "bic_kluster_mean",allresults.new$method)

### adding a new variable to be able to combine the four models
allresults.new$method0 = "BIC"
allresults.new$method0 = ifelse(allresults.new$method %in% c("PAM","pam_kluster_frq","pam_kluster_mean"), "PAM", allresults.new$method0)
allresults.new$method0 = ifelse(allresults.new$method %in% c("CAL","cal_kluster_mean","cal_kluster_frq"), "CAL", allresults.new$method0)
allresults.new$method0 = ifelse(allresults.new$method %in% c("AP","ap_kluster_mean","ap_kluster_frq"), "AP", allresults.new$method0)

##adding a new variable to distinguish time for kluster method
allresults.new$algorithm = "original algorithm"
allresults.new$algorithm = ifelse(allresults.new$method %in% c("ap_kluster_mean","ap_kluster_frq","cal_kluster_mean","cal_kluster_frq",
                                                       "pam_kluster_frq","pam_kluster_mean", "bic_kluster_frq", "bic_kluster_mean"), 
                              "kluster proc",
                              allresults.new$algorithm)

allresults.new$X <- seq.int(nrow(allresults.new))
allresults.new$data <- paste0(allresults.new$k_orig,"_",allresults.new$sepVal)
# allresults.new$err.absolute <- abs(allresults.new$err)
allresults.new$err.absolute.norm <- norm1.0(abs(allresults.new$err))

# for ranking, I have to create a new value because it ranks stuff by larger numbers
allresults.new$err.absolute.norm.rev <- 1-norm1.0(abs(allresults.new$err))

rm(output)
write.csv(allresults.new, "longresults1.csv")

##making the wide format as well
library(data.table)
setDT(allresults.new)
allresults.new.wide <- dcast(allresults.new, data ~ method, value.var="err.absolute.norm.rev")
write.csv(allresults.new.wide, "widesults1_errabsolutenormalrev.csv")
allresults.new.wide.err <- dcast(allresults.new, data ~ method, value.var="err")
write.csv(allresults.new.wide.err, "widesults1_err.csv")

###figure2
### let's look at time of processing for original algorithms
origs = subset(allresults.new,allresults.new$method %in% c("BIC","PAM","CAL","AP"))
## now project the processing time for data up to 2 million rows
p00 = ggplot(origs, aes(x=n, y=ptime/60))+
  stat_smooth(aes(color=as.factor(method0)),method="glm",
              formula = y ~ poly(x, 3),fullrange=TRUE, alpha=0.2)+
  # stat_smooth(aes(color=as.factor(method0)),method="loess",alpha=0.2,
  #             span = .7, fullrange = T, level = 0.95)+
  # geom_point(aes(color=as.factor(method0)),size = 0.6, alpha = 0.8) +
  # geom_point(aes(x=n, y=scale(ptime2)), color = "blue", size = 1, alpha = 0.4 ) +
  scale_x_sqrt(limits=c(0,100000))+
  scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Algorithm", nrow= 1)) + 
  geom_vline(xintercept=3000)+
  # theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "white"),
        #panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size=1, colour = "black"),
        #panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 10, face="plain"),
        axis.text.y=element_text(colour="black", size = 10, face="plain", angle = 90),
        legend.position="bottom") +
  coord_cartesian(xlim=c(600, 100000),ylim=c(0, 500))+
  xlab("SQRT Size of database (number of rows)") +
  ylab("processing time (minutes)")

ggsave(filename=paste0(getwd(),"fig2.PNG"), 
       plot=p00, dpi = 300, width = 8, height = 5)


#figure3
##density plot by method
fill <- "#4271AE"
lines <- "#1F3552"
p3 = ggplot(allresults.new, aes(x =err)) +
  geom_density(colour = "#4271AE", 
               size = 1,alpha = 0.9) + 
  scale_x_continuous(name = "Ratio of error to actual cluster number",
                     breaks = seq(-1, 1, .5),
                     limits=c(-1, 1)) +
  # scale_y_continuous(name = "Density") +
  coord_cartesian(ylim=c(0, 15))+
  # ggtitle("Density plot of cluster number identification error") +
  geom_vline(xintercept = 0.05, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  geom_vline(xintercept = -0.05, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  theme_bw() +
  theme(
    # axis.line = element_line(size=1, colour = "black"),
    #     panel.grid.major = element_line(colour = "#d3d3d3"),
    #     panel.grid.minor = element_blank(),
    #     panel.border = element_blank(), panel.background = element_blank(),
    plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
    text=element_text(family="Tahoma", face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9)) +
  ylab("probability") +
  facet_wrap(~ method , ncol = 3)

ggsave(filename=paste0(getwd(),"fig3.PNG"), 
       plot=p3, dpi = 300, width = 8, height = 4)


###pariwise stuff

library('scmamp')
# The goal is analyzing all the pair-wise comparisons.
# First, we check the differences using the Iman and Davenport omnibus test.
imanDavenportTest(as.matrix(allresults.new.wide.err[,-1]))
# Iman Davenport's correction of Friedman's rank sum test
# 
# data:  as.matrix(allresults.new.wide[, -1])
# Corrected Friedman's chi-squared = 42.185, df1 = 11, df2 = 990, p-value < 0.00000000000000022
# The p-value shown above denotes that there is at least one algorithm that performs differently than
# the rest and, therefore, we can proceed with the post-hoc analysis of the results.


##this wilcox test matrix is great. Compare is with SPSS output and potentially put it as an appendix.
wilcox.1 <- pairwise.wilcox.test(allresults.new$err,allresults.new$method,paired = T)

write.csv(wilcox.1$p.value, paste0(getwd(),"wilcox.csv"))

wilcox.2 <- pairwise.wilcox.test(allresults.new$err,allresults.new$method,paired = T,p.adjust.method = "none")

write.csv(wilcox.2$p.value, paste0(getwd(),"wilcox-nonadjust.csv"))
# wilcox.2 <- wilcox.test(allresults.new.wide$ap_kluster_frq,allresults.new.wide$ap_kluster_mean,alternative = "two.sided")
# wilcox.2$p.value


# the Nemenyi test, which is not a recommended choice in practice
nm <- nemenyiTest(as.matrix(allresults.new.wide[,-1]))
nm
nm$diff.matrix
# nemenyi <- plotCD(results.matrix = as.matrix(allresults.new.wide[,-1]), alpha = 0.05)
png(filename=paste0(getwd(),"fig4-new.PNG"), 
    units="in", 
    width=8, 
    height=5, 
    pointsize=10, 
    res=300)
nemenyi <- plotCD(results.matrix = as.matrix(allresults.new.wide[,-1]), alpha = 0.05)
dev.off()

# In this plot each algorithm is placed on an axis according to its
# average ranking. Then, those algorithms that show no significant differences are grouped together
# using a horizontal line. The plot also shows the size of the critical difference required for considering
# two algorithm as significantly different.


## we then proceed with the top performing algorithms for the Friedman post-hoc test with Bergmann and Hommel’s correction
allresults.new.wide.top <- dplyr::select(allresults.new.wide, data, CAL, BIC, bic_kluster_frq, bic_kluster_mean, pam_kluster_frq,pam_kluster_mean, PAM, AP)
test.res <- postHocTest(data = as.matrix(allresults.new.wide.top[,-1]), test = 'friedman', correct = 'bergmann')
average.ranking <- colMeans(rankMatrix(as.matrix(allresults.new.wide.top[,-1])))
drawAlgorithmGraph(pvalue.matrix = test.res$corrected.pval,mean.value = average.ranking)
write.csv(test.res$corrected.pval, paste0(getwd(),"friedman.csv"))



# res <- postHocTest(data = allresults.new.wide, algorithms = 2:13, test = 'friedman',correct = 'finner' , control = 'data')
# 
# 
# res <- postHocTest(data = allresults.new, group.by = 'method', test = 'wilcoxon',correct = 'finner', control = 'e')


## looking at good and bad results by sepVal and k_orig
agg1 <- aggregate(err ~ method+sepVal, data = allresults.new, mean)
ggplot(agg1,aes(x=as.factor(sepVal),y=err, group=method,color=method)) + geom_path() + geom_point()
agg2 <- aggregate(err ~ method+k_orig, data = allresults.new, mean)

goods = subset(allresults.new, err <= 0.05 & err >= -0.05)
sepVal_good = aggregate(err ~ method+sepVal, data = goods, length)
k_orig_good = aggregate(err ~ method+k_orig, data = goods, length)
ggplot(k_orig_good, aes(x= as.numeric(reorder(k_orig,err)), y= err)) + geom_point(aes(color = method)) + scale_x_continuous(limits=c(0,30))

## looking at good performance by 
### number of clusters
# k_orig_good2 = as.data.frame.matrix(table(goods$method,goods$k_orig))
k_orig_good = data.frame(table(goods$method,goods$k_orig))
k_orig_good$ratio = round(as.numeric(k_orig_good$Freq/7),4)
p1 = ggplot(k_orig_good, aes(Var2, reorder(Var1, ratio))) + 
  geom_tile(aes(fill = ratio),colour = "gray") + 
  scale_fill_gradient(low = "white", high = "darkblue") + 
  labs(fill='Prediction ratio with 95%+ accuracy') + #+ geom_label(label = k_orig_good$ratio)
  xlab("Number of Clusters") +
  ylab("") +
  theme(legend.position="bottom",
        axis.text=element_text(face="bold", size = 12)) +
  ggtitle("")


ggsave(filename=paste0(getwd(),"fig4-1.PNG"), 
       plot=p1, dpi = 300, width = 7.5, height = 7)


### separation value
# table(goods$method,goods$sepVal)
sepVal_good = data.frame(table(goods$method,goods$sepVal))
sepVal_good$ratio = round(as.numeric(sepVal_good$Freq/13),4)
p2 = ggplot(sepVal_good, aes(Var2, reorder(Var1, ratio))) + 
  geom_tile(aes(fill = ratio),colour = "gray") + 
  scale_fill_gradient(limits= c(0,1), low = "white", high = "darkblue") + 
  labs(fill='') + #+ geom_label(label = k_orig_good$ratio)
  xlab("Separation Value") +
  ylab("") +
  theme(legend.position="bottom",
        axis.text=element_text(face="bold", size = 12)) +
  ggtitle("")

ggsave(filename=paste0(getwd(),"fig4-2.PNG"), 
       plot=p2, dpi = 300, width = 5.5, height = 7)



### processing time over size of data
p5 = ggplot(allresults.new, aes(x=n, y=ptime))+
  # geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  geom_line(aes(color=as.factor(algorithm)),alpha = 0.7) +
  geom_point(size = 1, alpha = 0.4 ) +
  # geom_boxplot(alpha = 0.7)+
  geom_point(aes(x=n, y=scale(ptime)), color = "black", size = 1, alpha = 0.4 ) +
  # scale_y_continuous(limits=c(-2,2))+
  scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Procedure", nrow= 1)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 7, face="plain"),
        axis.text.y=element_text(colour="black", size = 7, face="plain"),
        legend.position="bottom") +
  xlab("Size of database (number of rows)") +
  ylab("processing time (seconds)") +
  facet_wrap(~method0, ncol = 2)

ggsave(filename=paste0(getwd(),"fig5.PNG"), 
       plot=p5, dpi = 300, width = 7, height = 7)
