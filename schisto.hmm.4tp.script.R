#### Schisto Script Paper 1 ####
# 4 time points with random walk
# 250 children 
# KK new data and CCA new data files


require(plyr)
require(tidyverse)
require(runjags)
require(fitdistrplus)
require(coda)
require(stats)
require(reshape2)
require(zoo)
require(egg)
require(ggridges)

#### Schisto modelling ####
source("schisto.hmm.all.funcs.R")

# load data normally just to look at if necessary #
kk.data <- read.csv("KK.4tp.csv")

cca.data <- read.csv("cca.clean.csv")

#### Load data ####

kk <- loadKKdata("KK.4tp.csv") 

KKIDs <- getKKChildIDs("KK.4tp.csv") # Get KKIDs

cca <- loadCCAdata("cca.clean.csv") # load and format CCA data
CCAIDs <- getCCAChildIDs("cca.clean.csv") # Get CCAIDs

CID <- intersect(KKIDs,CCAIDs) 
kk <- kk[match(CID,KKIDs),,]
cca <- cca[match(CID,CCAIDs),]

ccagscore <- loadgscoredata("cca.clean.csv")
gscoreCCAIDS <- getCCAChildIDsgscore("cca.clean.csv") 

gscoreCID <- intersect(KKIDs,gscoreCCAIDS) 
ccagscore <- ccagscore[match(CID,gscoreCID),]
#### KK only ####

results.kk <- do.KK.run(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk) # run jags model

prev.kk <- density(c(results.kk$mcmc[[1]][,"prev"],results.kk$mcmc[[2]][,"prev"]))

clear2KK <- density(c(results.kk$mcmc[[1]][,"clearance[2]"],results.kk$mcmc[[2]][,"clearance[2]"]))
clear3KK <- density(c(results.kk$mcmc[[1]][,"clearance[3]"],results.kk$mcmc[[2]][,"clearance[3]"]))
clear4KK <- density(c(results.kk$mcmc[[1]][,"clearance[4]"],results.kk$mcmc[[2]][,"clearance[4]"]))

reinf3KK <- density(c(results.kk$mcmc[[1]][,"reinfec[3]"],results.kk$mcmc[[2]][,"reinfec[3]"]))
reinf4KK <- density(c(results.kk$mcmc[[1]][,"reinfec[4]"],results.kk$mcmc[[2]][,"reinfec[4]"]))

rtnb.kk <- density(c(results.kk$mcmc[[1]][,"rtnb"], results.kk$mcmc[[2]][,"rtnb"]))

sh.kk <- density(c(results.kk$mcmc[[1]][,"sh"], results.kk$mcmc[[2]][,"sh"]))

rt.kk <- density(c(results.kk$mcmc[[1]][,"rt"], results.kk$mcmc[[2]][,"rt"]))

status.output <- as.data.frame(as.matrix(as.mcmc(results.kk))) 

variables <- status.output[,1:10]

status.output <- status.output[,-c(1:10)]
  
kk.time <- time.steps(status.output)
kk.time$model <- as.factor("Kato-Katz")

#### CCA and KK ####

results.both <- do.run(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk,cca) 

prev.both <- density(c(results.both$mcmc[[1]][,"prev"],results.both$mcmc[[2]][,"prev"]))

clear2.both <- density(c(results.both$mcmc[[1]][,"clearance[2]"],results.both$mcmc[[2]][,"clearance[2]"]))
clear3.both <- density(c(results.both$mcmc[[1]][,"clearance[3]"],results.both$mcmc[[2]][,"clearance[3]"]))
clear4.both <- density(c(results.both$mcmc[[1]][,"clearance[4]"],results.both$mcmc[[2]][,"clearance[4]"]))

reinf3.both <- density(c(results.both$mcmc[[1]][,"reinfec[3]"],results.both$mcmc[[2]][,"reinfec[3]"]))
reinf4.both <- density(c(results.both$mcmc[[1]][,"reinfec[4]"],results.both$mcmc[[2]][,"reinfec[4]"]))

#rtnb.both <- density(c(results.both$mcmc[[1]][,"rtnb"], results.both$mcmc[[2]][,"rtnb"]))

#sh.both <- density(c(results.both$mcmc[[1]][,"sh"], results.both$mcmc[[2]][,"sh"]))

#rt.both <- density(c(results.both$mcmc[[1]][,"rt"], results.both$mcmc[[2]][,"rt"]))

#tpp <- density(c(results.both$mcmc[[1]][,"TracPosProb"], results.both$mcmc[[2]][,"TracPosProb"]))

intercept.est <- density(c(results.both$mcmc[[1]][,"intercept"], results.both$mcmc[[2]][,"intercept"]))

k.est <- density(c(results.both$mcmc[[1]][,"k"], results.both$mcmc[[2]][,"k"]))

status.both <- as.data.frame(as.matrix(as.mcmc(results.both))) 

varsKKCCA <- status.both[,1:13]

status.both <- status.both[,-c(1:13)]

both.time <- time.steps(status.both)
both.time$model <- as.factor("Plus(+)")

#### KK and GScore ####

gscore.KK <- do.10.run(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk,ccagscore) 
prev.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"prev"],gscore.KK$mcmc[[2]][,"prev"]))

clear2.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"clearance[2]"], gscore.KK$mcmc[[2]][,"clearance[2]"]))
clear3.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"clearance[3]"], gscore.KK$mcmc[[2]][,"clearance[3]"]))
clear4.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"clearance[4]"], gscore.KK$mcmc[[2]][,"clearance[4]"]))

reinf3.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"reinfec[3]"], gscore.KK$mcmc[[2]][,"reinfec[3]"]))
reinf4.gscore.KK <- density(c( gscore.KK$mcmc[[1]][,"reinfec[4]"], gscore.KK$mcmc[[2]][,"reinfec[4]"]))

status.gscore.KK <- as.data.frame(as.matrix(as.mcmc(gscore.KK)))

varskkgscore <- status.gscore.KK[,1:13]

status.gscore.KK <- status.gscore.KK[,-c(1:13)]

gscore.time <- time.steps(status.gscore.KK)

intercept.kkgscore <- density(c(gscore.KK$mcmc[[1]][,"intercept"], gscore.KK$mcmc[[2]][,"intercept"]))

k.gscorekk <- density(c(gscore.KK$mcmc[[1]][,"k"], gscore.KK$mcmc[[2]][,"k"]))
gscore.time$model <- as.factor("G-Score")


#### FIGURES ####

# g score model intercept 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
pdf("intercept.kkgscore8.6.20.unif.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(intercept.kkgscore$x,intercept.kkgscore$y/max(intercept.kkgscore$y),xlim=c(0,100),ylim=c(0,1),type = "l", ylab = "",xlab="", col="#1f78b4")
#polygon(c(rev(intercept.kkgscore$x), intercept.kkgscore$x), c(rev(intercept.kkgscore$y/max(intercept.kkgscore$y)), rep(0,length(intercept.kkgscore$y))),
        #col = adjustcolor("#b2df8a", alpha=.2), border = NA)
mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("logistic function intercept ",side=1,cex=1,line=1.2)
dev.off()
graphics.off()

# logistic function plots 

logcurve <- function(x,k,intercept){
  y <- 1 / (1 + exp(-k*(x-intercept)))
  return(y)
}


x=seq(from=0, to=100, by =1)
k=mean(k.est$x)
intercept=mean(intercept.est$x)
ykkcca=logcurve(x,k,intercept)

k=mean(k.gscorekk$x)
intercept=mean(intercept.kkgscore$x)
ykkgscore=logcurve(x,k,intercept)

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,100), ylim = c(0,1), xlab="KK Count", ylab="f(x)", cex.axis=1)
lines(ykkcca, col="orange", lwd=2, lty = 2) # high dark red
lines(ykkgscore, col="darkred", lwd=1) # high dark red
add_legend("top", legend=c(c("KK & CCA (normal)", "KK & G-Score")), lty=1, lwd=2,
           col=c("orange","darkred"),
           horiz=TRUE, bty='n', cex=1)
dev.copy(pdf, "logcurve.11.6.20.calcprec.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# Number infected at each time point #

# need to take the outputs and make them into a dataframe with counts, can turn into props if need be later

# just KK model
kk.tp.counts <- as.data.frame(matrix(nrow=nrow(variables), ncol = 15))
colnames(kk.tp.counts) <- c("starting_number", "cl_BL_3wk", "No_Cl_BL_3wk", "total_inf_at3wk", "cl_3wk_9wk", "reinf_3wk_9wk", "No_cl_3to9", "No_reinf_3to9",
                            "total_inf_at9wks", "No_still_inf_frmBL", "Cl_9w_6m", "reinf_9w_6m", "No_cl_9wkto6m", "No_reinf_9wkto6m", "total_inf_ar6m")
kk.tp.counts[,1] <- 210*variables$prev
kk.tp.counts[,2] <- variables[,2]
kk.tp.counts[,3] <- kk.tp.counts[,1]*kk.tp.counts[,2]
kk.tp.counts[,4] <- kk.tp.counts[,1]-kk.tp.counts[,3]
kk.tp.counts[,5] <- variables[,3]
kk.tp.counts[,6] <- variables[,6]
kk.tp.counts[,7] <- kk.tp.counts[,4]*kk.tp.counts[,5]
kk.tp.counts[,8] <- kk.tp.counts[,3]*kk.tp.counts[,6]
kk.tp.counts[,9] <- (kk.tp.counts[,4]-kk.tp.counts[,7])+kk.tp.counts[,8]
kk.tp.counts[,10] <- kk.tp.counts[,4]-kk.tp.counts[,7]
kk.tp.counts[,11] <- variables[,4]
kk.tp.counts[,12] <- variables[,7]
kk.tp.counts[,13] <- kk.tp.counts[,10]*kk.tp.counts[,11]
kk.tp.counts[,14] <- kk.tp.counts[,7]*kk.tp.counts[,12]
kk.tp.counts[,15] <- (kk.tp.counts[,9]-kk.tp.counts[,13])+kk.tp.counts[,14]

kk.tp.inf <- as.data.frame(cbind(kk.tp.counts$starting_number, kk.tp.counts$total_inf_at3wk, kk.tp.counts$total_inf_at9wks, kk.tp.counts$total_inf_ar6m))
colnames(kk.tp.inf) <- c("pre-treatment", "3 weeks", "9 weeks", "6 months")
kk.tp.inf <- melt(kk.tp.inf)
colnames(kk.tp.inf) <- c("time", "number_infected")
kk.tp.inf$model <- as.factor("Kato-Katz")
kk.tp.inf$time <- as.factor(kk.tp.inf$time)

kkclearNo <- kk.tp.counts %>%
  dplyr::select(No_Cl_BL_3wk, No_cl_3to9, No_cl_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz", dynamic="Clearance")

kkreinfNo <- kk.tp.counts %>%
  dplyr::select(No_reinf_3to9, No_reinf_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz", dynamic="Reinfection")


# CCA and KK model 
ccakk.tp.counts <- as.data.frame(matrix(nrow=nrow(varsKKCCA), ncol = 15))
colnames(ccakk.tp.counts) <- c("starting_number", "cl_BL_3wk", "No_Cl_BL_3wk", "total_inf_at3wk", "cl_3wk_9wk", "reinf_3wk_9wk", "No_cl_3to9", "No_reinf_3to9",
                            "total_inf_at9wks", "No_still_inf_frmBL", "Cl_9w_6m", "reinf_9w_6m", "No_cl_9wkto6m", "No_reinf_9wkto6m", "total_inf_ar6m")
ccakk.tp.counts[,1] <- 210*varsKKCCA$prev
ccakk.tp.counts[,2] <- varsKKCCA[,2]
ccakk.tp.counts[,3] <- ccakk.tp.counts[,1]*ccakk.tp.counts[,2]
ccakk.tp.counts[,4] <- ccakk.tp.counts[,1]-ccakk.tp.counts[,3]
ccakk.tp.counts[,5] <- varsKKCCA[,3]
ccakk.tp.counts[,6] <- varsKKCCA[,6]
ccakk.tp.counts[,7] <- ccakk.tp.counts[,4]*ccakk.tp.counts[,5]
ccakk.tp.counts[,8] <- ccakk.tp.counts[,3]*ccakk.tp.counts[,6]
ccakk.tp.counts[,9] <- (ccakk.tp.counts[,4]-ccakk.tp.counts[,7])+ccakk.tp.counts[,8]
ccakk.tp.counts[,10] <- ccakk.tp.counts[,4]-ccakk.tp.counts[,7]
ccakk.tp.counts[,11] <- varsKKCCA[,4]
ccakk.tp.counts[,12] <- varsKKCCA[,7]
ccakk.tp.counts[,13] <- ccakk.tp.counts[,10]*ccakk.tp.counts[,11]
ccakk.tp.counts[,14] <- ccakk.tp.counts[,7]*ccakk.tp.counts[,12]
ccakk.tp.counts[,15] <- (ccakk.tp.counts[,9]-ccakk.tp.counts[,13])+ccakk.tp.counts[,14]

ccakk.tp.inf <- as.data.frame(cbind(ccakk.tp.counts$starting_number, ccakk.tp.counts$total_inf_at3wk, ccakk.tp.counts$total_inf_at9wks, ccakk.tp.counts$total_inf_ar6m))
colnames(ccakk.tp.inf) <- c("pre-treatment" ,"3 weeks", "9 weeks", "6 months")
ccakk.tp.inf <- melt(ccakk.tp.inf)
colnames(ccakk.tp.inf) <- c("time", "number_infected")
ccakk.tp.inf$model <- as.factor("CCA Plus (+)")
ccakk.tp.inf$time <- as.factor(ccakk.tp.inf$time)

kkccaclearNo <- ccakk.tp.counts %>%
  dplyr::select(No_Cl_BL_3wk, No_cl_3to9, No_cl_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz & CCA +", dynamic="Clearance")

kkccareinfNo <- ccakk.tp.counts %>%
  dplyr::select(No_reinf_3to9, No_reinf_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz & CCA +", dynamic="Reinfection")

# Gscore and KK model 

gscorekk.tp.counts <- as.data.frame(matrix(nrow=nrow(varskkgscore), ncol = 15))
colnames(gscorekk.tp.counts) <- c("starting_number", "cl_BL_3wk", "No_Cl_BL_3wk", "total_inf_at3wk", "cl_3wk_9wk", "reinf_3wk_9wk", "No_cl_3to9", "No_reinf_3to9",
                               "total_inf_at9wks", "No_still_inf_frmBL", "Cl_9w_6m", "reinf_9w_6m", "No_cl_9wkto6m", "No_reinf_9wkto6m", "total_inf_ar6m")
gscorekk.tp.counts[,1] <- 210*varskkgscore$prev
gscorekk.tp.counts[,2] <- varskkgscore[,2]
gscorekk.tp.counts[,3] <- gscorekk.tp.counts[,1]*gscorekk.tp.counts[,2]
gscorekk.tp.counts[,4] <- gscorekk.tp.counts[,1]-gscorekk.tp.counts[,3]
gscorekk.tp.counts[,5] <- varskkgscore[,3]
gscorekk.tp.counts[,6] <- varskkgscore[,6]
gscorekk.tp.counts[,7] <- gscorekk.tp.counts[,4]*gscorekk.tp.counts[,5]
gscorekk.tp.counts[,8] <- gscorekk.tp.counts[,3]*gscorekk.tp.counts[,6]
gscorekk.tp.counts[,9] <- (gscorekk.tp.counts[,4]-gscorekk.tp.counts[,7])+gscorekk.tp.counts[,8]
gscorekk.tp.counts[,10] <- gscorekk.tp.counts[,4]-gscorekk.tp.counts[,7]
gscorekk.tp.counts[,11] <- varskkgscore[,4]
gscorekk.tp.counts[,12] <- varskkgscore[,7]
gscorekk.tp.counts[,13] <- gscorekk.tp.counts[,10]*gscorekk.tp.counts[,11]
gscorekk.tp.counts[,14] <- gscorekk.tp.counts[,7]*gscorekk.tp.counts[,12]
gscorekk.tp.counts[,15] <- (gscorekk.tp.counts[,9]-gscorekk.tp.counts[,13])+gscorekk.tp.counts[,14]

gscorekk.tp.inf <- as.data.frame(cbind(gscorekk.tp.counts$starting_number , gscorekk.tp.counts$total_inf_at3wk, gscorekk.tp.counts$total_inf_at9wks, gscorekk.tp.counts$total_inf_ar6m))
colnames(gscorekk.tp.inf) <- c("pre-treatment" ,"3 weeks", "9 weeks", "6 months")
gscorekk.tp.inf <- melt(gscorekk.tp.inf)
colnames(gscorekk.tp.inf) <- c("time", "number_infected")
gscorekk.tp.inf$model <- as.factor("CCA G-Score")
gscorekk.tp.inf$time <- as.factor(gscorekk.tp.inf$time)

kkgscoreclearNo <- gscorekk.tp.counts %>%
  dplyr::select(No_Cl_BL_3wk, No_cl_3to9, No_cl_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz & G-Score", dynamic="Clearance")

kkgscorereinfNo <- gscorekk.tp.counts %>%
  dplyr::select(No_reinf_3to9, No_reinf_9wkto6m)%>%
  mutate(total=rowSums(.), model="Kato-Katz & G-Score", dynamic="Reinfection")

clearnums <- bind_rows(kkgscoreclearNo, kkccaclearNo, kkclearNo)%>%
  mutate_if(is.character, as.factor)

clearnums2 <- clearnums
colnames(clearnums2) <- c("Blto3wk", "3to9wk", "9wkto6mth", "total", "model", "dynamic" )
clearnums2 <- clearnums2 %>%
  pivot_longer(cols=c("Blto3wk", "3to9wk", "9wkto6mth", "total"),names_to = c("time"), values_to="count")


reinfnums <- bind_rows(kkgscorereinfNo, kkccareinfNo, kkreinfNo)
reinfnums2 <- reinfnums
colnames(reinfnums2) <- c("3to9wk", "9wkto6mth", "total", "model", "dynamic" )
reinfnums2 <- reinfnums2 %>% 
  pivot_longer(cols=c("3to9wk", "9wkto6mth","total"),names_to = c("time"), values_to="count")

crfacetdata <- bind_rows(reinfnums2, clearnums2)%>%mutate(time=factor(time, levels=c("9wkto6mth", "3to9wk", "Blto3wk","total")))

total_clear_reinf <- ggplot(crfacetdata, aes(x = count, y = time)) +
       geom_density_ridges(aes(fill=model)) +
       coord_cartesian(clip = "off")+
       facet_grid(dynamic~model,drop = TRUE, scales = "free_y")+
       theme_bw()+
       theme(axis.title.x = element_text(size = 16),
             axis.text.x = element_text(size = 14),
             axis.text.y = element_text(size = 14),
             axis.title.y = element_text(size = 16), 
             strip.text = element_text(size=12))+
       ylab("Time") + xlab("Number of Children")+
       scale_y_discrete(expand=expansion(add = c(0.2, 1.3)))+
       scale_fill_manual(values=c("#1b9e77","#7570b3","#d95f02"))

total_clear_reinf <- tag_facet2(total_clear_reinf)+ggsave("total_clear_reinf.pdf")

clear.SE<-as.data.frame(as.list(aggregate(total ~ model, data = clearnums, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))

clearSE <- summarySE(clearnums, measurevar = "total", groupvars = "model")%>%
  mutate(dynamic="Clearance")
reinfSE <- summarySE(reinfnums, measurevar = "total", groupvars = "model")%>%
  mutate(dynamic="Reinfection")

crmeans<- bind_rows(clearSE, reinfSE)%>%
  mutate_if(is.character, as.factor)

clear <- clearnums %>% dplyr::select(total, model, dynamic) 
reinf <- reinfnums %>% dplyr::select(total, model, dynamic) 

crnums <- bind_rows(clear, reinf)

# the proportion infected at each time point figure #

kkSE <- summarySE(kk.tp.inf, measurevar = "number_infected", groupvars = c("time", "model"))
kkccaSE <- summarySE(ccakk.tp.inf, measurevar = "number_infected", groupvars = c("time", "model"))
kkgscoreSE <- summarySE(gscorekk.tp.inf, measurevar = "number_infected", groupvars = c("time", "model"))

#00a8cc
# yellow = kkgscore
# green = kkCCA
# grey = KK
count.plot <- ggplot()+
  geom_point(data=kkSE, aes(x=time, y=mean/210), size=6, colour="#1b9e77")+
  geom_line(data=kkSE, aes(x=time, y=mean/210, group=model),  colour="#1b9e77")+
  geom_jitter(data=kk.tp.inf, aes(x=time, y=number_infected/210), alpha=0.05, size=0.01, colour="#1b9e77", position = position_jitter(width = .1))+
  geom_point(data=kkgscoreSE, aes(x=time, y=mean/210), size=4,  colour="#d95f02")+
  geom_line(data=kkgscoreSE, aes(x=time, y=mean/210, group=model),  colour="#d95f02")+
  geom_jitter(data=gscorekk.tp.inf, aes(x=time, y=number_infected/210), alpha=0.05, size=0.01, colour="#d95f02", position = position_jitter(width = .1))+
  geom_point(data=kkccaSE, aes(x=time, y=mean/210), size=4, colour="#7570b3")+
  geom_line(data=kkccaSE, aes(x=time, y=mean/210, group=model),  colour="#7570b3")+
  geom_jitter(data=ccakk.tp.inf, aes(x=time, y=number_infected/210), alpha=0.05, size=0.01, colour="#7570b3", position = position_jitter(width = .1))+
  ylab("Proportion of Cohort")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))+
  ggsave("count.plot.prop.20.7.20.png")

# the total numnber cleared and reinfected by each model figure 
# 95% confidence interval 
total.mean.plot <- ggplot()+
  geom_point(data=crnums, aes(x=model, y=mean, colour=dynamic), size=3)+
  geom_errorbar(data=crnums, aes(x=model, ymin=mean-1.96*se, ymax=mean+1.96*se))+
  ylab("Number of Children")
  theme_bw()
  
splitcounts <- ggplot(crnums, aes(x=model, y=total, fill = dynamic)) +
    geom_split_violin(draw_quantiles = c(0.025, 0.5, 0.975))+
    geom_point()
  theme_bw()+
  scale_fill_manual(values=c("#1f78b4", "#b2df8a"))+
  ggsave("clearreinfcounts.pdf")
  
  quantile(ccakkclearNo$total/ ccakk.tp.counts$starting_number, probs=c(0.025,0.5, 0.975))
  quantile(kkgscoreclearNo$total/ gscorekk.tp.counts$starting_number, probs=c(0.025,0.5, 0.975))
  quantile(kkccaclearNo$total/ kkcca.tp.counts$starting_number, probs=c(0.025,0.5, 0.975))
# Calculating the proportion cleared by reinfected by 9wks. 

kkcca9wkcr <- ccakk.tp.counts %>% dplyr::select(No_reinf_3to9, No_Cl_BL_3wk, No_cl_3to9, Cl_9w_6m)%>%
  mutate(totalclear=No_Cl_BL_3wk+ No_cl_3to9+ Cl_9w_6m,
         prop=No_reinf_3to9/totalclear,
         percent=prop*100)
confidence_interval(kkcca9wkcr$percent, 0.95)
median(kkcca9wkcr$percent)
sd(kkcca9wkcr$percent)/sqrt(nrow(kkcca9wkcr))

kkgscore9wkcr <- gscorekk.tp.counts %>% dplyr::select(No_reinf_3to9, No_Cl_BL_3wk, No_cl_3to9, Cl_9w_6m)%>%
  mutate(totalclear=No_Cl_BL_3wk+ No_cl_3to9+ Cl_9w_6m,
         prop=No_reinf_3to9/totalclear,
         percent=prop*100)
confidence_interval(kkgscore9wkcr$prop, 0.95)
median(kkgscore9wkcr$percent)
var(kkgscore9wkcr$percent)

kk9wkcr <- kk.tp.counts %>% dplyr::select(No_reinf_3to9, No_Cl_BL_3wk, No_cl_3to9, Cl_9w_6m)%>%
  mutate(totalclear=No_Cl_BL_3wk+ No_cl_3to9+ Cl_9w_6m,
         prop=No_reinf_3to9/totalclear,
         percent=prop*100)
confidence_interval(kk9wkcr$prop, 0.95)
median(kk9wkcr$percent)
sd(kk9wkcr$percent)/sqrt(nrow(kk9wkcr))

# Probability infected in each time point #

# This is a df for the geom_text with the co-ordinates of the labels 

anno <- data.frame(c("Baseline", "ThreeWeeks", "NineWeeks", "SixMonths"))
colnames(anno) <- "time"
anno$xcoord <- 0
anno$ycoord <-200
anno$label <- as.factor(c("A", "B", "C", "D"))

pdf("histogram.KK.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.kk <- qplot(value, data = kk.time, geom = "histogram", binwidth=0.02) +
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected (Kato-Katz)")

prob.inf.kk
dev.off()

pdf("histogram.KKCCA.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.both <- qplot(value, data = both.time, geom = "histogram", binwidth=0.02) +
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected (Plus (+))")

prob.inf.both
dev.off()

pdf("histogram.KKgscore.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.gscore <- qplot(value, data = gscore.time , geom = "histogram",  binwidth=0.02) +
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected(G-Score)")

prob.inf.gscore
dev.off()

# Posterior distributions from the models with no reinfection between BL - 3wk #

#prev
dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Prevalence", ylab="Scaled Density", cex.axis=1)
lines(prev.kk$x,prev.kk$y/max(prev.kk$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(prev.kk$x), prev.kk$x), c(rev(prev.kk$y/max(prev.kk$y)), rep(0,length(prev.kk$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(prev.gscore.KK$x,prev.gscore.KK$y/max(prev.gscore.KK$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(prev.gscore.KK$x), prev.gscore.KK$x), c(rev(prev.gscore.KK$y/max(prev.gscore.KK$y)), rep(0,length(prev.gscore.KK$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(prev.both$x,prev.both$y/max(prev.both$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(prev.both$x), prev.both$x), c(rev(prev.both$y/max(prev.both$y)), rep(0,length(prev.both$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topleft", c("Kato-Katz", "KK & G-Score", "KK & Plus (+)"),
       col=c("#1b9e77","#d95f02", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "prevalence.20.7.20.pdf", height = 6, width = 6)
dev.off()
graphics.off()


# clearance 
pdf("clearance.20.7.20.pdf",width=12,height=4)

par(font=2, cex.axis=0.75, lwd=2, mar=c(3,2.5,1.2,0)+0.1,mgp=c(3,0.4,0), bty = "n")
par(mfrow=c(1,3))

plot(clear2.gscore.KK$x, clear2.gscore.KK$y/max(clear2.gscore.KK$y), xlim = c(0,1), ylim=c(0,1), type= "l", ylab="", xlab="", col = "#d95f02")
lines(clear2.both$x, clear2.both$y/max(clear2.both$y), col = "#7570b3")
lines(clear2KK$x, clear2KK$y/max(clear2KK$y), col = "#1b9e77")
mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Clearance Baseline to 3 Weeks",side=1,cex=0.8,line=1.8)


plot(clear3.gscore.KK$x, clear3.gscore.KK$y/max(clear3.gscore.KK$y), xlim = c(0,1), ylim=c(0,1), type= "l", ylab="", xlab="", col = "#d95f02")
lines(clear3.both$x, clear3.both$y/max(clear3.both$y), col = "#7570b3")
lines(clear3KK$x, clear3KK$y/max(clear3KK$y), col = "#1b9e77", lty=1)
mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Clearance 3 Weeks to 9 Weeks",side=1,cex=0.8,line=1.8)


plot(clear4.gscore.KK$x, clear4.gscore.KK$y/max(clear4.gscore.KK$y), xlim = c(0,1), ylim=c(0,1), type= "l", ylab="", xlab="", col = "#d95f02")
lines(clear4.both$x, clear4.both$y/max(clear4.both$y), col = "#7570b3")
lines(clear4KK$x, clear4KK$y/max(clear4KK$y), col = "#1b9e77")

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Clearance 9 Weeks to 6 Months",side=1,cex=0.8,line=1.8)

legend("topright", c("Kato-Katz", "KK & G-Score", "KK & Plus (+)"),
       col=c("#1b9e77","#d95f02", "#7570b3"), 
       lty=c(1), bty='n',cex=0.75)
dev.off()

# reinfection 
pdf("reinfection.20.7.20.pdf",width=8,height=4)

par(font=2, cex.axis=0.75, lwd=2, mar=c(3,2.5,1.2,0)+0.1,mgp=c(3,0.4,0), bty = "n")
par(mfrow=c(1,2))

plot(reinf3.gscore.KK$x, reinf3.gscore.KK$y/max(reinf3.gscore.KK$y), xlim = c(0,1), ylim=c(0,1), type= "l", ylab="", xlab="", col = "#d95f02")
lines(reinf3.both$x, reinf3.both$y/max(reinf3.both$y), col = "#7570b3")
lines(reinf3KK$x, reinf3KK$y/max(reinf3KK$y), col = "#1b9e77", lty=1)
mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Reinfection 3 Weeks to 9 Weeks",side=1,cex=0.8,line=1.8)


plot(reinf4.gscore.KK$x, reinf4.gscore.KK$y/max(reinf4.gscore.KK$y), xlim = c(0,1), ylim=c(0,1), type= "l", ylab="", xlab="", col = "#d95f02")
lines(reinf4.both$x, reinf4.both$y/max(reinf4.both$y), col = "#7570b3")
lines(reinf4KK$x, reinf4KK$y/max(reinf4KK$y), col = "#1b9e77")

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Reinfection 9 Weeks to 6 Months",side=1,cex=0.8,line=1.8)

legend("topleft", c("Kato-Katz", "KK & G-Score", "KK & Plus (+)"),
       col=c("#1b9e77","#d95f02", "#7570b3"), 
       lty=c(1), bty='n',cex=0.75)
dev.off()


#### Just looking at the raw data ####

# restructure for categorical
raw.cca.data <- cca.data
raw.cca.data$Poppy <- as.character(raw.cca.data$Poppy)
raw.cca.data$Poppy[which(raw.cca.data$Poppy=="T")] <- "Trace"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==3)]<-"+++"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==2)]<- "++"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==1)]<- "+"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==0.5)]<- "Trace"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==0)]<- "Negative (0)"
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy)
raw.cca.data <- raw.cca.data %>% filter(Poppy!="-")
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy)
raw.cca.data <- raw.cca.data %>% filter(Poppy!="")
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy, levels=c("Negative (0)", "Trace", "+", "++","+++"))
raw.cca.data$gscore <- as.integer(raw.cca.data$gscore)

raw.cca.data$dateN <- NA
raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="25/09/2017" | raw.cca.data$date_on_tube=="26/09/2017" | raw.cca.data$date_on_tube=="27/09/2017"
                           | raw.cca.data$date_on_tube=="28/09/2017" | raw.cca.data$date_on_tube=="29/09/2017" | raw.cca.data$date_on_tube=="02/10/2017")] <- "Pre-T"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="23/10/2017" | raw.cca.data$date_on_tube=="24/10/2017" | raw.cca.data$date_on_tube=="25/10/2017"
                           | raw.cca.data$date_on_tube=="26/10/2017" | raw.cca.data$date_on_tube=="27/10/2017")] <- "3 weeks"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="04/12/2017" | raw.cca.data$date_on_tube=="05/12/2017")] <- "9 weeks"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="01/03/2018" | raw.cca.data$date_on_tube=="05/03/2018" | raw.cca.data$date_on_tube=="06/03/2018"
                           | raw.cca.data$date_on_tube=="07/03/2018" | raw.cca.data$date_on_tube=="08/03/2018" | raw.cca.data$date_on_tube=="09/03/2018")] <- "6 months"

# CCA normal by gscore 

ccaxgscore <- ggplot(data=raw.cca.data, aes(y=gscore, x = Poppy))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, colour = "#7570b3")+
  xlab("CCA original scale")+
  scale_y_continuous(breaks=1:10)+
  coord_flip()+
  ggsave("ccaxgscore.pdf")
ccaxgscore

ccadates <- raw.cca.data[,c("CID", "Poppy", "gscore", "dateN")]
ccadates$dateN <- as.factor(ccadates$dateN)

# sort out KK data for plot 

kkxgscore.data <- kk.data
kkxgscore.data$dateN <- NA
kkxgscore.data$dateN[which(kkxgscore.data$date=="25/09/2017" | kkxgscore.data$date=="26/09/2017" | kkxgscore.data$date=="27/09/2017"
               | kkxgscore.data$date=="28/09/2017" | kkxgscore.data$date=="29/09/2017" | kkxgscore.data$date=="02/10/2017")] <- "Pre-T"

kkxgscore.data$dateN[which(kkxgscore.data$date=="23/10/2017" | kkxgscore.data$date=="24/10/2017" | kkxgscore.data$date=="25/10/2017"
               | kkxgscore.data$date=="26/10/2017" | kkxgscore.data$date=="27/10/2017")] <- "3 weeks"

kkxgscore.data$dateN[which(kkxgscore.data$date=="04/12/2017" | kkxgscore.data$date=="05/12/2017")] <- "9 weeks"

kkxgscore.data$dateN[which(kkxgscore.data$date=="01/03/2018" | kkxgscore.data$date=="05/03/2018" | kkxgscore.data$date=="06/03/2018"
               | kkxgscore.data$date=="07/03/2018" | kkxgscore.data$date=="08/03/2018" | kkxgscore.data$date=="09/03/2018")] <- "6 months"

kkxgscore.data <- melt(kkxgscore.data,id.vars = c("child_id", "date", "dateN"), measure.vars = c("Sm_A", "Sm_B"))
kkxgscore.data$dateN <- as.factor(kkxgscore.data$dateN)
kkxgscore.data <- kkxgscore.data[-which(is.na(kkxgscore.data$value)==T),]

kkxgscore.data$child_id <- as.factor(kkxgscore.data$child_id)

kkxgscore.data <- kkxgscore.data %>% 
  group_by(dateN, child_id)%>%
  summarise(mean_count=mean(value, na.rm=T))

colnames(kkxgscore.data) <- c("dateN", "CID", "mean_count")

# merge data 

mergekkcca <- merge(ccadates, kkxgscore.data, by=c("CID", "dateN"))
mergekkcca$gscore <- as.factor(mergekkcca$gscore)
mergekkcca$dateN <- factor(mergekkcca$dateN, levels=c("Pre-T", "3 weeks", "9 weeks", "6 months"))
mergekkcca$logcount <- log(mergekkcca$mean_count*24)
mergekkcca[mergekkcca=="-Inf"] <- 0
# kato katz by normal CCA

kkxcca <- ggplot(data=mergekkcca, aes(y=mean_count, x = Poppy))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(group=dateN, colour=dateN))+
  xlab("CCA original scale")+
  ggsave("kkxcca.pdf")
kkxcca

# kato katz by gscore 

kkxgscore.plot <- ggplot(data=mergekkcca, aes(y=log(mean_count*24), x = gscore))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter( aes(group=dateN, colour=dateN))+
  xlab("G-Score")+
  ylab("Mean Kato-Katz")+
  ggsave("kkxgscorelog.pdf")
kkxgscore.plot

# kk by gscore with orignal cca score dots

kkxgscorexcca <- ggplot(data=mergekkcca, aes(y=logcount, x = gscore))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter( aes(group=Poppy, colour=Poppy, shape=dateN), size=3)+
  xlab("G-Score")+
  ylab("Mean Kato-Katz EPG on Log Scale")+
  labs(colour = "CCA Normal Scores", shape="Sampling Period")+
  
  ggsave("kkxgscorexccalogepg.pdf")
kkxgscorexcca

# see how many G-SCore categories there are for each KK count 
# gscore and kk
countcats <- mergekkcca
countcats$mean_count <- as.factor(countcats$mean_count)
countcats <- countcats %>% 
  group_by(mean_count, gscore) %>% 
  summarise(combos=n())%>%
  count(mean_count)

count.var <- var(countcats$n)
sdev <- sqrt(count.var)
dev <- sdev/2
prec <- 1/(dev^2)

# cca and kk
countcca <- mergekkcca %>%
  mutate_at(vars(mean_count),as.factor)%>%
  group_by(mean_count, Poppy) %>%
  summarise(ccakk_combos=n()) %>%
  count(mean_count)

kkccavar <- var(countcca$n)
sdevkkcca <- sqrt(kkccavar)
dev2 <- sdevkkcca/2
prec2 <- 1/dev2
