# Loading the libraries
library(tidyverse)

# Loading the data
ExpGrowth <- read.csv("../data/ExpGrowth.csv")

# Best Results
modelResults <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/modelResultsRun2.csv")

# Original Results
originalModel <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/allDataModelResultsOriginal.csv")

# Original Params 
params <- read.csv("~/Documents/Project/tempDepCode/data/model_enzyme_params_new_tagged.csv")

# Best Params
BestParams <- read.csv("~/Documents/Project/finalResults/bestRunYet/data/BestParamsTopt9.csv")

bestNow <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/allDataModelResultsBsetTillNow.csv")


######## Cleaning the data #############

#### Growth Data #######
growth <- ExpGrowth %>% filter(ID==0)
rate <- growth$r

############# Actual plotting ############
xlabel <- expression("Temperature " ( degree~C))
ylabel <- expression('Growth Rate'(h^-1))
co2 <- expression((CO_2))

########## Best fit #########
bestFit <- modelResults %>% filter(R2 == max(R2))
growth <- ExpGrowth %>% filter(ID==0)
rate <- growth$r
toPlot <- bestFit
toPlot <- toPlot %>% mutate(exp=rate)

svg('../graphs/goodfit.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=growth, color='growth'), size =1)
p <- p + geom_point(aes(y=exp, color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=20),
      axis.title=element_text(size=22), 
      axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                     size=20),
          axis.text.y = element_text( color="black", 
                                     size=20))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#000000", "#0072B2"),
                              name  =NULL,
                              breaks=c("exp", "growth"),
                              labels=c("Empirical Data", "Model Prediction"))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

#ggsave(filename = '../graphs/goodfit.svg', p + annotate("text", x= 22, y=1.3,label = 'R^2==0.8641298', parse=T, size=4.5)+ annotate("text", x= 22.5, y=1.2,label = 'MSE==0.02091063', parse=T, size=4.5))

p + annotate("text", x= 22, y=1.3,label = 'R^2==0.9432632', parse=T, size=4.5)+ annotate("text", x= 22.5, y=1.2,label = 'MSE==0.008731873', parse=T, size=4.5)
dev.off()


################ Plot for fluxes ##################
svg('../graphs/goodfitFlux.svg',  width = 8, height = 6)


p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=o2flux, color='o2flux'), size =1)
p <- p + geom_line(aes(y=gluFlux, color='gluFlux'), size=1) 
p <- p + geom_line(aes(y=-co2flux, color='co2flux'), size=1) 
p <- p + geom_line(aes(y=-acetateFlux, color='acetateFlux'), size=1) 


p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab('Reaction Flux')
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#D55E00", "#0072B2", "#CC79A7", 
                                     '#009E73'),
                            name  =NULL,
                            breaks=c("o2flux", "gluFlux", 
                                     'co2flux', 'acetateFlux'),
                            labels=c("Oxygen Uptake", 
                                     "Glucose Uptake", 'CO2 production',
                                     'Acetate Production'))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

p

dev.off()

############### Next Acutal Plot ################
bestFit <- modelResults %>% filter(R2 == max(R2))
worstFit <- modelResults %>% filter(R2 == min(R2))
growth <- ExpGrowth %>% filter(ID==0)
rate <- growth$r

toPlotShade <- bestFit %>% select(temp, growth)
toPlotShade <- toPlotShade %>% mutate(badGrowth= worstFit$growth)
toPlotShade <- toPlotShade %>% mutate(exp=rate)

svg('../graphs/fitRange.svg', width = 8, height = 6)
q <- ggplot(toPlotShade, aes(x=temp))
q <- q +  geom_line(aes(y=growth, color='growth'), size =1)
q <- q +  geom_line(aes(y=badGrowth, color='badGrowth'), size =1)
q <- q + geom_ribbon(aes(x = temp,ymin = growth,ymax = badGrowth), 
                     fill = "gray",alpha = 0.4)
q <- q + theme_minimal()
q <- q + xlab(xlabel) + ylab(ylabel)
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
q <- q + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
q <- q + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
q <- q + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
q <- q + scale_color_manual(values=c("#D55E00", "#0072B2"),
                            name  =NULL,
                            breaks=c("badGrowth", "growth"),
                            labels=c("Worst Prediction", "Best Prediction"))
q <- q + theme(legend.justification=c(0,1), legend.position=c(0,1))

q
dev.off()

############ Histograms ##################

toPlotHist <- BestParams %>% select(Topt) %>% 
    `colnames<-`(c("Temp"))
toPlotHistTm <- BestParams %>% select(Tm)%>% 
    `colnames<-`(c("Temp"))
toPlotHist <- toPlotHist %>% mutate(Tag=rep('T_opt', length(Temp)))
toPlotHistTm <- toPlotHistTm %>% mutate(Tag=rep('T_m', length(Temp)))

toPlotHist <- rbind(toPlotHist, toPlotHistTm)

svg('../graphs/finalHist.svg', width = 8, height = 6)

t <- ggplot(toPlotHist, aes(x = Temp, fill = Tag)) +   
    geom_histogram(position = "identity", alpha = 0.8, bins = 30)
t <- t + theme_minimal()
t <- t + xlab(xlabel) 
t <- t + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
t <- t + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
t <- t + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
t <- t + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 13),)
t <- t + scale_fill_manual(values=c("#D55E00", "#0072B2"),name  =NULL,
                            breaks=c("T_m", "T_opt"),
                            labels=c("Melting Temp", "Optimal Temp"))
t <- t + theme(legend.justification=c(1,1), legend.position=c(1,1))
t

dev.off()

############### Low Temp plot #################
bestFit <- bestNow %>% filter(iter==1)
growth <- ExpGrowth %>% filter(ID==1)
rate <- growth$r
toPlot <- bestFit
toPlot <- toPlot %>% mutate(exp=rate)


svg('../graphs/goodfitLow.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=growth, color='growth'), size =1)
p <- p + geom_point(aes(y=exp, color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#000000", "#0072B2"),
                            name  =NULL,
                            breaks=c("exp", "growth"),
                            labels=c("Empirical Data", "Model Prediction"))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

p + annotate("text", x= 12, y=1.3,label = 'R^2==0.9299357', parse=T, size=4.5)+ annotate("text", x= 12.5, y=1.2,label = 'MSE==0.01329684', parse=T, size=4.5)

dev.off()


############# Original model Plots #########
punchingBagSet <- originalModel %>% filter(iter==0)
growth <- ExpGrowth %>% filter(ID==0)
rate <- growth$r
toPlot <- punchingBagSet
toPlot <- toPlot %>% mutate(exp=rate)


svg('../graphs/original.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=growth, color='growth'), size =1) + 
    geom_point(aes(y=growth, color='growth'), size=2)
p <- p + geom_point(aes(y=exp, color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#000000", "#0072B2"),
                            name  =NULL,
                            breaks=c("exp", "growth"),
                            labels=c("Empirical Data", "Model Prediction"))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

#ggsave(filename = '../graphs/goodfit.svg', p + annotate("text", x= 22, y=1.3,label = 'R^2==0.8641298', parse=T, size=4.5)+ annotate("text", x= 22.5, y=1.2,label = 'MSE==0.02091063', parse=T, size=4.5))

#p + annotate("text", x= 12, y=1.3,label = 'R^2==0.9299357', parse=T, size=4.5)+ annotate("text", x= 12.5, y=1.2,label = 'MSE==0.01735985', parse=T, size=4.5)
p

dev.off()


################ Basal Salt Original ###################

basalSaltSet <- originalModel %>% filter(iter==22)
growth <- ExpGrowth %>% filter(ID==22)
rate <- growth$r
toPlot <- basalSaltSet
toPlot <- toPlot %>% mutate(exp=rate)


svg('../graphs/originalBS.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=growth, color='growth'), size =1) + 
    geom_point(aes(y=growth, color='growth'), size=2)
p <- p + geom_point(aes(y=exp, color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#000000", "#0072B2"),
                            name  =NULL,
                            breaks=c("exp", "growth"),
                            labels=c("Empirical Data", "Model Prediction"))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

#ggsave(filename = '../graphs/goodfit.svg', p + annotate("text", x= 22, y=1.3,label = 'R^2==0.8641298', parse=T, size=4.5)+ annotate("text", x= 22.5, y=1.2,label = 'MSE==0.02091063', parse=T, size=4.5))

#p + annotate("text", x= 12, y=1.3,label = 'R^2==0.9299357', parse=T, size=4.5)+ annotate("text", x= 12.5, y=1.2,label = 'MSE==0.01735985', parse=T, size=4.5)
p

dev.off()

################# Low Temp ############# 

LowTempSet <- originalModel %>% filter(iter==1)
growth <- ExpGrowth %>% filter(ID==1)
rate <- growth$r
toPlot <- LowTempSet
toPlot <- toPlot %>% mutate(exp=rate)


svg('../graphs/originalLow.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=growth, color='growth'), size =1) + 
    geom_point(aes(y=growth, color='growth'), size=2)
p <- p + geom_point(aes(y=exp, color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
p <- p + scale_color_manual(values=c("#000000", "#0072B2"),
                            name  =NULL,
                            breaks=c("exp", "growth"),
                            labels=c("Empirical Data", "Model Prediction"))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

#ggsave(filename = '../graphs/goodfit.svg', p + annotate("text", x= 22, y=1.3,label = 'R^2==0.8641298', parse=T, size=4.5)+ annotate("text", x= 22.5, y=1.2,label = 'MSE==0.02091063', parse=T, size=4.5))

#p + annotate("text", x= 12, y=1.3,label = 'R^2==0.9299357', parse=T, size=4.5)+ annotate("text", x= 12.5, y=1.2,label = 'MSE==0.01735985', parse=T, size=4.5)
p

dev.off()

########### Original Hist #############
toPlotHist <- params %>% filter(topt_source!='BullShit') %>% 
    select(Topt) %>% `colnames<-`(c("Temp"))
toPlotHistTm <- params %>% filter(TmTag=='Exp') %>% 
    select(Tm)%>% `colnames<-`(c("Temp"))
toPlotHist <- toPlotHist %>% mutate(Tag=rep('T_opt', length(Temp)))
toPlotHistTm <- toPlotHistTm %>% mutate(Tag=rep('T_m', length(Temp)))

toPlotHist <- rbind(toPlotHist, toPlotHistTm)

svg('../graphs/originalHist.svg', width = 8, height = 6)

t <- ggplot(toPlotHist, aes(x = Temp, fill = Tag)) +   
    geom_histogram(position = "identity", alpha = 0.8, bins = 30)
t <- t + theme_minimal()
t <- t + xlab(xlabel) 
t <- t + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
t <- t + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_text( color="black", 
                                           size=14))
t <- t + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
t <- t + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 13),)
t <- t + scale_fill_manual(values=c("#D55E00", "#0072B2"),name  =NULL,
                           breaks=c("T_m", "T_opt"),
                           labels=c("Melting Temp", "Optimal Temp"))
t <- t + theme(legend.justification=c(1,1), legend.position=c(1,1))
t

dev.off()

