# Loading the required librarires
library(tidyverse)

# loading the data
onlyKcat <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/onlyKcat.csv")

onlyNGAM <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/onlyNGAM.csv")

onlyTm <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/onlyTm.csv")

########### Axes Labels ############
xlabel <- expression("Temperature " ( degree~C))
ylabel <- expression('Growth Rate'(h^-1))
topt <- expression(T_opt)
tm <- expression(T_m)
cp <- expression(C_p)

################ Actual Plot ############### 
kcatGrowthMin <- onlyKcat %>% group_by(temp) %>% 
    summarise(minGrowthTopt=min(growth),.groups = 'drop')
kcatGrowthMax <- onlyKcat %>% group_by(temp) %>%
    summarise(maxGrowthTopt=max(growth),.groups = 'drop')
kcatGrowthMean <- onlyKcat %>% group_by(temp) %>% 
    summarise(meanGrowthTopt=mean(growth),.groups = 'drop')

NGAMGrowthMin <- onlyNGAM %>% group_by(temp) %>% summarise(minGrowthNGAM=min(growth),.groups = 'drop')
NGAMGrowthMax <- onlyNGAM %>% group_by(temp) %>% summarise(maxGrowthNGAM=max(growth),.groups = 'drop')
NGAMGrowthMean <- onlyNGAM %>% group_by(temp) %>% 
    summarise(meanGrowthNGAM=mean(growth),.groups = 'drop')

TmGrowthMin <- onlyTm %>% group_by(temp) %>% summarise(minGrowthTm=min(growth),
                                                       .groups = 'drop')
TmGrowthMax <- onlyTm %>% group_by(temp) %>% summarise(maxGrowthTm=max(growth),
                                                      .groups = 'drop')
TmGrowthMean <- onlyTm %>% group_by(temp) %>% 
    summarise(meanGrowthTm=mean(growth),.groups = 'drop')

toPlot <- merge(kcatGrowthMean, kcatGrowthMax, by = 'temp')
toPlot <- merge(toPlot, kcatGrowthMin, by = 'temp')
toPlot <- merge(toPlot, NGAMGrowthMax, by = 'temp')
toPlot <- merge(toPlot, NGAMGrowthMean, by = 'temp')
toPlot <- merge(toPlot, NGAMGrowthMin, by = 'temp')
toPlot <- merge(toPlot, TmGrowthMax, by = 'temp')
toPlot <- merge(toPlot, TmGrowthMean, by = 'temp')
toPlot <- merge(toPlot, TmGrowthMin, by = 'temp')
toPlot <- toPlot %>% mutate(ecGEM=0.6782999999999908)


svg('../graphs/ecAdditions.svg', width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p + geom_ribbon(aes(x = temp,ymin = minGrowthTm ,ymax = maxGrowthTm), 
                     fill = "#0072B2",alpha = 0.2)
p <- p + geom_ribbon(aes(x = temp,ymin = minGrowthNGAM ,ymax = maxGrowthNGAM), 
                     fill = "#CC79A7",alpha = 0.2)
p <- p + geom_ribbon(aes(x = temp,ymin = minGrowthTopt ,ymax = maxGrowthTopt), 
                     fill = "#D55E00",alpha = 0.2)
p <- p + geom_line(aes(y=meanGrowthTopt, colour='meanGrowthTopt'), size=1)
p <- p + geom_line(aes(y=meanGrowthTm, colour='meanGrowthTm'), size=1)
p <- p + geom_line(aes(y=meanGrowthNGAM, colour='meanGrowthNGAM'), size=1)
p <- p + geom_line(aes(y=ecGEM, colour='ecGEM'), size=1.2)

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel) + ylim(c(0, 0.82))
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
                            breaks=c("meanGrowthTopt", "meanGrowthTm",
                                     'meanGrowthNGAM', 'ecGEM'),
                            labels=c("ecGEM + Enzyme Activity", 
                                     "ecGEM + Enzyme Denaturation",
                                     "ecGEM + Maintenance",'ecGEM'))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

p
dev.off()


############# Loading more data #######################

wigCp <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/wigCp.csv")

wigTm <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/wigTm.csv")

wigTopt <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/results/wigTopt.csv")

originalModel <- read.csv("~/Documents/Project/tempDepCode/results/desktopOnes/allDataModelResultsOriginal.csv")

############ Cleaning data ############
smol <- originalModel %>% filter(iter==0)

Tm <- wigTm %>% group_by(temp) %>% 
    summarise(minTm=min(growth), maxTm=max(growth), 
              meanTm=mean(growth), .groups = 'drop')

Topt <- wigTopt %>% group_by(temp) %>% 
    summarise(minTopt=min(growth), maxTopt=max(growth), 
              meanTopt=mean(growth), .groups = 'drop')

Cp <- wigCp %>% group_by(temp) %>% 
    summarise(minCp=min(growth), maxCp=max(growth), 
              meanCp=mean(growth), .groups = 'drop')

toPlot <- merge(Tm, Topt, by='temp')
toPlot <- merge(toPlot, Cp, by='temp')
toPlot <- toPlot %>% mutate(etcGEM=smol$growth)

svg('../graphs/wiggle.svg', width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p + geom_ribbon(aes(x = temp,ymin = minTm ,ymax = maxTm), 
                     fill = "#0072B2",alpha = 0.2)
p <- p + geom_ribbon(aes(x = temp,ymin = minCp ,ymax = maxCp), 
                     fill = "#CC79A7",alpha = 0.2)
p <- p + geom_ribbon(aes(x = temp,ymin = minTopt ,ymax = maxTopt), 
                     fill = "#D55E00",alpha = 0.2)
p <- p + geom_line(aes(y=meanTopt, colour='meanTopt'), size=1)
p <- p + geom_line(aes(y=meanTm, colour='meanTm'), size=1)
p <- p + geom_line(aes(y=meanCp, colour='meanCp'), size=1)
p <- p + geom_line(aes(y=etcGEM), size=1.2)

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel) + ylim(c(0, 0.82))
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
                            breaks=c("meanTopt", "meanTm",
                                     'meanCp'),
                            labels=c(topt,tm,cp))
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))

p

dev.off()
