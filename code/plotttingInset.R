# Clearing the workspace
rm(list = ls())

# Loading the libraries
library(tidyverse)
library(cowplot)

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

######### 
bestFit <- bestNow %>% filter(iter==12)
growth <- ExpGrowth %>% filter(ID==12)
rate <- growth$r
toPlot <- bestFit
toPlot <- toPlot %>% mutate(exp=rate)


#svg('../graphs/wonkyMedia.svg',  width = 8, height = 6)

p <- ggplot(toPlot, aes(x=temp))
p <- p +  geom_line(aes(y=log(growth), color='growth'), size =1)
p <- p + geom_point(aes(y=log(exp), color='exp'), size=2.5) 

p <- p + theme_minimal()
p <- p + xlab(xlabel) + ylab(ylabel)
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_blank(), 
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
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1)) + theme(legend.position = "none") 

#p + annotate("text", x= 12, y=1.3,label = 'R^2==0.4532037', parse=T, size=4.5)+ annotate("text", x= 12.5, y=1.2,label = 'MSE==0.3573476', parse=T, size=4.5)

#dev.off()

q <- ggplot(toPlot, aes(x=temp))
q <- q +  geom_line(aes(y=growth, color='growth'), size =1)
q <- q + geom_point(aes(y=exp, color='exp'), size=2.5) 

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
q <- q + scale_color_manual(values=c("#000000", "#0072B2"),
                            name  =NULL,
                            breaks=c("exp", "growth"),
                            labels=c("Empirical Data", "Model Prediction"))
q <- q + theme(legend.justification=c(1,0), legend.position=c(1,0))

svg('../graphs/wonkyMediaInset.svg',  width = 8, height = 6)
plot.with.inset <-
    ggdraw() +
    draw_plot(q) +
    draw_plot(p, x = 0.07, y = .65, width = .3, height = .3)
plot.with.inset
dev.off()