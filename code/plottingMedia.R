# Clearing the workspace
rm(list = ls())

# Loading the libraries
library(tidyverse)

# Loading the data
ExpGrowth <- read.csv("../data/ExpGrowth.csv")

###### Possible Labels ####
xlabel <- expression("Temperature " ( degree~C))
ylabel <- expression('Growth Rate'(h^-1))
co2 <- expression((CO_2))

scatterXlabel <- expression(Log('Observed Growth Rate'(h^-1)))
scatterYlabel <- expression(Log('Predicted Growth Rate'(h^-1)))

##### Loading model predictions #####

allDataModelResultsBsetTillNow <- read.csv("~/Documents/Project/tempDepCode/results/allDataModelResultsBsetTillNow.csv")

mmConsB <- read.csv("~/Documents/Project/tempDepCode/results/mmConsB.csv")

mmOnly <- read.csv("~/Documents/Project/tempDepCode/results/mmOnly.csv")

onlyConsB <- read.csv("~/Documents/Project/tempDepCode/results/onlyConsB.csv")

lbOnly <- read.csv("~/Documents/Project/tempDepCode/results/lbOnly.csv")

lbRelaxB <- read.csv("~/Documents/Project/tempDepCode/results/lbRelaxB.csv")

onlyRelaxB <- read.csv("~/Documents/Project/tempDepCode/results/onlyRelaxB.csv")

########## Mutating the data #######
expMM <- ExpGrowth %>% filter(ID==25) %>% select(Ts,r)
expLB <- ExpGrowth %>% filter(ID==12) %>% select(Ts,r)

mmNormB <- allDataModelResultsBsetTillNow %>% 
    filter(iter==25) %>% 
    mutate(type='mmNormB', expGrowth=expMM$r) %>% 
    select(temp, growth, type, expGrowth)

lbNormB <- allDataModelResultsBsetTillNow %>% 
    filter(iter==12) %>% 
    mutate(type='lbNormB', expGrowth=expLB$r) %>% 
    select(temp, growth, type, expGrowth)

mmConsB <- mmConsB %>% mutate(type='mmConsB', expGrowth=expMM$r) %>% 
    select(temp, growth, type, expGrowth)

mmOnly <- mmOnly %>% mutate(type='mmOnly', expGrowth=expMM$r) %>% 
    select(temp, growth, type, expGrowth)

onlyConsB <- onlyConsB %>% mutate(type='onlyConsB', expGrowth=expMM$r) %>% 
    select(temp, growth, type, expGrowth)

lbRelaxB <- lbRelaxB %>% mutate(type='lbRelaxB', expGrowth=expLB$r) %>% 
    select(temp, growth, type, expGrowth)

lbOnly <- lbOnly %>% mutate(type='lbOnly', expGrowth=expLB$r)%>% 
    select(temp, growth, type, expGrowth)

onlyRelaxB <- onlyRelaxB %>% mutate(type='onlyRelaxB', expGrowth=expLB$r)%>% 
    select(temp, growth, type, expGrowth)

###### Making to plot frames #####

mmToPlot <- rbind(mmNormB, mmConsB, mmOnly, onlyConsB)
lbToPlot <- rbind(lbNormB, lbRelaxB, lbOnly, onlyRelaxB)

######## Actual Minimal Media plotting ########

svg('../graphs/mmScatter.svg',  width = 8, height = 6)

p <- ggplot(mmToPlot, aes(x=log(expGrowth)))
p <- p + geom_point((aes(y=log(growth), color=type)), size=3)
p <- p + geom_abline(slope = 1, intercept = 0)

p <- p + theme_minimal()
p <- p + xlab(scatterXlabel) + ylab(scatterYlabel) + xlim(c(-6, 0.5)) + 
    ylim(c(-6, 0.5))
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
p <- p + scale_color_manual(values=c("#D55E00", "#0072B2", "#CC79A7", 
                                     '#009E73'),
                            name  =NULL,
                            breaks=c("mmConsB", "mmNormB",
                                     "mmOnly", "onlyConsB"),
                            labels=c("Minimal Media Constrained Bounds", 
                                     "Normal Media Normal Bounds",
                                     "Minimal Media Normal Bounds",
                                     "Normal Media Constrained Bounds"))
p <- p + theme(legend.justification=c(1,0), legend.position=c(1,0))
p

dev.off()

########## Actual LB plotting ###########

svg('../graphs/lbScatter.svg',  width = 8, height = 6)

p <- ggplot(lbToPlot, aes(x=log(expGrowth)))
p <- p + geom_point((aes(y=log(growth), color=type)), size=3)
p <- p + geom_abline(slope = 1, intercept = 0)

p <- p + theme_minimal()
p <- p + xlab(scatterXlabel) + ylab(scatterYlabel)+ 
    xlim(c(-4.5, 1)) + ylim(c(-4.5, 1))
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
p <- p + scale_color_manual(values=c("#D55E00", "#0072B2", "#CC79A7", 
                                     '#009E73'),
                            name  =NULL,
                            breaks=c("lbNormB", "lbOnly",
                                     "lbRelaxB", "onlyRelaxB"),
                            labels=c("Regular Media Normal Bounds", 
                                     "LB Media Normal Bounds",
                                     "LB Media Relaxed Bounds",
                                     "Regular Media Relaxed Bounds"))
p <- p + theme(legend.justification=c(1,0), legend.position=c(1,0))
p

dev.off()

###### Single MM Curve ########
mmToPlotSing <- mmToPlot %>% filter(type=='mmConsB')

svg('../graphs/mmSing.svg',  width = 8, height = 6)

p <- ggplot(mmToPlotSing, aes(x=temp))
p <- p + geom_line((aes(y=growth, color='growth')), size=1)
p <- p + geom_point((aes(y=expGrowth)), size=2.5)


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
p <- p + scale_color_manual(values=c("#0072B2"),
                            name  =NULL,
                            breaks=c("growth"),
                            labels=c("Minimal Media Constrained Bounds"))
p <- p + theme(legend.justification=c(1,0), legend.position=c(1,0))
p

dev.off()

## Combined one

q <- ggplot(mmToPlotSing, aes(x=temp))
q <- q + geom_line((aes(y=log(growth), color='growth')), size=1)
q <- q + geom_point((aes(y=log(expGrowth))), size=2.5)


q <- q + theme_minimal()
q <- q + xlab(xlabel) + ylab(ylabel) 
q <- q + theme(axis.text=element_text(size=20),
               axis.title=element_blank(), 
               axis.line = element_line(colour = "black"))
q <- q + theme(axis.text.x = element_text(color="black", 
                                          size=20),
               axis.text.y = element_text( color="black", 
                                           size=20))
q <- q + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())
q <- q + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
q <- q + scale_color_manual(values=c("#0072B2"),
                            name  =NULL,
                            breaks=c("growth"),
                            labels=c("Minimal Media Constrained Bounds"))
q <- q + theme(legend.position = "none") 
q

svg('../graphs/mmSingInset.svg',  width = 8, height = 6)
plot.with.inset <-
    ggdraw() +
    draw_plot(p) +
    draw_plot(q, x = 0.15, y = .65, width = .3, height = .3)
plot.with.inset
dev.off()

###### Single LB Curve ########

lbToPlotSing <- lbToPlot %>% filter(type=='lbRelaxB')

svg('../graphs/lbSing.svg',  width = 8, height = 6)

p <- ggplot(lbToPlotSing, aes(x=temp))
p <- p + geom_line((aes(y=growth, color='growth')), size=1)
p <- p + geom_point((aes(y=expGrowth)), size=2.5)


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
p <- p + scale_color_manual(values=c("#0072B2"),
                            name  =NULL,
                            breaks=c("growth"),
                            labels=c("LB Media Relaxed Bounds"))
p <- p + theme(legend.justification=c(1,0), legend.position=c(1,0))
p

dev.off()

## Combined one

q <- ggplot(lbToPlotSing, aes(x=temp))
q <- q + geom_line((aes(y=log(growth), color='growth')), size=1)
q <- q + geom_point((aes(y=log(expGrowth))), size=2.5)


q <- q + theme_minimal()
q <- q + xlab(xlabel) + ylab(ylabel) 
q <- q + theme(axis.text=element_text(size=20),
               axis.title=element_blank(), 
               axis.line = element_line(colour = "black"))
q <- q + theme(axis.text.x = element_text(color="black", 
                                          size=20),
               axis.text.y = element_text( color="black", 
                                           size=20))
q <- q + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
q <- q + theme(legend.title = element_text(colour = "black", size = 12.5),
               legend.text = element_text(colour = "black", size = 15),)
q <- q + scale_color_manual(values=c("#0072B2"),
                            name  =NULL,
                            breaks=c("growth"),
                            labels=c("Minimal Media Constrained Bounds"))
q <- q + theme(legend.position = "none") 
q

svg('../graphs/lbSingInset.svg',  width = 8, height = 6)
plot.with.inset <-
    ggdraw() +
    draw_plot(p) +
    draw_plot(q, x = 0.1, y = .65, width = .3, height = .3)
plot.with.inset
dev.off()
