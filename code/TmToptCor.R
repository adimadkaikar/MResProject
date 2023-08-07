# Cleaning the environment
rm(list=ls())

# Loading the required packages
require(tidyverse)

# Loading the data
enzyme_params <- read_csv("../data/model_enzyme_params_new_tagged.csv")

# Filtering the data to keep only experimental values
small <- enzyme_params %>% filter(TmTag=='Exp', topt_source!='BullShit')

# Plotting
p <- ggplot(small, aes(x = Tm, y = Topt)) + 
    geom_point() + geom_smooth(method = "lm")
p <- p + theme_bw()
p

# Storing the observed correlation coefficient
o_cor <- cor(small$Tm, small$Topt)

# Setting the number of resampling itterations
n <- 10000
# Creating a vector to store the correlation coefficient for each resampling
cor_vec <- vector(, n)

# Resampling and calculating the correlation coefficient for each resampling 
# and storing it in the vector the loopy way
for(i in 1:n){
    tmp <-  sample(small$Topt, replace = F)
    cor_vec[i] <- cor(small$Tm, tmp)
}

# Plotting the histogram for the simulated cor and the overseved cor as an abline
g <- ggplot(data = as.data.frame(cor_vec), aes(x=cor_vec))+
    geom_histogram(aes(y=..density..), bins = 50, fill = I('#505F90')) +
    geom_density(colour = 'red')+
    theme_bw() 
g<- g + geom_vline(aes(xintercept=o_cor), colour="red", lwd = 1.1) +
    labs(title = "Histogram of the permutation analysis", 
         x = "Correlation coefficient") + 
    theme(plot.title = element_text(hjust = 0.5))
g

# Printing the probability of observed coefficient 
print(paste("The probabilty that the observed coefficient is by chance is:", sum(cor_vec > o_cor)/n))
