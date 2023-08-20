# load in ggplot
library(ggplot2)
library(nls.multstart)
library(rTPC)

# subset for the first TPC curve
data('chlorella_tpc')
d <- subset(chlorella_tpc, curve_id == 1)

# get start values and fit model
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
# fit model
mod <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref = 20),
                     data = d,
                     iter = c(3,3,3,3),
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                     upper = get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)

# look at model fit
summary(mod)

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

# plot
xlabel <- expression("Temperature " ( degree~C))

svg('../graphs/TPC.svg', width = 8, height = 6)
p <- ggplot(preds)  + 
    geom_rect(aes(xmin=31,xmax=35,ymin=-Inf,ymax=Inf),
              color="#CCCCCC",alpha=0.005) +
    geom_rect(aes(xmin=39.66667,xmax=43.66667,ymin=-Inf,ymax=Inf),
              color="#CCCCCC",alpha=0.005)
    
p <- p + geom_line(aes(temp, .fitted), col = '#0072B2', size=1) +
    theme_minimal()
p <- p + xlab(xlabel) + ylab('Physiological or Biochemical Rate')
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=18), 
               axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", 
                                          size=14),
               axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p + geom_point(aes(x=33, y = 1.0307816244), col='#D55E00', size=3.5)
p <- p + geom_point(aes(x=41.66667, 
                        y = 1.8126641716), col='#D55E00', size=3.5)
p <- p + geom_vline(xintercept=c(33,41.66667), linetype="dotted")
p
dev.off()
