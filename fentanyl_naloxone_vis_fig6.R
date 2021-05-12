# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 6 CODE -- PopPK: Smoking

#### PACKAGES ####

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')

root=getwd()

# Set themeing
my_theme <- theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))


################################################################################
#------------------------------ POP SMOKE -------------------------------------# 

# read in data
setwd(paste0(root,'/Data/pop_smoke'))
smoke_mOR <- as.data.frame(readMat("smokers_mOR.mat"))
smoke_AUC <- as.data.frame(readMat("popsmoke_AUC.mat"))
smoke_surv <- as.data.frame(readMat("popsmoke_surv.mat"))

smokeData <- cbind(smoke_mOR, smoke_AUC, smoke_surv)
names(smokeData) <- c("Available_mOR", "Output_AUC", "Survival")


# plot data
smokePlotAUC <- ggplot(data = NULL) +
          geom_line(data = smokeData, aes(x = Available_mOR, y = Output_AUC), size = 1.5) +
          xlab('Available mOR (nmol)') +
          ylab('AUC') +
          labs(title = 'Population Pharmacokinetics: The Effect of Smoking on AUC') +
          my_theme

smokePlotSURV <- ggplot(data = NULL) +
    geom_line(data = smokeData, aes(x = Available_mOR, y =Survival), size = 1.5) +
    xlab('Available mOR (nmol)') +
    ylab('Survival Rate') +
    labs(title = 'Population Pharmacokinetics: The Effect of Smoking on Survival') +
    my_theme

smokePlot <- (smokePlotAUC) / (smokePlotSURV) + plot_layout(guides = 'collect')

dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave(filename = 'Fig6.png', plot = smokePlot, width = 12, height = 8)
setwd(paste0(root))

