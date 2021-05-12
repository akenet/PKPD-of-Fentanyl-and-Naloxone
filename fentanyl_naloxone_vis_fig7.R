# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 7 CODE -- PopPK: Naloxone Absorption

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
#------------------------------ POP NOSE ------------------------------------# 

# read in data
setwd(paste0(root,'/Data/pop_nose'))
nose_abs <- as.data.frame(readMat("nose_abs.mat"))
nose_AUC <- as.data.frame(readMat("popnose_AUC.mat"))
nose_surv <- as.data.frame(readMat("popnose_surv.mat"))

noseData <- cbind(nose_abs, nose_AUC, nose_surv)
names(noseData) <- c("nose_abs", "Output_AUC", "Survival")


# plot data

nosePlotAUC <- ggplot(data = NULL) +
  geom_line(data = noseData, aes(x = nose_abs, y = Output_AUC), size = 1.5) +
  xlab('Naloxone Bioavailability (%)') +
  ylab('AUC') +
  labs(title = 'Population Pharmacokinetics: The Effect of Nasal Congestion on AUC') +
  my_theme


nosePlotSURV <- ggplot(data = NULL) +
  geom_line(data = noseData, aes(x = nose_abs, y =Survival), size = 1.5) +
  xlab('Naloxone Bioavailability (%)') +
  ylab('Survival Rate') +
  labs(title = 'Population Pharmacokinetics: The Effect of Nasal Congestion on Survival') +
  my_theme


nosePlot <- (nosePlotAUC) / (nosePlotSURV) + plot_layout(guides = 'collect')

dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave(filename = 'Fig7.png', plot = nosePlot, width = 12, height = 8)
setwd(paste0(root))