# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 9 CODE -- Sensitivity

#### PACKAGES ####

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')

root=getwd()

################################################################################
#------------------------------ SENSITIVITY -----------------------------------# 

# Set themeing
heatmap_theme <- theme_classic() + 
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.key.size=unit(0.5,'inch'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))

setwd(paste0(root,'/Data/sensitivity'))
sensitivityData <- as.data.frame(readMat("sensitivity.mat"))

drugs <- c("Fentanyl", "Naloxone")
parameters <- c('Vd_brain', 'amt_mOR', 'plasma_f', 
                'dose_f', 'dose_n', 'abs_n', 
                'k12_f', 'k12_n', 
                'k21_f', 'k21_n', 
                'k13_f', 'k13_n', 
                'k31_f', 'k31_n', 
                'kcl_f', 'kcl_n', 
                'kon_f', 'kon_n', 
                'koff_f', 'koff_n', 
                'Vd_f1', 'Vd_f2', 
                'Vd_n1', 'Vd_n2')

sensitivityData <- cbind(parameters, sensitivityData)
names(sensitivityData) <- c("Parameter", drugs)

sensitivityData <- sensitivityData[, 1:2]

# Reformat data so all of the output is in one column
sensitivityData <- pivot_longer(sensitivityData, !Parameter, names_to = 'Drug', values_to = 'Sens')

# Change parameter column to factor and specify order
sensitivityData$Parameter <- factor(sensitivityData$Parameter, levels = rev(parameters))

# plot heatmap
localSensPlot <- ggplot(data = sensitivityData, aes(x = Drug, y = Parameter, fill = Sens)) +
  geom_tile() +
  scale_fill_distiller(name = 'Sensitivity', palette = 'RdYlBu') + 
  labs(title = 'Parameter Local Sensitivity') +
  heatmap_theme 

dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave(filename = 'Fig9.png', plot = localSensPlot, width = 6, height = 6)
setwd(paste0(root))