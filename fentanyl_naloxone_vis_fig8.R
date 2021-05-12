# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 8 CODE -- PopPK: Smoking and Naloxone Absorption

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
#----------------------------- POP SMOSE ------------------------------------# 

# read in data
setwd(paste0(root,'/Data/pop_smose'))
smose_values <- as.data.frame(readMat("smose_values.mat"))
smoke_AUC <- as.data.frame(readMat("smose_AUC.mat"))
smoke_surv <- as.data.frame(readMat("smose_surv.mat"))

smoseData_AUC <- cbind(smose_values, smoke_AUC)
names(smoseData_AUC) <- c("Available_mOR", "nose_abs", "Output_AUC")
smoseData_surv <- cbind(smose_values, smoke_surv)
names(smoseData_surv) <- c("Available_mOR", "nose_abs", "Survival")


# plot data

# Set themeing
heatmap_theme <- theme_minimal() + 
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        # axis.ticks.x = element_blank(),
        # panel.grid = element_blank(),
        legend.key.size=unit(0.5,'inch'),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))

keeps <- c(1, 4, 6, 7)

smosePlotAUC <- ggplot(data = smoseData_AUC, aes(x = factor(nose_abs), y = factor(Available_mOR), fill = Output_AUC)) +
  geom_tile() +
  scale_fill_distiller(name = 'AUC', palette = 'Blues', direction = 1) + 
  labs(title = 'Population Pharmacokinetics: \n The Effect of Smoking and Nasal Congestion on AUC') +
  heatmap_theme +
  # scale_x_discrete(labels = NULL, breaks = NULL) +
  # scale_y_discrete(labels = NULL, breaks = NULL) +
  ylab('Available mOR (nmol)') +
  xlab('Bioavailability of Naloxone (%)') +
  scale_x_discrete(
    breaks = levels(smoseData_AUC$nose_abs)[keeps],
    labels = levels(smoseData_AUC$nose_abs)[keeps]) + 
  scale_y_discrete(
    breaks = levels(smoseData_AUC$Available_mOR)[keeps],
    labels = levels(smoseData_AUC$Available_mOR)[keeps]) # + theme_gray(base_size = 14)


smosePlotSurvival <- ggplot(data = smoseData_surv, aes(x = factor(nose_abs), y = factor(Available_mOR), fill = Survival)) +
  geom_tile() +
  scale_fill_distiller(name = 'Survival', palette = 'Blues', direction = 1) + 
  labs(title = 'Population Pharmacokinetics: \n The Effect of Smoking and Nasal Congestion on Survival') +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  ylab('Available mOR (nmol)') +
  xlab('Bioavailability of Naloxone (%)') +
  heatmap_theme # + theme_gray(base_size = 14)

smosePlot <- (smosePlotAUC) / (smosePlotSurvival) + plot_layout(guides = 'collect')

dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave(filename = 'Fig8.png', plot = smosePlot, width = 14, height = 15)
setwd(paste0(root))

