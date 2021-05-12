# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 4 CODE -- Missed Dose

#### PACKAGES ####

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')

################################################################################
root=getwd()

for(i in 1:4){
  # Load the data that was saved from MATLAB

  # base case
  Tbase <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_base/T_(missed_dose_base).mat')))
  Ybase <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_base/Y_(missed_dose_base).mat')))
  md_base <- cbind(Tbase, Ybase)
  
  # case 1
  T1 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_1/T_(missed_dose_1).mat')))
  Y1 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_1/Y_(missed_dose_1).mat')))
  md1 <- cbind(T1, Y1)
  
  # case 2
  T2 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_2/T_(missed_dose_2).mat')))
  Y2 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_2/Y_(missed_dose_2).mat')))
  md2 <- cbind(T2, Y2)
  
  # case 3
  T3 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_3/T_(missed_dose_3).mat')))
  Y3 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_3/Y_(missed_dose_3).mat')))
  md3 <- cbind(T3, Y3)
  
  # case 4
  T4 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_4/T_(missed_dose_4).mat')))
  Y4 <- as.data.frame(readMat(paste0(root,'/Data/missed_dose_4/Y_(missed_dose_4).mat')))
  md4 <- cbind(T4, Y4)
  
  # Change column names to be more informative
  cols <- c('Time', 'Free_Fentanyl_Blood', 'Free_Fentanyl_Body', 'Free_Fentanyl_Brain',
            'Free_Naloxone_Blood', 'Free_Naloxone_Body', 'Free_Naloxone_Brain',
            'Free_mOR_Blood', 'mOR_Fentanyl_Brain', 'mOR_Naloxone_Brain',
            'Cleared_Fentanyl', 'Cleared_Naloxone')
  
  names(md1) <- cols
  names(md2) <- cols
  names(md3) <- cols
  names(md4) <- cols
  names(md_base) <- cols
  
  # add mOR receptor occupancy
  mOR_conc <- 7.9365
  md1 <- mutate(md1, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  md2 <- mutate(md2, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  md3 <- mutate(md3, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  md4 <- mutate(md4, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  md_base <- mutate(md_base, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)

  # Set themeing
  my_theme <- theme_classic() +
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(margin = margin(10,0,0,0)),
          axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.text.x = element_text(margin = margin(5,0,0,0)),
          axis.text.y = element_text(margin = margin(0,5,0,0)))
  
  # Create plot
  cols <- c('Normal'='#000000','20% Late'='#f8766d','40% Late'='#b79f00','60% Late'='#00ba38','80% Late'='#00bfc4')
  
  if(i==1){
    plot1 <- ggplot(data = NULL) +
      geom_line(data = md_base, aes(x = Time, y = Free_Fentanyl_Blood, color='Normal'), size = 1.5) +
      geom_line(data = md1, aes(x = Time, y = Free_Fentanyl_Blood, color='20% Late'), size = 1.5) +
      geom_line(data = md2, aes(x = Time, y = Free_Fentanyl_Blood, color='40% Late'), size = 1.5) +
      geom_line(data = md3, aes(x = Time, y = Free_Fentanyl_Blood, color='60% Late'), size = 1.5) +
      geom_line(data = md4, aes(x = Time, y = Free_Fentanyl_Blood, color='80% Late'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,180,30)) +
      xlab('Time (min)') +
      ylab('[Fentanyl] (nM)') +
      labs(title = 'Blood', color = 'Case') +
      my_theme
  }
  else if(i==2){
    plot2 <- ggplot(data = NULL) +
      geom_line(data = md_base, aes(x = Time, y = Free_Fentanyl_Body, color='Normal'), size = 1.5) +
      geom_line(data = md1, aes(x = Time, y = Free_Fentanyl_Body, color='20% Late'), size = 1.5) +
      geom_line(data = md2, aes(x = Time, y = Free_Fentanyl_Body, color='40% Late'), size = 1.5) +
      geom_line(data = md3, aes(x = Time, y = Free_Fentanyl_Body, color='60% Late'), size = 1.5) +
      geom_line(data = md4, aes(x = Time, y = Free_Fentanyl_Body, color='80% Late'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,180,30)) +
      xlab('Time (min)') +
      ylab('[Fentanyl] (nM)') +
      labs(title = 'Body', color = 'Case') +
      my_theme
  }
  else if(i==3){
    plot3 <- ggplot(data = NULL) +
      geom_line(data = md_base, aes(x = Time, y = Free_Fentanyl_Brain, color='Normal'), size = 1.5) +
      geom_line(data = md1, aes(x = Time, y = Free_Fentanyl_Brain, color='20% Late'), size = 1.5) +
      geom_line(data = md2, aes(x = Time, y = Free_Fentanyl_Brain, color='40% Late'), size = 1.5) +
      geom_line(data = md3, aes(x = Time, y = Free_Fentanyl_Brain, color='60% Late'), size = 1.5) +
      geom_line(data = md4, aes(x = Time, y = Free_Fentanyl_Brain, color='80% Late'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,180,30)) +
      xlab('Time (min)') +
      ylab('[Free Fentanyl] (nM)') +
      labs(title = 'Brain', color = 'Case') +
      my_theme
  }
  else if(i==4){
    plot4 <- ggplot(data = NULL) +
      geom_line(data = md_base, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='Normal'), size = 1.5) +
      geom_line(data = md1, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='20% Late'), size = 1.5) +
      geom_line(data = md2, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='40% Late'), size = 1.5) +
      geom_line(data = md3, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='60% Late'), size = 1.5) +
      geom_line(data = md4, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='80% Late'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,180,30)) +
      xlab('Time (min)') +
      ylab('Fraction of \u03BCOR Bound to Fentanyl') +
      labs(title = '\u03BCOR Binding', color = 'Case') +
      my_theme
  }
}
myPlot <- (plot1 | plot2) / (plot3 | plot4) + plot_layout(guides = 'collect')


# Save plot
dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave('Fig4.png', plot = myPlot, width = 16, height = 12)
setwd(root)