# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 5 CODE -- Survival

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
  Tbase <- as.data.frame(readMat(paste0(root,'/Data/s2_add_naloxone_base/T_(s2_add_naloxone_base).mat')))
  Ybase <- as.data.frame(readMat(paste0(root,'/Data/s2_add_naloxone_base/Y_(s2_add_naloxone_base).mat')))
  s_base <- cbind(Tbase, Ybase)
  
  # case 1
  T1 <- as.data.frame(readMat(paste0(root,'/Data/s2_add_naloxone/T_(s2_add_naloxone).mat')))
  Y1 <- as.data.frame(readMat(paste0(root,'/Data/s2_add_naloxone/Y_(s2_add_naloxone).mat')))
  s1 <- cbind(T1, Y1)
  
  # case 2
  T2 <- as.data.frame(readMat(paste0(root,'/Data/s3_inc_delay_half_naloxone/T_(s3_inc_delay_half_naloxone).mat')))
  Y2 <- as.data.frame(readMat(paste0(root,'/Data/s3_inc_delay_half_naloxone/Y_(s3_inc_delay_half_naloxone).mat')))
  s2 <- cbind(T2, Y2)
  
  # case 3
  T3 <- as.data.frame(readMat(paste0(root,'/Data/s4_double_naloxone/T_(s4_double_naloxone).mat')))
  Y3 <- as.data.frame(readMat(paste0(root,'/Data/s4_double_naloxone/Y_(s4_double_naloxone).mat')))
  s3 <- cbind(T3, Y3)
  
  # Change column names to be more informative
  cols <- c('Time', 'Free_Fentanyl_Blood', 'Free_Fentanyl_Body', 'Free_Fentanyl_Brain',
            'Free_Naloxone_Blood', 'Free_Naloxone_Body', 'Free_Naloxone_Brain',
            'Free_mOR_Blood', 'mOR_Fentanyl_Brain', 'mOR_Naloxone_Brain',
            'Cleared_Fentanyl', 'Cleared_Naloxone')
  
  names(s1) <- cols
  names(s2) <- cols
  names(s3) <- cols
  names(s_base) <- cols
  
  # add mOR receptor occupancy
  mOR_conc <- 7.9365
  s1 <- mutate(s1, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  s2 <- mutate(s2, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  s3 <- mutate(s3, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  s_base <- mutate(s_base, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  
  # Set themeing
  my_theme <- theme_classic() +
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(margin = margin(10,0,0,0)),
          axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.text.x = element_text(margin = margin(5,0,0,0)),
          axis.text.y = element_text(margin = margin(0,5,0,0)))
  
  # Create plot
  cols <- c('No Naloxone'='#000000','Simulation 1'='#f8766d','Simulation 2'='#b79f00','Simulation 3'='#00ba38','80% Late'='#00bfc4')
  
  if(i==1){
    plot1 <- ggplot(data = NULL) +
      geom_line(data = s_base, aes(x = Time, y = Free_Fentanyl_Blood, color='No Naloxone'), size = 1.5) +
      geom_line(data = s1, aes(x = Time, y = Free_Fentanyl_Blood, color='Simulation 1'), size = 1.5) +
      geom_line(data = s2, aes(x = Time, y = Free_Fentanyl_Blood, color='Simulation 2'), size = 1.5) +
      geom_line(data = s3, aes(x = Time, y = Free_Fentanyl_Blood, color='Simulation 3'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,30,5)) +
      xlab('Time (min)') +
      ylab('[Fentanyl] (nM)') +
      labs(title = 'Blood', color = 'Legend') +
      my_theme
  }
  else if(i==2){
    plot2 <- ggplot(data = NULL) +
      geom_line(data = s_base, aes(x = Time, y = Free_Fentanyl_Body, color='No Naloxone'), size = 1.5) +
      geom_line(data = s1, aes(x = Time, y = Free_Fentanyl_Body, color='Simulation 1'), size = 1.5) +
      geom_line(data = s2, aes(x = Time, y = Free_Fentanyl_Body, color='Simulation 2'), size = 1.5) +
      geom_line(data = s3, aes(x = Time, y = Free_Fentanyl_Body, color='Simulation 3'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,30,5)) +
      xlab('Time (min)') +
      ylab('[Fentanyl] (nM)') +
      labs(title = 'Body', color = 'Legend') +
      my_theme
  }
  else if(i==3){
    plot3 <- ggplot(data = NULL) +
      geom_line(data = s_base, aes(x = Time, y = Free_Fentanyl_Brain, color='No Naloxone'), size = 1.5) +
      geom_line(data = s1, aes(x = Time, y = Free_Fentanyl_Brain, color='Simulation 1'), size = 1.5) +
      geom_line(data = s2, aes(x = Time, y = Free_Fentanyl_Brain, color='Simulation 2'), size = 1.5) +
      geom_line(data = s3, aes(x = Time, y = Free_Fentanyl_Brain, color='Simulation 3'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,30,5)) +
      xlab('Time (min)') +
      ylab('[Free Fentanyl] (nM)') +
      labs(title = 'Brain', color = 'Legend') +
      my_theme
  }
  else if(i==4){
    plot4 <- ggplot(data = NULL) +
      geom_line(data = s_base, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='No Naloxone'), size = 1.5) +
      geom_line(data = s1, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='Simulation 1'), size = 1.5) +
      geom_line(data = s2, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='Simulation 2'), size = 1.5) +
      geom_line(data = s3, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='Simulation 3'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,30,5)) +
      xlab('Time (min)') +
      ylab('Fraction of \u03BCOR Bound to Fentanyl') +
      labs(title = '\u03BCOR Binding', color = 'Legend') +
      my_theme
  }
}
myPlot <- (plot1 | plot2) / (plot3 | plot4) + plot_layout(guides = 'collect')


# Save plot
dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave('Fig5.png', plot = myPlot, width = 16, height = 12)
setwd(root)