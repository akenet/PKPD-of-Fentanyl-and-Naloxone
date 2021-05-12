# Systems Pharmacology and Personalized Medicine
# Amy van Ee, Lydia Fozo, Adam Kenet, Shiker Nair
# May 2021
# Final Project

# FIGURE 3 CODE -- Fentanyl and Naloxone

#### PACKAGES ####

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')

################################################################################
root=getwd()

for(i in 1:4){
  # Load the data that was saved from MATLAB
  setwd(paste0(root,'/Data/f5n3'))
  y_data <- readMat('Y_(f5n3).mat')
  t_data <- readMat('T_(f5n3).mat')
  
  # Convert the list from readMat() to a data frame
  y_data <- as.data.frame(y_data)
  t_data <- as.data.frame(t_data)
  
  # Combine times into one data frame
  dat <- cbind(t_data, y_data)
  
  # Change column names to be more informative
  names(dat) <- c('Time', 'Free_Fentanyl_Blood', 'Free_Fentanyl_Body', 'Free_Fentanyl_Brain',
                  'Free_Naloxone_Blood', 'Free_Naloxone_Body', 'Free_Naloxone_Brain',
                  'Free_mOR_Blood', 'mOR_Fentanyl_Brain', 'mOR_Naloxone_Brain',
                  'Cleared_Fentanyl', 'Cleared_Naloxone')
  
  # add mOR receptor occupancy
  mOR_conc <- 7.9365
  dat<- mutate(dat, Fentanyl_Receptor_Occupancy = mOR_Fentanyl_Brain/mOR_conc)
  dat<- mutate(dat, Naloxone_Receptor_Occupancy = mOR_Naloxone_Brain/mOR_conc)
  
  # Set themeing
  my_theme <- theme_classic() +
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(margin = margin(10,0,0,0)),
          axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.text.x = element_text(margin = margin(5,0,0,0)),
          axis.text.y = element_text(margin = margin(0,5,0,0)))
  
  # Create plot
  cols <- c('Naloxone'='#00bfc4','Fentanyl'='#f8766d')
  
  if(i==1){
    plot1 <- ggplot(data = NULL) +
      geom_line(data = dat, aes(x = Time, y = Free_Fentanyl_Blood, color='Fentanyl'), size = 1.5) +
      geom_line(data = dat, aes(x = Time, y = Free_Naloxone_Blood, color='Naloxone'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,200,30)) +
      xlab('Time (min)') +
      ylab('[Drug] (nM)') +
      labs(title = 'Blood', color = 'Drug') +
      my_theme
  }
  else if(i==2){
    plot2 <- ggplot(data = NULL) +
      geom_line(data = dat, aes(x = Time, y = Free_Fentanyl_Body, color='Fentanyl'), size = 1.5) +
      geom_line(data = dat, aes(x = Time, y = Free_Naloxone_Body, color='Naloxone'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,200,30)) +
      xlab('Time (min)') +
      ylab('[Drug] (nM)') +
      labs(title = 'Body', color = 'Drug') +
      my_theme
  }
  else if(i==3){
    plot3 <- ggplot(data = NULL) +
      geom_line(data = dat, aes(x = Time, y = Free_Fentanyl_Brain, color='Fentanyl'), size = 1.5) +
      geom_line(data = dat, aes(x = Time, y = Free_Naloxone_Brain, color='Naloxone'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,200,30)) +
      xlab('Time (min)') +
      ylab('[Free Drug] (nM)') +
      labs(title = 'Brain', color = 'Drug') +
      my_theme
  }
  else if(i==4){
    plot4 <- ggplot(data = NULL) +
      geom_line(data = dat, aes(x = Time, y = Fentanyl_Receptor_Occupancy, color='Fentanyl'), size = 1.5) +
      geom_line(data = dat, aes(x = Time, y = Naloxone_Receptor_Occupancy, color='Naloxone'), size = 1.5) +
      scale_x_continuous(breaks = seq(0,200,30)) +
      xlab('Time (min)') +
      ylab('Fraction of \u03BCOR Bound to Drug') +
      labs(title = '\u03BCOR Binding', color = 'Drug') +
      my_theme
  }
}
myPlot <- (plot1 | plot2) / (plot3 | plot4) + plot_layout(guides = 'collect')


# Save plot
dir.create(paste0(root,'/Plots'), showWarnings = FALSE)
setwd(paste0(root,'/Plots'))
ggsave('Fig3.png', plot = myPlot, width = 16, height = 12)
setwd(root)