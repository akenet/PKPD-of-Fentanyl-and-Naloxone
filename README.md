# PKPD-of-Fentanyl-and-Naloxone


# List of Files (17)
* fentanyl_naloxone_driver.m
* fentanyl_naloxone_fxnx.m
* fentanyl_naloxone_missed.m
* fentanyl_naloxone_survive.m
* fentanyl_naloxone_main.m
* fentanyl_naloxone_eqns.m
* fentanyl_naloxone_vis_fig2.R
* fentanyl_naloxone_vis_fig3.R
* fentanyl_naloxone_vis_fig4.R
* fentanyl_naloxone_vis_fig5.R
* fentanyl_naloxone_vis_fig6.R
* fentanyl_naloxone_vis_fig7.R
* fentanyl_naloxone_vis_fig8.R
* fentanyl_naloxone_vis_fig9.R
* fentanyl_naloxone_interactive.Rmd
* interactive_survival.m
* interactive_missed_dose.m


# Order to Run Files
### 1. Generate Data:
- fentanyl_naloxone_driver.m
### 2. Create Plots from Data:
- fentanyl_naloxone_vis_fig2.R
- fentanyl_naloxone_vis_fig3.R
- fentanyl_naloxone_vis_fig4.R
- fentanyl_naloxone_vis_fig5.R
- fentanyl_naloxone_vis_fig6.R
- fentanyl_naloxone_vis_fig7.R
- fentanyl_naloxone_vis_fig8.R
- fentanyl_naloxone_vis_fig9.R
### 3. Interactive Visualization:
- fentanyl_naloxone_interactive.Rmd






# File Structure

### 1. Generate Data:
* **(RUN)** _fentanyl_naloxone_driver.m_
    + fentanyl_naloxone_fxnx.m
         + fentanyl_naloxone_main.m
             + fentanyl_naloxone_eqns.m
    + fentanyl_naloxone_missed.m
         + fentanyl_naloxone_main.m
             + fentanyl_naloxone_eqns.m
    + fentanyl_naloxone_survive.m
         + fentanyl_naloxone_main.m
             + fentanyl_naloxone_eqns.m

### 2. Create Plots from Data:
* **(RUN)** _fentanyl_naloxone_vis_fig2.R_
* **(RUN)** _fentanyl_naloxone_vis_fig3.R_
* **(RUN)** _fentanyl_naloxone_vis_fig4.R_
* **(RUN)** _fentanyl_naloxone_vis_fig5.R_
* **(RUN)** _fentanyl_naloxone_vis_fig6.R_
* **(RUN)** _fentanyl_naloxone_vis_fig7.R_
* **(RUN)** _fentanyl_naloxone_vis_fig8.R_
* **(RUN)** _fentanyl_naloxone_vis_fig9.R_


### 3. Interactive Visualization:
* **(RUN)** _fentanyl_naloxone_interactive.Rmd_
    + interactive_survival.m
        + fentanyl_naloxone_survive.m
            + fentanyl_naloxone_main.m
                + fentanyl_naloxone_eqns.m
    + interactive_missed_dose.m
        + fentanyl_naloxone_missed.m
            + fentanyl_naloxone_main.m
                + fentanyl_naloxone_eqns.m 




# Raw Data Locations
* Data generated from _fentanyl_naloxone_driver.m_ will be saved in the directory _root/Data/scenario_, where "root" is the current working directory and "scenario" is which simulation is running.

* Plots generated from the _fentanyl_naloxone_vis_figX.R_ files will be saved in the directory _root/Plots/scenario_, where "root" is the current working directory and "scenario" is which simulation is running.

* Data generated from _fentanyl_naloxone_interactive.Rmd_ will be saved in the current working directory

# Software Requiremnts
* MATLAB
* R
    - Packages needed include: flexdashboard, tidyverse, dplyr, plotly, shiny, R.matlab, patchwork, matlabr
    
