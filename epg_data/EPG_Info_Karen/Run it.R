# Define working directory
setwd("./epg_data/Output/") #where I want the output

library("sciplot")

# Define directory of the script
source("C:\\RProjects\\stressbiome\\epg_data\\EPG_Info_Karen\\epgrun_0.332.R")

# Run
epgrun(loc=c("C:\\RProjects\\stressbiome\\epg_data\\Raw_data\\Riv_Ctr",
             "C:\\RProjects\\stressbiome\\epg_data\\Raw_data\\Riv_SA"),
epgstart=0, epgstop=8, epgbin=TRUE, epgname="Rivera")




