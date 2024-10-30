#NOTES: Things to discuss
# -- Code assumes no seasons so will not work correctly with a seasonal model do 
#    we need this capacity for any species we may want to include?
#     - Gulf assessment models do not yet include seasons
#
# -- No recruitment variability included at the moment. Is this important to add or 
#    just more noise to confuse the issue. Should we wait for reviewer feedback and 
#    just add it if they as for it? I don't think it will have any practical impact.
#     - Let's stick with simple first, and build in this complexity if/as needed
#
# -- Current approach uses fixed future F values for the simulation with approximates
#    a best case scenario where F_target stays the same and the fishery OFL is constantly 
#    being updated to achieve this (such as through accurate interim assessments). I 
#    think adding catch based removals will probably just exaggerate the impacts and 
#    distract from the red tide effect by allowing claims that future assessments
#    would correct for the impacts.
#
# -- Which species should we try to implement this for (Red Grouper, Gag, ???)
#     - Red and gag grouper only ones where we estimate red tide mortality
#
# -- What outputs do we want to record (Landings, SSB, ???)
#     - Added total biomass, red tide kills (biomass), exploitation rate, recruitment, stock status

library(r4ss)
library(RColorBrewer)
#Source setup file that should be named local.setup so that it will be 
#ignored by github tracking. 
# options(max.print = 100000000) #remove limit which messed up forecast file in a few base runs
#
DIR<-"C://Users/skyler.sagarese/Desktop/RT/RGR/" #Specify RGR or GAG
local.setup.location <- paste0(DIR,"local.setup.txt")
#local.setup.location <- "C:/Users/skyler.sagarese/Desktop/RT/RGR/local.setup.txt"
source(local.setup.location)

#Source in the SEFSC projections function
source(projection_script)

# #Copy files to the simulation folder
# if(dir.exists(file.path(working_dir))){
#   unlink(file.path(working_dir), recursive = TRUE)
# }
# dir.create(file.path(working_dir))
# dir.create(file.path(working_dir,"Base"))
# temp.files <- list.files(path=file.path(assessment_dir))
# file.copy(from = file.path(assessment_dir,temp.files), to = file.path(working_dir,"Base",temp.files))

#Read forecast file to get N_forecast years and base file for overwriting other runs
forecast_base <- r4ss::SS_readforecast(file=file.path(working_dir,"Base","forecast.ss"))
base_output <- r4ss::SS_output(file.path(working_dir,"Base"),covar = FALSE)

#Projection red tide values
rt_proj_ave <- sort(seq(0,0.1,0.01))

#True red tide averages
rt_mean <- sort(c(0.01,0.03,0.06)) #RGR range
#rt_mean <- sort(c(0.01,0.08,0.16)) #GAG range

#How many random red tide replicates to run
n_rand_reps <- 500

#How many years of known catch are entered in the projections
#Red tide mortality will not be added in these years for projections or 
#randomization.
N_fixed_years <- 0 #RGR does not error out, so no reason to subset years
#N_fixed_years <- 4 #GAG popn crashes in interim years, need to turn RTM off

#Set seed to allow replication of results
global.seed <- 1234
set.seed(global.seed)
rand_offset <- 0 #offset to avoid using same seed as previous runs
rand_seed <- floor(runif((n_rand_reps+rand_offset),100000,9999999))[(rand_offset+1):(n_rand_reps+rand_offset)]

#Identify the fleet associated with red tide
rt_fleet <- 5 #Red Grouper
#rt_fleet <- 6 #Gag Grouper

#Identify fleets to include in landings calculations
landings_fleets <- 1:4 #Red Grouper
#landings_fleets <- 1:5 #Gag Grouper

fleet_landings_cols <- grep("retain(B)",colnames(base_output$timeseries),fixed=TRUE)[landings_fleets]
fleet_dead_cols <- grep("dead(B)",colnames(base_output$timeseries),fixed=TRUE)[landings_fleets]

#Setup the random red tide mortality vector details
#Set the range for the number of red tide events in the projection period of 100 years
n_rt_events_min <- 5 #The minimum number of red tide events during the projection period
n_rt_events_max <- 20 #The maximum number of red tide events during the projection period

#Set the relative range for red tide in a single year these values will be rescaled 
#in each simulation so the total red tide mortality is always sums to the target mean
rt_min <- 0.1 #Relative value of the minimum red tide in a single year 
rt_max <- 0.4 #Relative value of the maximum red tide in a single year

#Set up output matrices for storing values of interest
#Data frame to track the iteration settings for each row of the results for indexing
results_setting <- data.frame(rt_projected=c(sort(rep(rt_proj_ave,length(rt_mean)*n_rand_reps))),
                              rt_mean=c(rep(sort(rep(rt_mean,n_rand_reps)),length(rt_proj_ave))),
                              replicate=c(rep(1:n_rand_reps,length(rt_mean)*length(rt_proj_ave))))
#Achieved OFL landings
results_landings <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Achieved SSB
results_SSB <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Target SPR
results_SPR <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Target depletion
results_dep <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Achieved Recruitment
results_recr <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Achieved Total Biomass
results_tbio <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Achieved Red Tide Kill (biomass)
results_RTkillbio <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))
#Achieved Exploitation Rate
results_Fexp <- matrix(data=NA,nrow=(length(rt_proj_ave)*length(rt_mean)*n_rand_reps),ncol=(forecast_base$Nforecastyrs+base_output$endyr-base_output$startyr+1))

#Index to track row for filling results data
index_row <- 1

# #First loop over the red tide rate included in projections as this will only
# #need the projections to be calculated once.
# for(i in seq_along(rt_proj_ave)){
#   #remove exisiting base folders if found
#   proj_dir <- file.path(working_dir,paste0("rtproj_",i))
#   if(dir.exists(proj_dir)){
#     unlink(proj_dir, recursive = TRUE)
#   }
#   #Create new base folders and copy over the original model files
#   dir.create(proj_dir)
#   dir.create(file.path(proj_dir,"Base"))
#   temp.files <- list.files(path=file.path(working_dir,"Base"))
#   file.copy(from = file.path(working_dir,"Base",temp.files), to = file.path(proj_dir,"Base",temp.files))
# 
#   #Adjust the redtide values for the base projection and rerun the projections
#   #to estimate OFL.
#   #For now I'm leaving out ABC as I think it distracts from the intent of the
#   #simulation because ABC is supposed to account for unknown uncertainty not
#   #offset an avoidable bias such as this.
#   #This uses the average red tide rate in every projection year
#   forecast_base$ForeCatch[forecast_base$ForeCatch$Fleet==rt_fleet & forecast_base$ForeCatch$Year>(base_output$endyr+N_fixed_years),4] <- rep(rt_proj_ave[i],(forecast_base$Nforecastyrs-N_fixed_years))
# 
#   r4ss::SS_writeforecast(mylist=forecast_base,dir=file.path(proj_dir,"Base"),overwrite=TRUE)
# 
#   base_proj <- run.projections(file.path(proj_dir,"Base"))
# 
#   #Loop over all mean red tide level scenarios
#   for(j in seq_along(rt_mean)){
#     #remove exisiting mean red tide folders if found
#     rt_dir <- file.path(proj_dir,paste0("rt_mean_",j))
#     if(dir.exists(rt_dir)){
#       unlink(rt_dir, recursive = TRUE)
#     }
#     #Create new folder for each mean red tide level and copy over the original model files
#     dir.create(rt_dir)
# 
#     #Loop over all random red tide sequences
#     for(k in 1:n_rand_reps){
#       #reset random seed for each random replicate seeds will be replicated across
#       #projected red tide levels and mean red tide levels
# 
#       set.seed(rand_seed[k])
#       #Create folders for each random sequence
#       dir.create(file.path(rt_dir,k))
#       temp.files <- list.files(path=file.path(proj_dir,"Base","OFL_target"))
#       file.copy(from = file.path(proj_dir,"Base","OFL_target",temp.files), to = file.path(rt_dir,k,temp.files))
# 
#       #Calculate a random red tide mortality vector based on specified mean and frequency
#       #draw a random number of red tide events from a uniform distribution between min and max number specified
#       n_rt_events <- sample(n_rt_events_min:n_rt_events_max,1) #
#       #calculate the total red tide mortality expected from the specified mean and number of projection years
#       rt_total <- rt_mean[j]*(forecast_base$Nforecastyrs-N_fixed_years)
#       #calculate random mortality rates from each event from a uniform distribution between min and max number specified
#       rt_mags <- runif(n_rt_events,rt_min,rt_max)
#       #rescale the red tide magnitudes so that they sum to the expected total mortality
#       rt_mags <- rt_mags*(rt_total/sum(rt_mags))
#       #create a zero mortality vector for all years
#       rand_red_tide <- rep(0,(forecast_base$Nforecastyrs))
#       #randomly select years for the red tide mortality to occur and replace zero's with random mortality rates
#       rand_red_tide[sample((N_fixed_years+1):forecast_base$Nforecastyrs,n_rt_events)] <- rt_mags
# 
# 
#       #Modify forecast file to include random red tide mortality sequence
#       forecast_rt <- r4ss::SS_readforecast(file=file.path(rt_dir,k,"forecast.ss"))
#       forecast_rt$ForeCatch[forecast_rt$ForeCatch$Fleet==rt_fleet & forecast_rt$ForeCatch$Year>(base_output$endyr+N_fixed_years),4] <- rand_red_tide[(N_fixed_years+1):forecast_base$Nforecastyrs]
#       #Write out the new forecast file and run model with new random mortality vector
#       r4ss::SS_writeforecast(mylist=forecast_rt,dir=file.path(rt_dir,k),overwrite=TRUE)
#       shell(paste("cd /d ",file.path(rt_dir,k)," && ss -nohess",sep=""))
# 
#       #Read in results and save values of interest for analysis
#       run_output <- r4ss::SS_output(dir=file.path(rt_dir,k),covar = FALSE)
# 
#       spr_series <- run_output$sprseries
# 
#       time_series <- run_output$timeseries
#       time_series_virg <- time_series[time_series$Era=="VIRG",]
#       time_series <- time_series[time_series$Era!="VIRG" & time_series$Era!="INIT",]
#       years <- unique(time_series$Yr)
#       for(i in seq_along(years)){
#         time_series_sub <- time_series[time_series$Yr==years[i],,drop=FALSE]
#         spr_series_sub <- spr_series[spr_series$Yr==years[i],,drop=FALSE]
#         results_landings[index_row,i] <- sum(time_series_sub[,fleet_landings_cols])
#         results_SSB[index_row,i] <- sum(time_series_sub[,'SpawnBio'])
#         results_SPR[index_row,i] <- sum(spr_series_sub[,'SPR'])
#         results_dep[index_row,i] <- sum(spr_series_sub[,'Deplete'])
#         results_recr[index_row,i] <- sum(time_series_sub[,'Recruit_0'])
#         results_tbio[index_row,i] <- sum(spr_series_sub[,'Bio_all.1'])
#         results_RTkillbio[index_row,i] <- sum(time_series_sub[,'dead(B):_5']) #RGR
# #        results_RTkillbio[index_row,i] <- sum(time_series_sub[,'dead(B):_6']) #GAG
#         results_Fexp[index_row,i] <- sum(spr_series_sub[,'F_report']) #note: this includes red tide mortality historically, but not in projections
#       }
#       index_row <- index_row+1
#     }
#   }
# }
# 
# all_results <- list()
# all_results[[1]] <- results_landings
# all_results[[2]] <- results_SSB
# all_results[[3]] <- results_SPR
# all_results[[4]] <- results_dep
# all_results[[5]] <- results_recr
# all_results[[6]] <- results_tbio
# all_results[[7]] <- results_RTkillbio
# all_results[[8]] <- results_Fexp
# 
# save(all_results,file=save_file)
load("results")
head(all_results)

results_landings<-all_results[[1]]
results_SSB<- all_results[[2]]
results_SPR<- all_results[[3]]
results_dep <- all_results[[4]]
results_recr <- all_results[[5]]
results_tbio <- all_results[[6]]
results_RTkillbio <- all_results[[7]]
results_Fexp <- all_results[[8]] 
# Note: exploitation in report file excludes red tide in forecast, but includes in timeseries
# We want to remove red tide from the exploitation time series for consistency
# 2005 (column 20)
RT2005<-base_output$timeseries[base_output$timeseries$Yr==2005,]
RT2005spr<-base_output$sprseries[base_output$sprseries$Yr==2005,]
FishKillBio2005<-sum(RT2005[,fleet_dead_cols])/RT2005spr$Bio_Smry.1
results_Fexp[,20]<-FishKillBio2005

# 2014 (column 29)
RT2014<-base_output$timeseries[base_output$timeseries$Yr==2014,]
RT2014spr<-base_output$sprseries[base_output$sprseries$Yr==2014,]
FishKillBio2014<-sum(RT2014[,fleet_dead_cols])/RT2014spr$Bio_Smry.1
results_Fexp[,29]<-FishKillBio2014

#Summarize results for display
summary_index <- 1

results_summary_setup <- data.frame(rt_projected=c(sort(rep(rt_proj_ave,length(rt_mean)))),
                            rt_mean=c(rep(sort(rt_mean),length(rt_proj_ave))))

results_landings_summary_mean <- results_landings[1:(length(rt_proj_ave)*length(rt_mean)),]
results_landings_summary_sd <- results_landings[1:(length(rt_proj_ave)*length(rt_mean)),]

results_SSB_summary_mean <- results_SSB[1:(length(rt_proj_ave)*length(rt_mean)),]
results_SSB_summary_sd <- results_SSB[1:(length(rt_proj_ave)*length(rt_mean)),]

results_SPR_summary_mean <- results_SPR[1:(length(rt_proj_ave)*length(rt_mean)),]
results_SPR_summary_sd <- results_SPR[1:(length(rt_proj_ave)*length(rt_mean)),]

results_dep_summary_mean <- results_dep[1:(length(rt_proj_ave)*length(rt_mean)),]
results_dep_summary_sd <- results_dep[1:(length(rt_proj_ave)*length(rt_mean)),]

results_recr_summary_mean <- results_recr[1:(length(rt_proj_ave)*length(rt_mean)),]
results_recr_summary_sd <- results_recr[1:(length(rt_proj_ave)*length(rt_mean)),]

results_tbio_summary_mean <- results_tbio[1:(length(rt_proj_ave)*length(rt_mean)),]
results_tbio_summary_sd <- results_tbio[1:(length(rt_proj_ave)*length(rt_mean)),]

results_RTkillbio_summary_mean <- results_RTkillbio[1:(length(rt_proj_ave)*length(rt_mean)),]
results_RTkillbio_summary_sd <- results_RTkillbio[1:(length(rt_proj_ave)*length(rt_mean)),]

results_Fexp_summary_mean <- results_Fexp[1:(length(rt_proj_ave)*length(rt_mean)),]
results_Fexp_summary_sd <- results_Fexp[1:(length(rt_proj_ave)*length(rt_mean)),]

for(i in seq_along(rt_proj_ave)){
  for(j in seq_along(rt_mean)){
    rows <- which(results_setting[,"rt_projected"]==rt_proj_ave[i] & results_setting[,"rt_mean"]==rt_mean[j])
    
    results_landings_summary_mean[summary_index,] <- apply(results_landings[rows,],2,mean) 
    results_landings_summary_sd[summary_index,] <- apply(results_landings[rows,],2,sd) 
    
    results_SSB_summary_mean[summary_index,] <- apply(results_SSB[rows,],2,mean) 
    results_SSB_summary_sd[summary_index,] <- apply(results_SSB[rows,],2,sd) 
    
    results_SPR_summary_mean[summary_index,] <- apply(results_SPR[rows,],2,mean) 
    results_SPR_summary_sd[summary_index,] <- apply(results_SPR[rows,],2,sd) 
    
    results_dep_summary_mean[summary_index,] <- apply(results_dep[rows,],2,mean) 
    results_dep_summary_sd[summary_index,] <- apply(results_dep[rows,],2,sd) 

    results_recr_summary_mean[summary_index,] <- apply(results_recr[rows,],2,mean) 
    results_recr_summary_sd[summary_index,] <- apply(results_recr[rows,],2,sd) 

    results_tbio_summary_mean[summary_index,] <- apply(results_tbio[rows,],2,mean) 
    results_tbio_summary_sd[summary_index,] <- apply(results_tbio[rows,],2,sd) 
    
    results_RTkillbio_summary_mean[summary_index,] <- apply(results_RTkillbio[rows,],2,mean) 
    results_RTkillbio_summary_sd[summary_index,] <- apply(results_RTkillbio[rows,],2,sd) 

    results_Fexp_summary_mean[summary_index,] <- apply(results_Fexp[rows,],2,mean) 
    results_Fexp_summary_sd[summary_index,] <- apply(results_Fexp[rows,],2,sd) 

    summary_index <- summary_index + 1
  }
} 


############
# Landings #
############

years<-c(base_output$startyr:(base_output$endyr+100))
cols<-brewer.pal(11, "RdYlBu") 
#Want this repeated for baseline simulations across true values
#Reverse colors so blue is good and red is bad ()
cols<-rev(c(rep(cols[1],3),rep(cols[2],3),rep(cols[3],3),rep(cols[4],3),rep(cols[5],3),
        rep(cols[6],3),rep(cols[7],3),rep(cols[8],3),rep(cols[9],3),rep(cols[10],3),
        rep(cols[11],3)))

jpeg(paste0(DIR,"Landings.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     ylab="",xlab="",las=1,xaxt="n",cex.axis=1.5)
text(2050,max(results_landings)/1000*0.9,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,xaxt="n",cex.axis=1.5)
text(2050,max(results_landings)/1000*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,max(results_landings)/1000*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
mtext("Landings (1000s metric tons)",side=2,outer=T)
dev.off()


##########################
# Spawning Stock Biomass #
##########################

jpeg(paste0(DIR,"SSB.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     ylab="",xlab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_SSB)/1000000*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_SSB)/1000000*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,max(results_SSB)/1000000*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}
mtext("Spawning Stock Biomass (Relative Number of Eggs)",side=2,outer=T,line=0.7)
dev.off()


#############
# SPR Ratio #
#############

jpeg(paste0(DIR,"SPR Ratio.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     ylab="SPR Ratio",xlab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_SPR),"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_SPR),"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,max(results_SPR),"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}
mtext("SPR Ratio",side=2,outer=T,line=0.7)
dev.off()


#############
# SSB Ratio #
#############

jpeg(paste0(DIR,"SSBratio.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="",xlab="",las=1,cex.axis=1.5,xaxt="n")
abline(h=0.3)
text(2050,max(results_dep)*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
abline(h=0.3)
text(2050,max(results_dep)*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1,cex.axis=1.5)
abline(h=0.3)
text(2050,max(results_dep)*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}
mtext("SSB Ratio (SSB/SSBunfished)",side=2,outer=T,line=0.7)
dev.off()


###############
# Recruitment #
###############

jpeg(paste0(DIR,"Recr.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     ylab="",xlab="",cex.axis=1.5,xaxt="n",las=1)
text(2050,max(results_recr)/1000*0.9,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",cex.axis=1.5,xaxt="n",las=1)
text(2050,max(results_recr)/1000*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",cex.axis=1.5,las=1)
text(2050,max(results_recr)/1000*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}
mtext("Recruitment (Millions of Fish)",side=2,outer=T,line=0.7)
dev.off()


#####################
# Exploitation Rate #
#####################

jpeg(paste0(DIR,"F.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     ylab="",xlab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,0.4*0.9,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,0.4*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,0.4*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}
mtext("Exploitation Rate (biomass)",side=2,outer=T,line=0.7)
dev.off()


#################
# Total Biomass #
#################

jpeg(paste0(DIR,"tbio.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     ylab="",xlab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_tbio)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,max(results_tbio)/1000*0.9,"True Future RTM mean=0.03",cex=1.5)
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,max(results_tbio)/1000*0.9,"True Future RTM mean=0.06",cex=1.5)
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}
mtext("Total Biomass (1000s metric tons)",side=2,outer=T,line=0.7)
dev.off()


#########################
# Red Tide Kill Biomass #
#########################

jpeg(paste0(DIR,"RTkillbio.jpeg"),res=300,height=2400,width=1600)
par(mfrow=c(3,1),mar=c(0.2,2,0.2,0.9),oma=c(2,2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     #plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     ylab="",xlab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,7,"True Future RTM mean=0.01",cex=1.5)
#text(2050,max(results_RTkillbio)/1000*0.9,"True Future RTM mean=0.01",cex=1.5,col="green")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

#plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     xlab="",ylab="",las=1,cex.axis=1.5,xaxt="n")
text(2050,7,"True Future RTM mean=0.03",cex=1.5)
#text(2050,max(results_RTkillbio)/1000*0.9,"True Future RTM mean=0.03",cex=1.5,col="green")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

#plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     xlab="",ylab="",las=1,cex.axis=1.5)
text(2050,7,"True Future RTM mean=0.06",cex=1.5)
#text(2050,max(results_RTkillbio)/1000*0.9,"True Future RTM mean=0.06",cex=1.5,col="green")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}
mtext("Red Tide Kill (1000s metric tons)",side=2,outer=T,line=0.7)
dev.off()


############################
# Results Combined Summary #
############################

jpeg(paste0(DIR,"Summary_Final.jpeg"),res=300,height=2700,width=2400)
par(mfrow=c(4,3),mar=c(2,4,0.5,1),oma=c(0.1,0.1,0.1,0.1))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     ylab="Landings (1000s metric tons)",xlab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.01",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.03",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.06",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
text(2060,5.25,"Baseline RTM")
legend(2020,5.2,legend=c("0","0.01","0.02","0.03","0.04","0.05"),lty=1,col=cols[seq(1,31,3)[1:6]],bty="n")
legend(2060,5.2,legend=c("0.06","0.07","0.08","0.09","0.1"),lty=1,col=cols[seq(1,31,3)[7:11]],bty="n")

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="SSB Ratio (SSB/SSBunfished)",xlab="",las=1,xaxt="n")
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1,xaxt="n")
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1,xaxt="n")
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     ylab="Total Biomass (1000s metric tons)",xlab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,7),
     ylab="Red Tide Kill (Biomass, 1000s metric tons)",xlab="",las=1)
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,7),
     xlab="",ylab="",las=1)
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,7),
     xlab="",ylab="",las=1)
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}
dev.off()


jpeg(paste0(DIR,"Summary.jpeg"),res=300,height=2400,width=2000)
par(mfrow=c(4,3),mar=c(2,4,2,1),oma=c(0.1,0.1,0.1,0.1))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     ylab="Landings (1000s metric tons)",xlab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.01",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.03",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_landings)/1000*0.9,"True Future RTM mean=0.06",cex=1)
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_landings_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     ylab="Spawning Stock Biomass (Relative # Eggs)",xlab="",las=1,xaxt="n")
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_SSB_summary_mean[i,]/1000000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     ylab="SPR Ratio",xlab="",las=1,xaxt="n")
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_SPR_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="SSB Ratio (SSB/SSBunfished)",xlab="",las=1)
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1)
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="",ylab="",las=1)
abline(h=0.3)
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_dep_summary_mean[i,],col=cols[i])
  }
}
dev.off()

jpeg(paste0(DIR,"Summary_2.jpeg"),res=300,height=2700,width=2000)
par(mfrow=c(4,3),mar=c(2,4,2,1),oma=c(0.1,0.1,0.1,0.1))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     ylab="Recruitment (Millions of Fish)",xlab="",las=1,xaxt="n")
text(2053,max(results_recr)/1000*0.9,"True Future RTM mean=0.01",cex=1)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_recr)/1000*0.9,"True Future RTM mean=0.03",cex=1)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2053,max(results_recr)/1000*0.9,"True Future RTM mean=0.06",cex=1)
for(i in seq_along(results_recr_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_recr_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     ylab="Total Biomass (1000s metric tons)",xlab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_tbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_tbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     ylab="Red Tide Kill (Biomass, 1000s metric tons)",xlab="",las=1,xaxt="n")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,8),
     xlab="",ylab="",las=1,xaxt="n")
for(i in seq_along(results_RTkillbio_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_RTkillbio_summary_mean[i,]/1000,col=cols[i])
  }
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     ylab="Exploitation Rate (biomass)",xlab="",las=1)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.01){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="",ylab="",las=1)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.03){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="",ylab="",las=1)
for(i in seq_along(results_Fexp_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
      lines(x=years,y=results_Fexp_summary_mean[i,],col=cols[i])
  }
}
dev.off()


###############################
# Individual Simulation Plots #
###############################

jpeg(paste0(DIR,"Sims_all.jpeg"),res=300,height=2400,width=2000)
par(mfrow=c(4,2),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
for(i in seq_along(results_landings[,1]/1000))
{
  lines(x=years,y=results_landings[i,]/1000)
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="Spawning Stock Biomass (Relative # Eggs)",las=1)
for(i in seq_along(results_SSB[,1]/1000000))
{
  lines(x=years,y=results_SSB[i,]/1000000)
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="SPR Ratio",las=1,xaxt="n")
for(i in seq_along(results_SPR[,1]))
{
  lines(x=years,y=results_SPR[i,])
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1,xaxt="n")
for(i in seq_along(results_dep[,1]))
{
  lines(x=years,y=results_dep[i,])
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="Year",ylab="Recruitment (Millions of Fish)",las=1,xaxt="n")
for(i in seq_along(results_recr[,1]))
{
  lines(x=years,y=results_recr[i,]/1000)
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
for(i in seq_along(results_tbio[,1]))
{
  lines(x=years,y=results_tbio[i,]/1000)
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
for(i in seq_along(results_RTkillbio[,1]))
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="Year",ylab="Exploitation Rate (biomass)",las=1)
for(i in seq_along(results_Fexp[,1]))
{
  lines(x=years,y=results_Fexp[i,])
}
dev.off()

####################################################################################
# Define simulation sequence (which runs correspond to mean = 0.01, 0.03, 0.06, etc)
## Folders run through Baseline 0 for mean 0.01, mean 0.03, and mean 0.06, then to Baseline 0.01

Start<-seq(1,16500,500)
End<-seq(500,16500,500)

Sims<-as.data.frame(cbind(Start,End))
Sims$Mean<-rep(seq(1:3),11)
Sims$Comb<-paste0(Start,":",End)
Sims$Comb[Sims$Mean==1]
Mean1<-c(1:500,1501:2000,3001:3500,4501:5000,6001:6500,7501:8000,9001:9500,10501:11000,
         12001:12500,13501:14000,15001:15500)
Mean1_sub<-c(1:5,1501:1505,3001:3005,4501:4505,6001:6005,7501:7505,9001:9005,10501:10505,
             12001:12005,13501:13505,15001:15005)
Mean2<-Mean1+500
Mean2_sub<-Mean1_sub+500
Mean3<-Mean2+500
Mean3_sub<-Mean2_sub+500

cols<-brewer.pal(4, "RdYlBu")
jpeg(paste0(DIR,"Deterministic.jpeg"),res=300,height=2800,width=2000)
par(mfrow=c(4,1),mar=c(1,4,0.1,0.5),oma=c(2,0.2,0.2,0.2))
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
#for(i in seq_along(results_landings[,1]/1000))
#text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
  lines(x=years,y=results_landings[501,]/1000,col=cols[1],lwd=2)
  lines(x=years,y=results_landings[2001,]/1000,col=cols[2],lwd=2)
  lines(x=years,y=results_landings[5001,]/1000,col=cols[3],lwd=2)
  lines(x=years,y=results_landings[9501,]/1000,col=cols[4],lwd=2)
abline(v=base_output$endyr)
text(2003,6.25,"Assessment Period",cex=1.1)
text(2028,6.2,"Projection Period",cex=1.1)
legend("topright",legend=c("Red tide mortality = 0",
                           "Red tide mortality = 0.01",
                           "Red tide mortality = 0.03",
                           "Red tide mortality = 0.06"),
       col=cols,bty="n",lty=c(1,1,1,1),lwd=c(2,2,2,2),cex=1.2)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),0.6),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1,xaxt="n")
  lines(x=years,y=results_dep[501,],col=cols[1],lwd=2)
  lines(x=years,y=results_dep[2001,],col=cols[2],lwd=2)
  lines(x=years,y=results_dep[5001,],col=cols[3],lwd=2)
  lines(x=years,y=results_dep[9501,],col=cols[4],lwd=2)
  abline(v=base_output$endyr)
  abline(h=0.3,lwd=2)
text(2004,0.27,"Target",cex=1.1)
abline(h=0.15,lwd=2,lty=3)
text(2020,0.12,"Minimum Stock Size Threshold",cex=1.1)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,35),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
lines(x=years,y=results_tbio[501,]/1000,col=cols[1],lwd=2)
lines(x=years,y=results_tbio[2001,]/1000,col=cols[2],lwd=2)
lines(x=years,y=results_tbio[5001,]/1000,col=cols[3],lwd=2)
lines(x=years,y=results_tbio[9501,]/1000,col=cols[4],lwd=2)
abline(v=base_output$endyr)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,7),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
lines(x=years,y=results_RTkillbio[9501,]/1000,col=cols[4],lwd=2)
lines(x=years,y=results_RTkillbio[5001,]/1000,col=cols[3],lwd=2)
lines(x=years,y=results_RTkillbio[2001,]/1000,col=cols[2],lwd=2)
lines(x=years,y=results_RTkillbio[501,]/1000,col=cols[1],lwd=2)
abline(v=base_output$endyr)
dev.off()

jpeg(paste0(DIR,"Deterministic_RT.jpeg"),res=300,height=2000,width=2000)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,7),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
lines(x=years,y=results_RTkillbio[9501,]/1000,col=cols[4],lwd=2)
lines(x=years,y=results_RTkillbio[5001,]/1000,col=cols[3],lwd=2)
lines(x=years,y=results_RTkillbio[2001,]/1000,col=cols[2],lwd=2)
lines(x=years,y=results_RTkillbio[501,]/1000,col=cols[1],lwd=2)
abline(v=base_output$endyr)
dev.off()

jpeg(paste0(DIR,"Sims_byTrueMean_1.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
#for(i in seq_along(results_landings[,1]/1000))
for(i in Mean1)
{
  lines(x=years,y=results_landings[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2)
{
  lines(x=years,y=results_landings[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3)
{
  lines(x=years,y=results_landings[i,]/1000)
}

#SSB
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="Spawning Stock Biomass (Relative # Eggs)",las=1)
for(i in Mean1)
{
  lines(x=years,y=results_SSB[i,]/1000000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
for(i in Mean2)
{
  lines(x=years,y=results_SSB[i,]/1000000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
for(i in Mean3)
{
  lines(x=years,y=results_SSB[i,]/1000000)
}

#SPR
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="SPR Ratio",las=1,xaxt="n")
for(i in Mean1)
{
  lines(x=years,y=results_SPR[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean2)
{
  lines(x=years,y=results_SPR[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean3)
{
  lines(x=years,y=results_SPR[i,])
}

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1)
for(i in Mean1)
{
  lines(x=years,y=results_dep[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1)
for(i in Mean2)
{
  lines(x=years,y=results_dep[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1)
for(i in Mean3)
{
  lines(x=years,y=results_dep[i,])
}
dev.off()

cols<-brewer.pal(11, "RdYlBu") 
#Want this repeated for baseline simulations across true values
#Reverse colors so blue is good and red is bad ()
cols<-rev(c(rep(cols[1],1500),rep(cols[2],1500),rep(cols[3],1500),rep(cols[4],1500),rep(cols[5],1500),
            rep(cols[6],1500),rep(cols[7],1500),rep(cols[8],1500),rep(cols[9],1500),rep(cols[10],1500),
            rep(cols[11],1500)))

jpeg(paste0(DIR,"Sims_byTrueMean_Subset_1.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in Mean1_sub) {
    lines(x=years,y=results_landings[i,]/1000,col=cols[i])
  }
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2_sub)
{
  lines(x=years,y=results_landings[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3_sub)
{
  lines(x=years,y=results_landings[i,]/1000,col=cols[i])
}

#SSB
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="Spawning Stock Biomass (Relative # Eggs)",las=1)
for(i in Mean1_sub)
{
  lines(x=years,y=results_SSB[i,]/1000000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
for(i in Mean2_sub)
{
  lines(x=years,y=results_SSB[i,]/1000000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
for(i in Mean3_sub)
{
  lines(x=years,y=results_SSB[i,]/1000000,col=cols[i])
}

#SPR
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="SPR Ratio",las=1,xaxt="n")
for(i in Mean1_sub)
{
  lines(x=years,y=results_SPR[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean2_sub)
{
  lines(x=years,y=results_SPR[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean3_sub)
{
  lines(x=years,y=results_SPR[i,],col=cols[i])
}

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1)
for(i in Mean1_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1)
for(i in Mean2_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1)
for(i in Mean3_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}
dev.off()

rc <- brewer.pal(2,"Dark2")
jpeg(paste0(DIR,"SimsPerc_byTrueMean_1.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.5)))
text(2080,5,"95th Percentile",col=rc[1])
text(2080,4,"50th Percentile",col=rc[2])

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.5)))

#SSB
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="Spawning Stock Biomass (Relative # Eggs)",las=1)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean1,]/1000000,2,quantile,probs=c(0.975)),
                              rev(apply(results_SSB[Mean1,]/1000000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean1,]/1000000,2,quantile,probs=c(0.75)),
                              rev(apply(results_SSB[Mean1,]/1000000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SSB[Mean1,]/1000000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean2,]/1000000,2,quantile,probs=c(0.975)),
                              rev(apply(results_SSB[Mean2,]/1000000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean2,]/1000000,2,quantile,probs=c(0.75)),
                              rev(apply(results_SSB[Mean2,]/1000000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SSB[Mean2,]/1000000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB)/1000000,max(results_SSB)/1000000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean3,]/1000000,2,quantile,probs=c(0.975)),
                              rev(apply(results_SSB[Mean3,]/1000000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SSB[Mean3,]/1000000,2,quantile,probs=c(0.75)),
                              rev(apply(results_SSB[Mean3,]/1000000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SSB[Mean3,]/1000000,2,quantile,probs=c(0.5)))

#SPR
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xaxt="n",ylab="SPR Ratio",las=1)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean1,],2,quantile,probs=c(0.975)),
                              rev(apply(results_SPR[Mean1,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean1,],2,quantile,probs=c(0.75)),
                              rev(apply(results_SPR[Mean1,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SPR[Mean1,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean2,],2,quantile,probs=c(0.975)),
                              rev(apply(results_SPR[Mean2,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean2,],2,quantile,probs=c(0.75)),
                              rev(apply(results_SPR[Mean2,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SPR[Mean2,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean3,],2,quantile,probs=c(0.975)),
                              rev(apply(results_SPR[Mean3,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_SPR[Mean3,],2,quantile,probs=c(0.75)),
                              rev(apply(results_SPR[Mean3,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_SPR[Mean3,],2,quantile,probs=c(0.5)))

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="SSB Ratio (SSB/SSBunfished)",las=1)
polygon(c(years,rev(years)),c(apply(results_dep[Mean1,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean1,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean1,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean1,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean1,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_dep[Mean2,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean2,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean2,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean2,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean2,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_dep[Mean3,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean3,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean3,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean3,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean3,],2,quantile,probs=c(0.5)))
dev.off()


jpeg(paste0(DIR,"Sims_byTrueMean_2.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
# Recruitment
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="Year",ylab="Recruitment (Millions of Fish)",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in Mean1)
{
  lines(x=years,y=results_recr[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2)
{
  lines(x=years,y=results_recr[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3)
{
  lines(x=years,y=results_recr[i,]/1000)
}

# Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
for(i in Mean1)
{
  lines(x=years,y=results_tbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2)
{
  lines(x=years,y=results_tbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3)
{
  lines(x=years,y=results_tbio[i,]/1000)
}

# Red Tide Kill (Biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1,xaxt="n")
for(i in Mean1)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}

#Exploitation Rate
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="Year",ylab="Exploitation Rate (biomass)",las=1)
for(i in Mean1)
{
  lines(x=years,y=results_Fexp[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="",ylab="",las=1)
for(i in Mean2)
{
  lines(x=years,y=results_Fexp[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     xlab="",ylab="",las=1)
for(i in Mean3)
{
  lines(x=years,y=results_Fexp[i,])
}
dev.off()

jpeg(paste0(DIR,"Sims_byTrueMean_Subset_2.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
# Recruitment
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="Year",ylab="Recruitment (Millions of Fish)",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
for(i in Mean1_sub)
{
  lines(x=years,y=results_recr[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2_sub)
{
  lines(x=years,y=results_recr[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xlab="",ylab="",las=1,xaxt="n")
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3_sub)
{
  lines(x=years,y=results_recr[i,]/1000,col=cols[i])
}

# Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
for(i in Mean1_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}

# Red Tide Kill (Biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1,xaxt="n")
for(i in Mean1_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}

#Exploitation Rate
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     xlab="Year",ylab="Exploitation Rate (biomass)",las=1)
for(i in Mean1_sub)
{
  lines(x=years,y=results_Fexp[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     xlab="",ylab="",las=1)
for(i in Mean2_sub)
{
  lines(x=years,y=results_Fexp[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),0.4),
     xlab="",ylab="",las=1)
for(i in Mean3_sub)
{
  lines(x=years,y=results_Fexp[i,],col=cols[i])
}
dev.off()

jpeg(paste0(DIR,"SimsPerc_byTrueMean_2.jpeg"),res=300,height=2400,width=3000)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Recruitment
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xaxt="n",ylab="Recruitment (Millions of Fish)",las=1)
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_recr[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_recr[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_recr[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_recr[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_recr[Mean1,]/1000,2,quantile,probs=c(0.5)))
text(2080,100,"95th Percentile",col=rc[1])
text(2080,90,"50th Percentile",col=rc[2])

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_recr[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_recr[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_recr[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_recr[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_recr[Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_recr)/1000,max(results_recr)/1000),
     xaxt="n",ylab="",las=1)
text(2050,max(results_recr)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_recr[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_recr[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_recr[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_recr[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_recr[Mean3,]/1000,2,quantile,probs=c(0.5)))

#Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="Total Biomass (1000s metric tons)",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.5)))

#Red Tide Kill (biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     xaxt="n",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean1,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean3,]/1000,2,quantile,probs=c(0.5)))

#Exploitation Rate 
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     ylab="Exploitation Rate (biomass)",las=1)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean1,],2,quantile,probs=c(0.975)),
                              rev(apply(results_Fexp[Mean1,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean1,],2,quantile,probs=c(0.75)),
                              rev(apply(results_Fexp[Mean1,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_Fexp[Mean1,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean2,],2,quantile,probs=c(0.975)),
                              rev(apply(results_Fexp[Mean2,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean2,],2,quantile,probs=c(0.75)),
                              rev(apply(results_Fexp[Mean2,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_Fexp[Mean2,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_Fexp),max(results_Fexp)),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean3,],2,quantile,probs=c(0.975)),
                              rev(apply(results_Fexp[Mean3,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_Fexp[Mean3,],2,quantile,probs=c(0.75)),
                              rev(apply(results_Fexp[Mean3,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_Fexp[Mean3,],2,quantile,probs=c(0.5)))
dev.off()

jpeg(paste0(DIR,"Sims_byTrueMean_Final.jpeg"),res=300,height=3000,width=2800)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
#for(i in seq_along(results_landings[,1]/1000))
for(i in Mean1)
{
  lines(x=years,y=results_landings[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2)
{
  lines(x=years,y=results_landings[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3)
{
  lines(x=years,y=results_landings[i,]/1000)
}

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1,xaxt="n")
for(i in Mean1)
{
  lines(x=years,y=results_dep[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean2)
{
  lines(x=years,y=results_dep[i,])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean3)
{
  lines(x=years,y=results_dep[i,])
}

# Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
for(i in Mean1)
{
  lines(x=years,y=results_tbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2)
{
  lines(x=years,y=results_tbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3)
{
  lines(x=years,y=results_tbio[i,]/1000)
}

# Red Tide Kill (Biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
for(i in Mean1)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1)
for(i in Mean2)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1)
for(i in Mean3)
{
  lines(x=years,y=results_RTkillbio[i,]/1000)
}
dev.off()

jpeg(paste0(DIR,"Sims_byTrueMean_Subset_Final.jpeg"),res=300,height=3000,width=2800)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
#for(i in seq_along(results_landings[,1]/1000))
for(i in Mean1_sub)
{
  lines(x=years,y=results_landings[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
for(i in Mean2_sub)
{
  lines(x=years,y=results_landings[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
for(i in Mean3_sub)
{
  lines(x=years,y=results_landings[i,]/1000,col=cols[i])
}

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="SSB Ratio (SSB/SSBunfished)",las=1,xaxt="n")
for(i in Mean1_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean2_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     xlab="Year",ylab="",las=1,xaxt="n")
for(i in Mean3_sub)
{
  lines(x=years,y=results_dep[i,],col=cols[i])
}

# Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="Year",ylab="Total Biomass (1000s metric tons)",las=1,xaxt="n")
for(i in Mean1_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean2_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xlab="",ylab="",las=1,xaxt="n")
for(i in Mean3_sub)
{
  lines(x=years,y=results_tbio[i,]/1000,col=cols[i])
}

# Red Tide Kill (Biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="Year",ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
for(i in Mean1_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1)
for(i in Mean2_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio)/1000,max(results_RTkillbio)/1000),
     xlab="",ylab="",las=1)
for(i in Mean3_sub)
{
  lines(x=years,y=results_RTkillbio[i,]/1000,col=cols[i])
}
dev.off()

rc <- brewer.pal(2,"Dark2")
jpeg(paste0(DIR,"SimsPerc_byTrueMean_Final.jpeg"),res=300,height=3000,width=2800)
par(mfrow=c(4,3),mar=c(2,4,1,1),oma=c(0.2,0.2,0.2,0.2))
#Landings
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="Landings (1000s metric tons)",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.01",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean1,]/1000,2,quantile,probs=c(0.5)))
text(2080,5,"95th Percentile",col=rc[1])
text(2080,4.5,"50th Percentile",col=rc[2])

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.03",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings)/1000,max(results_landings)/1000),
     xaxt="n",ylab="",las=1)
text(2053,max(results_landings)/1000*0.95,"True Future RTM mean=0.06",cex=1.5)
polygon(c(years,rev(years)),c(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_landings[Mean3,]/1000,2,quantile,probs=c(0.5)))

#SSB Ratio
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="SSB Ratio (SSB/SSBunfished)",las=1,xaxt="n")
polygon(c(years,rev(years)),c(apply(results_dep[Mean1,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean1,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean1,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean1,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean1,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="",las=1,xaxt="n")
polygon(c(years,rev(years)),c(apply(results_dep[Mean2,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean2,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean2,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean2,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean2,],2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)),
     ylab="",las=1,xaxt="n")
polygon(c(years,rev(years)),c(apply(results_dep[Mean3,],2,quantile,probs=c(0.975)),
                              rev(apply(results_dep[Mean3,],2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_dep[Mean3,],2,quantile,probs=c(0.75)),
                              rev(apply(results_dep[Mean3,],2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_dep[Mean3,],2,quantile,probs=c(0.5)))

#Total Biomass
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="Total Biomass (1000s metric tons)",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean1,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_tbio)/1000,max(results_tbio)/1000),
     xaxt="n",ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_tbio[Mean3,]/1000,2,quantile,probs=c(0.5)))

#Red Tide Kill (biomass)
plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     ylab="Red Tide Kill (Biomass, 1000s metric tons)",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean1,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean1,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean2,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean2,]/1000,2,quantile,probs=c(0.5)))

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_RTkillbio/1000),10),
     ylab="",las=1)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.975)),
                              rev(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.025)))),
        col=adjustcolor(rc[1],alpha.f=0.2),border=NA)
polygon(c(years,rev(years)),c(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.75)),
                              rev(apply(results_RTkillbio[Mean3,]/1000,2,quantile,probs=c(0.25)))),
        col=adjustcolor(rc[2],alpha.f=0.2),border=NA)
lines(x=years,y=apply(results_RTkillbio [Mean3,]/1000,2,quantile,probs=c(0.5)))
dev.off()




##############################################
# Calculate Stock Status for each simulation #
##############################################

Bench<-read.csv("Bench_targets.csv",header=T)

MSSTfrac<-0.5

# Current SSB (2017) / SSBMSY (2117)
#Mean 0.01
Sims[Sims$Mean==1,]
SSB_1<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(SSB_1)<-c("Mean1_0","Mean1_0.01","Mean1_0.02","Mean1_0.03","Mean1_0.04",
                   "Mean1_0.05","Mean1_0.06","Mean1_0.07","Mean1_0.08","Mean1_0.09","Mean1_0.10")
SSB_1$Mean1_0<-results_SSB[1:500,32]/Bench$GeomMeanSSB[1] #bench from original projection with baseline
SSB_1$Mean1_0.01<-results_SSB[1501:2000,32]/Bench$GeomMeanSSB[2]
SSB_1$Mean1_0.02<-results_SSB[3001:3500,32]/Bench$GeomMeanSSB[3]
SSB_1$Mean1_0.03<-results_SSB[4501:5000,32]/Bench$GeomMeanSSB[4]
SSB_1$Mean1_0.04<-results_SSB[6001:6500,32]/Bench$GeomMeanSSB[5]
SSB_1$Mean1_0.05<-results_SSB[7501:8000,32]/Bench$GeomMeanSSB[6]
SSB_1$Mean1_0.06<-results_SSB[9001:9500,32]/Bench$GeomMeanSSB[7]
SSB_1$Mean1_0.07<-results_SSB[10501:11000,32]/Bench$GeomMeanSSB[8]
SSB_1$Mean1_0.08<-results_SSB[12001:12500,32]/Bench$GeomMeanSSB[9]
SSB_1$Mean1_0.09<-results_SSB[13501:14000,32]/Bench$GeomMeanSSB[10]
SSB_1$Mean1_0.10<-results_SSB[15001:15500,32]/Bench$GeomMeanSSB[11]

#Mean 0.03
Sims[Sims$Mean==2,]
SSB_2<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(SSB_2)<-c("Mean2_0","Mean2_0.01","Mean2_0.02","Mean2_0.03","Mean2_0.04",
                   "Mean2_0.05","Mean2_0.06","Mean2_0.07","Mean2_0.08","Mean2_0.09","Mean2_0.10")
SSB_2$Mean2_0<-results_SSB[501:1000,32]/Bench$GeomMeanSSB[1] #bench from original projection with baseline
SSB_2$Mean2_0.01<-results_SSB[2001:2500,32]/Bench$GeomMeanSSB[2]
SSB_2$Mean2_0.02<-results_SSB[3501:4000,32]/Bench$GeomMeanSSB[3]
SSB_2$Mean2_0.03<-results_SSB[5001:5500,32]/Bench$GeomMeanSSB[4]
SSB_2$Mean2_0.04<-results_SSB[6501:7000,32]/Bench$GeomMeanSSB[5]
SSB_2$Mean2_0.05<-results_SSB[8001:8500,32]/Bench$GeomMeanSSB[6]
SSB_2$Mean2_0.06<-results_SSB[9501:10000,32]/Bench$GeomMeanSSB[7]
SSB_2$Mean2_0.07<-results_SSB[11001:11500,32]/Bench$GeomMeanSSB[8]
SSB_2$Mean2_0.08<-results_SSB[12501:13000,32]/Bench$GeomMeanSSB[9]
SSB_2$Mean2_0.09<-results_SSB[14001:14500,32]/Bench$GeomMeanSSB[10]
SSB_2$Mean2_0.10<-results_SSB[15501:16000,32]/Bench$GeomMeanSSB[11]

#Mean 0.06
Sims[Sims$Mean==3,]
SSB_3<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(SSB_3)<-c("Mean3_0","Mean3_0.01","Mean3_0.02","Mean3_0.03","Mean3_0.04",
                   "Mean3_0.05","Mean3_0.06","Mean3_0.07","Mean3_0.08","Mean3_0.09","Mean3_0.10")
SSB_3$Mean3_0<-results_SSB[1001:1500,32]/Bench$GeomMeanSSB[1] #bench from original projection with baseline
SSB_3$Mean3_0.01<-results_SSB[2501:3000,32]/Bench$GeomMeanSSB[2]
SSB_3$Mean3_0.02<-results_SSB[4001:4500,32]/Bench$GeomMeanSSB[3]
SSB_3$Mean3_0.03<-results_SSB[5501:6000,32]/Bench$GeomMeanSSB[4]
SSB_3$Mean3_0.04<-results_SSB[7001:7500,32]/Bench$GeomMeanSSB[5]
SSB_3$Mean3_0.05<-results_SSB[8501:9000,32]/Bench$GeomMeanSSB[6]
SSB_3$Mean3_0.06<-results_SSB[10001:10500,32]/Bench$GeomMeanSSB[7]
SSB_3$Mean3_0.07<-results_SSB[11501:12000,32]/Bench$GeomMeanSSB[8]
SSB_3$Mean3_0.08<-results_SSB[13001:13500,32]/Bench$GeomMeanSSB[9]
SSB_3$Mean3_0.09<-results_SSB[14501:15000,32]/Bench$GeomMeanSSB[10]
SSB_3$Mean3_0.10<-results_SSB[16001:16500,32]/Bench$GeomMeanSSB[11]

# Current F (geom mean 2015-2017) / FMSY (2117)
#Mean 0.01
Sims[Sims$Mean==1,]
Fexp_1<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(Fexp_1)<-c("Mean1_0","Mean1_0.01","Mean1_0.02","Mean1_0.03","Mean1_0.04",
                    "Mean1_0.05","Mean1_0.06","Mean1_0.07","Mean1_0.08","Mean1_0.09","Mean1_0.10")
#values are identical across rows, so just using geom mean
Fexp_1$Mean1_0<-exp(mean(log((results_Fexp[1:500,30:32]))))/Bench$GeomMeanF[1] #bench from original projection with baseline
Fexp_1$Mean1_0.01<-exp(mean(log((results_Fexp[1501:200,30:32]))))/Bench$GeomMeanF[2]
Fexp_1$Mean1_0.02<-exp(mean(log((results_Fexp[3001:3500,30:32]))))/Bench$GeomMeanF[3]
Fexp_1$Mean1_0.03<-exp(mean(log((results_Fexp[4501:5000,30:32]))))/Bench$GeomMeanF[4]
Fexp_1$Mean1_0.04<-exp(mean(log((results_Fexp[6001:6500,30:32]))))/Bench$GeomMeanF[5]
Fexp_1$Mean1_0.05<-exp(mean(log((results_Fexp[7501:8000,30:32]))))/Bench$GeomMeanF[6]
Fexp_1$Mean1_0.06<-exp(mean(log((results_Fexp[9001:9500,30:32]))))/Bench$GeomMeanF[7]
Fexp_1$Mean1_0.07<-exp(mean(log((results_Fexp[10501:11000,30:32]))))/Bench$GeomMeanF[8]
Fexp_1$Mean1_0.08<-exp(mean(log((results_Fexp[12001:12500,30:32]))))/Bench$GeomMeanF[9]
Fexp_1$Mean1_0.09<-exp(mean(log((results_Fexp[13501:14000,30:32]))))/Bench$GeomMeanF[10]
Fexp_1$Mean1_0.10<-exp(mean(log((results_Fexp[15001:15500,30:32]))))/Bench$GeomMeanF[11]

#Mean 0.03
Sims[Sims$Mean==2,]
Fexp_2<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(Fexp_2)<-c("Mean2_0","Mean2_0.01","Mean2_0.02","Mean2_0.03","Mean2_0.04",
                    "Mean2_0.05","Mean2_0.06","Mean2_0.07","Mean2_0.08","Mean2_0.09","Mean2_0.10")
#values are identical across rows, so just using geom mean
Fexp_2$Mean2_0<-exp(mean(log((results_Fexp[501:1000,30:32]))))/Bench$GeomMeanF[1] #bench from original projection with baseline
Fexp_2$Mean2_0.01<-exp(mean(log((results_Fexp[2001:2500,30:32]))))/Bench$GeomMeanF[2]
Fexp_2$Mean2_0.02<-exp(mean(log((results_Fexp[3501:4000,30:32]))))/Bench$GeomMeanF[3]
Fexp_2$Mean2_0.03<-exp(mean(log((results_Fexp[1500:5500,30:32]))))/Bench$GeomMeanF[4]
Fexp_2$Mean2_0.04<-exp(mean(log((results_Fexp[6501:7000,30:32]))))/Bench$GeomMeanF[5]
Fexp_2$Mean2_0.05<-exp(mean(log((results_Fexp[8001:8500,30:32]))))/Bench$GeomMeanF[6]
Fexp_2$Mean2_0.06<-exp(mean(log((results_Fexp[9501:10000,30:32]))))/Bench$GeomMeanF[7]
Fexp_2$Mean2_0.07<-exp(mean(log((results_Fexp[11001:11500,30:32]))))/Bench$GeomMeanF[8]
Fexp_2$Mean2_0.08<-exp(mean(log((results_Fexp[12501:13000,30:32]))))/Bench$GeomMeanF[9]
Fexp_2$Mean2_0.09<-exp(mean(log((results_Fexp[14001:14500,30:32]))))/Bench$GeomMeanF[10]
Fexp_2$Mean2_0.10<-exp(mean(log((results_Fexp[15501:16000,30:32]))))/Bench$GeomMeanF[11]

#Mean 0.06
Sims[Sims$Mean==3,]
Fexp_3<-as.data.frame(matrix(ncol=11,nrow=500)) #500 sims vs 11 scenarios (base RTM)
colnames(Fexp_3)<-c("Mean3_0","Mean3_0.01","Mean3_0.02","Mean3_0.03","Mean3_0.04",
                    "Mean3_0.05","Mean3_0.06","Mean3_0.07","Mean3_0.08","Mean3_0.09","Mean3_0.10")
#values are identical across rows, so just using geom mean
Fexp_3$Mean3_0<-exp(mean(log((results_Fexp[1001:1500,30:32]))))/Bench$GeomMeanF[1] #bench from original projection with baseline
Fexp_3$Mean3_0.01<-exp(mean(log((results_Fexp[2501:3000,30:32]))))/Bench$GeomMeanF[2]
Fexp_3$Mean3_0.02<-exp(mean(log((results_Fexp[4001:4500,30:32]))))/Bench$GeomMeanF[3]
Fexp_3$Mean3_0.03<-exp(mean(log((results_Fexp[5501:6000,30:32]))))/Bench$GeomMeanF[4]
Fexp_3$Mean3_0.04<-exp(mean(log((results_Fexp[7001:7500,30:32]))))/Bench$GeomMeanF[5]
Fexp_3$Mean3_0.05<-exp(mean(log((results_Fexp[8501:9000,30:32]))))/Bench$GeomMeanF[6]
Fexp_3$Mean3_0.06<-exp(mean(log((results_Fexp[10001:10500,30:32]))))/Bench$GeomMeanF[7]
Fexp_3$Mean3_0.07<-exp(mean(log((results_Fexp[11501:12000,30:32]))))/Bench$GeomMeanF[8]
Fexp_3$Mean3_0.08<-exp(mean(log((results_Fexp[13001:13500,30:32]))))/Bench$GeomMeanF[9]
Fexp_3$Mean3_0.09<-exp(mean(log((results_Fexp[14501:15000,30:32]))))/Bench$GeomMeanF[10]
Fexp_3$Mean3_0.10<-exp(mean(log((results_Fexp[16001:16500,30:32]))))/Bench$GeomMeanF[11]

#All the SSB & F refs are the same
tail(SSB_1)
tail(SSB_2)
tail(SSB_3)
tail(Fexp_1)
tail(Fexp_2)
tail(Fexp_3)

Ref<-rbind(as.numeric(SSB_1[500,]),as.numeric(SSB_2[500,]),as.numeric(SSB_3[500,]),
           as.numeric(Fexp_1[500,]),as.numeric(Fexp_2[500,]),as.numeric(Fexp_3[500,]))
colnames(Ref)<-c("Base0","Base0.01","Base0.02","Base0.03","Base0.04","Base0.05",
                 "Base0.06","Base0.07","Base0.08","Base0.09","Base0.1")
rownames(Ref)<-c("SSB1","SSB2","SSB3","Fexp1","Fexp2","Fexp3")
write.csv(Ref,"StockStatus.csv")
