#NOTES: Things to discuss
# -- Code assumes no seasons so will not work correctly with a seasonal model do 
#    we need this capacity for any species we may want to include?
#
# -- No recruitment variability included at the moment. Is this important to add or 
#    just more noise to confuse the issue. Should we wait for reviewer feedback and 
#    just add it if they as for it? I don't think it will have any practical impact.
#
# -- Current approach uses fixed future F values for the simulation with approximates
#    a best case scenario where F_target stays the same and the fishery OFL is constantly 
#    being updated to achieve this (such as through accurate interim assessments). I 
#    think adding catch based removals will probably just exaggerate the impacts and 
#    distract from the red tide effect by allowing claims that future assessments
#    would correct for the impacts.
#
# -- Which species should we try to implement this for (Red Grouper, Gag, ???)
#
# -- What outputs do we want to record (Landings, SSB, ???)

library(r4ss)
#Source setup file that should be named local.setup so that it will be 
#ignored by github tracking. 
#
local.setup.location <- "C:/Users/Nathan/Documents/GitHub/Red_tide_benchmarks/local.setup"
source(local.setup.location)

#Source in the SEFSC projections function
source(projection_script)

#Copy files to the simulation folder 
if(dir.exists(file.path(working_dir))){
  unlink(file.path(working_dir), recursive = TRUE)
}
dir.create(file.path(working_dir))
dir.create(file.path(working_dir,"Base"))
temp.files <- list.files(path=file.path(assessment_dir))
file.copy(from = file.path(assessment_dir,temp.files), to = file.path(working_dir,"Base",temp.files))

#Read forecast file to get N_forecast years and base file for overwriting other runs
forecast_base <- r4ss::SS_readforecast(file=file.path(working_dir,"Base","forecast.ss"))
base_output <- r4ss::SS_output(file.path(working_dir,"Base"),covar = FALSE)

#Projection red tide values
rt_proj_ave <- sort(seq(0,0.1,0.01))

#True red tide averages
rt_mean <- sort(c(0.01,0.03,0.06))

#How many random red tide replicates to run
n_rand_reps <- 500

#Set seed to allow replication of results
global.seed <- 1234
set.seed(global.seed)
rand_offset <- 0 #offset to avoid using same seed as previous runs
rand_seed <- floor(runif((n_rand_reps+rand_offset),100000,9999999))[(rand_offset+1):(n_rand_reps+rand_offset)]

#Identify the fleet associated with red tide
rt_fleet <- 5

#Identify fleets to include in landings calculations
landings_fleets <- 1:4

fleet_landings_cols <- grep("retain(B)",colnames(base_output$timeseries),fixed=TRUE)[landings_fleets]

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

#Index to track row for filling results data
index_row <- 1


#First loop over the red tide rate included in projections as this will only 
#need the projections to be calculated once.
for(i in seq_along(rt_proj_ave)){
  #remove exisiting base folders if found
  proj_dir <- file.path(working_dir,paste0("rtproj_",i))
  if(dir.exists(proj_dir)){
    unlink(proj_dir, recursive = TRUE)
  }
  #Create new base folders and copy over the original model files
  dir.create(proj_dir)
  dir.create(file.path(proj_dir,"Base"))
  temp.files <- list.files(path=file.path(working_dir,"Base"))
  file.copy(from = file.path(working_dir,"Base",temp.files), to = file.path(proj_dir,"Base",temp.files))
  
  #Adjust the redtide values for the base projection and rerun the projections
  #to estimate OFL.
  #For now I'm leaving out ABC as I think it distracts from the intent of the 
  #simulation because ABC is supposed to account for unknown uncertainty not 
  #offset an avoidable bias such as this.
  #This uses the average red tide rate in every projection year
  forecast_base$ForeCatch[forecast_base$ForeCatch$Fleet==rt_fleet,4] <- rep(rt_proj_ave[i],forecast_base$Nforecastyrs)
  
  r4ss::SS_writeforecast(mylist=forecast_base,dir=file.path(proj_dir,"Base"),overwrite=TRUE)

  base_proj <- run.projections(file.path(proj_dir,"Base")) 
  
  #Loop over all mean red tide level scenarios 
  for(j in seq_along(rt_mean)){
    #remove exisiting mean red tide folders if found
    rt_dir <- file.path(proj_dir,paste0("rt_mean_",j))
    if(dir.exists(rt_dir)){
      unlink(rt_dir, recursive = TRUE)
    }
    #Create new folder for each mean red tide level and copy over the original model files
    dir.create(rt_dir)
    
    #Loop over all random red tide sequences 
    for(k in 1:n_rand_reps){
      #reset random seed for each random replicate seeds will be replicated across
      #projected red tide levels and mean red tide levels
      
      set.seed(rand_seed[k])
      #Create folders for each random sequence
      dir.create(file.path(rt_dir,k))
      temp.files <- list.files(path=file.path(proj_dir,"Base","OFL_target"))
      file.copy(from = file.path(proj_dir,"Base","OFL_target",temp.files), to = file.path(rt_dir,k,temp.files))
      
      #Calculate a random red tide mortality vector based on specified mean and frequency
      #draw a random number of red tide events from a uniform distribution between min and max number specified
      n_rt_events <- sample(n_rt_events_min:n_rt_events_max,1) #
      #calculate the total red tide mortality expected from the specified mean and number of projection years
      rt_total <- rt_mean[j]*forecast_base$Nforecastyrs
      #calculate random mortality rates from each event from a uniform distribution between min and max number specified
      rt_mags <- runif(n_rt_events,rt_min,rt_max)
      #rescale the red tide magnitudes so that they sum to the expected total mortality
      rt_mags <- rt_mags*(rt_total/sum(rt_mags))
      #create a zero mortality vector for all years
      rand_red_tide <- rep(0,forecast_base$Nforecastyrs)
      #randomly select years for the red tide mortality to occur and replace zero's with random mortality rates
      rand_red_tide[sample(1:forecast_base$Nforecastyrs,n_rt_events)] <- rt_mags
      
      
      #Modify forecast file to include random red tide mortality sequence
      forecast_rt <- r4ss::SS_readforecast(file=file.path(rt_dir,k,"forecast.ss")) 
      forecast_rt$ForeCatch[forecast_rt$ForeCatch$Fleet==rt_fleet,4] <- rand_red_tide
      #Write out the new forecast file and run model with new random mortality vector
      r4ss::SS_writeforecast(mylist=forecast_rt,dir=file.path(rt_dir,k),overwrite=TRUE)
      shell(paste("cd /d ",file.path(rt_dir,k)," && ss -nohess",sep=""))
      
      #Read in results and save values of interest for analysis
      run_output <- r4ss::SS_output(dir=file.path(rt_dir,k),covar = FALSE)
      
      spr_series <- run_output$sprseries
      
      time_series <- run_output$timeseries
      time_series_virg <- time_series[time_series$Era=="VIRG",]
      time_series <- time_series[time_series$Era!="VIRG" & time_series$Era!="INIT",]
      years <- unique(time_series$Yr)
      for(i in seq_along(years)){
        time_series_sub <- time_series[time_series$Yr==years[i],,drop=FALSE]
        spr_series_sub <- spr_series[spr_series$Yr==years[i],,drop=FALSE]
        results_landings[index_row,i] <- sum(time_series_sub[,fleet_landings_cols])
        results_SSB[index_row,i] <- sum(time_series_sub[,'SpawnBio'])
        results_SPR[index_row,i] <- sum(spr_series_sub[,'SPR'])
        results_dep[index_row,i] <- sum(spr_series_sub[,'Deplete'])
      }
      index_row <- index_row+1
    }
  }
}

all_results <- list()
all_results[[1]] <- results_landings
all_results[[2]] <- results_SSB
all_results[[3]] <- results_SPR
all_results[[4]] <- results_dep

save(all_results,file=save_file)

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
    
    summary_index <- summary_index + 1
  }
} 


plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings),max(results_landings)))
for(i in seq_along(results_landings_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
    if(results_summary_setup[i,"rt_projected"]<results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_landings_summary_mean[i,],col="red")
    }else if(results_summary_setup[i,"rt_projected"]>results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_landings_summary_mean[i,],col="blue")
    }else{
      lines(x=years,y=results_landings_summary_mean[i,],col="green")
    }
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB),max(results_SSB)))
for(i in seq_along(results_SSB_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
    if(results_summary_setup[i,"rt_projected"]<results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_SSB_summary_mean[i,],col="red")
    }else if(results_summary_setup[i,"rt_projected"]>results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_SSB_summary_mean[i,],col="blue")
    }else{
      lines(x=years,y=results_SSB_summary_mean[i,],col="green")
    }
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)))
for(i in seq_along(results_SPR_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
    if(results_summary_setup[i,"rt_projected"]<results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_SPR_summary_mean[i,],col="red")
    }else if(results_summary_setup[i,"rt_projected"]>results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_SPR_summary_mean[i,],col="blue")
    }else{
      lines(x=years,y=results_SPR_summary_mean[i,],col="green")
    }
  }
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)))
for(i in seq_along(results_dep_summary_mean[,1]))
{
  if(results_summary_setup[i,"rt_mean"]==0.06){
    if(results_summary_setup[i,"rt_projected"]<results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_dep_summary_mean[i,],col="red")
    }else if(results_summary_setup[i,"rt_projected"]>results_summary_setup[i,"rt_mean"]){
      lines(x=years,y=results_dep_summary_mean[i,],col="blue")
    }else{
      lines(x=years,y=results_dep_summary_mean[i,],col="green")
    }
  }
}




plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_landings),max(results_landings)))
for(i in seq_along(results_landings[,1]))
{
  lines(x=years,y=results_landings[i,])
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SSB),max(results_SSB)))
for(i in seq_along(results_SSB[,1]))
{
  lines(x=years,y=results_SSB[i,])
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_SPR),max(results_SPR)))
for(i in seq_along(results_SPR[,1]))
{
  lines(x=years,y=results_SPR[i,])
}

plot(x=NA,y=NA,xlim=c(min(years),max(years)),ylim=c(min(results_dep),max(results_dep)))
for(i in seq_along(results_dep[,1]))
{
  lines(x=years,y=results_dep[i,])
}

