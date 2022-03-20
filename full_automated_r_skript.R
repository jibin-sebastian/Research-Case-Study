###### Loading library ###### 

library(plyr)
library(dostats)
library(stringr)
library(ggplot2)
library(plotly)
library(XML)
library(methods)

# We set this parameter this way, because we work with micro seconds
options(digits.secs = 6)
###### Defining path to data ###### 

path_to_smartphone_sct <- "C:/Users/jibin/OneDrive/Documents/MS Germany/Academic/DS 3rd Sem/RCS/smartphon_sct"
path_to_smartphone_6mwt <- "C:/Users/jibin/OneDrive/Documents/MS Germany/Academic/DS 3rd Sem/RCS/smartphone_6mwt"
path_to_fitbit <- "C:/Users/jibin/OneDrive/Documents/MS Germany/Academic/DS 3rd Sem/RCS/fitbit/fitbit"
path_to_accelerometer <- "C:/Users/jibin/OneDrive/Documents/MS Germany/Academic/DS 3rd Sem/RCS/accelerometer/accelerometer"


# Create folder structure for the outputs
# In my case I created a folder with the name Graphs_Plots to store all plots
maindir <- "C:/Users/jibin/OneDrive/Documents/MS Germany/Academic/DS 3rd Sem/RCS/Graphs_Plots"

# Create Toplevel folders 
dir.create(file.path(paste0(maindir,"/", "smartphone")), showWarnings = F) # smartphone folder
dir.create(file.path(paste0(maindir,"/", "fitbit")),showWarnings = F)    # fitbit folder 
dir.create(file.path(paste0(maindir,"/", "accelerometer")),showWarnings = F) # accelerometer


subdir_smartphone_sct <- "smartphone/sct"
subdir_smartphone_6mwt <- "smartphone/6mwt"
subdir_fitbit_sct <- "fitbit/sct"
subdir_fitbit_6mwt <- "fitbit/6mwt"
subdir_accelerometer_sct <- "accelerometer/sct"
subdir_accelerometer_6mwt <- "accelerometer/6mwt"


dir.create(file.path(maindir, subdir_smartphone_sct),showWarnings = F)
dir.create(file.path(maindir, subdir_smartphone_6mwt),showWarnings = F)
dir.create(file.path(maindir, subdir_fitbit_sct),showWarnings = F)
dir.create(file.path(maindir, subdir_fitbit_6mwt),showWarnings = F)
dir.create(file.path(maindir, subdir_accelerometer_sct),showWarnings = F)
dir.create(file.path(maindir, subdir_accelerometer_6mwt),showWarnings = F)


###### Identifiying the timeframe of interest ###### 

#' Because of the difference of the document with the name title "near end of" 
#' and "start of" we will exclude them from this procedure, and observe them
#' later on. 
#' 
#' we also exclude the data 2b9ffa95 from an automized searching , because 
#' it does not have enough data points/ behave structural different from the rest
#' 

smartphone_file_to_sct <- list.files(path_to_smartphone_sct
                                     ,pattern = ".csv",
                                     full.names = TRUE)

smartphone_file_to_6mwt <- list.files(path_to_smartphone_6mwt, 
                                      pattern = ".csv",
                                      full.names =  TRUE)


# Remove the unwanted data that would be loader otherwise 

smartphone_file_to_sct <- smartphone_file_to_sct[c(seq(1,7,1),
                                                   seq(9,41,1))]

smartphone_file_to_6mwt <- smartphone_file_to_6mwt[c(seq(1,7,1),
                                                     seq(9,41,1))]


smartphone_file_to_6mwt[1]

# Fitbit files
fitbit_file <- paste0(path_to_fitbit,
                      "/",
                      substring(smartphone_file_to_6mwt,
                                first = (nchar(path_to_smartphone_6mwt)+2)))


#Accelrometer files 

accelerometer_file <- paste0(path_to_accelerometer,
                             "/",
                             substring(smartphone_file_to_6mwt,
                                       first = (nchar(path_to_smartphone_6mwt)+2), 
                                       last = (nchar(smartphone_file_to_6mwt)-4)),
                             "/")


# Create data frame to store values 

df <- data.frame()
df <- rbind(df,  c("id",
                   'Minimum_of_smartphone_6mwt',
                   'Maximum_of_smartphone_6mwt',
                   'Deviation_of_smartphone_6mwt',
                   'Minimum_of_smartphone_sct',
                   'Maximum_of_smartphone_sct',
                   'Deviation_of_smartphone_sct',
                   'Minimum_of_fitbit_6mwt',
                   'Maximum_of_fitbit_6mwt',
                   'Deviation_of_fitbit_6mwt',
                   'amount_of_confidence_zero_in_6mwt',
                   'amount_of_confidence_one_in_6mwt',
                   'amount_of_confidence_two_in_6mwt',
                   'amount_of_confidence_three_in_6mwt',
                   'Minimum_of_fitbit_sct',
                   'Maximum_of_fitbit_sct',
                   'Deviation_of_fitbit_sct',
                   'amount_of_confidence_zero_in_sct',
                   'amount_of_confidence_one_in_sct',
                   'amount_of_confidence_two_in_sct',
                   'amount_of_confidence_three_in_sct',
                   'Maximum_of_accelerometer_6mwt_split1_v1',
                   'Minimum_of_accelerometer_6mwt_split1_v1',
                   'Deviation_of_accelerometer_6mwt_split1_v1',
                   'Maximum_of_accelerometer_6mwt_split2_v1',
                   'Minimum_of_accelerometer_6mwt_split2_v1',
                   'Deviation_of_accelerometer_6mwt_split2_v1',
                   'Maximum_of_accelerometer_6mwt_split3_v1',
                   'Minimum_of_accelerometer_6mwt_split3_v1',
                   'Deviation_of_accelerometer_6mwt_split3_v1',
                   'Maximum_of_accelerometer_6mwt_split4_v1',
                   'Minimum_of_accelerometer_6mwt_split4_v1',
                   'Deviation_of_accelerometer_6mwt_split4_v1',
                   'Maximum_of_accelerometer_6mwt_split1_v2',
                   'Minimum_of_accelerometer_6mwt_split1_v2',
                   'Deviation_of_accelerometer_6mwt_split1_v2',
                   'Maximum_of_accelerometer_6mwt_split2_v2',
                   'Minimum_of_accelerometer_6mwt_split2_v2',
                   'Deviation_of_accelerometer_6mwt_split2_v2',
                   'Maximum_of_accelerometer_6mwt_split3_v2',
                   'Minimum_of_accelerometer_6mwt_split3_v2',
                   'Deviation_of_accelerometer_6mwt_split3_v2',
                   'Maximum_of_accelerometer_6mwt_split4_v2',
                   'Minimum_of_accelerometer_6mwt_split4_v2',
                   'Deviation_of_accelerometer_6mwt_split4_v2',
                   'Maximum_of_accelerometer_6mwt_split1_v3',
                   'Minimum_of_accelerometer_6mwt_split1_v3',
                   'Deviation_of_accelerometer_6mwt_split1_v3',
                   'Maximum_of_accelerometer_6mwt_split2_v3',
                   'Minimum_of_accelerometer_6mwt_split2_v3',
                   'Deviation_of_accelerometer_6mwt_split2_v3',
                   'Maximum_of_accelerometer_6mwt_split3_v3',
                   'Minimum_of_accelerometer_6mwt_split3_v3',
                   'Deviation_of_accelerometer_6mwt_split3_v3',
                   'Maximum_of_accelerometer_6mwt_split4_v3',
                   'Minimum_of_accelerometer_6mwt_split4_v3',
                   'Deviation_of_accelerometer_6mwt_split4_v3',
                   'Maximum_of_accelerometer_sct_split1_v1',
                   'Minimum_of_accelerometer_sct_split1_v1',
                   'Deviation_of_accerometer_sct_split1_V1',
                   'Maximum_of_accelerometer_sct_split2_v1',
                   'Minimum_of_accelerometer_sct_split2_v1',
                   'Deviation_of_accerometer_sct_split2_V1',
                   'Maximum_of_accelerometer_sct_split1_v2',
                   'Minimum_of_accelerometer_sct_split1_v2',
                   'Deviation_of_accerometer_sct_split1_V2',
                   'Maximum_of_accelerometer_sct_split2_v2',
                   'Minimum_of_accelerometer_sct_split2_v2',
                   'Deviation_of_accerometer_sct_split2_V2',
                   'Maximum_of_accelerometer_sct_split1_v3',
                   'Minimum_of_accelerometer_sct_split1_v3',
                   'Deviation_of_accerometer_sct_split1_V3',
                   'Maximum_of_accelerometer_sct_split2_v3',
                   'Minimum_of_accelerometer_sct_split2_v3',
                   'Deviation_of_accerometer_sct_split2_V3'
                   
))


for(j in 1:length(smartphone_file_to_sct)){
  
  
  ##### Load the file of the sct data #####
  
  tmp <-  read.csv(smartphone_file_to_sct[j], sep = ";")
  
  
  
  
  # Filter the file such that we have only look at the data of interest
  tmp <- tmp%>%filter(type_profile == "basic" & type_category == "acceleration")
  
  # Adjust the time problem (R is not good when it comes to micro seconds, so 
  # we need to help us out here with this work around)
  date_and_time <- data.frame(str_sub(tmp$time, end=-4))
  
  miliseconds <- as.numeric(str_sub(tmp$time, -3,-1))/1000
  
  time <- as.POSIXct(as.numeric(date_and_time$str_sub.tmp.time..end....4.),
                     origin="1970-01-01",
                     tz = "UTC")
  time <- time + miliseconds
  
  
  # Bind time and value to one data.frame
  data_of_interest <- data.frame(time,tmp$value)
  
  
  # To delete oulier values that are just way off from our normal distribution 
  # This would be a big problem, because we identify our 
  # 99.7% of the data that would stem from the same distribution
  data_of_interest <- data_of_interest[data_of_interest$tmp.value
                                       < (mean(data_of_interest$tmp.value)
                                          +3*sd(data_of_interest$tmp.value)),]
  
  
  
  frame_to_iterate <- which(data_of_interest$tmp.value > mean(data_of_interest$tmp.value))
  
  start <- list()
  end  <- list()
  Result <- list()
  
  for(i in frame_to_iterate){
    if((i+59) <= dim(data_of_interest)[1]){
      
      start[[i]] <- i
      end[[i]] <- i+59
      Result[[i]] <-sum(data_of_interest[i:(59+i),2]>mean(data_of_interest$tmp.value)) 
    }
    else{
      break
    }
    
    
  }
  start <- do.call(rbind.data.frame,
                   start)
  
  end <- do.call(rbind.data.frame,
                 end)
  
  Result <- do.call(rbind.data.frame,
                    Result)
  
  # Bind values together to one data.frame
  intermediate_Result <- cbind(start,
                               end,
                               Result)
  
  colnames(intermediate_Result) <- c("start",
                                     "end",
                                     "Result")
  
  
  
  
  # Find the interval that has has the most values above the mean
  
  Maximum <- max(intermediate_Result$Result, 
                 na.rm = T)
  
  helper <- na.omit(intermediate_Result[intermediate_Result$Result == Maximum, ])
  
  
  iter_for_sct <- helper[ceiling(dim(helper)[1]/2),1]
  
  Final_data_sct <- data_of_interest[iter_for_sct:(iter_for_sct+59),]
  
  # We add additional time otherwise our Graph get's hard to read
  
  
  
  Data_for_plot <- data_of_interest[(if(iter_for_sct -30 <=0){
    1
  }else{
    iter_for_sct-30
  }):((iter_for_sct+59)+30),]
  p <- ggplot(data = Data_for_plot, aes(x = time,
                                        y = tmp.value))+
    geom_point()+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_sct[c(1, dim(Final_data_sct)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_sct[j],
                           start = (nchar(path_to_smartphone_sct)+2),
                           end = -5)))+
    ylab("Value")+
    xlab("Time")
  
  # Define path where we save our output to
  setwd(paste0(maindir,
               "/",
               subdir_smartphone_sct))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_sct[j],
                                start = (nchar(path_to_smartphone_sct)+2)),
                        ".png"),
         width = 14,
         height = 10,
         units = "cm")
  
  
  ##### Load the file of the 6mwt data  ######
  
  tmp <-  read.csv(smartphone_file_to_6mwt[j], sep = ";")
  
  # Filter the file such that we have only look at the data of interest
  tmp <- tmp%>%filter(type_profile == "basic" & type_category == "acceleration")
  
  # Adjust the time problem (R is not good when it comes to micro seconds, so 
  # we need to help us out here with this work around)
  date_and_time <- data.frame(str_sub(tmp$time, end=-4))
  
  miliseconds <- as.numeric(str_sub(tmp$time, -3,-1))/1000
  
  time <- as.POSIXct(as.numeric(date_and_time$str_sub.tmp.time..end....4.),
                     origin="1970-01-01",
                     tz = "UTC")
  time <- time + miliseconds
  
  
  # Bind time and value to one data.frame
  data_of_interest <- data.frame(time,tmp$value)
  
  
  # To delete oulier values that are just way off from our normal distribution 
  # This would be a big problem, because we identify our 
  # 99.7% of the data that would stem from the same distribution
  data_of_interest <- data_of_interest[data_of_interest$tmp.value
                                       < (mean(data_of_interest$tmp.value)
                                          +3*sd(data_of_interest$tmp.value)),]
  
  
  
  frame_to_iterate <- which(data_of_interest$tmp.value > mean(data_of_interest$tmp.value))
  
  start <- list()
  end  <- list()
  Result <- list()
  
  for(i in frame_to_iterate){
    if((i+359) <= dim(data_of_interest)[1]){
      
      start[[i]] <- i
      end[[i]] <- i+359
      Result[[i]] <-sum(data_of_interest[i:(359+i),2]>mean(data_of_interest$tmp.value)) 
    }
    else{
      break
    }
    
    
  }
  start <- do.call(rbind.data.frame,
                   start)
  
  end <- do.call(rbind.data.frame,
                 end)
  
  Result <- do.call(rbind.data.frame,
                    Result)
  
  # Bind values together to one data.frame
  intermediate_Result <- cbind(start,
                               end,
                               Result)
  
  colnames(intermediate_Result) <- c("start",
                                     "end",
                                     "Result")
  
  
  
  
  # Find the interval that has has the most values above the mean
  
  Maximum <- max(intermediate_Result$Result,
                 na.rm = T)
  
  helper <- na.omit(intermediate_Result[intermediate_Result$Result == Maximum, ])
  
  
  iter_for_6mwt <- helper[ceiling(dim(helper)[1]/2),1]
  
  Final_data_6mwt <- data_of_interest[iter_for_6mwt:(iter_for_6mwt+359),]
  
  # We add additional time otherwise our Graph get's hard to read
  
  Data_for_plot <- data_of_interest[(if(iter_for_6mwt -30 <=0){
    1
  }else{
    iter_for_6mwt-30
  }):((iter_for_6mwt+359)+30),]
  
  p <- ggplot(data = Data_for_plot,
              aes(x = time,
                  y = tmp.value))+
    geom_point()+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_6mwt[c(1, dim(Final_data_6mwt)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_6mwt[j],
                           start = (nchar(path_to_smartphone_6mwt)+2),
                           end = -5)))+
    ylab("Value")+
    xlab("Time")
  
  # Define path where we save our output to
  setwd(paste0(maindir, "/", subdir_smartphone_6mwt))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_6mwt[j], 
                                start = (nchar(path_to_smartphone_6mwt)+2)),
                        ".png"),
         width = 14,
         height = 10,
         units = "cm")
  
  
  
  
  
  
  
  
  ##### Fitbit Data #####
  
  fitbit_data <- read.csv(fitbit_file[j], sep = ";")
  fitbit_data$confidence <- as.factor(fitbit_data$confidence)
  
  fitbit_data$dateTime <- as.POSIXct(as.character(fitbit_data$dateTime), 
                                     format="%m/%d/%y %H:%M:%S",
                                     tz = "UTC")
  
  
  
  fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= Final_data_6mwt[1,1] 
                             & fitbit_data$dateTime <= Final_data_6mwt[360,1],]
  
  
  for(i in seq(0,5, by =0.1)){
    if(as.numeric(difftime(fitbit_6mwt[dim(fitbit_6mwt)[1],1],fitbit_6mwt[1,1])) < 6){
      fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= (Final_data_6mwt[1,1]-i)
                                 & fitbit_data$dateTime <= (Final_data_6mwt[360,1]+i),]
      
    }else{
      break
    }
  }
  
  
  Data_to_plot <- fitbit_data[fitbit_data$dateTime >= (Final_data_6mwt[1,1]-20)
                              & fitbit_data$dateTime <= (Final_data_6mwt[360,1]+20),]
  
  
  p <- ggplot(data = Data_to_plot,
              aes(x = dateTime,
                  y = bpm, 
                  color = confidence))+
    geom_point()+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_6mwt[c(1, dim(Final_data_6mwt)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_6mwt[j],
                           start = (nchar(path_to_smartphone_6mwt)+2), 
                           end = -5)))+
    xlab("Time")
  
  # Define path where we save our output to
  setwd(paste0(maindir, "/", subdir_fitbit_6mwt))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_6mwt[j],
                                start = (nchar(path_to_smartphone_6mwt)+2)),
                        ".png"), 
         width = 14,
         height = 10,
         units = "cm")
  
  
  
  
  
  fitbit_sct <- fitbit_data[fitbit_data$dateTime >= Final_data_sct[1,1] 
                            & fitbit_data$dateTime <= Final_data_sct[60,1],]
  
  
  for(i in seq(0,5, by =0.1)){
    if(as.numeric(difftime(fitbit_sct[dim(fitbit_sct)[1],1],
                           fitbit_sct[1,1],
                           units = "mins")) < 1){
      fitbit_sct <- fitbit_data[fitbit_data$dateTime >= (Final_data_sct[1,1]-i) 
                                & fitbit_data$dateTime <= (Final_data_sct[60,1]+i),]
      
    }else{
      break
    }
  }
  
  
  Data_to_plot <- fitbit_data[fitbit_data$dateTime >= (Final_data_sct[1,1]-20) 
                              & fitbit_data$dateTime <= (Final_data_sct[60,1]+20),]
  
  
  p <- ggplot(data = Data_to_plot,
              aes(x = dateTime,
                  y = bpm, 
                  color = confidence))+
    geom_point()+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_sct[c(1, dim(Final_data_sct)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_6mwt[j],
                           start = (nchar(path_to_smartphone_6mwt)+2),
                           end = -5)))+
    xlab("Time")
  
  # Define path where we save our output to
  setwd(paste0(maindir, "/", subdir_fitbit_sct))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_6mwt[j], 
                                start = (nchar(path_to_smartphone_6mwt)+2)),
                        ".png"),
         width = 14,
         height = 10,
         units = "cm")
  
  
  
  
  
  
  
  ##### Accelerometer Data #####
  
  
  accelerometer_data <- read.csv(paste0(accelerometer_file[j],
                                        "acc.csv"),
                                 sep = ";",
                                 header = F)
  
  
  xmlfile <- xmlParse(file = paste0(accelerometer_file[j],
                                    "unisens.xml"))
  
  
  rootnode <- xmlRoot(xmlfile)
  helperlist <- xmlToList(rootnode)
  
  
  starting_point <- as.POSIXct(gsub("T",
                                    " ",
                                    helperlist$.attrs[[4]]))
  #Subtract the two hours from the CEST format to get to UTC
  starting_point <- starting_point -7200
  
  
  starting_point <- as.POSIXct(as.character(starting_point),
                               tz = "UTC")
  
  time <- seq(from = starting_point,
              to =(starting_point +
                     (dim(accelerometer_data)[1]/64)-1),
              length.out = dim(accelerometer_data)[1] )
  
  acc_with_time <- cbind(time, accelerometer_data)
  
  
  # 6mwt
  
  accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= Final_data_6mwt[1,1]
                                           & acc_with_time$time <= Final_data_6mwt[360,1],]
  
  
  for(i in seq(0,5, by =0.1)){
    if(as.numeric(difftime(accelerometer_data_6mwt[dim(accelerometer_data_6mwt)[1],1],
                           accelerometer_data_6mwt[1,1])) < 6)
    {
      accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= (Final_data_6mwt[1,1]-i) 
                                               & acc_with_time$time <= (Final_data_6mwt[360,1]+i),]
      
    }else{
      break
    }
  }
  
  
  Data_to_plot <- acc_with_time[acc_with_time$time >= (Final_data_6mwt[1,1]-20) 
                                & acc_with_time$time <= (Final_data_6mwt[360,1]+20),]
  
  
  p <- ggplot(data = Data_to_plot,
              aes(x = time))+
    geom_line(aes(y = V1 ,
                  color = "V1"))+
    geom_line(aes(y = V2,
                  color = "V2"))+
    geom_line(aes(y = V3, 
                  color = "V3"))+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_6mwt[c(1, dim(Final_data_6mwt)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_6mwt[j],
                           start = (nchar(path_to_smartphone_6mwt)+2),
                           end = -5)))+
    scale_color_manual(name = "Group", 
                       values = c("V1" = "blue",
                                  "V2" = "green",
                                  "V3" = "purple"),
                       labels = c("V1",
                                  "V2",
                                  "V3"))+
    xlab("Time")+
    ylab("Value")
  
  
  
  # Define path where we save our output to
  setwd(paste0(maindir, "/", subdir_accelerometer_6mwt))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_6mwt[j], 
                                start = (nchar(path_to_smartphone_6mwt)+2)),
                        ".png"),
         width = 14,
         height = 10, 
         units = "cm")
  
  
  # sct
  
  accelerometer_data_sct <- acc_with_time[acc_with_time$time >= Final_data_sct[1,1]
                                          & acc_with_time$time <= Final_data_sct[60,1],]
  
  
  for(i in seq(0,5, by =0.1)){
    if(as.numeric(difftime(accelerometer_data_sct[dim(accelerometer_data_sct)[1],1],
                           accelerometer_data_sct[1,1])) < 1)
    {
      accelerometer_data_sct <- acc_with_time[acc_with_time$time >= (Final_data_sct[1,1]-i)
                                              & acc_with_time$time <= (Final_data_sct[60,1]+i),]
      
    }else{
      break
    }
  }
  
  
  
  Data_to_plot <- acc_with_time[acc_with_time$time >= (Final_data_sct[1,1]-20) & acc_with_time$time <= (Final_data_sct[60,1]+20),]
  
  
  
  p <- ggplot(data = Data_to_plot,
              aes(x = time))+
    geom_line(aes(y = V1 ,
                  color = "V1"))+
    geom_line(aes(y = V2, 
                  color = "V2"))+
    geom_line(aes(y = V3,
                  color = "V3"))+
    theme_bw()+
    geom_vline(xintercept = as.numeric(Final_data_sct[c(1, dim(Final_data_sct)[1]),1]),
               color = "red")+
    ggtitle(paste0("Plot of data file ", 
                   str_sub(smartphone_file_to_6mwt[j],
                           start = (nchar(path_to_smartphone_6mwt)+2),
                           end = -5)))+
    scale_color_manual(name = "Group", 
                       values = c("V1" = "blue",
                                  "V2" = "green", 
                                  "V3" = "purple"),
                       labels = c("V1",
                                  "V2",
                                  "V3"))+
    xlab("Time")+
    ylab("Value")
  
  
  # Define path where we save our output to
  setwd(paste0(maindir, "/", subdir_accelerometer_sct))
  ggplotly(p)
  ggsave(p, file=paste0(str_sub(smartphone_file_to_6mwt[j],
                                start = (nchar(path_to_smartphone_6mwt)+2)),
                        ".png"),
         width = 14,
         height = 10,
         units = "cm")
  
  
  
  
  
  Minimum_of_smartphone_6mwt <- min(Final_data_6mwt$tmp.value)
  Maximum_of_smartphone_6mwt <-max(Final_data_6mwt$tmp.value)
  Deviation_of_smartphone_6mwt <-sd(Final_data_6mwt$tmp.value)
  
  Minimum_of_smartphone_sct <-min(Final_data_sct$tmp.value)
  Maximum_of_smartphone_sct <-max(Final_data_sct$tmp.value)
  Deviation_of_smartphone_sct <-sd(Final_data_sct$tmp.value)
  
  
  Minimum_of_fitbit_6mwt <-min(fitbit_6mwt$bpm)
  Maximum_of_fitbit_6mwt <-max(fitbit_6mwt$bpm)
  Deviation_of_fitbit_6mwt <-sd(fitbit_6mwt$bpm)
  
  distribution <- table(fitbit_6mwt$confidence)/dim(fitbit_6mwt)[1]
  
  amount_of_confidence_zero_in_6mwt <- distribution[1][[1]]
  amount_of_confidence_one_in_6mwt <- distribution[2][[1]]
  amount_of_confidence_two_in_6mwt <- distribution[3][[1]]
  amount_of_confidence_three_in_6mwt <- distribution[4][[1]]
  
  
  Minimum_of_fitbit_sct <- min(fitbit_sct$bpm)
  Maximum_of_fitbit_sct <- max(fitbit_sct$bpm)
  Deviation_of_fitbit_sct <- sd(fitbit_sct$bpm)
  
  distribution <- table(fitbit_sct$confidence)/dim(fitbit_sct)[1]
  
  amount_of_confidence_zero_in_sct <- distribution[1][[1]]
  amount_of_confidence_one_in_sct <- distribution[2][[1]]
  amount_of_confidence_two_in_sct <- distribution[3][[1]]
  amount_of_confidence_three_in_sct <- distribution[4][[1]]
  
  
  splitter <- split(accelerometer_data_6mwt,
                    (seq(nrow(accelerometer_data_6mwt))-1) %/% (nrow(accelerometer_data_6mwt)/4))
  
  Maximum_of_accelerometer_6mwt_split1_v1 <- max(splitter$`0`[,2])
  Minimum_of_accelerometer_6mwt_split1_v1 <- min(splitter$`0`[,2])
  Deviation_of_accelerometer_6mwt_split1_v1 <- sd(splitter$`0`[,2])
  
  Maximum_of_accelerometer_6mwt_split2_v1 <- max(splitter$`1`[,2])
  Minimum_of_accelerometer_6mwt_split2_v1 <- min(splitter$`1`[,2])
  Deviation_of_accelerometer_6mwt_split2_v1 <- sd(splitter$`1`[,2])
  
  Maximum_of_accelerometer_6mwt_split3_v1 <- max(splitter$`2`[,2])
  Minimum_of_accelerometer_6mwt_split3_v1 <- min(splitter$`2`[,2])
  Deviation_of_accelerometer_6mwt_split3_v1 <- sd(splitter$`2`[,2])
  
  Maximum_of_accelerometer_6mwt_split4_v1 <- max(splitter$`3`[,2])
  Minimum_of_accelerometer_6mwt_split4_v1 <- min(splitter$`3`[,2])
  Deviation_of_accelerometer_6mwt_split4_v1 <- sd(splitter$`3`[,2])
  
  
  # V2
  Maximum_of_accelerometer_6mwt_split1_v2 <- max(splitter$`0`[,3])
  Minimum_of_accelerometer_6mwt_split1_v2 <- min(splitter$`0`[,3])
  Deviation_of_accelerometer_6mwt_split1_v2 <- sd(splitter$`0`[,3])
  
  Maximum_of_accelerometer_6mwt_split2_v2 <- max(splitter$`1`[,3])
  Minimum_of_accelerometer_6mwt_split2_v2 <- min(splitter$`1`[,3])
  Deviation_of_accelerometer_6mwt_split2_v2 <- sd(splitter$`1`[,3])
  
  Maximum_of_accelerometer_6mwt_split3_v2 <- max(splitter$`2`[,3])
  Minimum_of_accelerometer_6mwt_split3_v2 <- min(splitter$`2`[,3])
  Deviation_of_accelerometer_6mwt_split3_v2 <- sd(splitter$`2`[,3])
  
  Maximum_of_accelerometer_6mwt_split4_v2 <- max(splitter$`3`[,3])
  Minimum_of_accelerometer_6mwt_split4_v2 <- min(splitter$`3`[,3])
  Deviation_of_accelerometer_6mwt_split4_v2 <- sd(splitter$`3`[,3])
  
  # V3
  
  Maximum_of_accelerometer_6mwt_split1_v3 <- max(splitter$`0`[,4])
  Minimum_of_accelerometer_6mwt_split1_v3 <- min(splitter$`0`[,4])
  Deviation_of_accelerometer_6mwt_split1_v3 <- sd(splitter$`0`[,4])
  
  Maximum_of_accelerometer_6mwt_split2_v3 <- max(splitter$`1`[,4])
  Minimum_of_accelerometer_6mwt_split2_v3 <- min(splitter$`1`[,4])
  Deviation_of_accelerometer_6mwt_split2_v3 <- sd(splitter$`1`[,4])
  
  Maximum_of_accelerometer_6mwt_split3_v3 <- max(splitter$`2`[,4])
  Minimum_of_accelerometer_6mwt_split3_v3 <- min(splitter$`2`[,4])
  Deviation_of_accelerometer_6mwt_split3_v3 <- sd(splitter$`2`[,4])
  
  Maximum_of_accelerometer_6mwt_split4_v3 <- max(splitter$`3`[,4])
  Minimum_of_accelerometer_6mwt_split4_v3 <- min(splitter$`3`[,4])
  Deviation_of_accelerometer_6mwt_split4_v3 <- sd(splitter$`3`[,4])
  
  
  splitter <- split(accelerometer_data_sct,
                    (seq(nrow(accelerometer_data_sct))-1) %/% (nrow(accelerometer_data_sct)/2))
  
  Maximum_of_accelerometer_sct_split1_v1 <- max(splitter$`0`[,2])
  Minimum_of_accelerometer_sct_split1_v1 <- min(splitter$`0`[,2])
  Deviation_of_accerometer_sct_split1_V1 <- sd(splitter$`0`[,2])
  
  Maximum_of_accelerometer_sct_split2_v1 <- max(splitter$`1`[,2])
  Minimum_of_accelerometer_sct_split2_v1 <- min(splitter$`1`[,2])
  Deviation_of_accerometer_sct_split2_V1 <- sd(splitter$`1`[,2])
  
  #V2
  Maximum_of_accelerometer_sct_split1_v2 <- max(splitter$`0`[,3])
  Minimum_of_accelerometer_sct_split1_v2 <- min(splitter$`0`[,3])
  Deviation_of_accerometer_sct_split1_V2 <- sd(splitter$`0`[,3])
  
  Maximum_of_accelerometer_sct_split2_v2 <- max(splitter$`1`[,3])
  Minimum_of_accelerometer_sct_split2_v2 <- min(splitter$`1`[,3])
  Deviation_of_accerometer_sct_split2_V2 <- sd(splitter$`1`[,3])
  
  # V3
  
  Maximum_of_accelerometer_sct_split1_v3 <- max(splitter$`0`[,4])
  Minimum_of_accelerometer_sct_split1_v3 <- min(splitter$`0`[,4])
  Deviation_of_accerometer_sct_split1_V3 <- sd(splitter$`0`[,4])
  
  Maximum_of_accelerometer_sct_split2_v3 <- max(splitter$`1`[,4])
  Minimum_of_accelerometer_sct_split2_v3 <- min(splitter$`1`[,4])
  Deviation_of_accerometer_sct_split2_V3 <- sd(splitter$`1`[,4])
  
  
  
  
  output <- c(str_sub(smartphone_file_to_6mwt[j],
                      start = (nchar(path_to_smartphone_6mwt)+2),
                      end = -5),
              Minimum_of_smartphone_6mwt, 
              Maximum_of_smartphone_6mwt, 
              Deviation_of_smartphone_6mwt, 
              Minimum_of_smartphone_sct, 
              Maximum_of_smartphone_sct, 
              Deviation_of_smartphone_sct, 
              Minimum_of_fitbit_6mwt, 
              Maximum_of_fitbit_6mwt, 
              Deviation_of_fitbit_6mwt, 
              amount_of_confidence_zero_in_6mwt, 
              amount_of_confidence_one_in_6mwt, 
              amount_of_confidence_two_in_6mwt, 
              amount_of_confidence_three_in_6mwt, 
              Minimum_of_fitbit_sct, 
              Maximum_of_fitbit_sct,
              Deviation_of_fitbit_sct, 
              amount_of_confidence_zero_in_sct, 
              amount_of_confidence_one_in_sct, 
              amount_of_confidence_two_in_sct, 
              amount_of_confidence_three_in_sct, 
              Maximum_of_accelerometer_6mwt_split1_v1,
              Minimum_of_accelerometer_6mwt_split1_v1,
              Deviation_of_accelerometer_6mwt_split1_v1,
              Maximum_of_accelerometer_6mwt_split2_v1,
              Minimum_of_accelerometer_6mwt_split2_v1,
              Deviation_of_accelerometer_6mwt_split2_v1,
              Maximum_of_accelerometer_6mwt_split3_v1,
              Minimum_of_accelerometer_6mwt_split3_v1,
              Deviation_of_accelerometer_6mwt_split3_v1,
              Maximum_of_accelerometer_6mwt_split4_v1,
              Minimum_of_accelerometer_6mwt_split4_v1,
              Deviation_of_accelerometer_6mwt_split4_v1,
              Maximum_of_accelerometer_6mwt_split1_v2,
              Minimum_of_accelerometer_6mwt_split1_v2,
              Deviation_of_accelerometer_6mwt_split1_v2,
              Maximum_of_accelerometer_6mwt_split2_v2,
              Minimum_of_accelerometer_6mwt_split2_v2,
              Deviation_of_accelerometer_6mwt_split2_v2,
              Maximum_of_accelerometer_6mwt_split3_v2,
              Minimum_of_accelerometer_6mwt_split3_v2,
              Deviation_of_accelerometer_6mwt_split3_v2,
              Maximum_of_accelerometer_6mwt_split4_v2,
              Minimum_of_accelerometer_6mwt_split4_v2,
              Deviation_of_accelerometer_6mwt_split4_v2,
              Maximum_of_accelerometer_6mwt_split1_v3,
              Minimum_of_accelerometer_6mwt_split1_v3,
              Deviation_of_accelerometer_6mwt_split1_v3,
              Maximum_of_accelerometer_6mwt_split2_v3,
              Minimum_of_accelerometer_6mwt_split2_v3,
              Deviation_of_accelerometer_6mwt_split2_v3,
              Maximum_of_accelerometer_6mwt_split3_v3,
              Minimum_of_accelerometer_6mwt_split3_v3,
              Deviation_of_accelerometer_6mwt_split3_v3,
              Maximum_of_accelerometer_6mwt_split4_v3,
              Minimum_of_accelerometer_6mwt_split4_v3,
              Deviation_of_accelerometer_6mwt_split4_v3,
              Maximum_of_accelerometer_sct_split1_v1,
              Minimum_of_accelerometer_sct_split1_v1,
              Deviation_of_accerometer_sct_split1_V1,
              Maximum_of_accelerometer_sct_split2_v1,
              Minimum_of_accelerometer_sct_split2_v1,
              Deviation_of_accerometer_sct_split2_V1,
              Maximum_of_accelerometer_sct_split1_v2,
              Minimum_of_accelerometer_sct_split1_v2,
              Deviation_of_accerometer_sct_split1_V2,
              Maximum_of_accelerometer_sct_split2_v2,
              Minimum_of_accelerometer_sct_split2_v2,
              Deviation_of_accerometer_sct_split2_V2,
              Maximum_of_accelerometer_sct_split1_v3,
              Minimum_of_accelerometer_sct_split1_v3,
              Deviation_of_accerometer_sct_split1_V3,
              Maximum_of_accelerometer_sct_split2_v3,
              Minimum_of_accelerometer_sct_split2_v3,
              Deviation_of_accerometer_sct_split2_V3
              
              
  )
  
  df <- rbind(df, output)
}


# For the case that 6mwt and sct wher in one file in the smartphone data




smartphone_file_to_sct <- list.files(path_to_smartphone_sct
                                     ,pattern = ".csv",
                                     full.names = TRUE)

smartphone_file_to_6mwt <- list.files(path_to_smartphone_6mwt, 
                                      pattern = ".csv",
                                      full.names =  TRUE)


# Remove the unwanted data that would be loader otherwise 

smartphone_file_to_sct <- smartphone_file_to_sct[c(42)]

smartphone_file_to_6mwt <- smartphone_file_to_6mwt[c(42)]




# Fitbit files
fitbit_file <- paste0(path_to_fitbit,
                      "/",
                      substring(smartphone_file_to_6mwt,
                                first = (nchar(path_to_smartphone_6mwt)+11)))


#Accelrometer files 

accelerometer_file <- paste0(path_to_accelerometer,
                             "/",
                             substring(smartphone_file_to_6mwt,
                                       first = (nchar(path_to_smartphone_6mwt)+11), 
                                       last = (nchar(smartphone_file_to_6mwt)-4)),
                             "/")



tmp <-  read.csv(smartphone_file_to_sct, sep = ";")




# Filter the file such that we have only look at the data of interest
tmp <- tmp%>%filter(type_profile == "basic" & type_category == "acceleration")

# Adjust the time problem (R is not good when it comes to micro seconds, so 
# we need to help us out here with this work around)
date_and_time <- data.frame(str_sub(tmp$time, end=-4))

miliseconds <- as.numeric(str_sub(tmp$time, -3,-1))/1000

time <- as.POSIXct(as.numeric(date_and_time$str_sub.tmp.time..end....4.),
                   origin="1970-01-01",
                   tz = "UTC")
time <- time + miliseconds


# Bind time and value to one data.frame
data_of_interest <- data.frame(time,tmp$value)

data_of_interest <- data_of_interest%>%filter(time > as.POSIXct("2021-06-16 16:00:00", tz = "UTC") & time < as.POSIXct("2021-06-16 17:00:00", tz = "UTC") )


# To delete oulier values that are just way off from our normal distribution 
# This would be a big problem, because we identify our 
# 99.7% of the data that would stem from the same distribution
data_of_interest <- data_of_interest[data_of_interest$tmp.value
                                     < (mean(data_of_interest$tmp.value)
                                        +3*sd(data_of_interest$tmp.value)),]



frame_to_iterate <- which(data_of_interest$tmp.value > mean(data_of_interest$tmp.value))

start <- list()
end  <- list()
Result <- list()

for(i in frame_to_iterate){
  if((i+59) <= dim(data_of_interest)[1]){
    
    start[[i]] <- i
    end[[i]] <- i+59
    Result[[i]] <-sum(data_of_interest[i:(59+i),2]>mean(data_of_interest$tmp.value)) 
  }
  else{
    break
  }
  
  
}
start <- do.call(rbind.data.frame,
                 start)

end <- do.call(rbind.data.frame,
               end)

Result <- do.call(rbind.data.frame,
                  Result)

# Bind values together to one data.frame
intermediate_Result <- cbind(start,
                             end,
                             Result)

colnames(intermediate_Result) <- c("start",
                                   "end",
                                   "Result")




# Find the interval that has has the most values above the mean

Maximum <- max(intermediate_Result$Result, 
               na.rm = T)

helper <- na.omit(intermediate_Result[intermediate_Result$Result == Maximum, ])


iter_for_sct <- helper[ceiling(dim(helper)[1]/2),1]

Final_data_sct <- data_of_interest[iter_for_sct:(iter_for_sct+59),]




tmp <-  read.csv(smartphone_file_to_6mwt, sep = ";")

# Filter the file such that we have only look at the data of interest
tmp <- tmp%>%filter(type_profile == "basic" & type_category == "acceleration")

# Adjust the time problem (R is not good when it comes to micro seconds, so 
# we need to help us out here with this work around)
date_and_time <- data.frame(str_sub(tmp$time, end=-4))

miliseconds <- as.numeric(str_sub(tmp$time, -3,-1))/1000

time <- as.POSIXct(as.numeric(date_and_time$str_sub.tmp.time..end....4.),
                   origin="1970-01-01",
                   tz = "UTC")
time <- time + miliseconds


# Bind time and value to one data.frame
data_of_interest <- data.frame(time,tmp$value)

data_of_interest <- data_of_interest%>%filter(time < as.POSIXct("2021-06-16 16:00:00", tz = "UTC") )



# To delete oulier values that are just way off from our normal distribution 
# This would be a big problem, because we identify our 
# 99.7% of the data that would stem from the same distribution
data_of_interest <- data_of_interest[data_of_interest$tmp.value
                                     < (mean(data_of_interest$tmp.value)
                                        +3*sd(data_of_interest$tmp.value)),]



frame_to_iterate <- which(data_of_interest$tmp.value > mean(data_of_interest$tmp.value))

start <- list()
end  <- list()
Result <- list()

for(i in frame_to_iterate){
  if((i+359) <= dim(data_of_interest)[1]){
    
    start[[i]] <- i
    end[[i]] <- i+359
    Result[[i]] <-sum(data_of_interest[i:(359+i),2]>mean(data_of_interest$tmp.value)) 
  }
  else{
    break
  }
  
  
}
start <- do.call(rbind.data.frame,
                 start)

end <- do.call(rbind.data.frame,
               end)

Result <- do.call(rbind.data.frame,
                  Result)

# Bind values together to one data.frame
intermediate_Result <- cbind(start,
                             end,
                             Result)

colnames(intermediate_Result) <- c("start",
                                   "end",
                                   "Result")




# Find the interval that has has the most values above the mean

Maximum <- max(intermediate_Result$Result,
               na.rm = T)

helper <- na.omit(intermediate_Result[intermediate_Result$Result == Maximum, ])


iter_for_6mwt <- helper[ceiling(dim(helper)[1]/2),1]

Final_data_6mwt <- data_of_interest[iter_for_6mwt:(iter_for_6mwt+359),]



##### Fitbit Data #####

fitbit_data <- read.csv(fitbit_file, sep = ";")
fitbit_data$confidence <- as.factor(fitbit_data$confidence)

fitbit_data$dateTime <- as.POSIXct(as.character(fitbit_data$dateTime), 
                                   format="%m/%d/%y %H:%M:%S",
                                   tz = "UTC")



fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= Final_data_6mwt[1,1] 
                           & fitbit_data$dateTime <= Final_data_6mwt[360,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(fitbit_6mwt[dim(fitbit_6mwt)[1],1],fitbit_6mwt[1,1])) < 6){
    fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= (Final_data_6mwt[1,1]-i)
                               & fitbit_data$dateTime <= (Final_data_6mwt[360,1]+i),]
    
  }else{
    break
  }
}



fitbit_sct <- fitbit_data[fitbit_data$dateTime >= Final_data_sct[1,1] 
                          & fitbit_data$dateTime <= Final_data_sct[60,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(fitbit_sct[dim(fitbit_sct)[1],1],
                         fitbit_sct[1,1],
                         units = "mins")) < 1){
    fitbit_sct <- fitbit_data[fitbit_data$dateTime >= (Final_data_sct[1,1]-i) 
                              & fitbit_data$dateTime <= (Final_data_sct[60,1]+i),]
    
  }else{
    break
  }
}



##### Accelerometer Data #####


accelerometer_data <- read.csv(paste0(accelerometer_file,
                                      "acc.csv"),
                               sep = ";",
                               header = F)


xmlfile <- xmlParse(file = paste0(accelerometer_file,
                                  "unisens.xml"))


rootnode <- xmlRoot(xmlfile)
helperlist <- xmlToList(rootnode)


starting_point <- as.POSIXct(gsub("T",
                                  " ",
                                  helperlist$.attrs[[4]]))
#Subtract the two hours from the CEST format to get to UTC
starting_point <- starting_point -7200


starting_point <- as.POSIXct(as.character(starting_point),
                             tz = "UTC")

time <- seq(from = starting_point,
            to =(starting_point +
                   (dim(accelerometer_data)[1]/64)-1),
            length.out = dim(accelerometer_data)[1] )

acc_with_time <- cbind(time, accelerometer_data)


# 6mwt

accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= Final_data_6mwt[1,1]
                                         & acc_with_time$time <= Final_data_6mwt[360,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(accelerometer_data_6mwt[dim(accelerometer_data_6mwt)[1],1],
                         accelerometer_data_6mwt[1,1])) < 6)
  {
    accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= (Final_data_6mwt[1,1]-i) 
                                             & acc_with_time$time <= (Final_data_6mwt[360,1]+i),]
    
  }else{
    break
  }
}



# sct

accelerometer_data_sct <- acc_with_time[acc_with_time$time >= Final_data_sct[1,1]
                                        & acc_with_time$time <= Final_data_sct[60,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(accelerometer_data_sct[dim(accelerometer_data_sct)[1],1],
                         accelerometer_data_sct[1,1])) < 1)
  {
    accelerometer_data_sct <- acc_with_time[acc_with_time$time >= (Final_data_sct[1,1]-i)
                                            & acc_with_time$time <= (Final_data_sct[60,1]+i),]
    
  }else{
    break
  }
}

Minimum_of_smartphone_6mwt <- min(Final_data_6mwt$tmp.value)
Maximum_of_smartphone_6mwt <-max(Final_data_6mwt$tmp.value)
Deviation_of_smartphone_6mwt <-sd(Final_data_6mwt$tmp.value)

Minimum_of_smartphone_sct <-min(Final_data_sct$tmp.value)
Maximum_of_smartphone_sct <-max(Final_data_sct$tmp.value)
Deviation_of_smartphone_sct <-sd(Final_data_sct$tmp.value)


Minimum_of_fitbit_6mwt <-min(fitbit_6mwt$bpm)
Maximum_of_fitbit_6mwt <-max(fitbit_6mwt$bpm)
Deviation_of_fitbit_6mwt <-sd(fitbit_6mwt$bpm)

distribution <- table(fitbit_6mwt$confidence)/dim(fitbit_6mwt)[1]

amount_of_confidence_zero_in_6mwt <- distribution[1][[1]]
amount_of_confidence_one_in_6mwt <- distribution[2][[1]]
amount_of_confidence_two_in_6mwt <- distribution[3][[1]]
amount_of_confidence_three_in_6mwt <- distribution[4][[1]]


Minimum_of_fitbit_sct <- min(fitbit_sct$bpm)
Maximum_of_fitbit_sct <- max(fitbit_sct$bpm)
Deviation_of_fitbit_sct <- sd(fitbit_sct$bpm)

distribution <- table(fitbit_sct$confidence)/dim(fitbit_sct)[1]

amount_of_confidence_zero_in_sct <- distribution[1][[1]]
amount_of_confidence_one_in_sct <- distribution[2][[1]]
amount_of_confidence_two_in_sct <- distribution[3][[1]]
amount_of_confidence_three_in_sct <- distribution[4][[1]]


splitter <- split(accelerometer_data_6mwt,
                  (seq(nrow(accelerometer_data_6mwt))-1) %/% (nrow(accelerometer_data_6mwt)/4))

Maximum_of_accelerometer_6mwt_split1_v1 <- max(splitter$`0`[,2])
Minimum_of_accelerometer_6mwt_split1_v1 <- min(splitter$`0`[,2])
Deviation_of_accelerometer_6mwt_split1_v1 <- sd(splitter$`0`[,2])

Maximum_of_accelerometer_6mwt_split2_v1 <- max(splitter$`1`[,2])
Minimum_of_accelerometer_6mwt_split2_v1 <- min(splitter$`1`[,2])
Deviation_of_accelerometer_6mwt_split2_v1 <- sd(splitter$`1`[,2])

Maximum_of_accelerometer_6mwt_split3_v1 <- max(splitter$`2`[,2])
Minimum_of_accelerometer_6mwt_split3_v1 <- min(splitter$`2`[,2])
Deviation_of_accelerometer_6mwt_split3_v1 <- sd(splitter$`2`[,2])

Maximum_of_accelerometer_6mwt_split4_v1 <- max(splitter$`3`[,2])
Minimum_of_accelerometer_6mwt_split4_v1 <- min(splitter$`3`[,2])
Deviation_of_accelerometer_6mwt_split4_v1 <- sd(splitter$`3`[,2])


# V2
Maximum_of_accelerometer_6mwt_split1_v2 <- max(splitter$`0`[,3])
Minimum_of_accelerometer_6mwt_split1_v2 <- min(splitter$`0`[,3])
Deviation_of_accelerometer_6mwt_split1_v2 <- sd(splitter$`0`[,3])

Maximum_of_accelerometer_6mwt_split2_v2 <- max(splitter$`1`[,3])
Minimum_of_accelerometer_6mwt_split2_v2 <- min(splitter$`1`[,3])
Deviation_of_accelerometer_6mwt_split2_v2 <- sd(splitter$`1`[,3])

Maximum_of_accelerometer_6mwt_split3_v2 <- max(splitter$`2`[,3])
Minimum_of_accelerometer_6mwt_split3_v2 <- min(splitter$`2`[,3])
Deviation_of_accelerometer_6mwt_split3_v2 <- sd(splitter$`2`[,3])

Maximum_of_accelerometer_6mwt_split4_v2 <- max(splitter$`3`[,3])
Minimum_of_accelerometer_6mwt_split4_v2 <- min(splitter$`3`[,3])
Deviation_of_accelerometer_6mwt_split4_v2 <- sd(splitter$`3`[,3])

# V3

Maximum_of_accelerometer_6mwt_split1_v3 <- max(splitter$`0`[,4])
Minimum_of_accelerometer_6mwt_split1_v3 <- min(splitter$`0`[,4])
Deviation_of_accelerometer_6mwt_split1_v3 <- sd(splitter$`0`[,4])

Maximum_of_accelerometer_6mwt_split2_v3 <- max(splitter$`1`[,4])
Minimum_of_accelerometer_6mwt_split2_v3 <- min(splitter$`1`[,4])
Deviation_of_accelerometer_6mwt_split2_v3 <- sd(splitter$`1`[,4])

Maximum_of_accelerometer_6mwt_split3_v3 <- max(splitter$`2`[,4])
Minimum_of_accelerometer_6mwt_split3_v3 <- min(splitter$`2`[,4])
Deviation_of_accelerometer_6mwt_split3_v3 <- sd(splitter$`2`[,4])

Maximum_of_accelerometer_6mwt_split4_v3 <- max(splitter$`3`[,4])
Minimum_of_accelerometer_6mwt_split4_v3 <- min(splitter$`3`[,4])
Deviation_of_accelerometer_6mwt_split4_v3 <- sd(splitter$`3`[,4])


splitter <- split(accelerometer_data_sct,
                  (seq(nrow(accelerometer_data_sct))-1) %/% (nrow(accelerometer_data_sct)/2))

Maximum_of_accelerometer_sct_split1_v1 <- max(splitter$`0`[,2])
Minimum_of_accelerometer_sct_split1_v1 <- min(splitter$`0`[,2])
Deviation_of_accerometer_sct_split1_V1 <- sd(splitter$`0`[,2])

Maximum_of_accelerometer_sct_split2_v1 <- max(splitter$`1`[,2])
Minimum_of_accelerometer_sct_split2_v1 <- min(splitter$`1`[,2])
Deviation_of_accerometer_sct_split2_V1 <- sd(splitter$`1`[,2])

#V2
Maximum_of_accelerometer_sct_split1_v2 <- max(splitter$`0`[,3])
Minimum_of_accelerometer_sct_split1_v2 <- min(splitter$`0`[,3])
Deviation_of_accerometer_sct_split1_V2 <- sd(splitter$`0`[,3])

Maximum_of_accelerometer_sct_split2_v2 <- max(splitter$`1`[,3])
Minimum_of_accelerometer_sct_split2_v2 <- min(splitter$`1`[,3])
Deviation_of_accerometer_sct_split2_V2 <- sd(splitter$`1`[,3])

# V3

Maximum_of_accelerometer_sct_split1_v3 <- max(splitter$`0`[,4])
Minimum_of_accelerometer_sct_split1_v3 <- min(splitter$`0`[,4])
Deviation_of_accerometer_sct_split1_V3 <- sd(splitter$`0`[,4])

Maximum_of_accelerometer_sct_split2_v3 <- max(splitter$`1`[,4])
Minimum_of_accelerometer_sct_split2_v3 <- min(splitter$`1`[,4])
Deviation_of_accerometer_sct_split2_V3 <- sd(splitter$`1`[,4])




output <- c(str_sub(smartphone_file_to_6mwt,
                    start = (nchar(path_to_smartphone_6mwt)+11),
                    end = -5),
            Minimum_of_smartphone_6mwt, 
            Maximum_of_smartphone_6mwt, 
            Deviation_of_smartphone_6mwt, 
            Minimum_of_smartphone_sct, 
            Maximum_of_smartphone_sct, 
            Deviation_of_smartphone_sct, 
            Minimum_of_fitbit_6mwt, 
            Maximum_of_fitbit_6mwt, 
            Deviation_of_fitbit_6mwt, 
            amount_of_confidence_zero_in_6mwt, 
            amount_of_confidence_one_in_6mwt, 
            amount_of_confidence_two_in_6mwt, 
            amount_of_confidence_three_in_6mwt, 
            Minimum_of_fitbit_sct, 
            Maximum_of_fitbit_sct, 
            Deviation_of_fitbit_sct, 
            amount_of_confidence_zero_in_sct, 
            amount_of_confidence_one_in_sct, 
            amount_of_confidence_two_in_sct, 
            amount_of_confidence_three_in_sct, 
            Maximum_of_accelerometer_6mwt_split1_v1,
            Minimum_of_accelerometer_6mwt_split1_v1,
            Deviation_of_accelerometer_6mwt_split1_v1,
            Maximum_of_accelerometer_6mwt_split2_v1,
            Minimum_of_accelerometer_6mwt_split2_v1,
            Deviation_of_accelerometer_6mwt_split2_v1,
            Maximum_of_accelerometer_6mwt_split3_v1,
            Minimum_of_accelerometer_6mwt_split3_v1,
            Deviation_of_accelerometer_6mwt_split3_v1,
            Maximum_of_accelerometer_6mwt_split4_v1,
            Minimum_of_accelerometer_6mwt_split4_v1,
            Deviation_of_accelerometer_6mwt_split4_v1,
            Maximum_of_accelerometer_6mwt_split1_v2,
            Minimum_of_accelerometer_6mwt_split1_v2,
            Deviation_of_accelerometer_6mwt_split1_v2,
            Maximum_of_accelerometer_6mwt_split2_v2,
            Minimum_of_accelerometer_6mwt_split2_v2,
            Deviation_of_accelerometer_6mwt_split2_v2,
            Maximum_of_accelerometer_6mwt_split3_v2,
            Minimum_of_accelerometer_6mwt_split3_v2,
            Deviation_of_accelerometer_6mwt_split3_v2,
            Maximum_of_accelerometer_6mwt_split4_v2,
            Minimum_of_accelerometer_6mwt_split4_v2,
            Deviation_of_accelerometer_6mwt_split4_v2,
            Maximum_of_accelerometer_6mwt_split1_v3,
            Minimum_of_accelerometer_6mwt_split1_v3,
            Deviation_of_accelerometer_6mwt_split1_v3,
            Maximum_of_accelerometer_6mwt_split2_v3,
            Minimum_of_accelerometer_6mwt_split2_v3,
            Deviation_of_accelerometer_6mwt_split2_v3,
            Maximum_of_accelerometer_6mwt_split3_v3,
            Minimum_of_accelerometer_6mwt_split3_v3,
            Deviation_of_accelerometer_6mwt_split3_v3,
            Maximum_of_accelerometer_6mwt_split4_v3,
            Minimum_of_accelerometer_6mwt_split4_v3,
            Deviation_of_accelerometer_6mwt_split4_v3,
            Maximum_of_accelerometer_sct_split1_v1,
            Minimum_of_accelerometer_sct_split1_v1,
            Deviation_of_accerometer_sct_split1_V1,
            Maximum_of_accelerometer_sct_split2_v1,
            Minimum_of_accelerometer_sct_split2_v1,
            Deviation_of_accerometer_sct_split2_V1,
            Maximum_of_accelerometer_sct_split1_v2,
            Minimum_of_accelerometer_sct_split1_v2,
            Deviation_of_accerometer_sct_split1_V2,
            Maximum_of_accelerometer_sct_split2_v2,
            Minimum_of_accelerometer_sct_split2_v2,
            Deviation_of_accerometer_sct_split2_V2,
            Maximum_of_accelerometer_sct_split1_v3,
            Minimum_of_accelerometer_sct_split1_v3,
            Deviation_of_accerometer_sct_split1_V3,
            Maximum_of_accelerometer_sct_split2_v3,
            Minimum_of_accelerometer_sct_split2_v3,
            Deviation_of_accerometer_sct_split2_V3
            
            
)

df <- rbind(df, output)



# In the case of 2b9ffa95 the sct smartphone file has some information loss, 
# so we have close to none information about it 



smartphone_file_to_sct <- list.files(path_to_smartphone_sct
                                     ,pattern = ".csv",
                                     full.names = TRUE)

smartphone_file_to_6mwt <- list.files(path_to_smartphone_6mwt, 
                                      pattern = ".csv",
                                      full.names =  TRUE)


# Remove the unwanted data that would be loader otherwise 

smartphone_file_to_sct <- smartphone_file_to_sct[c(8)]

smartphone_file_to_6mwt <- smartphone_file_to_6mwt[c(8)]




# Fitbit files
fitbit_file <- paste0(path_to_fitbit,
                      "/",
                      substring(smartphone_file_to_6mwt,
                                first = (nchar(path_to_smartphone_6mwt)+2)))


#Accelrometer files 

accelerometer_file <- paste0(path_to_accelerometer,
                             "/",
                             substring(smartphone_file_to_6mwt,
                                       first = (nchar(path_to_smartphone_6mwt)+2), 
                                       last = (nchar(smartphone_file_to_6mwt)-4)),
                             "/")



tmp <-  read.csv(smartphone_file_to_6mwt, sep = ";")



# Filter the file such that we have only look at the data of interest
tmp <- tmp%>%filter(type_profile == "basic" & type_category == "acceleration")

# Adjust the time problem (R is not good when it comes to micro seconds, so 
# we need to help us out here with this work around)
date_and_time <- data.frame(str_sub(tmp$time, end=-4))

miliseconds <- as.numeric(str_sub(tmp$time, -3,-1))/1000

time <- as.POSIXct(as.numeric(date_and_time$str_sub.tmp.time..end....4.),
                   origin="1970-01-01",
                   tz = "UTC")
time <- time + miliseconds


# Bind time and value to one data.frame
data_of_interest <- data.frame(time,tmp$value)




# To delete oulier values that are just way off from our normal distribution 
# This would be a big problem, because we identify our 
# 99.7% of the data that would stem from the same distribution
data_of_interest <- data_of_interest[data_of_interest$tmp.value
                                     < (mean(data_of_interest$tmp.value)
                                        +3*sd(data_of_interest$tmp.value)),]



frame_to_iterate <- which(data_of_interest$tmp.value > mean(data_of_interest$tmp.value))

start <- list()
end  <- list()
Result <- list()

for(i in frame_to_iterate){
  if((i+359) <= dim(data_of_interest)[1]){
    
    start[[i]] <- i
    end[[i]] <- i+359
    Result[[i]] <-sum(data_of_interest[i:(359+i),2]>mean(data_of_interest$tmp.value)) 
  }
  else{
    break
  }
  
  
}
start <- do.call(rbind.data.frame,
                 start)

end <- do.call(rbind.data.frame,
               end)

Result <- do.call(rbind.data.frame,
                  Result)

# Bind values together to one data.frame
intermediate_Result <- cbind(start,
                             end,
                             Result)

colnames(intermediate_Result) <- c("start",
                                   "end",
                                   "Result")




# Find the interval that has has the most values above the mean

Maximum <- max(intermediate_Result$Result,
               na.rm = T)

helper <- na.omit(intermediate_Result[intermediate_Result$Result == Maximum, ])


iter_for_6mwt <- helper[ceiling(dim(helper)[1]/2),1]

Final_data_6mwt <- data_of_interest[iter_for_6mwt:(iter_for_6mwt+359),]



##### Fitbit Data #####

fitbit_data <- read.csv(fitbit_file, sep = ";")
fitbit_data$confidence <- as.factor(fitbit_data$confidence)

fitbit_data$dateTime <- as.POSIXct(as.character(fitbit_data$dateTime), 
                                   format="%m/%d/%y %H:%M:%S",
                                   tz = "UTC")



fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= Final_data_6mwt[1,1] 
                           & fitbit_data$dateTime <= Final_data_6mwt[360,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(fitbit_6mwt[dim(fitbit_6mwt)[1],1],fitbit_6mwt[1,1])) < 6){
    fitbit_6mwt <- fitbit_data[fitbit_data$dateTime >= (Final_data_6mwt[1,1]-i)
                               & fitbit_data$dateTime <= (Final_data_6mwt[360,1]+i),]
    
  }else{
    break
  }
}


##### Accelerometer Data #####


accelerometer_data <- read.csv(paste0(accelerometer_file,
                                      "acc.csv"),
                               sep = ";",
                               header = F)


xmlfile <- xmlParse(file = paste0(accelerometer_file,
                                  "unisens.xml"))


rootnode <- xmlRoot(xmlfile)
helperlist <- xmlToList(rootnode)


starting_point <- as.POSIXct(gsub("T",
                                  " ",
                                  helperlist$.attrs[[4]]))
#Subtract the two hours from the CEST format to get to UTC
starting_point <- starting_point -7200


starting_point <- as.POSIXct(as.character(starting_point),
                             tz = "UTC")

time <- seq(from = starting_point,
            to =(starting_point +
                   (dim(accelerometer_data)[1]/64)-1),
            length.out = dim(accelerometer_data)[1] )

acc_with_time <- cbind(time, accelerometer_data)


# 6mwt

accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= Final_data_6mwt[1,1]
                                         & acc_with_time$time <= Final_data_6mwt[360,1],]


for(i in seq(0,5, by =0.1)){
  if(as.numeric(difftime(accelerometer_data_6mwt[dim(accelerometer_data_6mwt)[1],1],
                         accelerometer_data_6mwt[1,1])) < 6)
  {
    accelerometer_data_6mwt <- acc_with_time[acc_with_time$time >= (Final_data_6mwt[1,1]-i) 
                                             & acc_with_time$time <= (Final_data_6mwt[360,1]+i),]
    
  }else{
    break
  }
}

#Because our method does work for the sct data, we eyeball them. From the fitbit 
# plot we can guess that the time interval has to be in between 
# 10:42:36 &  10:43:38

accelerometer_data_sct <- acc_with_time %>% filter(time >= as.POSIXct("2021-07-13 10:42:36", tz = "UTC") & time <=as.POSIXct("2021-07-13 10:43:38", tz = "UTC") )

fitbit_sct <- fitbit_data %>% filter(dateTime >= as.POSIXct("2021-07-13 10:42:36", tz = "UTC") & dateTime <=as.POSIXct("2021-07-13 10:43:38", tz = "UTC"))




Minimum_of_smartphone_6mwt <- min(Final_data_6mwt$tmp.value)
Maximum_of_smartphone_6mwt <-max(Final_data_6mwt$tmp.value)
Deviation_of_smartphone_6mwt <-sd(Final_data_6mwt$tmp.value)

Minimum_of_smartphone_sct <-NA
Maximum_of_smartphone_sct <-NA
Deviation_of_smartphone_sct <-NA


Minimum_of_fitbit_6mwt <-min(fitbit_6mwt$bpm)
Maximum_of_fitbit_6mwt <-max(fitbit_6mwt$bpm)
Deviation_of_fitbit_6mwt <-sd(fitbit_6mwt$bpm)

distribution <- table(fitbit_6mwt$confidence)/dim(fitbit_6mwt)[1]

amount_of_confidence_zero_in_6mwt <- distribution[1][[1]]
amount_of_confidence_one_in_6mwt <- distribution[2][[1]]
amount_of_confidence_two_in_6mwt <- distribution[3][[1]]
amount_of_confidence_three_in_6mwt <- distribution[4][[1]]


Minimum_of_fitbit_sct <- min(fitbit_sct$bpm)
Maximum_of_fitbit_sct <- max(fitbit_sct$bpm)
Deviation_of_fitbit_sct <- sd(fitbit_sct$bpm)

distribution <- table(fitbit_sct$confidence)/dim(fitbit_sct)[1]

amount_of_confidence_zero_in_sct <- distribution[1][[1]]
amount_of_confidence_one_in_sct <- distribution[2][[1]]
amount_of_confidence_two_in_sct <- distribution[3][[1]]
amount_of_confidence_three_in_sct <- distribution[4][[1]]


splitter <- split(accelerometer_data_6mwt,
                  (seq(nrow(accelerometer_data_6mwt))-1) %/% (nrow(accelerometer_data_6mwt)/4))

Maximum_of_accelerometer_6mwt_split1_v1 <- max(splitter$`0`[,2])
Minimum_of_accelerometer_6mwt_split1_v1 <- min(splitter$`0`[,2])
Deviation_of_accelerometer_6mwt_split1_v1 <- sd(splitter$`0`[,2])

Maximum_of_accelerometer_6mwt_split2_v1 <- max(splitter$`1`[,2])
Minimum_of_accelerometer_6mwt_split2_v1 <- min(splitter$`1`[,2])
Deviation_of_accelerometer_6mwt_split2_v1 <- sd(splitter$`1`[,2])

Maximum_of_accelerometer_6mwt_split3_v1 <- max(splitter$`2`[,2])
Minimum_of_accelerometer_6mwt_split3_v1 <- min(splitter$`2`[,2])
Deviation_of_accelerometer_6mwt_split3_v1 <- sd(splitter$`2`[,2])

Maximum_of_accelerometer_6mwt_split4_v1 <- max(splitter$`3`[,2])
Minimum_of_accelerometer_6mwt_split4_v1 <- min(splitter$`3`[,2])
Deviation_of_accelerometer_6mwt_split4_v1 <- sd(splitter$`3`[,2])


# V2
Maximum_of_accelerometer_6mwt_split1_v2 <- max(splitter$`0`[,3])
Minimum_of_accelerometer_6mwt_split1_v2 <- min(splitter$`0`[,3])
Deviation_of_accelerometer_6mwt_split1_v2 <- sd(splitter$`0`[,3])

Maximum_of_accelerometer_6mwt_split2_v2 <- max(splitter$`1`[,3])
Minimum_of_accelerometer_6mwt_split2_v2 <- min(splitter$`1`[,3])
Deviation_of_accelerometer_6mwt_split2_v2 <- sd(splitter$`1`[,3])

Maximum_of_accelerometer_6mwt_split3_v2 <- max(splitter$`2`[,3])
Minimum_of_accelerometer_6mwt_split3_v2 <- min(splitter$`2`[,3])
Deviation_of_accelerometer_6mwt_split3_v2 <- sd(splitter$`2`[,3])

Maximum_of_accelerometer_6mwt_split4_v2 <- max(splitter$`3`[,3])
Minimum_of_accelerometer_6mwt_split4_v2 <- min(splitter$`3`[,3])
Deviation_of_accelerometer_6mwt_split4_v2 <- sd(splitter$`3`[,3])

# V3

Maximum_of_accelerometer_6mwt_split1_v3 <- max(splitter$`0`[,4])
Minimum_of_accelerometer_6mwt_split1_v3 <- min(splitter$`0`[,4])
Deviation_of_accelerometer_6mwt_split1_v3 <- sd(splitter$`0`[,4])

Maximum_of_accelerometer_6mwt_split2_v3 <- max(splitter$`1`[,4])
Minimum_of_accelerometer_6mwt_split2_v3 <- min(splitter$`1`[,4])
Deviation_of_accelerometer_6mwt_split2_v3 <- sd(splitter$`1`[,4])

Maximum_of_accelerometer_6mwt_split3_v3 <- max(splitter$`2`[,4])
Minimum_of_accelerometer_6mwt_split3_v3 <- min(splitter$`2`[,4])
Deviation_of_accelerometer_6mwt_split3_v3 <- sd(splitter$`2`[,4])

Maximum_of_accelerometer_6mwt_split4_v3 <- max(splitter$`3`[,4])
Minimum_of_accelerometer_6mwt_split4_v3 <- min(splitter$`3`[,4])
Deviation_of_accelerometer_6mwt_split4_v3 <- sd(splitter$`3`[,4])


splitter <- split(accelerometer_data_sct,
                  (seq(nrow(accelerometer_data_sct))-1) %/% (nrow(accelerometer_data_sct)/2))

Maximum_of_accelerometer_sct_split1_v1 <- max(splitter$`0`[,2])
Minimum_of_accelerometer_sct_split1_v1 <- min(splitter$`0`[,2])
Deviation_of_accerometer_sct_split1_V1 <- sd(splitter$`0`[,2])

Maximum_of_accelerometer_sct_split2_v1 <- max(splitter$`1`[,2])
Minimum_of_accelerometer_sct_split2_v1 <- min(splitter$`1`[,2])
Deviation_of_accerometer_sct_split2_V1 <- sd(splitter$`1`[,2])

#V2
Maximum_of_accelerometer_sct_split1_v2 <- max(splitter$`0`[,3])
Minimum_of_accelerometer_sct_split1_v2 <- min(splitter$`0`[,3])
Deviation_of_accerometer_sct_split1_V2 <- sd(splitter$`0`[,3])

Maximum_of_accelerometer_sct_split2_v2 <- max(splitter$`1`[,3])
Minimum_of_accelerometer_sct_split2_v2 <- min(splitter$`1`[,3])
Deviation_of_accerometer_sct_split2_V2 <- sd(splitter$`1`[,3])

# V3

Maximum_of_accelerometer_sct_split1_v3 <- max(splitter$`0`[,4])
Minimum_of_accelerometer_sct_split1_v3 <- min(splitter$`0`[,4])
Deviation_of_accerometer_sct_split1_V3 <- sd(splitter$`0`[,4])

Maximum_of_accelerometer_sct_split2_v3 <- max(splitter$`1`[,4])
Minimum_of_accelerometer_sct_split2_v3 <- min(splitter$`1`[,4])
Deviation_of_accerometer_sct_split2_V3 <- sd(splitter$`1`[,4])




output <- c(str_sub(smartphone_file_to_6mwt,
                    start = (nchar(path_to_smartphone_6mwt)+2),
                    end = -5),
            Minimum_of_smartphone_6mwt, 
            Maximum_of_smartphone_6mwt, 
            Deviation_of_smartphone_6mwt, 
            Minimum_of_smartphone_sct, 
            Maximum_of_smartphone_sct, 
            Deviation_of_smartphone_sct, 
            Minimum_of_fitbit_6mwt, 
            Maximum_of_fitbit_6mwt, 
            Deviation_of_fitbit_6mwt, 
            amount_of_confidence_zero_in_6mwt, 
            amount_of_confidence_one_in_6mwt, 
            amount_of_confidence_two_in_6mwt, 
            amount_of_confidence_three_in_6mwt, 
            Minimum_of_fitbit_sct, 
            Maximum_of_fitbit_sct, 
            Deviation_of_fitbit_sct, 
            amount_of_confidence_zero_in_sct, 
            amount_of_confidence_one_in_sct, 
            amount_of_confidence_two_in_sct, 
            amount_of_confidence_three_in_sct, 
            Maximum_of_accelerometer_6mwt_split1_v1,
            Minimum_of_accelerometer_6mwt_split1_v1,
            Deviation_of_accelerometer_6mwt_split1_v1,
            Maximum_of_accelerometer_6mwt_split2_v1,
            Minimum_of_accelerometer_6mwt_split2_v1,
            Deviation_of_accelerometer_6mwt_split2_v1,
            Maximum_of_accelerometer_6mwt_split3_v1,
            Minimum_of_accelerometer_6mwt_split3_v1,
            Deviation_of_accelerometer_6mwt_split3_v1,
            Maximum_of_accelerometer_6mwt_split4_v1,
            Minimum_of_accelerometer_6mwt_split4_v1,
            Deviation_of_accelerometer_6mwt_split4_v1,
            Maximum_of_accelerometer_6mwt_split1_v2,
            Minimum_of_accelerometer_6mwt_split1_v2,
            Deviation_of_accelerometer_6mwt_split1_v2,
            Maximum_of_accelerometer_6mwt_split2_v2,
            Minimum_of_accelerometer_6mwt_split2_v2,
            Deviation_of_accelerometer_6mwt_split2_v2,
            Maximum_of_accelerometer_6mwt_split3_v2,
            Minimum_of_accelerometer_6mwt_split3_v2,
            Deviation_of_accelerometer_6mwt_split3_v2,
            Maximum_of_accelerometer_6mwt_split4_v2,
            Minimum_of_accelerometer_6mwt_split4_v2,
            Deviation_of_accelerometer_6mwt_split4_v2,
            Maximum_of_accelerometer_6mwt_split1_v3,
            Minimum_of_accelerometer_6mwt_split1_v3,
            Deviation_of_accelerometer_6mwt_split1_v3,
            Maximum_of_accelerometer_6mwt_split2_v3,
            Minimum_of_accelerometer_6mwt_split2_v3,
            Deviation_of_accelerometer_6mwt_split2_v3,
            Maximum_of_accelerometer_6mwt_split3_v3,
            Minimum_of_accelerometer_6mwt_split3_v3,
            Deviation_of_accelerometer_6mwt_split3_v3,
            Maximum_of_accelerometer_6mwt_split4_v3,
            Minimum_of_accelerometer_6mwt_split4_v3,
            Deviation_of_accelerometer_6mwt_split4_v3,
            Maximum_of_accelerometer_sct_split1_v1,
            Minimum_of_accelerometer_sct_split1_v1,
            Deviation_of_accerometer_sct_split1_V1,
            Maximum_of_accelerometer_sct_split2_v1,
            Minimum_of_accelerometer_sct_split2_v1,
            Deviation_of_accerometer_sct_split2_V1,
            Maximum_of_accelerometer_sct_split1_v2,
            Minimum_of_accelerometer_sct_split1_v2,
            Deviation_of_accerometer_sct_split1_V2,
            Maximum_of_accelerometer_sct_split2_v2,
            Minimum_of_accelerometer_sct_split2_v2,
            Deviation_of_accerometer_sct_split2_V2,
            Maximum_of_accelerometer_sct_split1_v3,
            Minimum_of_accelerometer_sct_split1_v3,
            Deviation_of_accerometer_sct_split1_V3,
            Maximum_of_accelerometer_sct_split2_v3,
            Minimum_of_accelerometer_sct_split2_v3,
            Deviation_of_accerometer_sct_split2_V3
            
            
)

df <- rbind(df, output)


colnames(df) <- df[1,]
df <- df[-1,]


# Fixing data error
df$amount_of_confidence_three_in_6mwt[1] <- 0
df$amount_of_confidence_three_in_sct[1] <- 0


write.csv(df,paste0(maindir,"/additional_data.csv"),
          sep = ";",
          dec = ".",
          col.names = TRUE, 
          row.names = FALSE)



