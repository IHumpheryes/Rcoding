
#---###########################################
#---  title: "RICT IV Helper files version 1"
#--- output: html_document
#--- Author: Dr K. M Muyeba aka Maybin
#--- #########################################

  ####################################################################################################################################
  # This R Script is a functional module that aims to define a series of small utility functions to help with mundane processes within
  # prediction and / classification functions. We factor out these functions to help readability of the R codes
  # All RICT IV models, GB, NI etc can use this script for allthe routine processes ,all in one place
  ####################################################################################################################################



library(dplyr)
library(data.table)

# 1. getLogAlkalinity
# This function gets Alkalinity from hardness, calcium, conduct when supplied, else return the alkalinity rovided as log10
# Returns log10 of the raw value. Use mapply as a "m"ultivariate version of sapply, applies to each element
# Can use "switch"" for better readability.
# Produce dataframe with value of alkalinity in bounds [1.2, 366] as warning, and a message "Succ" or "Fail" to remove record from dataset for prediction

getLogAlkalinity <- function (hardness, calcium, conduct, alkal, lower_b, upper_b) {
  #lower_b = 1.2, upper_b = 366, i.e. the log10 boundary values
    msg <- "Succ"
    ALK <- alkal
    # Case none provided
    if(is.null(hardness) & is.null(calcium) & is.null(conduct) & is.null(ALK)) { #MVP but.. or is.na(hardness) & is.na(calcium) etc etc
      # Means check alkalinity bounds and calculate log10 from alkal

      msg <- "Fail: one of hard, Calc, cond, alk is null, bounds exceeded"
      alk_frame <- data.frame(log10(0), msg)
      colnames(alk_frame) <- c("log","msg")
      return (alk_frame)
    }#
    #Case provided, start with calculations
    if(!is.null(hardness) & (!is.na(hardness))) { #' dont use is.na()
      ALK = 4.677 + 0.6393*hardness
      #print(!is.na(hardness))
      #print("hardness not null")
    }
    else
      if(!is.null(calcium) & (!is.na(calcium))) {
        ALK = 14.552 + 1.7606*calcium
        #print("calcium not null")
      }
    else
      if(!is.null(conduct) & (!is.na(conduct))) {
        ALK = 0.3201*conduct - 8.0593
        #print("Conduct not null")
      }

  #print(ALK)
  if(!is.null(ALK) & (!is.na(ALK))){
     ALK.LOG <- log10(ALK)

     #First check Fail conditions  Case 2
     if(alkal<0 ) { # Completely out of bounds, should fail with null value for "alt"
       msg <- "Fail: ALKALINITY bounds exceeded"
       alk_frame <- data.frame(ALK.LOG, msg) # A null value means nothing or zero
       colnames(alk_frame) <- c("log","msg")
       return (alk_frame)
     }
     # Do the normal things
     if( (alkal>=0 & alkal<lower_b) | (alkal >upper_b) ){
         msg <- "Warn: ALKALINITY bounds reached"
         alk_frame <- data.frame(ALK.LOG, msg)
         colnames(alk_frame) <- c("log","msg") # past0 to "log" for name of attribute/parameter
         return (alk_frame) # Return both message and log value
       }
       else {
         alk_frame <- data.frame(ALK.LOG, msg)
         colnames(alk_frame) <- c("log","msg")
         return(alk_frame)
       }
  }
}


#Usage: Apply the mapply multivariate function
    # copy_data$ALKALINITY <- with(copy_data, mapply(getAlkalinity, HARDNESS, CALCIUM, CONDUCTIVITY, ALKALINITY))
    # head(copy_data$ALKALINITY,length(copy_data$ALKALINITY))

# 2.The function getValidEnvInput() returns a dataframe of log10 value and a message, given limit bound [lower_b, upper_b]
#  This function is generic to these input environment variables:
#  MEAN_WIDTH, # MEAN_DEPTH, # SLOPE, # DIST_FROM_SOURCE, #ALTITUDE
#  Use "New GB (43 Group level)
#  Check if it is "Fail" or "Warn", substr(getValidEnvInput(0)[2][[1]], 1,4)
#  BUG: Added fail condition to Fail case 2 :| (env_var<lower_b) on 25-11-2019

getValidEnvInput <- function(env_var, lower_b, upper_b, envName) {
  # lower_b <- 1
  # upper_b <- 590
  msg <- "Succ"
  #First check Fail conditions  Case 1
  if(is.null(env_var) | lower_b<0 | is.null(upper_b) | (env_var<0) ) { #Should fail with null value for "alt"
    msg <- paste("Fail: ",envName," bounds exceeded")
    env_frame <- data.frame(-999, msg) # A null value means nothing or -999
    colnames(env_frame) <- c("log","msg")
    return (env_frame)
  }

  #First check Warn conditions  Case 2
  if( (0<env_var & env_var<lower_b) | (env_var>upper_b)) { # Completely out of bounds, should fail # Check MEAN_WIDTH FAILS
    msg <- paste("Warn: ",envName," bounds exceeded")
    env_frame <- data.frame(log10(env_var), msg) # A null value means nothing or zero
    colnames(env_frame) <- c("log","msg")
    return (env_frame)
  }
  else { # Do the normal things
            env_frame <- data.frame(log10(env_var), msg)
            colnames(env_frame) <- c("log","msg")
            return(env_frame)
      }
  # endif
}# end getValidEnvInput

# 3. Function getAltitude(), change to base 10, Use "New GB (43 Group level)
# This function is deprecated, use getEnvVariableTwoBounds()


getAltitude <- function(alt, lower_b, upper_b) { # CHANGE AZURE
  lower_lowb <- 0
  upper_highb <- 1345
  msg <- "Succ"
  if(is.null(alt) | is.na(alt)) { # Null cant be evaluated, return 0
    msg <- "Fail: Altitude bounds exceeded with null values"
    alt_frame <- data.frame(0, msg)
    colnames(alt_frame) <- c("log","msg")
    return (alt_frame)
  }
  else {##Fail case
      if(alt<lower_lowb | alt>upper_highb){ # Completely out of bounds, should fail with null value for "alt"
        msg <- "Fail: Altitude bounds exceeded"
        alt_frame <- data.frame(0 , msg)
        colnames(alt_frame) <- c("log","msg")
        return (alt_frame)
      }
      else {
            if((alt>=lower_lowb & alt<lower_b) | (alt>upper_b & alt<=upper_highb)){ # Deal with warn cases, otherwise normal cases
                msg <- "Warn: Altitude bounds reached"
                alt_frame <- data.frame(log10(alt), msg)
                colnames(alt_frame) <- c("log","msg")
                return (alt_frame)
            }
            else { # normal case
                  alt_frame <- data.frame(log10(alt), msg)  # msg = "Succ"
                  colnames(alt_frame) <- c("log","msg")
                  return (alt_frame)
            }# end else
    }# end else
   }# end first else
}#end getAltitude


# 3.1.1 getEnvVariable () - validates every environment variable with two bounds: #ALTITUDE,
# The bounds are different: Warn [1,590], Fail [0,1345]. Value is less than 0 or greater than 1345, in which case "Fail".

getEnvVariableTwoBounds <- function(envValue, lower_b, upper_b, lower_lowb, upper_highb, envName ) { # CHANGE AZURE
   msg <- "Succ"
  if(is.null(envValue) | is.na(envValue)) { # Null cant be evaluated, return 0
    msg <- paste("Fail: ",envName," bounds exceeded")
    env_frame <- data.frame(0, msg)
    colnames(env_frame) <- c("log","msg")
    return (env_frame)
  }
  else {##Fail case
    if(envValue<lower_lowb | envValue>upper_highb){ # Completely out of bounds, should fail with null value for "alt"
      msg <- "Fail: Altitude bounds exceeded"
      env_frame <- data.frame(0 , msg)
      colnames(env_frame) <- c("log","msg")
      return (env_frame)
    }
    else {
      if((envValue>=lower_lowb & envValue<lower_b) | (envValue>upper_b & envValue<=upper_highb)){ # Deal with warn cases, otherwise normal cases
        msg <- paste("Warn: ",envName," bounds reached")
        env_frame <- data.frame(log10(envValue), msg)
        colnames(env_frame) <- c("log","msg")
        return (env_frame)
      }
      else { # normal case
        env_frame <- data.frame(log10(envValue), msg)  # msg = "Succ"
        colnames(env_frame) <- c("log","msg")
        return (env_frame)
      }# end else
    }# end else
  }# end first else
}#end getAltitude



# 3.1 getLong - similar to getAltitude
# NOT USED

getLONG <- function(longit, lower_b, upper_b) {
  lower_lowb <- -11
  upper_highb <- 2
  msg <- "Succ"
  if(is.null(longit) | is.na(longit)) { # Null cant be evaluated, return log10(0)
    msg <- "Fail: LONGITITUDE bounds exceeded with null values"
    longit_frame <- data.frame(0, msg)
    return (longit_frame)
  }
  else {##Fail case
    if(longit<lower_lowb | longit>upper_highb){ # Completely out of bounds, should fail with null value for "alt"
      msg <- "Fail: LONGITUDE bounds exceeded"
      longit_frame <- data.frame(0, msg)
      colnames(longit_frame) <- c("_log","_msg")
      return (longit_frame)
    }
    else {
      if((longit>lower_lowb & longit<lower_b) | (longit>upper_b & longit<upper_highb)){ # Deal with warn cases, otherwise normal cases
        msg <- "Warn: LONGITUDE bounds reached"
        longit_frame <- data.frame(longit, msg)
        colnames(longit_frame) <- c("_log","_msg")
        return (longit_frame)
      }
      else { # normal case
        longit_frame <- data.frame(longit, msg)  # msg = "Succ"
        colnames(longit_frame) <- c("_log","_msg")
        return (longit_frame)
      }# end else
    }# end else
  }# end first else
}#end getLONG


# 3.2 getLAT - similar to getAltitude
# Input [54, 55.2]
# NOT USED
getLAT <- function(lat, lower_b, upper_b) {
  lower_lowb <- 49
  upper_highb <- 71
  msg <- "Succ"
  if(is.null(lat) | is.na(lat)) { # Null cant be evaluated, return log10(0)
    msg <- "Fail: LATITUDE bounds exceeded with null values"
    lat_frame <- data.frame(0, msg)
    return (lat_frame)
  }
  else {##Fail case
    if(lat<lower_lowb | lat>upper_highb){ # Completely out of bounds, should fail with null value for "alt"
      msg <- "Fail: LATITUDE bounds exceeded"
      lat_frame <- data.frame(0, msg)
      colnames(lat_frame) <- c("_log","_msg")
      return (lat_frame)
    }
    else {
      if((lat>lower_lowb & lat<lower_b) | (lat>upper_b & lat<upper_highb)){ # Deal with warn cases, otherwise normal cases
        msg <- "Warn: LATITUDE bounds reached"
        lat_frame <- data.frame(lat, msg)
        colnames(lat_frame) <- c("_log","_msg")
        return (lat_frame)
      }
      else { # normal case
        lat_frame <- data.frame(lat, msg)  # msg = "Succ"
        colnames(lat_frame) <- c("_log","_msg")
        return (lat_frame)
      }# end else
    }# end else
  }# end first else
}#end getLAT

# 4. Change "Distance from Source" to base 10
getDistanceFromSource <- function(dist_source) {
  return (log10(dist_source))
}

# 5. Change "mean width of river" to base 10, and check bounds [a,b]. Use "New GB (43 Group level)
getMeanWidth <- function(mean_width) {
  lower_b <- 0.4
  upper_b <- 117
  if(mean_width>=lower_b & mean_width<=upper_b) {
      return (log10(mean_width))
  }
  else
    return (log10(lower_b))  # could be log10(mean(lower_b, upper_b))
}

# 6. Change "mean depth of river" to base 10, check bounds [a,b]. Use "New GB (43 Group level)

getMeanDepth <- function(mean_depth) {
  lower_b <- 1.7
  upper_b <- 300
  if(mean_depth >=lower_b & mean_depth<=upper_b){
      return (log10(mean_depth))
  }
  else
    return (log10(lower_b))
}


# 7. Calculate Substrate , use mapply
# Check that the input files - bould_cob+pebbpes_gr+snd+silt_cl=[97,103], and MSUBST = [-8,8], the latter is hard coded
getSubstrate <- function(bould_cob, pebbles_gr, snd, silt_cl, lower_b, upper_b) {
  TOTSUB <- sum(bould_cob, pebbles_gr, snd, silt_cl)
  msg    <- "Succ"
  MSUBST <- -999 #Error that no calculation takes place
  if(TOTSUB<lower_b | TOTSUB > upper_b) { # Check bounds [97,103] for TOTSUB
      msg  <- "Fail: TOTSUB bounds exceeded for substrate, "
      msubst_frame <- data.frame(MSUBST, msg)
      colnames(msubst_frame) <- c("log","msg")
    return (msubst_frame)
  }
  else { #Check normal way
    if(TOTSUB>=lower_b & TOTSUB <= upper_b) {
      MSUBST <- ((-7.75 * bould_cob) - (3.25 * pebbles_gr) + (2 * snd) + (8 * silt_cl)) / TOTSUB
      #print(MSUBST)

      if(MSUBST< -8.0 | MSUBST >8.0) { # Check if bounds broken, #hard coded for now
          msg  <- "Fail: MSUBST bounds exceeded for substrate"
          msubst_frame <- data.frame(MSUBST, msg)
          colnames(msubst_frame) <- c("log","msg")
        return (msubst_frame )
      }
      else {# Warn case
        if(TOTSUB==lower_b | TOTSUB==upper_b) {
              msg  <- "Warn: bounds reached for substrate"
              msubst_frame <- data.frame(MSUBST, msg)
              colnames(msubst_frame) <- c("log","msg")
          return (msubst_frame )
        }
        else {# Succ case
            msubst_frame <- data.frame(MSUBST, msg)
            colnames(msubst_frame) <- c("log","msg")
          return (msubst_frame )
        }
      }#end else
   }# endif
  }# endelse
}# end getSubstrate

#Test getsubstrate()
# Usage:
    # copy_data$MEAN.SUBSTRATUM <- with(copy_data, mapply(getSubstrate_old, BOULDER_COBBLES, PEBBLES_GRAVEL, SAND, SILT_CLAY))

# 8. Get Discharge
getRawDischarge <- function(discharge) {

  return (discharge)
}

# 9. Get Raw Alkalinity
getRawAlkalinity <- function (alkalinity){

  return (alkalinity)
}

# 10.1 Get latLong, use package rnrfa for function osg_parse
#library(rnrfa)

getLatLong <- function (nat_grid_ref, easting, northing,coordsys_latlon) {
  lat_long    <- NA
  concat_bng  <- paste(nat_grid_ref, substr(easting,1,3), substr(northing,1,3), sep="")
  lat_long    <- osg_parse (concat_bng, coord_system = coordsys_latlon)
  return (lat_long)
}

#10.12

getLatLong_AZURE <- function (nat_grid_ref, easting, northing,coordsys_latlon) {
  lat_long    <- NA
  concat_bng  <- paste(nat_grid_ref, substr(easting,1,3), substr(northing,1,3), sep="")
  lat_long    <- osg_parse (concat_bng, CoordSystem = coordsys_latlon)
  return (lat_long)
}

# This function gets the Easting, NMorthiung coordinatres and returns Latitude/Longitude pairs.
# It uses package library "sf". sf is Simple Features.

# 10.2
getLatLong_NI <- function (x, y) {
  # Declare the lat_long
  lat_long_all <- NULL
  # It is assumed the East/North have 5 digits, just add a ZERO at the end
  xx <- as.numeric(paste0(x, 0))
  yy <- as.numeric(paste0(y, 0))

  # Loop throuhg the codes to extract the Easting and Northing

  for(i in 1:length(x)) {
    xy = data.frame(easting_x=xx[i], northing_y=yy[i]) # Edited, just to give site identifier
    # 1. create sf object in Irish National Grid (CRS 29903)
    irish.ng <- st_as_sf(xy, coords = c("easting_x", "northing_y"), crs = 29903)
    lat_long <- st_transform(irish.ng, crs = 4326)
    results <- c(lat_long$geometry[[1]][[2]], lat_long$geometry[[1]][1])
    lat_long_all  <- rbind(lat_long_all,results)
  }
  #Remove rownames
  rownames(lat_long_all) <- c()
  #names columns appropriately as "Easting",and "Northing"
  colnames(lat_long_all) <- c("Latitude","Longitude")
  return (as.data.frame(lat_long_all))
}


# 11.3. Get latLong, use package rnrfa for function osg_parse
#library(rnrfa)
# NOT USED

getLatLongNI <- function (nat_grid_ref, easting, northing,coordsys_latlon) {
  lat_long    <- NA
  concat_bng  <- paste(nat_grid_ref, substr(easting,1,3), substr(northing,1,3), sep="")
  print(concat_bng)
  # Add "I" for NOrthern Ireland
  lat_long    <- osg_parse (paste0("I",concat_bng), CoordSystem = coordsys_latlon)
  #print(lat_long)
  return (lat_long)
}

# NOte :
# Can view all Grid references with "SE..." from  allStations <- catalogue(), using the rNRFA package
# this <- filter(allStations, grepl("SE", gridReference, fixed = TRUE))$gridReference
# sort(this, decresing  = TRUE)
# this
# BNG <- osg_parse(concatenatedNGR, CoordSystem = "BNG")
# lat long is RICT explanatory variables.
# Easting and Northing used for temperature lookup
# lat.long <- osg_parse(concatenatedNGR, CoordSystem = "WGS84")

# 11.4 Get BNG

getBNG <- function (nat_grid_ref, easting, northing,coordsys_bng) {
  lat_long    <- NA
  concat_bng  <- paste(nat_grid_ref, substr(easting,1,3), substr(northing,1,3), sep="")
  new_bng    <- osg_parse (concat_bng, coord_system = coordsys_bng)
  return (new_bng)
}

# 11.45 Get BNG

# 11. Get BNG - This version works with AZURE

getBNG_AZURE <- function (nat_grid_ref, easting, northing,coordsys_bng) {
  lat_long    <- NA
  concat_bng  <- paste(nat_grid_ref, substr(easting,1,3), substr(northing,1,3), sep="")
  new_bng    <- osg_parse (concat_bng, CoordSystem = coordsys_bng)
  return (new_bng)
}
# 12. sPECIAL CASE OF discharge from VELOCITY
# Note the bounds for Discharge are (0,10) and [1,9], and bounds for velocity ar [1]
# Note that lower_b==0, upper_b==10, low_b==1, up_b==9

getLogDischarge <- function (depth, width, discharge, velocity, lower_b, upper_b, low_b, up_b) {
  vel_lower <- 1
  vel_upper_b <- 5
  msg <- "Succ"
  #check value of velocity to be ZERO or not supplied (NULL), if so assign 1
  if(is.na(velocity) | velocity==0)
    velocity <- 1

  if( is.na(discharge)) { # Case no discharge provided
    # Check velocity
    if(!is.na(velocity) & !is.na(depth) & !is.na(width) ){ #do calculations
      #Check bounds for non-null velocity
      if(velocity <vel_lower | velocity>vel_upper_b){#out of bounds
        msg <- "Fail: bounds exceeded velocity with null values"
        disch_frame <- data.frame(0, msg)
        colnames(disch_frame) <- c("log","msg")
        return (disch_frame)
      }
      else  #Warn case
        if (velocity ==vel_lower | velocity==vel_upper_b){
          msg <- "Warn: bounds reached velocity"
          disch <- depth/100 * width * velocity/100
          disch_frame <- data.frame(log10(disch), msg)
          colnames(disch_frame) <- c("log","msg")
          return (disch_frame)
        } else {# normal velocity, do calculations
            disch <- depth/100 * width * velocity/100
            disch_frame <- data.frame(log10(disch), msg)
            colnames(disch_frame) <- c("log","msg")
            return (disch_frame)
      }#end else
    }#end if
    else { # In case all are NA including depth, width and velocity
      msg <- "Fail: bounds exceeded with  everything null"
      disch_frame <- data.frame(0, msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
  }#
  else{ # Normal service resumes, discharge not null and provided
    #check bounds, then do a log if ok
    if(lower_b>discharge  | upper_b < discharge){
      msg <- "Fail: bounds exceeded discharge"
      print("Failing discharge")
      disch_frame <- data.frame(log10(discharge), msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
    else # Warnings bounds reached
      if ((discharge>lower_b & discharge<low_b)  | (discharge>up_b & discharge<=upper_b)){
        msg <- "Warn: bounds reached discharge"
        disch_frame <- data.frame(log10(discharge), msg)
        colnames(disch_frame) <- c("log","msg")
        return (disch_frame)
    }# Success
    else {# All normal, discharge provided
      disch_frame <- data.frame(log10(discharge), msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
  }#end else
}

# 12.1 sPECIAL CASE OF discharge from VELOCITY for NI model. We Can factor out parameters to make same as GB model
# Note the bounds for Discharge are [1,8], and bounds for velocity ar [1,10]
#Overall range [0,10] needs to be incorporated??? 30/11/2019

getLogDischarge_NI <- function (depth, width, discharge, velocity, lower_b, upper_b) {
  vel_lower <- 1
  vel_upper_b <- 10
  msg <- "Succ"

  if(is.na(velocity) | velocity==0)
    velocity <- 1

  if( is.na(discharge)) { # Case no discharge provided
    # Check velocity
    if(!is.na(velocity) & !is.na(depth) & !is.na(width) ){ #do calculations
      #Check bounds for non-null velocity
      if(velocity <vel_lower | velocity>vel_upper_b){#out of bounds
        msg <- "Fail: bounds exceeded velocity with null values"
        disch_frame <- data.frame(log10(0), msg)
        colnames(disch_frame) <- c("log","msg")
        return (disch_frame)
      }
      else  #Warn case
        if (velocity ==vel_lower | velocity==vel_upper_b){
          msg <- "Warn: bounds reached velocity"
          disch <- depth/100 * width * velocity/100
          disch_frame <- data.frame(log10(disch), msg)
          colnames(disch_frame) <- c("log","msg")
          return (disch_frame)
        } else {# normal velocity, do calculations
          disch <- depth/100 * width * velocity/100
          disch_frame <- data.frame(log10(disch), msg)
          colnames(disch_frame) <- c("log","msg")
          return (disch_frame)
        }#end else
    }#end if
    else { # In case all are NA including depth, width and velocity
      msg <- "Fail: bounds exceeded with  everything null"
      disch_frame <- data.frame(log10(0), msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
  }#
  else{ # Normal service resumes, discharge not null and provided
    #check bounds, then do a log if ok
    if(lower_b>discharge  | upper_b < discharge){
      msg <- "Fail: bounds exceeded discharge"
      disch_frame <- data.frame(log10(discharge), msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
    else # bounds reached
      if (lower_b==discharge  | upper_b == discharge){
        msg <- "Warn: bounds reached discharge"
        disch_frame <- data.frame(log10(discharge), msg)
        colnames(disch_frame) <- c("log","msg")
        return (disch_frame)
      }# Success
    else {# All normal, discharge provided
      disch_frame <- data.frame(log10(discharge), msg)
      colnames(disch_frame) <- c("log","msg")
      return (disch_frame)
    }
  }#end else
}

#13. Special writeTo File function - replaces old file if one exists aleardy
# 3.2. Write to a file using the pathname
writeToFile <- function (toWriteFrame, pathname, filename) {

  if( file.exists(paste0(file =pathname,filename))) {
    file.remove(file = paste0(pathname,filename))
    write.csv(toWriteFrame, file = paste0(pathname,filename))
  } else {
    write.csv(toWriteFrame, file = paste0(pathname,filename))
  }
}


#14. Function to return 5 character long Northing, Easting with leading ZEROS if short of 5 characters
# Alternative function used previously - this only replaces one zero, we need multiple application for more zero replacements
#         mycopyEasting <- apply(mycopyEasting, 1, function(x) ifelse(nchar(x)<5, paste0(0, x), x)) # Use "1" in apply() for row
#

getCorrectCodes <- function (coordinates) {

  localCopy <- coordinates
  for(i in 1:length(coordinates)){
    if(5-nchar(localCopy[i])==1)
      localCopy[i] <- paste0(0, localCopy[i])
    if(5-nchar(localCopy[i])==2)
      localCopy[i] <- paste0(00, localCopy[i])
    if(5-nchar(localCopy[i])==3)
      localCopy[i] <- paste0(000, localCopy[i])
    if(5-nchar(localCopy[i])==4)
      localCopy[i] <- paste0(0000, localCopy[i])
    if(5-nchar(localCopy[i])==5)
      localCopy[i] <- paste0(0000, localCopy[i])
  }

  return (localCopy)
}


#15. Validate the longitude

getLongitude <- function(longtde, lower_b, upper_b) { # CHANGE AZURE
  lower_lowb <-  -11.0
  upper_highb <- 2
  msg <- "Succ"
  if(is.null(longtde) | is.na(longtde)) { # Null cant be evaluated, return 0
    msg <- "Fail: lLongitude bounds exceeded with null values"
    #print(msg)
    longtde_frame <- data.frame(0, msg)
    colnames(longtde_frame) <- c("log","msg")
    return (longtde_frame)
  }
  else {##Fail case
    if(longtde<=lower_lowb | longtde>=upper_highb){ # Completely out of bounds, should fail with null value for "longtde"
      msg <- "Fail: Longitude bounds exceeded"
      longtde_frame <- data.frame(0, msg)
      colnames(longtde_frame) <- c("log","msg")
      #print(msg)
      return (longtde_frame)
    }
    else {
      if((longtde>lower_lowb & longtde<lower_b) | (longtde>upper_b & longtde<upper_highb)){ # Deal with warn cases, otherwise normal cases
        msg <- "Warn: Longitude bounds reached"
        longtde_frame <- data.frame( longtde, msg)
        colnames(longtde_frame) <- c("log","msg")
        #print(msg)
        return (longtde_frame)
      }
      else { # normal case
        longtde_frame <- data.frame(longtde, msg)  # msg = "Succ"
        colnames(longtde_frame) <- c("log","msg")
        #print(msg)
        return (longtde_frame)
      }# end else
    }# end else
  }# end first else
}#end getLongitude


#16. Validate the latitude

getLatitude <- function(latde, lower_b, upper_b) { # CHANGE AZURE
  lower_lowb <-  49
  upper_highb <- 71
  msg <- "Succ"
  if(is.null(latde) | is.na(latde)) { # Null cant be evaluated, return 0
    msg <- "Fail: Latitude bounds exceeded with null values"
    #print(msg)
    latde_frame <- data.frame(0, msg)
    colnames(latde_frame) <- c("log","msg")
    return (latde_frame)
  }
  else {##Fail case
    if(latde<=lower_lowb | latde>=upper_highb){ # Completely out of bounds, should fail with null value for "latde"
      msg <- "Fail: Latitude bounds exceeded"
      latde_frame <- data.frame(0, msg)
      colnames(latde_frame) <- c("log","msg")
      #print(msg)
      return (latde_frame)
    }
    else {
      if((latde>lower_lowb & latde<lower_b) | (latde>upper_b & latde<upper_highb)){ # Deal with warn cases, otherwise normal cases
        msg <- "Warn: Latitude bounds reached"
        latde_frame <- data.frame( latde, msg)
        colnames(latde_frame) <- c("log","msg")
        #print(msg)
        return (latde_frame)
      }
      else { # normal case
        latde_frame <- data.frame(latde, msg)  # msg = "Succ"
        colnames(latde_frame) <- c("log","msg")
        #print(msg)
        return (latde_frame)
      }# end else
    }# end else
  }# end first else
}#end getLatitude


# 17. Validates all the input Environment Variables by changing their out of scope values
# to preset valueds

# # Data validation and conversion
# 0. Validation of Various Env. variables before transformation NO NEED for this as its been done by Helpfulfunctions.
# ALKALINITY : GB model range: If the value is less than 1.2 or greater than 366 = "Warn"
# If Altitude = 0 then replace value with 1

# If Altitude = 0 then replace value with 1
# If Distance from source is <0.1 then replace value with 0.1
# If Width is <0.1 then replace value with 0.1
# If Depth is <1 then replace value with 1
# If Discharge category is 0 then replace value with 1
# If Alkalinity is <0.1 then replace value with 0.1
# If Slope is =0 then replace value with 0.1


validateEnvData <- function (data) {

  data$ALTITUDE <- ifelse(data$ALTITUDE==0,1,  data$ALTITUDE)

  # DIST_FROM_SOURCE
  data$DIST_FROM_SOURCE <- ifelse(data$DIST_FROM_SOURCE<0.1, 0.1, data$DIST_FROM_SOURCE)

  # MEAN_WIDTH
  data$MEAN_WIDTH <- ifelse(data$MEAN_WIDTH<0.1, 0.1, data$MEAN_WIDTH)

  # MEAN_DEPTH
  data$MEAN_DEPTH <- ifelse(data$MEAN_DEPTH<1.00, 1, data$MEAN_DEPTH)

  # DISCHARGE VALIDATION
  data$DISCHARGE <- ifelse(data$DISCHARGE==0,1.0, data$DISCHARGE)


  #ALKALINITY VALIDATION
  data$ALKALINITY <- ifelse(data$ALKALINITY <0.1, 0.1, data$ALKALINITY)

  #SLOPE validation
  data$SLOPE <- ifelse(data$SLOPE ==0, 0.1, data$SLOPE)

  return (data)
}


# 18 Reading big files with fread, ... is to add parametric needs e.g. skip=.., stringsAsFactors

readCSVFile <- function (path, ...) {
  return (fread(path, ...)) #

}


# EndGrpProb_Replacement - function inputs x endgroups, any of 1-43, and appends a 'p' to make p1,p2,..,p43
# It then selects the probabilities pi = p1,p2,p3, .. from any SITE i of the predicted values y[i,pi] = final.predictionsxxx dataframe[i,pi]
# usage : site1 <- EndGrpProb_Replacement (a$EndGrp_Probs,final.predictors_try2, 1) # for site k=1

EndGrpProb_Replacement <- function(x,final_data,k) {
   allprobss <- sapply(x, function (x)  ifelse(x==x,paste0("p",x), 0)) # apply to all values in EndGroup
   #allprobss <- noquote(allprobss)
   c <- as.double(unname(final_data[k,c(allprobss)]))
  return (c) # How do you select with duplicates included???

}


#
# groupSitesFunction - Calculates endgroup probabilities with index values per group of similar site, taxa level, season code and furse_Code
# returns a site with calculated site index, and binds to the whole site in main data frame when called

groupSitesFunction <- function (allSites, k, siteindex, b1) {
  perGroupSite  <- b1[allSites[k,1:5]$Season_Code==b1$Season_Code & allSites[k,1:5]$TL==b1$TL & allSites[k,1:5]$Furse_Code==b1$Furse_Code,]
  a             <-  data.frame(t(colSums(perGroupSite$mlist_endGrps[[1]]*perGroupSite[,11:18])))
  siteName      <- paste0("TST-", paste0(siteindex,"-R"))
  siteX         <- cbind(data.frame(siteName,perGroupSite[1,1:10]),a) # #Add this site to AllSites - for final putput
  return (siteX)
}
