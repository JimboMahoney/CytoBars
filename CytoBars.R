## Draws one (or up to three if certain conditions are met) bar graphs of the DCs from an FCS file
## Intended to quickly summarise DC intensities of each marker above 190 BCKG, showing 
## colour-coded high CVs (i.e. how much the intensities vary for that parameter) and 
## high (>20k solution modeDCs) intensity - e.g. possible contaminants or overstaining of Ir / Pt
## Also shows % of events that contained that parameter above the bar
## Ir staining should be ~ 1000 counts in event mode DC
## If two or three plots are produced, use the Next/Prev button on the Plots panel to view
## Optional conversion factor to convert from Event mode DCs -> Solution Mode DCs

# Data Import from file chosen by user

library(svDialogs)
# Get user input for file
testfile<-dlg_open()
# Convert to string value
testfile <- capture.output(testfile)[7]
{
  
  if ((testfile)=="character(0)")
    stop("File input cancelled")
  
  #Remove invalid characters from file input location
  testfile <- gsub("[\"]","",testfile)
  testfile<-substring (testfile,5)
  
  #Set file and directory
  filename <- basename (testfile)
  dir <- dirname (testfile)
  
  # Set working directory accoding to file chosen
  setwd(dir)
  
  library(flowCore)
  
  # this read.FCS() function imports the flow data:
  raw_fcs<-read.FCS(filename, alter.names = TRUE)
  
  
  # Preparation work for arcsinh transform (columns is also used later for naming changes)
  # Create list of parameters
  columns<-colnames(raw_fcs)
  # Remove "Time" column to avoid it being transformed
  columns<-setdiff(columns,"Time")
  # Remove "Cell_Length" and Gaussians column to avoid it being transformed
  columns<-setdiff(columns,"Event_length")
  columns<-setdiff(columns,"Cell_length")
  columns<-setdiff(columns,"Center")
  columns<-setdiff(columns,"Offset")
  columns<-setdiff(columns,"Width")
  columns<-setdiff(columns,"Residual")
  ## Remove FSC and SSC
  removefscssc<-grep("FSC|SSC",columns,value=TRUE)
  columns<-columns[! columns %in% removefscssc]
  
  
  
  # Read data into a data frame
  FCSDATA <- as.data.frame(exprs(raw_fcs))
  
  ## DIFFERENT FROM CytobankGraphs ##
  ## Estimate conversion factor for 20,000 DCs in 1 second solution mode -> Event mode
  ## Used later to determine if a marker is "positive"
  cutoffFactor <- 1000000/(mean(FCSDATA[,grep("length",colnames(FCSDATA))])*13)
  ## END OF DIFFERENCE
  
  
  ############ Optional Data Transform section
  
  #Remove comments from code lines to transform using asinh
  ## Automatically estimate the logicle transformation based on the data
  #lgcl <- estimateLogicle(raw_fcs, channels = c(columns))
  ## transform  parameters using the estimated logicle transformation
  #raw_fcs_trans <- transform(raw_fcs, lgcl)
  # Load into data frame
  #FCSDATA <- as.data.frame(exprs(raw_fcs_trans))
  
  ########### End of optional Data Transform section
  
  
  
  #Remove unnecessary parameter text
  names(FCSDATA)[-1] <- sub("Di", "", names(FCSDATA)[-1])
  names(FCSDATA)[-1] <- sub("Dd", "", names(FCSDATA)[-1])
  # Create list of channel / parameter descriptions 
  params<-parameters(raw_fcs)[["desc"]]
  # Replace parameters with descriptions, keeping things like Time, Event Length unchanged
  colnames(FCSDATA)[!is.na(params)] <- na.omit(params)
  
  # Determine whether data is CyTOF or Flow by presence of FSC
  # isflow will be 0 for a CyTOF or greater than 1 if flow
  isflow <-sum(grep("FSC",colnames(FCSDATA)))
  # Determine whether data is pre CyTOF 3 (Helios) by presence of "Cell_length", rather than "Event_length"
  isCyTOF2 <-sum(grep("Cell_length",colnames(FCSDATA)))
  
  ## Remove Time, Event_Length & Gaussian Parameters
  removecolumns <- c("Event_length", "Center", "Offset", "Width", "Residual", "Cell_length")
  FCSDATA <- FCSDATA[,!(names(FCSDATA) %in% removecolumns)]
  
  
  ## Remove FSC and SSC
  library(tidyverse) 
  FCSDATA <- FCSDATA %>% select(-contains("FSC"))
  FCSDATA <- FCSDATA %>% select(-contains("SSC"))
  
  
  
  # Get number of cell events (based on "193" - i.e. Iridium)
  if(isflow==0){
    cellevents<-as.data.frame(apply(FCSDATA, 2, function(c)sum(c!=0)))
    colnames(cellevents) <-c("Events")
    # Note that this only works correctly because "Time" has been removed by a previous step - otherwise the position would be wrong.
    irpos<-grep("193",columns)
    cellevents<-cellevents$Events[irpos]
    kcellevents <-round(cellevents/1000,0)
  }
  
  
  #For converting FCS time to mins - flow uses 10ms units, CyTOF uses ms
  if(isflow>0){
    div = (60*100)
  }else{
    div = (60*1000)
  }
  
  # Find total acquisition time
  maxtime<-round(max(FCSDATA$Time)/div,2)
  
  # Now that we have the total time, we can calculate the number of cell events/sec
  if(isflow==0){
    eventspersec <- round(cellevents/maxtime/60,0)
  }
  
  
  # Create number formatted list of intensity values and event counts
  # Changed EventsecList to total number of events
  # Needed to change "parameter" column to "Marker"
  Meanintensitylist <- c(format(c(round(colMeans(FCSDATA)),1),big.mark = ",",trim=TRUE))
  # For some reason, another item is added, so we need to remove that
  EventList <- c(format(c(round((colSums(FCSDATA !=0)),0),trim=TRUE)))
  # Remove the last row that is added by format
  Meanintensitylist<-Meanintensitylist[-length(Meanintensitylist)]
  EventList<-EventList[-length(EventList)]
  # Create data frame for labels to print mean intensity on plots
  datalabels <- data.frame(
    Meanintensity=c(Meanintensitylist),
    Marker = c(colnames(FCSDATA)),
    Events = c(EventList)
  )
  
  
  #Calculate size of dataset - NOT NEEDED FOR THIS?
#  DataSizeM <- (ncol(FCSDATA)*nrow(FCSDATA))/1000000
  #Subsample if greater than 10,000
#  if (DataSizeM>2.5){
    #using random 10% of original rows
    #FCSDATA <- FCSDATA[sample(nrow(FCSDATA),nrow(FCSDATA)/10),]
    #OR
    #Subsample using a number of random rows, where the number is defined by numrows
#    numrows <- 5000
#    FCSDATA <- FCSDATA[sample(nrow(FCSDATA),numrows),]
#  }
  
## Not all of the below is needed for this, but the code works in so far as adding the additional parameter names.  
  
  # Add a blank to the columns list to match its length to that of FCSDATA (i.e. the time row)
  columns<-columns<-append(columns,"Time",after=0)
  # Add back the original marker names
  datalabels[,"OrigMarkers"]<-columns
  # Remove Di / Dd
  datalabels$OrigMarkers <- sub("Di", "", datalabels$OrigMarkers)
  datalabels$OrigMarkers <- sub("Dd", "", datalabels$OrigMarkers)
  
  # This is needed for pre-Helios data to ensure we don't mess with the parameter names
  if (isCyTOF2>1){
    
    # Remove other symbols
    datalabels$OrigMarkers <- gsub("[[:punct:]]", "", datalabels$OrigMarkers)
    
    # Create a function to extract the last n characters from a string
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    # Extract only last 5 characters (i.e the element and mass) - this is clumsy and doesn't work well for flow data, which may have longer names for markers
    # But I can't figure out a way to remove duplicate text
    datalabels$OrigMarkers<-substrRight(datalabels$OrigMarkers,5)
    
    # Compare columns and keep only original markers if they are different
    datalabels$OrigMarkers<-ifelse(datalabels$Marker==datalabels$OrigMarkers,"",paste("/",datalabels$OrigMarkers))
    # Replace parameters column with orignal marker names / parameters
    datalabels[,"Marker"]<-paste(datalabels$Marker,datalabels$OrigMarkers)
  } #End of flow data name comparison / CyTOF paramater rename loop
  
  # Remove the OrigMarkers column as it's no longer needed
  datalabels<-datalabels[,-4]
  
  
  # Make sure the FCSDATA matches the datalabels
  colnames(FCSDATA)<-datalabels$Marker
  
  #Trim the trailing whitespace added by paste
  colnames(FCSDATA)<-trimws(colnames(FCSDATA),"r")
  datalabels$Marker<-trimws(datalabels$Marker,"r")
  
  
  # Remove Time from labels
  datalabels <- datalabels[!(rownames(datalabels) %in% "Time"),]
  # Change rownames to numeric 
  rownames(datalabels) <- 1:nrow(datalabels)
  # Change parameters to factors to control facet order
  datalabels$Marker<-as.factor(datalabels$Marker)

## DIFFERS FROM CytobankGraphs here

  # Calculate approx. cutoff based on 20k DC counts in solution mode
  # CutoffFactor is calculated on the initial FCSDATA further up
  cutoff<-20000/cutoffFactor
  
  # OR use all parameters, irrespective of whether they are named or not
  SDs <- sapply((FCSDATA),sd)
  meanDCs <- colMeans(FCSDATA)
  markers <- colnames(FCSDATA)
  
  
  
  # Put into a data frame
  alldata <- list(Marker=markers,MeanDC=meanDCs,SD=SDs)
  summarydata <- as.data.frame(alldata)
  
  
  # Remove Time  
  summarydata<-summarydata[-(grep("Time",summarydata$Marker)),]
  
  # Different from CytobankGraphs
  # Add # of events to this
  summarydata$Events<-as.numeric(as.character(datalabels$Events))
  # Convert to % of events
  summarydata$Events<-round(summarydata$Events/max(summarydata$Events)*100,0)
  
  # Change rownames to numeric - not essential I don't think. but I prefer it
  rownames(summarydata) <- 1:nrow(summarydata)
  
  
  # Get value of 190 BCKG channel - we don't care about anything if it's lower than this.
  meanBCKG <-mean(FCSDATA[["190BCKG"]])
  # If we can't find this marker, set to 0
  if (is.na(meanBCKG)==TRUE){
    meanBCKG=0
  }
  
  
  # Three colour ramp
  colfunc <- colorRampPalette(c( "black","purple4", "red"))
  
  library(ggplot2)
  
  ggplot(summarydata, aes(x=Marker,y=MeanDC, fill=SD/MeanDC)) +
    # Apply colour ramp
    scale_fill_gradientn(colours = colfunc(128)) +
    # Rename fill to "CV"
    labs(fill="CV") +
    # Use actual values from MeanDC column
    geom_bar(stat="identity") +
    # Rotate text 90 degrees
    theme(axis.text.x = element_text(angle = 90)) +
    # Don't reorder the Markers (default is alphanumeric)
    scale_x_discrete(limits=summarydata$Marker) +
    # Log scale for Y  show numbers instead of notation
    #scale_y_continuous(trans="log10", labels=scales::comma)+
    # Or linear
    scale_y_continuous(labels=scales::comma)+
    ylab("Intensity") +
    xlab("Channel") +
    # Add horizontal line at the DC value that concerns us
    geom_hline(yintercept=cutoff,color="red",linetype="dashed") +
    # Zoom plot to only the values of interest
    coord_cartesian(ylim=c(1, max(summarydata$MeanDC)))+
    ggtitle(paste(filename, " - All Markers","=",nrow(summarydata))) +
    # Rotate X Axis Text 45 degrees
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    # Add frequency of events as label
    geom_text(size=3,vjust=-1,aes(label=paste(Events,"%")))
  
  
} # End of File cancel loop  

## Create and plot only data less than cutoff  

lowmarkers <- subset(summarydata,MeanDC<cutoff)

# Don't draw second plot if all markers or data is flow or there are no lowmarkers or file cancelled.
if ((testfile)=="character(0)" || nrow(lowmarkers)/nrow(summarydata)==1 || nrow(lowmarkers)==0 || isflow>0){
  stop("Only drawing one plot")
}else{
  
  ggplot(lowmarkers, aes(x=Marker,y=MeanDC, fill=SD/MeanDC)) +
  # Apply colour ramp
  scale_fill_gradientn(colours = colfunc(128)) +
  # Rename fill to "CV"
  labs(fill="CV") +
  # Use actual values from MeanDC column
  geom_bar(stat="identity") +
  # Rotate text 90 degrees
  theme(axis.text.x = element_text(angle = 90)) +
  # Don't reorder the Markers (default is alphanumeric)
  scale_x_discrete(limits=lowmarkers$Marker) +
  # Log scale for Y  show numbers instead of notation
  #scale_y_continuous(trans="log10", labels=scales::comma)+
  ylab("Intensity") +
  xlab("Channel") +
  # Add horizontal line at the DC value that concerns us
  geom_hline(yintercept=cutoff,color="red",linetype="dashed") +
  # Zoom plot to only the values of interest
  coord_cartesian(ylim=c(0, max(lowmarkers$MeanDC)))+
  scale_y_continuous(labels=scales::comma)+
  ggtitle(paste(filename, "- Markers below cutoff (",round(cutoff,1),")","=",nrow(lowmarkers))) +
  # Rotate X Axis Text 45 degrees
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  # Add frequency of events as label
  geom_text(size=3,vjust=-1,aes(label=paste(Events,"%")))
}

## Create and plot only data greater than cutoff  

highmarkers <- subset(summarydata,MeanDC>cutoff)

# Don't draw second plot if all markers or data is flow or there are no highmarkers or file cancelled.
if ((testfile)=="character(0)" || nrow(highmarkers)/nrow(summarydata)==1 || nrow(highmarkers)==0 || isflow>0){
  stop("Only drawing one plot")
}else{
  
  ggplot(highmarkers, aes(x=Marker,y=MeanDC, fill=SD/MeanDC)) +
  # Apply colour ramp
  scale_fill_gradientn(colours = colfunc(128)) +
  # Rename fill to "CV"
  labs(fill="CV") +
  # Use actual values from MeanDC column
  geom_bar(stat="identity") +
  # Rotate text 90 degrees
  theme(axis.text.x = element_text(angle = 90)) +
  # Don't reorder the Markers (default is alphanumeric)
  scale_x_discrete(limits=highmarkers$Marker) +
  # Log scale for Y  show numbers instead of notation
  #scale_y_continuous(trans="log10", labels=scales::comma)+
  ylab("Intensity") +
  xlab("Channel") +
  # Add horizontal line at the DC value that concerns us
  geom_hline(yintercept=cutoff,color="red",linetype="dashed") +
  # Zoom plot to only the values of interest
  coord_cartesian(ylim=c(cutoff, max(highmarkers$MeanDC)))+
  scale_y_continuous(labels=scales::comma)+
  ggtitle(paste(filename, "- Markers above cutoff (",round(cutoff,1),")","=",nrow(highmarkers))) +
  # Rotate X Axis Text 45 degrees
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  # Add frequency of events as label
  geom_text(size=3,vjust=-1,aes(label=paste(Events,"%")))
}