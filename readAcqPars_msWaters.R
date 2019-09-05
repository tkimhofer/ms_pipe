# read in MS metadata from original Waters files (not CDFs)
# WatersDataDir='/Volumes/Torben_1/DiaObesity/DiaObesity_Lipid_POS.PRO/'

readWaters_acquisitonPars=function (WatersDataDir, plot=T) 
{
  require(plyr)
  require(ggplot2)
  require(reshape2)
  
  if(grepl('/$', WatersDataDir)){WatersDataDir=gsub('/$', '', WatersDataDir)}
  warnDef <- options("warn")$warn
  warnRead <- options(warn = -1)
  datapath <- WatersDataDir
  
  print(datapath)
  
  
  # searches for parameter files in header (data, time)
  hfile=list.files(path=datapath, pattern='_HEADER', full.names=T, recursive=T)
  
  # searches for inlet files in header (injection volume)
  ifile=list.files(path=datapath, pattern='_INLET', full.names=T, recursive=T)
  
  # searches for extern files (various parameter values)
  efile=list.files(path=datapath, pattern='_extern', full.names=T, recursive=T)
  
  print(length(efile))
  print(efile[1:10])
  
  if(length(hfile)!=length(ifile) | length(efile)!=length(ifile)){
    print('Number of paramter files do not match - corrupt files?'); return(NULL)
  }

  L <- length(hfile)

    out=lapply(rbind(1:L), function(i, hf=hfile, jf=ifile, ef=efile){
      print(ef[1:10])
      # extern files (various parameter values)
      con <- file(ef[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      aLine=gsub('\\(\xb0\\C)', '', aLine)
      myV <- strsplit(aLine, "\\t")
      close(con)
      myV=lapply(myV, function(x){
       x[which(x !='')]
      })
      myV=myV[sapply(myV, length)==2]
      parsSpec=as.data.frame(do.call(rbind, myV[1:length(myV)]), stringsAsFactors = F)

      # inlet files in header (injection volume)
      con=file(jf[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      close(con)
      tt=strsplit(aLine[grep('Injection Volume', aLine)], '   -   ')
      parsSpec=rbind(parsSpec, tt[[1]])
      
      # header files (acquisition status parameters)
      con <- file(hf[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      close(con)
      
      idx=grep('Acquired Name|Date|Time|Sample Description|Bottle Number|Solvent Delay|Cal Temperature|Inlet Method|MS Method|Tune Method', aLine)
      aLine=aLine[idx]
      strsplit(gsub('^\\$\\$ ', '', aLine), ': ')
      myV=strsplit(gsub('^\\$\\$ ', '', aLine), ': ')
      acquSpec=as.data.frame(do.call(rbind, myV), stringsAsFactors = F)
     
      return(rbind(acquSpec, parsSpec))
      #return(c(acquSpec, parsSpec))
    })
  
  out1=lapply(out,function(y){ 
    y=as.data.frame(t(y),stringsAsFactors=FALSE);
    colnames(y)=y[1,]; 
    y=y[-1,]})
  
  out1=rbind.fill(out1)
  
  out1$date=paste(out1$`Acquired Date`, out1$`Acquired Time`)
  head(out1$date)
  out1$date=as.POSIXlt(out1$date, tz = "", format,
             format = c("%d-%b-%Y %H:%M:%OS"),
             optional = FALSE)
  head(out1$date)
  
  # add rack position
  test=strsplit(out1$`Bottle Number`, ':')
  # test[which(sapply(test, length)<5)]
  df=as.data.frame(do.call(rbind, test))
  colnames(df)=c('RackInst', 'RackPosID')
  out1=cbind(out1, df)
  
  if(plot==T){ # only the status files an acquisition finished date
    # get the time of acquisition and run order
    out1$date=as.POSIXct(out1$date)
    out1$DateFinishedWeekday=weekdays(out1$date, abbreviate = T)
    out1$HoursFromStart=as.numeric(abs(min(out1$date)-out1$date)/60/60)
    out1$RunOrder=rank(as.numeric(out1$date-min(out1$date)), ties.method = 'min')
    
    out1$DateFinishedWeekday=factor(out1$DateFinishedWeekday, levels=c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'))
    g1=ggplot(out1, aes(date, RackPosID))+
      geom_point(aes(colour=`Sample Description`))+ #=='study QC serum'
      theme_bw()+facet_grid(~DateFinishedWeekday, scales='free')+
      coord_flip()
    plot(g1)
  }
  
  
  return(out1)
}
