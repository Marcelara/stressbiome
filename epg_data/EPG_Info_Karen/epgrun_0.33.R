#' Electrical Penetration Graph data conversion
#' version 0.33 updates:
#' Waveforms that are not performed by an aphid are treated as:
#' (1) zero (0) value for occ, num, sum, sing.e1
#' (2) missing value (NA) for mean, max and prop.e1
#' (3) end of recording (epgstopmin) for lat
#'
#' version 0.32 updates:
#' (1) Fix of bug in wilcox and kruskal tests (some data were occasionally classified as factor)
#' (2) Chi square test works now correct. N.B. In case of zero observations among groups (e.g. 100% and 100% occ.e1), no test can be conducted and NA values are produced for P value and test statistic.
#' (3) Pds (wave 3) variables are now included and are not mutually exclusive for C (wave 2), in order that C is not interrupted by pds.
#' (4) Pds do not interrupt Cs and Css anymore
#'
#' \code{epgrun} converts Electrical Penetration Graph (EPG) annotation data from
#' 'EPG Stylet +a' to EPG variables.
#'
#' @param loc Vector with directory name(s) where the .ANA output files are stored
#'   from 'EPG Stylet +a' (www.epgsystems.eu). ANA files belonging to different groups
#'   (e.g. plant lines or treatments) should be saved in separate directories.
#'   The folder name will be used as group name.
#' @param epgstart Number indicating the start time (in hours) of the EPG
#'   recordings. Should be > 0 if the first part of the recording should
#'   be excluded. By default 0.
#' @param epgstop Number indicating the end time (in hours) of the EPG recordings.
#'   Should be larger than epgstart. This is the total duration of the actual
#'   recording, or a shorter duration if the last part of the recording should
#'   be excluded. Maximum 16.
#' @param epgplot Logical scalar indicating whether plots should be made. TRUE
#'   or FALSE. By default TRUE.
#' @param epgtable Logical scalar indicating whether a table with mean values and
#'  standard error should be made and univariate statistical tests should be performed.
#'  TRUE or FALSE. By default TRUE.
#' @param epgbin Logical scalar indicating whether timebins should be calculated.
#'  TRUE or FALSE. By default TRUE.
#' @param epgname String with project name which will appear in output
#'   file names. By default a project name will be created with the first and
#'   the last group name and the total number of groups.
#'
#' @return Maximum 7 output files will be saved to the working directory:
#' \itemize{
#'      \item rawdata.csv (raw data with user-defined end time of recording)
#'      \item probes.csv (extended information per probe, including previous and next waveform)
#'      \item variables.csv (EPG summary variables per insect, 59 variables. See detailed description below)
#'      \item epgtable.csv (if epgtable = TRUE: table with mean plus minus standard error, and outcomes
#'      of univariate statistical tests on summary variables. Chi square tests are performed on
#'      occurrence data (e.g. occ.e1). Other variables are analyzed with Wilcoxon signed rank tests
#'      "wilcox.test()" (2 groups) or Kruskal-Wallis tests "kruskal.test()" (> 2 groups). All tests
#'      are performed in default settings. Samples or groups with missing values for a variable are not
#'      included in a test, tests are not executed in case < 2 groups have data, or - in case of a chi
#'      square test- if there are no values > 0. 'test_statistic' = value of test statistic
#'      (e.g. W for Wilcoxon, ), 'n_' = number of biological replicates per group for the respective
#'      variable, Pvalue is NA in case of missing or equal values, standard error is NA in case of only
#'      <= 1 biological replicate)
#'      \item boxplots.pdf (if epgplot = TRUE: boxplots of total recording)
#'      \item timebins.csv (if epgbin = TRUE: probes categorized in timebins of 1 hour.
#'      If probes are spanning several hours, they will be chopped up)
#'      \item timebinplots.pdf (if epgbin & epgplot = TRUE: lineplots)
#'}
#'
#'This function is meant to facilitate EPG data processing and analysis. Correct usage and
#'interpretation of the data and statistical tests are the responsibility of the user of this package.
#'
#' Description of waveforms:
#' \itemize{
#'     \item waveform 1 - Non probing
#'     \item waveform 2 - Pathway
#'     \item waveform 4 - Salivation in the phloem
#'     \item waveform 5 - Phloem feeding
#'     \item waveform 6 - Penetration difficulties
#'     \item waveform 7 - Xylem ingestion
#'     \item waveform 8 - Potential drops
#'     \item waveform 11 - User defined
#'     \item waveform 12 - User defined
#'}
#'
#' N.B. It is not required to have annotated all waveforms described above.
#' N.B.B. Waveforms are considered to be mutually exclusive, except waveform 8, which is not mutually exclusive from waveform 2.
#' This means, that waveform 2 is not considered to end when waveform 8 starts,
#' and that their duration may include waveform 8 events. In general, waveform 8 is
#' used to annotate potential drops.
#'
#' Description of EPG variables (adapted from www.epgsystems.eu):
#' \itemize{
#'     \item sum.np - Total duration non probing (min)
#'     \item lat.c - Latency to first C from start of recording (min)
#'     \item sum.c - Total duration C (min)
#'     \item num.c - Total number of C
#'     \item num.cs - Total number of C < 3 min ('short probes')
#'     \item num.css - Total number of C < 0.5 min
#'     \item mean.c - Mean duration C (min)
#'     \item occ.pd - Does insect perform pd? (1=yes, 0=no)
#'     \item num.pd - Total number of potential drops
#'     \item sum.pd - Total duration of potential drops (min)
#'     \item sum.pd - Mean duration of potential drops (min)
#'     \item rate.pd - Number of potential drops per min of C (num.pd divided by sum.c)
#'     \item occ.e1 - Does insect perform E1? (1=yes, 0=no)
#'     \item lat.e1 - Latency to first E1 from start of recording (min)
#'     \item lat.e1p - Latency to first E1 from start of first probe (min)
#'     \item num.e1 - Total number of E1
#'     \item sum.e1 - Total duration E1 (min)
#'     \item mean.e1 - Mean duration E1 (min)
#'     \item max.e1 - Maximum duration E1 (min)
#'     \item prop.e1 - Proportion of time spent on E1 in phloem phase (E1+E2, excluding single E1s) (%)
#'     \item occ.sing.e1 - Does insect perform single E1? (1=yes, 0=no)
#'     \item sing.e1 - Number of single E1s
#'     \item num.e1f - Number of E1 before E2 starts (excluding single E1s)
#'     \item num.e1m - Number of E1s in between E2
#'     \item num.e1l - Number of E1s at end of phloem feeding
#'     \item sum.e1f - Total duration E1 before starting E2 (min)
#'     \item sum.e1m - Total duration E1 in between E2 (min)
#'     \item sum.e1l - Total duration E1 when ending E2 (min)
#'     \item prop.e1f - Proportion of time spent on E1 before starting E2 in phloem phase (E1+E2) (%)
#'     \item prop.e1m - Proportion of time spent on E1 in between E2 in phloem phase (E1+E2) (%)
#'     \item prop.e1l - Proportion of time spent on E1 when ending E2 in phloem phase (E1+E2) (%)
#'     \item occ.e2 - Does insect perform E2? (1=yes, 0=no)
#'     \item occ.e2s - Does insect perform sustained E2 (>10 min)? (1=yes, 0=no)
#'     \item lat.e2 - Latency to first E2 from start of recording (min)
#'     \item lat.e2p - Latency to first E2 from start of first probe (min)
#'     \item lat.e2s - Latency to first sustained E2 (>10 min) from start of recording (min)
#'     \item lat.e2sp - Latency to first sustained E2 (>10 min) from start of first probe (min)
#'     \item num.e2 - Number of E2
#'     \item num.e2s - Number of sustained E2 (>10 min)
#'     \item sum.e2 - Total duration of E2 (min)
#'     \item sum.e2s - Total duration of sustained E2 (>10 min) (min)
#'     \item mean.e2 - Mean duration of E2 (min)
#'     \item mean.e2t - Mean total duration of E2 phase including potential e1m (min)
#'     \item mean.e2s - Mean duration of sustained E2 (>10 min) (min)
#'     \item mean.e2st - Mean total duration of E2s phase including potential e1m (min)
#'     \item max.e2 - Maximum duration E2 (min)
#'     \item occ.f - Does insect perform F? (1=yes, 0=no)
#'     \item num.f - Number of F
#'     \item sum.f - Total duration F (min)
#'     \item mean.f - Mean duration F (min)
#'     \item occ.g - Does insect perform G? (1=yes, 0=no)
#'     \item num.g - Number of G
#'     \item sum.g - Total duration G (min)
#'     \item mean.g - Mean duration G (min)
#'     \item occ.w11 - Does insect perform waveform 11? (1=yes, 0=no)
#'     \item lat.w11 - Latency to first waveform 11 from start of first probe (min)
#'     \item num.w11 - Number of waveform 11
#'     \item sum.w11 - Total duration waveform 11 (min)
#'     \item mean.w11 - Mean duration waveform 11 (min)
#'     \item occ.w12 - Does insect perform waveform 11? (1=yes, 0=no)
#'     \item lat.w12 - Latency to first waveform 12 from start of first probe (min)
#'     \item num.w12 - Number of waveform 12
#'     \item sum.w12 - Total duration waveform 12 (min)
#'     \item mean.w12 - Mean duration waveform 12 (min)
#'     \item N.B. If a waveform does not occur it will have a missing value (NA) for 'lat',
#'     'sum' and 'mean'.
#'     }
#'
#' @examples
#' # Define working directory
#' setwd("C:/Myworking/directory/")
#'
#' # Convert EPG data. In this case 2 plant lines and 8 h recordings
#' epgrun(loc = c("C:/folder/subfolder/plantline1", "C:/folder/subfolder/plantline2"),
#' epgstop = 8)
#'
#' # Output files are saved to the working directory

# Start function -----------------------------------------------

#' @export
epgrun <- function(loc, epgstart, epgstop, epgtable, epgplot, epgbin, epgname)
{

  # Preparations -------------------------------------------------

  ## User environment
  # Save old
  old_options <- getOption("stringsAsFactors")
  old_mfrow <- par("mfrow")
  old_mar <- par("mar")
  # Restore at end
  on.exit(options(stringsAsFactors = old_options), add = TRUE)
  on.exit(par(old_mfrow), add = TRUE)
  on.exit(par(old_mar), add = TRUE)
  # Adjust for read.table() and read.csv()
  options(stringsAsFactors = TRUE)
  # Define separator loc (directory .ANA files)
  if (length(grep("/", loc)) < 1) {
    locsep <- "\\"
  } else {
    locsep <- "/"
  }

  ## Set defaults
  # epgstart
  if (missing(epgstart)) {
    epgstart <- 0
  }

  # epgplot
  if (missing(epgplot)) {
    epgplot <- TRUE
  }

  # table
  if (missing(epgtable)) {
    epgtable <- TRUE
  }

  # epgbin
  if (missing(epgbin)) {
    epgbin <- TRUE
  }

  # epgname
  outp <- numeric(0)

  if (missing(epgname)) {
    for (i in 1:length(loc)) {
      df1 <- numeric(0)
      foldernames <- numeric(0)
      foldernames <- strsplit(loc[i], as.character(locsep), fixed = TRUE)
      df1 <- unlist(foldernames)[length(unlist(foldernames))]
      outp <- c(outp, df1)
    }
    if (length(outp) > 1) {
      epgname <- paste(outp[1], outp[length(outp)], length(loc), "groups", sep = "-")
    } else {
      epgname <- paste(outp[1])
    }}

  ## Warnings messages
  if (epgstop <= epgstart) {
    print("Warning: epgstop is equal to or smaller than epgstart")
  }


  if ((epgstop - epgstart) > 16) {
    print("Warning: Total EPG duration is more than 16 hours.
          Please check epgstart and epgstop (both in hours)")
  }


  # Start subfunctions----------------------------------------------------------------------

  # Raw-----------------------------------------------------------
  # Read .ANA files
  x <- numeric(0)

  for (i in 1:length(loc)) {
    anaobj <- sapply(dir(path = loc[i], pattern = ".ANA$", full.names = TRUE),
                     read.table, fileEncoding = "UTF-16")

    for (j in 1:dim(anaobj)[2]) {
      df1 <- numeric(0)
      df1 <- cbind(data.frame(unlist(rep(colnames(anaobj)[j], length(anaobj[[1, j]])))),
                   data.frame(unlist(anaobj[1, j])),
                   data.frame(unlist(anaobj[2, j]),
                              data.frame(unlist(anaobj[3, j]))))
      df1[,3] <- gsub(",", ".",df1[,3]) #in case of "," as decimal separator, change into "."
      df1 <- rbind(df1, c(levels(df1[,1]), 1, 0, NA))
      df1 <- rbind(df1, c(levels(df1[,1]), 0, epgstop*3600, NA))
      foldernames <- strsplit(loc[i], as.character(locsep), fixed = TRUE)
      df1$group <- rep(unlist(foldernames)[length(unlist(foldernames))], nrow(df1))
      x <- rbind(df1, x)
    }
  }

  colnames(x) <- c("file", "wave", "start", "volt", "group")

  # Add column with unique number to each recording
  x$id_rec <- as.numeric(x$file)

  # Order rows in ascending order of time
  x$start <- as.numeric(x$start)
  x <- x[order(x$id_rec, x$start),]

  # Order columns
  x <- x[, c(1, 6, 5, 2:4)]

  # Save as csv file
  write.table(x, paste(epgname, "_rawdata_", "0-", as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)

  # Probes------------------------------------------------------------
  outp <- numeric(0)

  # Convert sec or hours to min
  x$start <- as.numeric(x$start)/60
  epgstopmin <- epgstop*60
  epgstartmin <- epgstart*60

  # Make factors
  x$file <- factor(x$file)
  x$id_rec <- factor(x$id_rec)
  x$id_rec <- factor(x$id_rec)
  x$id_rec <- factor(x$id_rec)
  x$group <- factor(x$group)
  x$wave <- factor(x$wave)

    # Add column with unique id per probe per aphid
  id.i <- unique(x$id_rec)
  id_probe <- NULL
  id_probe_all <- NULL

  for (i in id.i) {
    id_probe <- NULL
    if(sum(x$id_rec == i) != 0) {  #for each individual rec
      id_probe = c(1:sum(x$id_rec == i)) #number of rows for each rec. vector 1:max probe number
      id_probe_all = c(id_probe_all, id_probe)
    }
  }
  x$id_probe <- id_probe_all
  x$id_probe <- as.numeric(x$id_probe)

  ## Adjust start time to epgstart
  # Temporary column
  x$lab <- x$start - epgstartmin
  df <- x
  outp <- matrix(nr=0, nc=1)

  for (i in id.i) {
    if (epgstartmin > 0) {
      df1 <- matrix(nr=0, nc=1)
      df2 <- matrix(nr=0, nc=1)
      df3 <- matrix(nr=0, nc=1)
      df1 <- subset(df, df$id_rec == i) # select each aphid
      df2 <- subset(df1, df1$lab < 0) # probes before epgstart
      startwave <- subset(df2, df2$lab == max(df2$lab)) # last probe before epgstart
      startwave[, 5] <- epgstartmin # adjust start time last probe to epgstart
      df3 <- subset(df1, df1$lab >= 0) # probes after epgstart
      df3 <- rbind(startwave, df3)
      outp <- rbind(outp, df3)
    }
  }

  for (i in 1) {
    if (epgstart > 0) {
      x <- outp[, 1:(ncol(outp)-1)]
    }
  }

  ## Remove probes starting after epgstop
  x$lab <- ifelse(x$start > epgstopmin, 1, 0)
  x <- subset(x, x$lab != 1)
  x <- x[,-ncol(x)]

 ## Previous and next waveforms
  # a. Define previous and next wave
  x$wave <- as.numeric(levels(x$wave))[x$wave]
  x$id_rec <- as.numeric(levels(x$id_rec)) [x$id_rec]

  df <- x
  outp <- matrix(nr=0, nc=1)
  outp2 <- matrix(nr=0, nc=1)
  df2 <- numeric(0)
  df3 <- numeric(0)
  df4 <- numeric(0)
  df5 <- numeric(0)
  probnum <- numeric(0)

  # make column next wave
  for (i in id.i) {
    df2 <- cbind(df[which(df$id_rec == i), ])
    probnum <- as.numeric(dim(df2)[1])
    if(probnum > 1) {
      df3 <- c(df2[2:probnum, 4], 0) #column 4=waveform; paste next waveform
      outp <- c(outp, df3)
    }
    else {
      if (probnum == 1) {
        outp <- c(outp, 0)
      }
    }
  }

  x$nextwave <- outp

  # make column previous wave
  for (i in id.i) {

    df4 <- cbind(df[which(df$id_rec == i), ])
    probnum <- as.numeric(dim(df4)[1])
    if(probnum > 1){
      df5 <- c(0, df4[1:(probnum-1), 4]) #column 4=waveform; paste previous waveform
      outp2 <- c(outp2, df5)
    }
    else {
      if (probnum == 1) {
        outp2 <- c(outp2, 0)
      }}}
  x$prevwave <- outp2

  ## determine end time of waveform
  df <- x
  outp = matrix(nr = 0, nc = 1)
  df2 = numeric(0)
  df3 = numeric(0)
  df4 = numeric(0)
  probnum = numeric(0)

  for (i in id.i) {
    df2 <- cbind(df[which(df$id_rec == i), ])
    probnum <- as.numeric(dim(df2)[1])
    if (probnum>1) {
      df3 <- c(df2[2:probnum, 5], epgstopmin) # 5=start time
      outp <- c(outp, df3)
    }
    else {
      if (probnum==1) {
        outp <- c(outp, epgstopmin)
      }}}

  x <- cbind(x, outp)
  colnames (x)[ncol(x)] <- "stop"

  #c. duration waveform
  x$dur <- x$stop-x$start
  # to remove potential double events for timepoint zero
  x <- subset(x,x$dur>0)

  # name waveform
  x$wavename <- ifelse(x$wave == 1, "non probing",
                       ifelse(x$wave == 2, "pathway",
                              ifelse(x$wave == 4, "salivation phloem",
                                     ifelse(x$wave == 5, "phloem feeding",
                                            ifelse(x$wave == 6, "penetration difficulties",
                                                   ifelse(x$wave == 7, "xylem ingestion",
                                                          ifelse(x$wave == 8, "potential drop",
                                                                 x$wave)))))))

  # Order columns
  x <- x[,c(1:3, 7, 4, 12, 6, 5, 10, 11, 8, 9)]

  # Save as csv file
  write.table(x, paste(epgname, "_probes_", as.character(epgstart), "-", as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)

  #Sum-----------------------------------------------------------------------
  # Summary variables--------------------------------------------------------------

  ### make dataframe with the recording id's and their respective group
  # file names and their respective numbers
  fnames <- data.frame(cbind(1:length(levels(x$file)),levels(x$file)))
  colnames(fnames) <- c("id_rec","file")
  # groups and their respective group numbers
  grnames <- data.frame(cbind(1:length(levels(x$group)),levels(x$group)))
  colnames(grnames) <- c("group_nr","group")
  # recordings and their respective group numbers
  recs <- tapply(x$group, x$id_rec, "[",1)
  x$id_rec <- factor(x$id_rec)
  recs_df <- data.frame(cbind(row.names(recs),recs))
  colnames(recs_df) <- c("id_rec","group_nr")
  # recordings and their respective group and file names
  recs_df <- merge(recs_df,grnames,by="group_nr",all.x=T,all.y=F)
  recs_df <- merge(recs_df,fnames,by="id_rec",all.x=T,all.y=F)
  recs_df <- recs_df[,c(4,1,3)]

  # Subset per waveform: Duration per event
  x1 <- ifelse(x$wave==1,x$dur,0)
  x2 <- ifelse(x$wave==2,x$dur,0)
  x2.lat <- ifelse(x$wave==2,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording
  x4 <- ifelse(x$wave==4,x$dur,0)
  x4.lat <- ifelse(x$wave==4,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording
  x4.sing <- ifelse(x$wave!=4,0,
                    ifelse(x$wave==4 & x$nextwave==5,0,
                           ifelse(x$wave==4 & x$prevwave==5,0,x$dur)))
  x5 <- ifelse(x$wave==5,x$dur,0)
  x5.lat <- ifelse(x$wave==5,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording
  x5s <- ifelse(x$wave==5 & x$dur>10,x$dur,0)
  x5s.lat <- ifelse(x$wave==5 & x$dur>10,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording
  x6 <- ifelse(x$wave==6,x$dur,0)
  x7 <- ifelse(x$wave==7,x$dur,0)
  x8 <- ifelse(x$wave==8,x$dur,0)
  x11 <- ifelse(x$wave==11,x$dur,0)
  x11.lat <- ifelse(x$wave==11,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording
  x12 <- ifelse(x$wave==12,x$dur,0)
  x12.lat <- ifelse(x$wave==12,x$start,epgstopmin) #10,000 for non-existing events#NEW: end of recording

  # Subset per waveform: Number of events per waveform
  x1.num <- ifelse(x$wave==1,1,0)
  x2.num <- ifelse(x$wave==2 & x$prevwave != 8,1,0) # Wave 2's separated by wave 8 (pd) should be considered as one event)
  df.x2s <- x
  df.x2s$lab <- ifelse(df.x2s$wave==2 & df.x2s$dur<3,1,0)
  x2s.num <- ifelse(df.x2s$lab==1 & df.x2s$prevwave!=8,1,0)
  #x2s.num <- ifelse(x$wave==2 & x$dur<3,1,0)
  df.x2s$lab2 <- ifelse(df.x2s$wave==2 & df.x2s$dur<0.5,1,0)
  x2ss.num <- ifelse(df.x2s$lab2==1 & df.x2s$prevwave!=8,1,0)
  #x2ss.num <- ifelse(x$wave==2 & x$dur<0.5,1,0)
  x4.num <- ifelse(x$wave==4,1,0)
  x5.num <- ifelse(x$wave==5,1,0)
  x5s.num <- ifelse(x$wave==5 & x$dur>10,1,0)
  x6.num <- ifelse(x$wave==6,1,0)
  x7.num <- ifelse(x$wave==7,1,0)
  x8.num <- ifelse(x$wave==8,1,0)
  x11.num <- ifelse(x$wave==11,1,0)
  x12.num <- ifelse(x$wave==12,1,0)

  # TOTAL DURATION
  # duration np (1)
  # duration C (2)
  # duration E1 (4)
  # duration E2 (5)
  # duration sustained E2 > 10 min (5)
  # duration xylem uptake G (7)

  sum.np <- tapply(x1,x$id_rec,sum)
  sum.c <- tapply(x2,x$id_rec,sum)
  sum.pd <- tapply(x8,x$id_rec,sum)
  sum.e1 <- tapply(x4,x$id_rec,sum)
  max.e1 <- tapply(x4,x$id_rec,max)
  sum.sing.e1 <- tapply(x4.sing,x$id_rec,sum)
  sum.e2 <- tapply(x5,x$id_rec,sum)
  sum.e2s <- tapply(x5s,x$id_rec,sum)
  max.e2 <- tapply(x5,x$id_rec,max)
  sum.f <- tapply(x6,x$id_rec,sum)
  sum.g <- tapply(x7,x$id_rec,sum)
  sum.w11 <- tapply(x11,x$id_rec,sum)
  sum.w12 <- tapply(x12,x$id_rec,sum)

  # NUMBER
  # C (2)
  # C < 3 min (2)
  # E1 (4)
  # E2 (5)
  # E2 > 10 min (5)
  num.c <- tapply(x2.num,x$id_rec,sum)
  num.cs <- tapply(x2s.num,x$id_rec,sum)
  num.css <- tapply(x2ss.num,x$id_rec,sum)
  num.pd <- tapply(x8.num,x$id_rec,sum)
  num.e1 <- tapply(x4.num,x$id_rec,sum)
  num.e2 <- tapply(x5.num,x$id_rec,sum)
  num.e2s <- tapply(x5s.num,x$id_rec,sum)
  num.f <- tapply(x6.num,x$id_rec,sum)
  num.g <- tapply(x7.num,x$id_rec,sum)
  num.w11 <- tapply(x11.num,x$id_rec,sum)
  num.w12 <- tapply(x12.num,x$id_rec,sum)
  occ.pd <- ifelse(num.pd<1,0,1)
  occ.e1 <- ifelse(num.e1<1,0,1)
  occ.e2 <- ifelse(num.e2<1,0,1)
  occ.e2s <- ifelse(num.e2s<1,0,1)
  occ.f <- ifelse(num.f<1,0,1)
  occ.g <- ifelse(num.g<1,0,1)
  occ.w11 <- ifelse(num.w11<1,0,1)
  occ.w12 <- ifelse(num.w12<1,0,1)
  rate.pd <- num.pd/sum.c

  ## number of single E1 (4)
  # definition single E1: an E1 event which is not preceded or followed by E2.
  x$sing.e1 <- ifelse(x$wave!=4,0,
                      ifelse(x$wave==4 & x$nextwave==5,0,
                             ifelse(x$wave==4 & x$prevwave==5,0,1)))
  sing.e1 <- tapply(x$sing.e1,x$id_rec,sum)
  occ.sing.e1 <- ifelse(sing.e1<1,0,1)

  # Proportion of time spend on salivation during phloem phase (excluding single e1s):
  prop.e1 <- 100*((sum.e1-sum.sing.e1)/(sum.e1+sum.e2))
  prop.e1 <- ifelse(prop.e1==100,NA,prop.e1) # only %E1 in E in case E2 takes place

  # latency from start recording to 1st E1 and E2
  lat.e1 <- tapply(x4.lat,x$id_rec,min)
  lat.e2 <- tapply(x5.lat,x$id_rec,min)
  lat.e2s <- tapply(x5s.lat,x$id_rec,min)

  # latency from start of first probe to 1st E1 and E2
  lat.c <- tapply(x2.lat,x$id_rec,min)
  lat.w11 <- tapply(x11.lat,x$id_rec,min)
  lat.w12 <- tapply(x12.lat,x$id_rec,min)
  lat.e1p <- lat.e1-lat.c
  lat.e2p <- lat.e2-lat.c
  lat.e2sp <- lat.e2s-lat.c
  lat.w11p <- lat.w11-lat.c
  lat.w12p <- lat.w12-lat.c

  ## E1 separate for:
  # E1 preceding E2 'first' (e1f)
  # E1 in between E2 'middle' (e1m)
  # E1 at the end of E2 'last' (e1l)

  # duration
  e1f <- ifelse(x$wave==4 & x$prevwave!=5 & x$nextwave==5,x$dur,0)
  e1m <- ifelse(x$wave==4 & x$prevwave==5 & x$nextwave==5,x$dur,0)
  e1l <- ifelse(x$wave==4 & x$prevwave==5 & x$nextwave!=5,x$dur,0)
  # per aphid
  sum.e1f <- tapply(e1f,x$id_rec,sum)
  sum.e1m <- tapply(e1m,x$id_rec,sum)
  sum.e1l <- tapply(e1l,x$id_rec,sum)
  # number
  e1f.num <- ifelse(x$wave==4 & x$prevwave!=5 & x$nextwave==5,1,0)
  e1m.num <- ifelse(x$wave==4 & x$prevwave==5 & x$nextwave==5,1,0)
  e1l.num <- ifelse(x$wave==4 & x$prevwave==5 & x$nextwave!=5,1,0)
  # per aphid
  num.e1f <- tapply(e1f.num,x$id_rec,sum)
  num.e1m <- tapply(e1m.num,x$id_rec,sum)
  num.e1l <- tapply(e1l.num,x$id_rec,sum)
  #prop
  prop.e1f <- 100*(sum.e1f/(sum.e1+sum.e2))
  prop.e1m <- 100*(sum.e1m/(sum.e1+sum.e2))
  prop.e1l <- 100*(sum.e1l/(sum.e1+sum.e2))

  # MEAN DURATION
  # mean duration C
  # mean duration E2
  # mean duration sustained E2
  # mean duration E1
  mean.c <- sum.c/num.c
  mean.pd <- sum.pd/num.pd
  mean.e1 <- sum.e1/num.e1
  mean.e2 <- sum.e2/num.e2
  mean.e2t <- (sum.e2+sum.e1m)/num.e2
  mean.e2s <- sum.e2s/num.e2s
  mean.e2st <- (sum.e2s+sum.e1m)/num.e2s
  mean.f <- sum.f/num.f
  mean.g <- sum.g/num.g
  mean.w11 <- sum.w11/num.w11
  mean.w12 <- sum.w12/num.w12

  ## Merge variables with independent variables
  sum.np_df <- data.frame(cbind(row.names(sum.np),sum.np))
  colnames(sum.np_df) <- c("id_rec","sum.np")


  y <- data.frame(cbind(
    sum.np_df,lat.c,
    num.c,num.cs,num.css,sum.c,mean.c,
    occ.pd,num.pd,sum.pd,mean.pd,rate.pd,
    occ.e1,lat.e1,lat.e1p,
    num.e1,sum.e1,mean.e1,max.e1,
    prop.e1,occ.sing.e1,sing.e1,
    num.e1f,num.e1m,num.e1l,sum.e1f,sum.e1m,sum.e1l,
    prop.e1f,prop.e1m,prop.e1l,
    occ.e2,occ.e2s,lat.e2,lat.e2p,lat.e2s,lat.e2sp,
    num.e2,num.e2s,sum.e2,sum.e2s,
    mean.e2,mean.e2t,mean.e2s,mean.e2st,max.e2,
    occ.f,num.f,sum.f,mean.f,
    occ.g,num.g,sum.g,mean.g,
    occ.w11,lat.w11p,num.w11,sum.w11,mean.w11,
    occ.w12,lat.w12p,num.w12,sum.w12,mean.w12))
  y <- merge(recs_df,y,by="id_rec",all.y=TRUE)


  # If aphids do not perform an event, then sum=NA (and not zero) and latency=NA ##NEW sum=0, lat=epgstopmin (see above)
  # if pd=0 dan NA (not annotated)
  #y$sum.e1 <- ifelse(y$occ.e1==0,NA,y$sum.e1)
  #y$sum.e2 <- ifelse(y$occ.e2==0,NA,y$sum.e2)
  #y$sum.e2s <- ifelse(y$occ.e2s==0,NA,y$sum.e2s)
  #y$sum.f <- ifelse(y$occ.f==0,NA,y$sum.f)
  #y$sum.g <- ifelse(y$occ.g==0,NA,y$sum.g)
  #y$sum.w11 <- ifelse(y$occ.w11==0,NA,y$sum.w11)
  #y$sum.w12 <- ifelse(y$occ.w12==0,NA,y$sum.w12)
  #y$lat.c <- ifelse(y$lat.c<=0,NA,y$lat.c)
  #y$lat.w11p <- ifelse(y$occ.w11==0,NA,y$lat.w11p)
  #y$lat.w12p <- ifelse(y$occ.w12==0,NA,y$lat.w12p)
  #y$lat.e1 <- ifelse(y$lat.e1<=0,NA,y$lat.e1)
  #y$lat.e1p <- ifelse(y$lat.e1p<=0,NA,y$lat.e1p)
  #y$lat.e2 <- ifelse(y$lat.e2<=0,NA,y$lat.e2)
  #y$lat.e2p <- ifelse(y$lat.e2p<=0,NA,y$lat.e2p)
  #y$lat.e2s <- ifelse(y$lat.e2s<=0,NA,y$lat.e2s)
  #y$lat.e2sp <- ifelse(y$lat.e2sp<=0,NA,y$lat.e2sp)
  y$max.e1 <- ifelse(y$max.e1 == 0, NA, y$max.e1)
  y$max.e2 <- ifelse(y$max.e2 == 0, NA, y$max.e2)
  y$num.pd <- ifelse(y$num.pd == 0, NA, y$num.pd)
  y$rate.pd <- ifelse(y$rate.pd == 0, NA, y$rate.pd)
  y$sum.pd <- ifelse(y$occ.pd == 0, NA, y$sum.pd)

  # Order rows according to group
  y <- y[order(as.numeric(as.character(y$id_rec))),]

  # Save as csv file
  write.table(y, paste(epgname, "_variables_", as.character(epgstart), "-", as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)

  # Epgstat-mean and standard error---------------------------------------------------------------------

  if (epgtable == "TRUE") {

  outp <- numeric(0)
  tabley <- numeric(0)
  for (i in 4:ncol(y)) { # for each variable:
      df1 <- numeric(0)
      df2 <- numeric(0)
      df1 <- y[!is.na(y[,i]),] # exclude samples with NA value
      df2 <- paste(round(tapply(as.numeric(as.character(df1[,i])), df1$group, mean), 1),
                   intToUtf8(177),
            round(tapply(as.numeric(as.character(df1[,i])), df1$group, sd)/sqrt(tapply(as.numeric(as.character(df1[,i])), df1$group, length)), 1),
            sep = " ")
      outp <- rbind(outp, df2)
  }

  # Add groups as column names
  suppressWarnings(outp <- data.frame(outp))
  colnames(outp) <- row.names(tapply(as.numeric(as.character(y[,i])), y$group, mean))
  # Add variable names as a first column
  outp$variable <- colnames(y)[4:ncol(y)]
  outp <- outp[,c(ncol(outp),1:(ncol(outp)-1))]
  tabley <- outp
  # Change class of data of all data columns
   for (i in 2:(ncol(tabley))) {
     tabley[,i] <- as.character(tabley[,i])
   }

  # occ data
  yocc <- y[,c(3,grep("occ",colnames(y)))]
  outp <- numeric(0)
  for (i in 2:ncol(yocc)) { # for each variable:
    df1 <- numeric(0)
    df1 <- paste(round(100*(tapply(as.numeric(as.character(yocc[,i])), yocc$group, sum)/
                              tapply(as.numeric(as.character(yocc[,i])), yocc$group, length)), 0), #percentage of insects performing it
                 "%", sep = "")
    outp <- rbind(outp, df1)
  }
  # Add groups as column names
  suppressWarnings(outp <- data.frame(outp))
  colnames(outp) <- row.names(tapply(as.numeric(as.character(y[,i])), y$group, mean))
  # Add variable names as a first column
  outp$variable <- colnames(yocc)[2:ncol(yocc)]
  outp <- outp[,c(ncol(outp),1:(ncol(outp)-1))]
  # Change class of data of all data columns
  for (i in 2:(ncol(outp))) {
    outp[,i] <- as.character(outp[,i])
  }

  # Replace occ data with percentage
  tabley[grep("occ", tabley$variable),] <- outp
  rownames(tabley) <- 1:nrow(tabley)

  # Epgstat-statistical tests-------------------------------------------------------------------

  # Two groups - Mann Whitney U pairwise comparisons-----
  df1 <- numeric(0)
  df2 <- numeric(0)
  varnam <- numeric(0)
  repls <- numeric(0)
  stat <- numeric(0)
  test <- numeric(0)


  if (length(levels(y$group)) == 2) {
    yy <- y[,-grep("occ",colnames(y))] # all variables except occurrence data (needs other test)
    for (i in 4:ncol(yy)) { # for each variable
      df1 <- yy[!is.na(yy[,i]),c(3,i)] # save variable (i) and group (col3) in new df, exclude samples with NA value
      if (nrow(df1) > 0 & length(unique(df1$group)) == 2) { # if there is data of 2 groups

        if (class(df1[,2]) == "factor") {
          suppressWarnings(w <- wilcox.test(as.numeric(levels(df1[,2]))[df1[,2]] ~ df1$group)) # Mann Whitney U test (column 2 is variable excl NA values)
        } else {
          suppressWarnings(w <- wilcox.test(as.numeric(df1[,2]) ~ df1$group)) # Mann Whitney U test (column 2 is variable excl NA values)
        }

        df2 <- rbind(df2,w$p.value)
        df2 <- data.frame(df2)
        df2$sign <- ifelse(df2$df2<0.001,"***",
                           ifelse(df2$df2<0.01,"**",
                                  ifelse(df2$df2<0.05,"*","ns")))
        varnam <- c(varnam,as.character(colnames(yy)[i])) # variable name
        repls <- rbind(repls, rbind(tapply(df1[,2],df1[,1],length))) # number of replicates per group
        test <- c(test, as.character(w$method))
        stat <- c(stat, as.numeric(w$statistic))
      }
    }

    yo <- y[,c(1:3,grep("occ",colnames(y)))]
    for (i in 4:ncol(yo)) {
      if (sum(yo[,i]) > 0) { # at least one group should have value > 1
        suppressWarnings(chi <- chisq.test(rbind(tapply(yo[,i], y$group, sum),
                                                 tapply(yo[,i], y$group, length)-tapply(yo[,i], y$group, sum))))
      df2 <- rbind(df2,chi$p.value)
      df2 <- data.frame(df2)
      df2$sign <- ifelse(df2$df2<0.001,"***",
                         ifelse(df2$df2<0.01,"**",
                                ifelse(df2$df2<0.05,"*","ns")))
      varnam <- c(varnam,as.character(colnames(yo)[i])) # variable name
      repls <- rbind(repls, rbind(tapply(yo[,i],yo[,3],length))) # number of replicates per group
      test <- c(test, as.character(chi$method))
      stat <- c(stat, as.numeric(chi$statistic))

      }
    }

    # Add and reorder columns
    df2$variable <- varnam
    df2 <- cbind(df2, test, stat, repls)
    df2 <- df2[,c(3,4,5,1,2,6:(5+length(levels(y$group))))]

    # names of columns
    v1 <- numeric(0)
    for (i in 1:length(levels(y$group))) {
      naml <- numeric(0)
      naml <- paste("n_", as.character(levels(y$group)[i]), sep="")
      v1 <- c(v1, naml)
    }
    colnames(df2) <- c("variable","test", "test_statistic", "Pvalue","sign", strsplit(v1, " "))

    # Round P values
    df2$Pvalue <- round(df2$Pvalue,4)

    # Merge with mean and standard errors, use temporary column 'sort' to order the rows
    tabley$sort <- rownames(tabley)
    tabley <- merge(tabley, df2, by = "variable", all.x = TRUE)
    tabley <- tabley[order(as.numeric(tabley$sort)),]
    tabley <- tabley[,-grep("sort", colnames(tabley))]

    # Save outcomes
    write.table(tabley, paste(epgname, "_epgtable_", as.character(epgstart), "-",
                           as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)
  }

  # >2 groups - Kruskal Wallis groupwise comparison
  df1 <- numeric(0)
  df2 <- numeric(0)
  varnam <- numeric(0)
  repls <- numeric(0)
  stat <- numeric(0)
  test <- numeric(0)


  if (length(levels(y$group)) > 2) {
    yy <- y[,-grep("occ",colnames(y))] # all variables except occurrence data (needs other test)
    for (i in 4:ncol(yy)) {
      df1 <- yy[!is.na(yy[,i]),c(3,i)] # exclude groups with NA values (column 3 is group, column i is variable)
      if (nrow(df1) > 0 & length(unique(df1$group)) >= 2) { # if there is data of at least 2 groups

        if (class(df1[,2]) == "factor") {
          suppressWarnings(w <- kruskal.test(as.numeric(levels(df1[,2]))[df1[,2]] ~ df1$group)) # Mann Whitney U test (column 2 is variable excl NA values)
        } else {
          suppressWarnings(w <- kruskal.test(as.numeric(df1[,2]) ~ df1$group)) # kruskal test (column 2 is variable excl NA values)
        }

        df2 <- rbind(df2,w$p.value)
        df2 <- data.frame(df2)
        df2$sign <- ifelse(df2$df2<0.001,"***",
                           ifelse(df2$df2<0.01,"**",
                                  ifelse(df2$df2<0.05,"*","ns")))
        varnam <- c(varnam,as.character(colnames(yy)[i])) # variable name
        repls <- rbind(repls, rbind(tapply(df1[,2],df1[,1],length))) # number of replicates per group
        test <- c(test, as.character(w$method))
        stat <- c(stat, as.numeric(w$statistic))
      }
    }

    yo <- y[,c(1:3,grep("occ",colnames(y)))]
    for (i in 4:ncol(yo)) {
      if (sum(yo[,i]) > 0) { # at least one group should have value > 1
        suppressWarnings(chi <- chisq.test(rbind(tapply(yo[,i], y$group, sum),
                                                 tapply(yo[,i], y$group, length)-tapply(yo[,i], y$group, sum))))
        df2 <- rbind(df2,chi$p.value)
        df2 <- data.frame(df2)
        df2$sign <- ifelse(df2$df2<0.001,"***",
                           ifelse(df2$df2<0.01,"**",
                                  ifelse(df2$df2<0.05,"*","ns")))
        varnam <- c(varnam,as.character(colnames(yo)[i])) # variable name
        repls <- rbind(repls, rbind(tapply(yo[,i],yo[,3],length))) # number of replicates per group
        test <- c(test, as.character(chi$method))
        stat <- c(stat, as.numeric(chi$statistic))

      }
    }

    # Add and reorder columns
    df2$variable <- varnam
    df2 <- cbind(df2, test, stat, repls)
    df2 <- df2[,c(3,4,5,1,2,6:(5+length(levels(y$group))))]

    # names of columns
    v1 <- numeric(0)
    for (i in 1:length(levels(y$group))) {
      naml <- numeric(0)
      naml <- paste("n_", as.character(levels(y$group)[i]), sep="")
      v1 <- c(v1, naml)
    }
    colnames(df2) <- c("variable","test", "test_statistic", "Pvalue","sign", strsplit(v1, " "))

    # Round P values
    df2$Pvalue <- round(df2$Pvalue,4)

    # Merge with mean and standard errors, use temporary column 'sort' to order the rows
    tabley$sort <- rownames(tabley)
    tabley <- merge(tabley, df2, by = "variable", all.x = TRUE)
    tabley <- tabley[order(tabley$sort),]
    tabley <- tabley[,-grep("sort", colnames(tabley))]

    # Save outcomes
    write.table(tabley, paste(epgname, "_epgtable_", as.character(epgstart), "-",
                              as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)
  }
}

  # Epgplots---------------------------------------------------------------------------

  # If argument epgplot is TRUE

  if (epgplot == "TRUE") {

    # factor to numeric
    for (i in 4:ncol(y)) {
      y[,i] <- as.numeric(as.character(y[,i]))
    }

    # Remove columns with only missing values (otherwise error message boxplot())
    y <- Filter(function(x)!all(is.na(x)), y)

    # Colors
    #epgcolors <- rep(c("blue", "brown1", "honeydew4", "gold"),5)
    epgcolors <- rep(c(rgb(0,0,1,0.8), rgb(1,0,0,0.8),rgb(0,0,1,0.4), rgb(1,0,0,0.4), rgb(0.3,0.8,0,0.8), rgb(0.8,0.8,0,0.8)),3)

    # Make plot
    pdf(paste(epgname, "_boxplots_", as.character(epgstart), "-", as.character(epgstop),
              "h", ".pdf",sep = ""),width=10+(2*(length(loc)-2)),height=30)
    par(mfrow=c(15,4), mar=c(2,2,2,1))

    for (i in 4:ncol(y)) {
      suppressWarnings(boxplot(y[,i]~y$group, main = colnames(y)[i], col=epgcolors[1:length(loc)]))
    }

    dev.off()

  }

  # Epgbin-------------------------------------------------------------------------------
  if (epgbin == "TRUE") {

    # Make time bins---------------------------------------------------

    # convert hour to min
    epgstartmin <- epgstart*60
    epgstopmin <- epgstop*60

    ## make 2 new columns to define the hour (levels 1 to 8) in which the event starts and stops
    x$start.hour<-0
    x$stop.hour<-0

    #Hours 1 to 8 in min:0-60, 60-120, 120-180, 180-240, 240-300, 300-360, 360-420, 420-480

    ## for each event (=row) define start and stop hour
    for (i in 1:nrow(x)){
      if (x$start[i]>=epgstartmin & x$start[i]<epgstartmin+60) {x$start.hour[i] <- round((epgstartmin+60)/60,1)}
      if (x$start[i]>=epgstartmin+60 & x$start[i]<epgstartmin+120) {x$start.hour[i] <- round((epgstartmin+120)/60,1)}
      if (x$start[i]>=epgstartmin+120 & x$start[i]<epgstartmin+180) {x$start.hour[i] <- round((epgstartmin+180)/60,1)}
      if (x$start[i]>=epgstartmin+180 & x$start[i]<epgstartmin+240) {x$start.hour[i] <- round((epgstartmin+240)/60,1)}
      if (x$start[i]>=epgstartmin+240 & x$start[i]<epgstartmin+300) {x$start.hour[i] <- round((epgstartmin+300)/60,1)}
      if (x$start[i]>=epgstartmin+300 & x$start[i]<epgstartmin+360) {x$start.hour[i] <- round((epgstartmin+360)/60,1)}
      if (x$start[i]>=epgstartmin+360 & x$start[i]<epgstartmin+420) {x$start.hour[i] <- round((epgstartmin+420)/60,1)}
      if (x$start[i]>=epgstartmin+420 & x$start[i]<epgstartmin+480) {x$start.hour[i] <- round((epgstartmin+480)/60,1)}
      if (x$start[i]>=epgstartmin+480 & x$start[i]<epgstartmin+540) {x$start.hour[i] <- round((epgstartmin+540)/60,1)}
      if (x$start[i]>=epgstartmin+540 & x$start[i]<epgstartmin+600) {x$start.hour[i] <- round((epgstartmin+600)/60,1)}
      if (x$start[i]>=epgstartmin+600 & x$start[i]<epgstartmin+660) {x$start.hour[i] <- round((epgstartmin+660)/60,1)}
      if (x$start[i]>=epgstartmin+660 & x$start[i]<epgstartmin+720) {x$start.hour[i] <- round((epgstartmin+720)/60,1)}
      if (x$start[i]>=epgstartmin+720 & x$start[i]<epgstartmin+780) {x$start.hour[i] <- round((epgstartmin+780)/60,1)}
      if (x$start[i]>=epgstartmin+780 & x$start[i]<epgstartmin+840) {x$start.hour[i] <- round((epgstartmin+840)/60,1)}
      if (x$start[i]>=epgstartmin+840 & x$start[i]<epgstartmin+900) {x$start.hour[i] <- round((epgstartmin+900)/60,1)}
      if (x$start[i]>=epgstartmin+900 & x$start[i]<epgstartmin+960) {x$start.hour[i] <- round((epgstartmin+960)/60,1)}


      if (x$stop[i]>epgstartmin & x$stop[i]<=epgstartmin+60) {x$stop.hour[i] <- round((epgstartmin+60)/60,1)}
      if (x$stop[i]>epgstartmin+60 & x$stop[i]<=epgstartmin+120) {x$stop.hour[i] <- round((epgstartmin+120)/60,1)}
      if (x$stop[i]>epgstartmin+120 & x$stop[i]<=epgstartmin+180) {x$stop.hour[i] <- round((epgstartmin+180)/60,1)}
      if (x$stop[i]>epgstartmin+180 & x$stop[i]<=epgstartmin+240) {x$stop.hour[i] <- round((epgstartmin+240)/60,1)}
      if (x$stop[i]>epgstartmin+240 & x$stop[i]<=epgstartmin+300) {x$stop.hour[i] <- round((epgstartmin+300)/60,1)}
      if (x$stop[i]>epgstartmin+300 & x$stop[i]<=epgstartmin+360) {x$stop.hour[i] <- round((epgstartmin+360)/60,1)}
      if (x$stop[i]>epgstartmin+360 & x$stop[i]<=epgstartmin+420) {x$stop.hour[i] <- round((epgstartmin+420)/60,1)}
      if (x$stop[i]>epgstartmin+420 & x$stop[i]<=epgstartmin+480) {x$stop.hour[i] <- round((epgstartmin+480)/60,1)}
      if (x$stop[i]>epgstartmin+480 & x$stop[i]<=epgstartmin+540) {x$stop.hour[i] <- round((epgstartmin+540)/60,1)}
      if (x$stop[i]>epgstartmin+540 & x$stop[i]<=epgstartmin+600) {x$stop.hour[i] <- round((epgstartmin+600)/60,1)}
      if (x$stop[i]>epgstartmin+600 & x$stop[i]<=epgstartmin+660) {x$stop.hour[i] <- round((epgstartmin+660)/60,1)}
      if (x$stop[i]>epgstartmin+660 & x$stop[i]<=epgstartmin+720) {x$stop.hour[i] <- round((epgstartmin+720)/60,1)}
      if (x$stop[i]>epgstartmin+720 & x$stop[i]<=epgstartmin+780) {x$stop.hour[i] <- round((epgstartmin+780)/60,1)}
      if (x$stop[i]>epgstartmin+780 & x$stop[i]<=epgstartmin+840) {x$stop.hour[i] <- round((epgstartmin+840)/60,1)}
      if (x$stop[i]>epgstartmin+840 & x$stop[i]<=epgstartmin+900) {x$stop.hour[i] <- round((epgstartmin+900)/60,1)}
      if (x$stop[i]>epgstartmin+900 & x$stop[i]<=epgstartmin+960) {x$stop.hour[i] <- round((epgstartmin+960)/60,1)}
    }

    ## For events > 1 hour, make new rows for each hour after the 1st hour
    first.nr <- nrow(x) #To compare after the loop if any row have been added
    nr <- first.nr #baseline for row number of newly added rows

    for (i in 1:first.nr) {
      if (x$start.hour[i]!=x$stop.hour[i]) { # if a probe spans more than 1 timebin
        A <- x$stop.hour[i] - x$start.hour[i] # the number of timebins the probe is occurring
        co <- -1 #counts the new start and stop hour of the extra timebins
        while (A > 0) { # for every hour the probe spans:
          x <- rbind(x,x[i,]) # add a new timebin (row)
          nr<- nr + 1 # calculate row number of the new timebin
          co <- co + 1 # calculate the extra hour(s) which will be added to the probe starthour
          x$start[nr] <- (x$start.hour[i] + co)*60 # new timebin: adjust start of probe to start of new timebin
          if (x$stop[nr] <= (x$start.hour[i] + co + 1)*60) { # new timebin: if the probe stops before the end of the bin
            x$stop[nr] <- x$stop[nr] # then don't adjust original stop time
          } else { # if probe continues after the end of the new timebin
            x$stop[nr] <- (x$start.hour[i] + co + 1)*60 # stop it at the end of this new timebin
          }
          x$start.hour[nr] <- x$stop[nr]/60 # new timebin: adjust starthour
          x$stop.hour[nr] <- x$start.hour[nr] # new timebin: adjust stop hour so that startbin=stopbin
          A <- A-1 # one round (one extra timebin) is finished. While loop repeats until every extra hour has its own timebin
        }
        x$stop[i] <- x$start.hour[i]*60 # the original starting timebin: adjust stop time of probe to the end of the bin
      }
    }


    ## Remove probes in last hour if last hour is incomplete
    x$lab <- ifelse(x$start.hour > epgstop,1,0)
    x <- subset(x,x$lab != 1)
    x <- x[, 1:(ncol(x)-1)]

    ## Calculate the percentage of time spent on the event in the specific hour
    x$dur <- x$stop - x$start
    x$dur_perc <- (x$dur/60)*100

    # Order rows according to recording id and start probe
    x$id_rec <- as.numeric(x$id_rec)
    x<-x[order(x$id_rec, x$start),]

    # Calculate percentage of time, number of events and mean duration of events per bin
    outp <- numeric(0)
    id_i <- unique(x$id_rec)
    id_hour <- min(x$start.hour):max(x$start.hour)
    id_wave <- unique(x$wave)

    for (i in id_i) {
      for (j in id_hour) {
        for (k in id_wave) {
          df1 <- numeric(0)
          df1 <- subset(x, x$id_rec == i)
          df2 <- numeric(0)
          df2 <- subset(df1, df1$start.hour == j)
          df3 <- numeric(0)
          df3 <- subset(df2, df2$wave == k)
          v1 <- numeric(0)
          v1 <- c(i, j, k, sum(df3$dur_perc), nrow(df3), (sum(df3$dur_perc)/nrow(df3)))
          outp <- rbind(outp, v1)
        }
      }
    }
    colnames(outp) <- c("id_rec", "hour", "wave", "perc", "num", "mean")
    namebins <- x[, 1:3]
    bins <- merge(unique(namebins), outp, by="id_rec", all.x = F, all.y = T)

    # name waveform
    bins$wavename <- ifelse(bins$wave == 1, "non probing",
                            ifelse(bins$wave == 2, "pathway",
                                   ifelse(bins$wave == 4, "salivation phloem",
                                          ifelse(bins$wave == 5, "phloem feeding",
                                                 ifelse(bins$wave == 6, "penetration difficulties",
                                                        ifelse(bins$wave == 7, "xylem ingestion",
                                                               ifelse(bins$wave == 8, "potential drop",
                                                                      bins$wave)))))))
    bins <- bins[,c(1:5,9,6:8)]

    write.table(bins, paste(epgname, "_timebins_", as.character(epgstart), "-",
                            as.character(epgstop), "h", ".csv",sep = ""), sep = ",", row.names = FALSE)

    # Timebinplots -------------------------------------------------------

    if (epgplot == "TRUE") {

      # Plot with {sciplot}
      pdf(paste(epgname, "_timebinplots_", as.character(epgstart), "-", as.character(epgstop),
                "h", ".pdf",sep = ""),width = 8, height = 3*(length(id_wave)))
      par(mfrow=c(length(id_wave),2))
      for (i in id_wave) {
        df1 <- numeric(0)
        df1 <- subset(bins, bins$wave == i)
        suppressWarnings(sciplot::lineplot.CI(df1$hour,df1$perc,group=df1$group,type="b", ylab="% Time",
                                     xlab="Time (h)",err.width=0.002,col= epgcolors,main=as.character(df1[1,6])))
        suppressWarnings(sciplot::lineplot.CI(df1$hour,df1$num,group=df1$group,type="b", ylab="# Events",
                                     xlab="Time (h)",err.width=0.002,col= epgcolors,main=as.character(df1[1,6])))
      }
      dev.off()

    }
  }
  #library(sciplot)



  # end function--------------------------------------------------
  }

