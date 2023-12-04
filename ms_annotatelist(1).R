library(xcms)
library(RAMClustR) 
library(CAMERA)
library(readxl)
library(tidyverse)
library(pheatmap)
library(Spectra)
library(MsBackendMsp)
library(MetaboAnnotation)



#First of all TRANSFORM the data with msConvert peakpicking 1-2 
#After do a excel with sample names (including .mzML example prem_va6098.mzML) and sample groups (like C , WS or whatever)

#Put the directory with the data
direccion <- "/home/acien/MS/MapiPuri/raw"
setwd(direccion)

#In this part it will pop up a csv in the directory we are, and we should 
#fill what we know and save. After we press enter in rstudio to continue. 
#Here It will open a csv archive but open it with txt editor!!. Full everything you want mandatory:
#name, LC-MS, ce2, mode and in mslevel its 2.
experiment <- defineExperiment(csv = TRUE) 


#Read data 
info <- read_excel("Infopos.xlsx")
pd2 <- paste0(file.path(direccion, info$sample_names))
raw_data <- readMSData(files = pd2, pdata = new("NAnnotatedDataFrame", info) ,mode = "onDisk")

#Starting the preprocessing, We implement IPO package for ajusting the parameters.In this script adjust for RP in Synapt
xcmsparam <- CentWaveParam(snthresh = 3,noise = 50,ppm = 25, mzdiff = 0.05, peakwidth = c(5,20),fitgauss = TRUE, verboseColumns = TRUE)
xdata <- findChromPeaks(raw_data, param = xcmsparam, msLevel = 1L )
xcmsparam <- CentWaveParam(snthresh = 3, ppm = 25, noise = 500 ,mzdiff = 0.05, peakwidth = c(5,20),fitgauss = TRUE, verboseColumns = TRUE)
xdata <- findChromPeaks(xdata, param = xcmsparam, msLevel = 2L, add=TRUE )
#Refine peaks could be implemented, just to remove noise, but I have a error so I have to fixed 
#mpp <- MergeNeighboringPeaksParam(expandRt = 2)
#xdata <- refineChromPeaks(xdata, mpp)
pdp <- PeakDensityParam(sampleGroups = xdata$sample_groups ,minFraction = max(0.8, 1), bw =2, binSize = 0.05, maxFeatures = 50)
xdata <- groupChromPeaks(xdata, param = pdp)
pgp <- PeakGroupsParam( minFraction = 0.5, family="symmetric", span = 1, extraPeaks = 0)
xdata<- adjustRtime(xdata, param = pgp)
xdata <- groupChromPeaks(xdata, param = pdp)
fpp <- FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 0)
xdata <- fillChromPeaks(xdata, param = fpp) 

#This is the part where ms1 and ms2 are connected.
rc<- ramclustR(xcmsObj = xdata, MStag = info$sample_names, idMSMStag = info$sample_names ,taglocation = "pheno" , ExpDes = experiment, usePheno = TRUE)

#This function will create a folder with 2 folder each have the posible compounds with the ms1, rt and ms2 associated, in two format.
RC <- do.findmain(ramclustObj = rc, mode = "positive") #Add positive or negative

#This part is neccesary for conecting the object for the spectra package which is the one matching to a database
write.Realmsp <- function(
    ramclustObj = NULL,
    one.file = FALSE
) {
  
  if(!is(ramclustObj, "hclus") & 
     ramclustObj$dist.method != "RAMClustR") {
    stop("this is not a RAMClustR object")
  }
  
  if(is.null(ramclustObj$precursor.mz)) {
    ms2 <- FALSE
  } else {
    ms2 <- TRUE
  }
  
  if(!dir.exists('spectra')) {
    dir.create('spectra')
  }
  
  if(!one.file) {
    if(!dir.exists('spectra/msp')) {
      dir.create('spectra/msp')
    }
  }
  ion.mode <- as.character(ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]]) == "msmode"),1])
  if(toupper(substring(ion.mode, 1, 1)) == "P") {
    ion.mode = "Positive"
  } else {
    ion.mode = "Negative"
  }
  
  out.list <- as.list(rep(NA, length(ramclustObj$cmpd)))
  
  for(i in 1:length(ramclustObj$cmpd)) {
    ions <- which(ramclustObj$featclus == i)
    
    if(ms2) {
      spectrum <- data.frame(
        'mz' = ramclustObj$fmz[ions],
        'int' = ramclustObj$msmsint[ions]
      )
    } else {
      spectrum <- data.frame(
        'mz' = ramclustObj$fmz[ions],
        'int' = ramclustObj$msint[ions]
      )
    }
    
    spectrum <- spectrum[order(spectrum[,"int"], decreasing = TRUE), ]
    
    out <- paste0(
      "NAME: ", ramclustObj$cmpd[i], '\n', 
      "IONMODE: ", ion.mode, '\n',
      "SPECTRUMTYPE:Centroid", '\n',
      "RETENTIONTIME: ", round(ramclustObj$clrt[i], 2), '\n'
    )
    
    if(any(names(ramclustObj)=="clri")) {paste0(
      out,
      "RETENTIONINDEX:", round(ramclustObj$clri[i], 2),  '\n'
    )
    }
    
    if(ms2) {
      out <- paste0(out, 
                    "PRECURSORMZ: ", ramclustObj[["M.ann"]][[i]][["mz"]][[which.max(ramclustObj[["M.ann"]][[i]][["int"]])]] ,'\n', 
                    "PRECURSORTYPE: ", ramclustObj$precursor.type[i] , '\n'
      )
    }
    out <- paste0(out,
                  "Num Peaks: ", nrow(spectrum), '\n'
                  # m/z intensity pair (tab, comma, space can be used as the delimiter.)
    )
    for(j in 1:nrow(spectrum)) {
      out <- paste0(out, 
                    spectrum[j,"mz"], 
                    " ", 
                    round(spectrum[j,"int"]),
                    '\n'
      )
    }
    out.list[[i]] <- out
  }
  
  if(one.file) {
    out <- vector(mode = "character")
    for(i in 1:length(out.list)) {
      out <- paste0(out, out.list[[i]], '\n')
    }
    exp.name <- ramclustObj$ExpDes[[1]][which(row.names(ramclustObj$ExpDes[[1]]) == "Experiment"),1]
    if(nchar(exp.name) == 0) {
      exp.name <- "spectra"
    }
    sink(paste0("spectra/", exp.name, ".msp"))
    cat(out)
    sink()
  } else {
    for(i in 1:length(out.list)) {
      sink(paste0("spectra/msp/", ramclustObj$cmpd[[i]], ".msp"))
      cat(out.list[[i]], '\n')
      sink()
    }
  }
  
}

write.Realmsp(ramclustObj = RC, one.file = TRUE)

#As this R object have all the information is important to keep if we want to reuse
saveRDS(RC, file="Mapipos.RDS")

#Matching part
#Adding the database whatever we are interested, we could create a database with our identification. With ms1 and rt is posible to do the matching too.
compuestos2 <- Spectra("Spectra/Mapipos.msp", source = MsBackendMsp())
basedatos <- Spectra("BioMSMS-Pos-PlaSMA.msp", source = MsBackendMsp())
#basedatos <- Spectra("MoNA-export-VF-NPL_QTOF.msp", source = MsBackendMsp(),mapping = spectraVariableMapping(MsBackendMsp(), "mona"))

basedatos2 <- Spectra("MSMS-Pos-Respect.msp", source = MsBackendMsp())
basedatos3 <- Spectra("MSMS_Public_EXP_Pos_VS17.msp", source = MsBackendMsp())

#Here Its the function for comparing to the database, we can change the ppm for the precursor (I think bigger than 5ppm maybe its to much)
csp <- CompareSpectraParam(requirePrecursor = TRUE, ppm = 5)
mtches <- matchSpectra(compuestos2, target = c(basedatos, basedatos2, basedatos3), param = csp)
#Two different matching formulas I like more second one (in the information say that Its similar to msdial one)
mp <- MatchForwardReverseParam(requirePrecursor = TRUE, ppm = 5)
mtches2 <- matchSpectra(compuestos2, target = c(basedatos,basedatos3, basedatos2), param = mp)

#Here we rescue the data which have a positive and take some info, if the database have other info like category like rna-seq for pathways
mtches_csp <- spectraData(mtches, columns = c("rtime","precursorMz", "name","adduct", "target_rtime", "target_precursorMz", "target_name", "score"))
mtches_mp <- spectraData(mtches2, columns = c("rtime","precursorMz", "name","adduct", "target_rtime", "target_precursorMz", "target_name", "score", "presence_ratio","matched_peaks_count" ))

#Just for export to excel
guardar <- as.data.frame(mtches_mp@listData)
guardar <- guardar %>% subset(guardar$target_name != "NA")
write.csv(guardar, "identificados_mp.csv")
guardar <- as.data.frame(mtches_csp@listData)
guardar <- guardar %>% subset(guardar$target_name != "NA")
write.csv(guardar, "identificados_csp.csv")

#Here we can export the intensity values, I have to check if I can do the differentreport in this dataframe, that would be nice
specab<- rc$SpecAbund
write.csv(t(specab), "intesidad.csv")



#This part its the camera part which is important for confirmation of the 
xfill2 <- filterMsLevel(xdata, msLevel = 1)
#We need to transform the object to the old one.
xfill2 <- as(xfill2, "xcmsSet")
sampclass(xfill2) <- pData(xdata)$sample_groups
sampnames(xfill2) <- pData(xdata)$sample_names
xfill2 <- fillPeaks(xfill2, BPPARAM = SerialParam())
  an<-xsAnnotate(xfill2)
#group according to retention time
an_groupFWHM<-groupFWHM(an)

#check grouping with EIC correlation, when indicated regroup
an_groupCorr<-groupCorr(an_groupFWHM)

#search for isotopes
an_findIsotopes<-findIsotopes(an_groupCorr)

#calculate possible adducts
an_findAdducts<-findAdducts(an_findIsotopes,polarity="positive")
#an_findAdducts<-findAdducts(an_findIsotopes,polarity="negative")

#get annotated peaklist
xset_an<-getPeaklist(an_findAdducts)

file='pruebacameraMapipos.csv'   
write.csv(xset_an,file=file)
dr <- diffreport(xfill2,filebase='anovacameraMapipos.csv')

#Now create a document with the anova info or other info of diffreport to camera informatition. For this we order both document for mzmed and add anova to the variable with aduct info.
dr <- dr[order(dr$mzmed),]
xset_an$anova <- dr$anova
#we can filter by anova value, in order to reduce the number of rows
xset_an2 <- xset_an %>% filter(xset_an$anova < 0.01)
write.csv(xset_an2, "Mapiposfileter.csv")






