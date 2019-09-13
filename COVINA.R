library(xcms)
#library(faahKO)
#library(RColorBrewer)
#library(pander)
#library(magrittr)
library(Kendall)

##Specify the path for the COVINA test dataset
mzXMLPath <- "/Adduct/COVINA Test Dataset"

## Get the full path to the mzXML files
mzXMLs <- dir(mzXMLPath, full.names = TRUE,
            recursive = TRUE)

##The raw data files are read in the order of asending CID energy levels
raw_data <- readMSData(files = mzXMLs[c(1,5:8,2:4)], 
                       mode = "onDisk") 

mzs <- mz(raw_data)
intensities <- intensity(raw_data)

mzs_by_file <- split(mzs, f = fromFile(raw_data))
intensities_by_file <- split(intensities, f = fromFile(raw_data))

Chrom <- function(sample=1,scans,mzlist,pool.method=sum) {
  
  Intensities.Selected.Scans <- intensities_by_file[[sample]][scans]
  mzs.Selected.Scans <- mzs_by_file[[sample]][scans]
  
  #int1 <- intensities_by_file[[1]]
  #mzs1 <- mzs_by_file[[1]]
  n.scan <- length(scans)
  n.mzs <- sapply(Intensities.Selected.Scans, length)
  scan.num <- rep(scans, n.mzs)
  int <- unlist(Intensities.Selected.Scans)
  mzs <- unlist(mzs.Selected.Scans)
  
  if(class(mzlist)=="matrix") {
    start <- mzlist[,1]
    end <- mzlist[,2]
  } else {
    start <- mzlist[1]
    end <- mzlist[2]
  }

  
  cond1 <- outer(mzs, start, ">=")
  cond2 <- outer(mzs, end, "<=")
  filter <- cond1 & cond2
  
  pool.int <- matrix(0, nrow = n.scan, ncol = ncol(filter))
  
  for(i in 1:ncol(filter)) {
    f <- filter[, i]
    int.f <- int[f]
    scan.f <- scan.num[f]
    int.by.scan <- split(int.f, scan.f)
    temp <- sapply(int.by.scan, pool.method)
    if(length(temp)!=0) {pool.int[(as.integer(names(temp))-scans[1]+1), i] <- temp}
  }
  return(pool.int)
}

FlipSum <- function (Number.Sequence) {
  N <- length(Number.Sequence)
  FS <- rep(0,N)
  for (i in 1:N) {
    FS[i] <- sum(c(rep(-1,(i-1)),rep(1,(N+1-i)))*Number.Sequence)
  }
  return(FS)
}

InCIDR <- function(Scan.Number=0, Base.mz,  Chromatogram.HalfWidth=20, Mass.Tolerance=20,
                   Number.Of.Samples=1,Intensity.Threshold=0, Correlation.Threshold=0.9,samplenames=NULL) {
  
  #Number.Of.Samples <- 4
  Base.Sample <- 1
  #Scan.Number <- 199
  #Intensity.Threshold <- 50000
  Mass.Tolerance <- Mass.Tolerance/10^6
  #Base.mz <- 130.087265
  #Chromatogram.HalfWidth <- 12
  #Correlation.Threshold <- 0.9

  Full.Scan.Range <- c(1:length(mzs_by_file[[Base.Sample]]))
  Full.Chromatograms <- Chrom(scans=Full.Scan.Range,mzlist=Base.mz*(1+c(-Mass.Tolerance,Mass.Tolerance)))
  
  if(Scan.Number==0) {
    Scan.Number<-which.max(Full.Chromatograms)
  }
  
  print(paste("Peak is found at Scan",Scan.Number))

  SS <- Scan.Number+c((-Chromatogram.HalfWidth*3):(Chromatogram.HalfWidth*3))
  PS <- sum(Full.Chromatograms[SS]>(Full.Chromatograms[Scan.Number]/10))
  print(paste("Peak width is",PS,"scans at 10% height."))
  
  Scan.Range <- Scan.Number+c(-Chromatogram.HalfWidth:Chromatogram.HalfWidth)
  
  Query.Spectrum <- cbind(mzs_by_file[[Base.Sample]][[Scan.Number]],intensities_by_file[[Base.Sample]][[Scan.Number]])
  Query.Spectrum.Trimmed <- Query.Spectrum[intensities_by_file[[Base.Sample]][[Scan.Number]]>Intensity.Threshold,]
  
  mz.Table <- matrix(c(Query.Spectrum.Trimmed[,1]*(1-Mass.Tolerance),Query.Spectrum.Trimmed[,1]*(1+Mass.Tolerance)),ncol=2)
  
  Query.Chromatograms <- Chrom(scans=Scan.Range,mzlist=Base.mz*(1+c(-Mass.Tolerance,Mass.Tolerance)))
  Results.Chromatograms <- matrix(0,ncol=Number.Of.Samples*nrow(mz.Table),nrow=length(Scan.Range))
  
  for (i in 1:Number.Of.Samples) {
    Results.Chromatograms[,(nrow(mz.Table)*(i-1)+1):(nrow(mz.Table)*i)] <- Chrom(sample=i,scans=Scan.Range,mzlist=mz.Table)
  }

  Spectra.Cor <- cor(Query.Chromatograms,Results.Chromatograms[,1:nrow(Query.Spectrum.Trimmed)])
  if(!sum(Spectra.Cor>Correlation.Threshold)>1) 
    return(print("No Co-Variant ion found!")) 
  else (print(paste(sum(Spectra.Cor>Correlation.Threshold)-1,"Co-variant ions found!")))
  Spectrum.Covariants <- cbind(Query.Spectrum.Trimmed[Spectra.Cor>Correlation.Threshold,1],
                               Spectra.Cor[Spectra.Cor>Correlation.Threshold])
  
  InCIDR.Results <- matrix(0,nrow=nrow(Spectrum.Covariants),ncol=(Number.Of.Samples+4))
  InCIDR.Results[,1:2] <- Spectrum.Covariants

  
  for (i in 1:Number.Of.Samples) {
    CoEluents.Intensities <- Results.Chromatograms
    InCIDR.Results[,(2+i)] <- colSums(CoEluents.Intensities[,(match(Spectrum.Covariants[,1],
                                                                    Query.Spectrum.Trimmed[,1])+(i-1)*nrow(mz.Table))])
  }
  
  if(Number.Of.Samples>1) {
    InCIDR.Results[,(ncol(InCIDR.Results)-1)] <- apply(InCIDR.Results[,3:(ncol(InCIDR.Results)-2)],MARGIN = 1,
                                                   FUN=function(x) Kendall(c(Number.Of.Samples:1),x)$tau)
    InCIDR.Results[,(ncol(InCIDR.Results))] <- FlipSum(InCIDR.Results[,(ncol(InCIDR.Results)-1)])
  }
  
  if(is.null(samplenames)) {
    samplenames <- paste("Sample",c(1:Number.Of.Samples),sep="")
  }
  colnames(InCIDR.Results) <- c("mz","Correlation",samplenames,"Tau","RankScore")
  
  InCIDR.Results
  
  return(InCIDR.Results)
}

Sample.Names <- c("CID_0eV","CID_2eV",
                              "CID_4eV","CID_6eV","CID_8eV","CID_10eV",
                              "CID_15eV","CID_20eV")

Results.Lactate <- InCIDR(Base.mz = 89.024361,Chromatogram.HalfWidth = 40,
                                      Number.Of.Samples = 8,samplenames=Sample.Names) # Lactate

Results.Pyruvate <- InCIDR(Base.mz = 87.008743,Chromatogram.HalfWidth = 40,
                                      Number.Of.Samples = 8) # Pyruvate

Results.Leucine <- InCIDR(Scan.Number = 429,Base.mz = 130.087265,
                                      Chromatogram.HalfWidth = 40,
                                      Number.Of.Samples = 8) # Leucine

Results.Isoleucine <- InCIDR(Scan.Number = 471,Base.mz = 130.087341,
                                      Chromatogram.HalfWidth = 40,
                                      Number.Of.Samples = 8) # Isoleucine

Results.Malate <- InCIDR(Base.mz = 133.014221,Chromatogram.HalfWidth = 40,
                                      Number.Of.Samples = 8) # Malate

Results.Glc6P <- InCIDR(Scan.Number = 1180,Base.mz = 259.02240,
                                     Chromatogram.HalfWidth = 40,
                                     Number.Of.Samples = 8) # Glc6P

Results.Fru6P <- InCIDR(Scan.Number = 1109,Base.mz = 259.02231,
                                  Chromatogram.HalfWidth = 40,
                                  Number.Of.Samples = 8) # Fru6P

Results.FBP <- InCIDR(Base.mz = 338.988770,
                                  Chromatogram.HalfWidth = 80,
                                  Intensity.Threshold = 10000,
                                  Number.Of.Samples = 8) # FBP

Results.NAD <- InCIDR(Base.mz = 662.101746,
                                  Chromatogram.HalfWidth = 40,
                                  Number.Of.Samples = 8) # NAD

Results.ATP <- InCIDR(Base.mz = 505.988586,
                                  Chromatogram.HalfWidth = 40,
                                  Number.Of.Samples = 8) # ATP



