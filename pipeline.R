# xcms-based peak-picking of LC-MS data based on mzXML files (V.2)
#################################################
# Torben Kimhofer, 01 Aug. 2019
# Murdoch University
# email: torben.kimhofer@murdoch.edu.au
#
# prerequisites:
# each individual assay and run requires peak picking parametrisation,
# this should have been performed prior to running this script
# a user-freindly peak picking-optimisation tool can be found under http://www.github.com/kimsche/msexplorer
# get in touch if you have any questions (see email address above)
#################################################

# load packages
library("BiocParallel")
library("xcms")
library("gridExtra")

### set up parallel processing to speed up computation time
# for windows
# bpinst=register(bpstart(SnowParam(7)))
# for mac
bpinst=MulticoreParam(progressbar = T, timeout=60L*60L*15L)

### sepecify filte path of mzXML files
path <-'/Volumes/ANPC_ext1/HeatSealFoilComparison_mzXML/'
files <- list.files(path, full.names = T, pattern="\\.mzXML")

### create dataframe of metadata information
# include analytical run order, group memberships (eg., case, control, QC)
# MODIFY THE FOLLOWING


df <- data.frame(ID=basename(files), 
                 fullN=files,
                 RunOrder=1:length(files),
                 #Group=rep(paste('Group', LETTERS[1:7]), each=4),
                 stringsAsFactors = F)

df$Group='Sample'
df$Group[grep('QC', files)]='QC'
df$Group[grep('dil', files, ignore.case = T)]='Dilution'


meta=readWaters_acquisitonPars(WatersDataDir = '/Volumes/Torben_1/DiaObesity/Diaobesity_Lipid_NEG.PRO/', plot=T)

idx=match(gsub('01$', '', gsub('\\.CDF', '', basename(df$ID))), meta$`Acquired Name`)
meta=meta[idx,]
df$RunOrder=rank(meta$date)

table(meta$`Inlet Method`)

idx=!grepl('END', meta$`Inlet Method`)
df=df[idx,]
meta=meta[idx,]





### read-in data and perform peak picking using centwave method
# PARAMETRISATION SHOULD HAVE BEEN DONE BEFORE RUNNING THIS SCRIPT (this is very important!)
raw_xcms <- xcmsSet(files, 
                    phenoData=AnnotatedDataFrame(data=df),
                    method='centWave', 
                    
                    # always optimise the following three parameters for individual analytical runs
                    ppm=5, 
                    peakwidth=c(3, 100),
                    prefilter=c(3,100),
                    
                    # instrument-specific parameters (these do not require optimisation for each run and are fixed for a single instrument)
                    mzdiff=-1,
                    noise=100,
                    snthresh=2,
                    
                    # global parameters (do not require changes at all)
                    profparam = list(profstep=0), 
                    mslevel=1, 
                    includeMSn=F, 
                    mzCenterFun='meanApex3',
                    integrate=2,
                    fitgauss=F,
                    
                    # parallel processing (remove if using only one worker)
                    BPPARAM = bpinst
)


### perform peak alignment/retention time correction (step is optional)
# only required for longer analytical runs (n>40)
# not suitable for runs showing batch effects

# 1. find well-aligned peaks across all samples and use these as reference points to align other peaks
# peak grouping: establishes which peaks are present across samples
# Minimum fraction (minfrac) filter for peak grouping: If a signal is present in a mimimum of x (samples/n samples) in at least one group (see column Group in variabel df, defined above)
# minfrac can rages between 0-1 (higher values)
xcms_anchor <- group.density(raw_xcms, minfrac=0.9)

# 2. perform retention time correction
xcms_rtcor <- retcor(xcms_anchor, method="obiwarp", plottype = "deviation", profStep=0.01, response=80, distFunc="cor_opt", gapInit=NULL, gapExtend=NULL, factorDiag=2, factorGap=1, localAlignment=0, initPenalty=0)

# 3. re-group peaks based on aligned peaks
xcms_group <- group.density(xcms_group, minfrac=0.8)

### count missing signals (extract not peak-filled data matrix, ie., where NA values indicate missing peaks)
ds=groupval(xcms_group, value="into")

# count missing signals
nas <- list(
  samples=apply(ds, 1, is.na), 
  features=apply(ds, 2, is.na))

### peak filling (replace NA values with noise)
xcms_nfil=fillPeaks(xcms_group, method='chrom')

### extract data for further pre-processing and normalisation
ds=t(groupval(xcms_nfil, value="into"))

# get information about MS feature (mz, rt)
f_id=data.frame(xcms_nfil@groups[,1:7], stringsAsFactors = F)
f_id$rt_m=round(out[,4]/60, 3)
f_id$id=paste(out[,4], out[,1], sep='_')


### rsd filter
rsd_fct=function(x){ sd(x, na.rm=T) / mean(x, na.rm=T) }

idx=which(df$Group='QC')
f_id$rsd=apply(ds[idx,], 2, rsd_fct)
# dilution filter

### dilution series filter
idx=which(df$Group='Dilution')
f_id$dilfilter=dil_fct(dsn[idx,], dilFs=c(1,2,4,8,16))

# normalisation using pqn
# library(MetaboMate)
dsn=pqn(ds, plot=T, add.DilF='dilF')
df$dilF=dilF

### PCA
# library(MetaboMate)
mod=pca(dsn, pc = 6)
plotscores(mod, an=list(df$dilF, df$Group))

#  some plotting fcts
g1=plotscores(mod, an=list(Group=paste(df$gr, df$cl1, sep = ': ')))
g2=coeff.influence(mod, labs='q', Curriculum = list('All Features'))
g3=ggplot(f_id, aes(rt_s,mz, colour=Intave))+
  geom_point()+
  scale_x_continuous(sec.axis = sec_axis(~./60, name='RT (min)'), name='RT (s)')+
  labs(size='Average Int', colour='PCA scores cluster')

grid.arrange(g1, g2, g3, ncol=3)

