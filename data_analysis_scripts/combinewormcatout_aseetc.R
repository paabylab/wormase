# Analyze wormcat outputs from many comparisons (each run using wormcat_givebackgroundset.R)
# by Avery Davis Bell, begun 2023.03.13. Updated 2024.02.29
require(data.table)
require(argparser)
require(ggplot2)

#### Functions ####
proctorunf<-function(nftorunf){
  # Breaks wormcat directory names from input file into useable chunks for downstream processing
  # In: Path describing comparisons that were run through wormcat (with nextflow workflow): --torun parameter.
  #         Key column is testname which describes output directory names in **ASSUMED** format
  #         <tested genes>_vs_<background genes><_exclhypdivbadcov if those excluded>_<STRAIN>
  # Out: data table with columns:
  #         testname, as in input
  #         Strain, worm strain
  #         testset, tested genes descriptor
  #         backgroundset, background geneset descriptor,
  #         exclhypdivbadcov, T or F, were hypdiv bad cov genes excluded (as flagged in file name)
  
  datdescrip<-fread(nftorunf, header = T)[, .(testname, exclhypdivbadcov = grepl("exclhypdivbadcov", testname))]
  spl1<-datdescrip[, tstrsplit(testname, split = "_vs_", fixed = F)] # V1 is testset
  spl2<-spl1[, tstrsplit(V2, "_")] # V1 is backgroundset
  
  datdescrip[, `:=`(Strain = sapply(strsplit(spl1$V2, split = "_"), function(x) x[length(x)]),
                    testset = spl1$V1,
                    backgroundset = spl2$V1)]
  setcolorder(datdescrip, c("testname", "Strain", "testset", "backgroundset", "exclhypdivbadcov"))
  return(datdescrip)
}

colcatdotplot<-function(pdat, xcol = "Strain", ycol = "Category", colcol = "foldEnrichmentTestinCateg",
                        xorder = NA, yorder = NA, myxlab = xcol, myylab = ycol, mycollab = NA){
  # Plots 2 factors against each other - dot if there's an observation for this pairing. Dot colored by colcol
  # In: pdat, data.table with all data needed for plot
  #     xcol, character name of column to plot on x axis
  #     ycol, character name of column to plot on y axis
  #     colcol, character name of continuous column to use for color shading
  #     xorder, optional character vector ordering all xcol observations
  #     yorder, optional character vector ordering all ycol observations
  #     myxlab, x axis label
  #     myylab, y axis label
  #     mycollab, label for color gradient legend
  # Out: ggplot2 of this plot
  
  # Set up data
  pdat<-copy(pdat) # going to monkey with it, copying in case it was internal to something
  pdat[, `:=`(myx=get(xcol), myy = get(ycol), mycolnums = get(colcol))]
  if(!is.na(xorder[1])){
    pdat[, myx:=factor(myx, levels = xorder)]
  }
  if(!is.na(yorder[1])){
    pdat[, myy:=factor(myy, levels = yorder)]
  }
  if(is.na(mycollab)){
    mycollab<-colcol
  }
  
  # make plot
  plt<-ggplot(pdat, aes(myx, myy)) + geom_point(aes(color = mycolnums)) +
    scale_color_gradient() + # can customize this more if desired
    xlab(myxlab) + ylab(myylab) + labs(color = mycollab) + theme_bw() + 
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust=1), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17),legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  # Return
  return(plt)
}

#### Arguments ####
p<-arg_parser("Analyze wormcat outputs from many comparisons (each run using wormcat_givebackgroundset.R)", 
              name = "combinewormcatout_aseetc.R", hide.opts = TRUE)

p<-add_argument(p, "--nftorun",
                help = "Path describing comparisons that were run through wormcat (with nextflow workflow): --torun parameter.
                Key column is testname which describes output directory names in **ASSUMED** format
                <tested genes>_vs_<background genes><_exclhypdivbadcov if those excluded>_<STRAIN>",
                type = "character")
p<-add_argument(p, "--wormcatparentdir",
                help = "Directory containing all wormcat output directories to process together. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Output directory path. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                default = "out")
p<-add_argument(p, "--strains",
                help = "Strains to summarize across - path to no header list or comma-separated (i.e., all processed strains excluding any tests done in metacategory.
                So number of strains in which comparison is significant will be of these.)",
                default = "JU1088,EG4348,CB4856,QX1211")

p<-parse_args(p)

# Output directory
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# Input directory
if(p$wormcatparentdir=="getwd()"){
  wormcatparentdir<-getwd()
}else{
  wormcatparentdir<-p$wormcatparentdir
}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

#### Read in & combine data, numerical summary ####
# Get in data
mycats<-paste0("cat", 1:3)

datinfo<-proctorunf(p$nftorun)
##  the following probably should've been a function!!
alldat<-lapply(1:3, function(catnum){ # For each level of wormcat analysis
  onedat<-rbindlist(lapply(1:nrow(datinfo), function(x){ # Read in data for each test of interest
    # significant categories
    dat<-fread(file.path(wormcatparentdir, datinfo[x, testname], paste0("rgs_fisher_cat", catnum, "_apv.csv")),
               header = T)
    # total annotated genes: read in all categories and sum
    forn<-fread(file.path(wormcatparentdir, datinfo[x, testname], paste0("rgs_fisher_cat", catnum, ".csv")),
                header = T)
    return(data.table(datinfo[x, ], dat,
                      nAnnotGenesTestset = forn[,sum(RGS)],
                      nAnnotGenesBackground = forn[,sum(AC)]))
  }))
  setkey(onedat, testname)
  # Compute percentages, ratios
  onedat[, `:=` (pTestInCategofAllTest = RGS/nAnnotGenesTestset,
                 pBackgroundInCategofAllBg = AC/nAnnotGenesBackground,
                 pTestInCategofCateg = RGS/AC,
                 foldEnrichmentTestinCateg = (RGS/nAnnotGenesTestset)/(AC/nAnnotGenesBackground))]
})
names(alldat)<-mycats

# Get numbers per test
nsumm<-lapply(alldat, function(onedat){
  onedat[, .(Strain = Strain[1], testset = testset[1], backgroundset = backgroundset[1], 
             exclhypdivbadcov = exclhypdivbadcov[1],
             nEnrichedCats = sum(!is.na(Category))), by = testname]
})

# Get numbers per test summarized across strains (for quick looks): for each category, how many strains have ANY significant
#     **NOT saying what these are, if they're shared yet
ncatstrain<-lapply(alldat, function(onedat){
  testdat<-onedat[Strain%in%strains, ]
  setkey(testdat, testset, backgroundset, exclhypdivbadcov)
  out<-testdat[, length(unique(Strain[!is.na(Category)])),
          by = c("testset", "backgroundset", "exclhypdivbadcov")]
  setnames(out, "V1", "nStrainsWithEnrichedCats")
  return(out)
})

# Save summaries & underlying data
invisible(lapply(mycats, function(mycat){
  write.table(alldat[[mycat]],
              file.path(p$outdir, paste0(p$baseoutname, "_combinedwormcatapv_", mycat, ".txt")),
              row.names = F, sep = "\t", quote = F)
  write.table(nsumm[[mycat]],
              file.path(p$outdir, paste0(p$baseoutname, "_nwormcatsigcats_", mycat, ".txt")),
              row.names = F, sep = "\t", quote = F)
  write.table(ncatstrain[[mycat]],
              file.path(p$outdir, paste0(p$baseoutname, "_percatnstrainsenriched_", mycat, ".txt")),
              row.names = F, sep = "\t", quote = F)
}))


#### Plots/combinatorial & comparative analyses ####
#--- Dot plots - dot where there's an enriched category
# Make allll the plots
allplts<-(lapply(mycats, function(mycat){
  plts<-lapply(c(T,F), function(hypdiv){ # Separately for including/excluding hypdivbadcov genes
    # Set up data
    dat<-alldat[[mycat]][exclhypdivbadcov==hypdiv, ]
    dat[,pname:=paste(testset, "vs\n", backgroundset)]
    ptitle<-ifelse(hypdiv, paste(mycat, "- hyperdivergent & bad coverage genes EXCLUDED"),
                   paste(mycat, "- hyperdivergent & bad coverage genes INCLUDED"))
    
    # Plot strains across x axis, separate facet for each gene set test
    bystr<-colcatdotplot(pdat = dat[Strain%in%strains, ], xcol = "Strain", ycol = "Category",
                         colcol = "foldEnrichmentTestinCateg", xorder = strains, myxlab = "Strain",
                         myylab = "Enriched category", mycollab = "Fold enrichment\nin test set") +
      facet_wrap(~pname) + ggtitle(ptitle) 
    ## NAs removed
    bystr.nona<-colcatdotplot(pdat = dat[Strain%in%strains & !is.na(Category), ], xcol = "Strain", ycol = "Category",
                         colcol = "foldEnrichmentTestinCateg", xorder = strains, myxlab = "Strain",
                         myylab = "Enriched category", mycollab = "Fold enrichment\nin test set") +
      facet_wrap(~pname) + ggtitle(ptitle, subtitle = "Tests with no enrichments exlcuded") 
    
    # Plot gene set tests across x axis, separate facet for each strain (?)
    bytest<-colcatdotplot(pdat = dat[Strain%in%strains, ], xcol = "pname", ycol = "Category",
                          colcol = "foldEnrichmentTestinCateg", myxlab = "",
                          myylab = "Enriched category",  mycollab = "Fold enrichment\nin test set") +
      facet_wrap(~factor(Strain, levels = strains), ) + ggtitle(ptitle)
    ## NAs removed (not a big deal in this orientation at least for cat1, but for completeness)
    bytest.nona<-colcatdotplot(pdat = dat[Strain%in%strains & !is.na(Category), ], xcol = "pname", ycol = "Category",
                          colcol = "foldEnrichmentTestinCateg", myxlab = "",
                          myylab = "Enriched category",  mycollab = "Fold enrichment\nin test set") +
      facet_wrap(~factor(Strain, levels = strains), ) + ggtitle(ptitle, subtitle = "Tests with no enrichments exlcuded")
    
    # Return plots
    return(list(bystr = bystr, bystr.nona = bystr.nona, bytest = bytest, bytest.nona = bytest.nona))
  })
  names(plts)<-c("exclhypdivbadcov", "inclhypdivbadcov")
  return(plts)
}))
names(allplts)<-mycats

# Save plots - one PDF per plot type, gene set (categories, NA incl vs not in same PDF)
invisible(lapply(c("exclhypdivbadcov", "inclhypdivbadcov"), function(hypdiv){
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_dotplotstrainxaxis_", hypdiv, ".pdf")), 16, 12)
  lapply(allplts, function(x) print(x[[hypdiv]]$bystr))
  lapply(allplts, function(x) print(x[[hypdiv]]$bystr.nona))
  invisible(dev.off())
  
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_dotplottestxaxis_", hypdiv, ".pdf")), 17, 12)
  lapply(allplts, function(x) print(x[[hypdiv]]$bytest))
  lapply(allplts, function(x) print(x[[hypdiv]]$bytest.nona))
  invisible(dev.off())
  
  return(NULL)
}))


# ?? Trend of # genes vs # categories enriched??

#### Script completion message & session information ####
cat("....combinewormcatout_aseetc.R processing complete! Session information:....\n")
sessionInfo()


