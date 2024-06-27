#! /usr/bin/env/ Rscript
# Examine cis-trans estimates for artifactual correlation (see Fraser 2019; Zhang & Emerson 2019), experiment with corrections/ways to think about this
# FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R and initial cis/trans classifications from ase_de_cistransclassifications.R
# by Avery Davis Bell, begun 2024.02.26
require(data.table)
require(DESeq2)
require(argparser)
require(ggplot2)

#### Functions ####
freadcombstrains<-function(exampinput, strains){
  # Reads in exampinput file for each strain in strains; combines the data.tables with added (first) column specifying strain
  # In: exampinput, example filepath to file to read in as data.table, with strain identified in filename. Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz
  #     strains, vector of strain IDs to substitute into exampinput and read in
  # Out: data.table with columns strain, all columns in input. Not keyed.
  
  out<-rbindlist(lapply(strains, function(x){
    cat(paste("Reading in file for strain", x, "; recommend checking that this in fact is for the right strain! Path is: \n",
              file.path(gsub("STRAIN", replacement = x, exampinput, fixed = T), "\n")))
    out1<-fread(file.path(gsub("STRAIN", replacement = x, exampinput, fixed = T)))
    out1[,strain:=x]
  }))
  
  setcolorder(out, "strain")
  return(out)
}

cistranscorr<-function(onedat, ciscol = "log2FoldChange.ASE", parentcol = "log2FoldChange.ParentVsParentN2",
                       transcol = NA, 
                       metadat = data.table()){
  # Computes several correlations, p-values, etc for cis vs trans per-gene estimates (trans = parental - cis)
  # In: onedat, data with columns with names provided in other arguments
  #     ciscol, character name of column where the cis (ASE) effect is in log2FC space
  #     parentcol, character name of column where the among-parent effect matching the ASE effect direction is in log2FC space
  #                   ***OR NA. only provide this OR transcol**
  #     transcol, character name of column for trans effect to directly correlate with ciscol. 
  #             ***OR NA. only provide this OR parentcol**
  #     metadat, one-row data.table of metadata to append to output
  # Out: one-row data table with columns:
  #     <any in metadat, with values from that option >
  #     r, Pearson's R
  #     r.squared, Pearson's R squared
  #      r.pval, p-value for correlation test vs 0 for Pearson's R
  #     rho, Spearman's rho
  #     rho.pval, Spearman's rho correlation test vs 0 p value
  #     tau, Kendall's tau
  #     tau.pval, kendall's tau correlation test vs 0 p value
  
  myx<-onedat[, get(ciscol)]
  if(is.na(transcol)){
    myy<-onedat[, get(parentcol) - get(ciscol)]
  }else if(is.na(parentcol)){
    myy<-onedat[, get(transcol)]
  }
  
  pearstest<-cor.test(myx, myy, method = "pearson")
  speartest<-cor.test(myx, myy, method = "spearman")
  kendtest<-cor.test(myx, myy, method = "kendall")
  
  out<-data.table(metadat, r = pearstest$estimate, r.squared = pearstest$estimate^2, r.pval = pearstest$p.value,
                  rho = speartest$estimate, rho.pval = speartest$p.value,
                  tau = kendtest$estimate, tau.pval = kendtest$p.value)
  
  return(out)
}

myscatter<-function(pdat, xname, yname, xlabel = "", ylabel = "", mytitle = "",
                    mysubt = "", eqaxes = F, giveaxes = c(NA, NA)){
  # Makes pretty scatter plot
  # In: pdat, data
  #     xname, character name of column for X values
  #     yname, character name of column for Y values
  #     xlabel, x axis label
  #     ylabel, y axis label
  #     mytitle, title for plot
  #     mysubt, subtitlefor plot
  #     eqaxes, make x and y axes equal and span data
  #     giveaxes, only works if eqaxes is F: manual axis bounds for both
  # Out: ggplot
  
  plt<-ggplot(pdat, aes(eval(as.name(xname)), eval(as.name(yname)))) + geom_point(alpha = 0.1) +
    xlab(xlabel) + ylab(ylabel) + ggtitle(mytitle, subtitle = mysubt) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 12),
                       strip.text.y = element_text(size = 15))
  
  
  if(eqaxes){
    # make axes equal, span biggest data
    myaxes<-c(min(c(pdat[,min(get(xname))], pdat[,min(get(yname))])), max(c(pdat[,max(get(xname))], pdat[,max(get(yname))])))
    plt<-plt + xlim(myaxes) + ylim(myaxes)
  }else{
    if(!is.na(giveaxes[1]) & !is.na(giveaxes[2])){
      plt<-plt + xlim(giveaxes) + ylim(giveaxes)
    }
  }
  
  return(plt)
}

asepersamp<-function(aseddsf, gstrains, zeropad = 0.1){
  # Generates ASE estimates (log2(ref/alt)) for each sample individually for each gene in gstrains - simply computes from unnormalized counts
  # In: aseddsf, path to R saved object containing dds.ase, the DESeq2 ASE testing object made by ase_de_annotategenes_deseq2_fromemaseout.R
  #     gstrains, two-column data.table: strain, gene_id. For all the genes we want to do this testing in. 
  #     zeropad, what to pad 0s with to get an approx LFC from ratio when 0 is denominator [decisions like these are why I let the experts handle this...]
  # Out: data.table with one row per input gene/strain combination in gstrains. Columns:
  #       strain, strain these estimates are for
  #       gene_id, gene ID - gene in this strain estimates are for
  #       samp_<1-n>_aselfc, one per sample in sdat; suffix number is number of sample when sorted by name order.
  #             Contains log2(allele 1 counts / allele 2 counts) for this sample and gene
  
  # Subfunctions
  onease<-function(gdat, sdat, zeropad = 0.1){
    # computes ASE (log2(Allele1/ Allele1)) for all genes for all samples
    # In: gdat, gene counts data to compute. 
    #     sdat, sample data information. row names are the column names in grow. Used to find the columns to compare/divide
    #     zeropad, what to pad 0s with to get an approx LFC from ratio when 0 is denominator
    # Out: data.table with one row per input gene in gdat. Columns:
    #       gene_id, gene ID
    #       samp_<1-n>_aselfc, one per sample in sdat; suffix number is number of sample when sorted by name order.
    #             Contains log2(allele 1 counts / allele 2 counts) for this sample and gene
    
    # Do per sample - sorted, mind you!
    out<-lapply(sort(unique(sdat$SampleID)), function(samp){
      sone<-sdat[sdat$SampleID==samp, ]
      al1<-rownames(sone[sone$Allele==unique(sone$Allele1), ])
      al2<-rownames(sone[sone$Allele==unique(sone$Allele2), ])
      
      mygdat<-data.table(al1 = as.numeric(gdat[, al1]),
                         al2 = as.numeric(gdat[, al2]))
      # Zero pad
      mygdat[al1==0, al1:=zeropad]
      mygdat[al2==0, al2:=zeropad]
      
      # Test & return
      return(mygdat[, log2(al1/al2)])
    })
    
    outdt<-data.table(gene_id = row.names(gdat), do.call(cbind, out))
    setnames(outdt, c("gene_id", paste("samp", 1:length(out), "aselfc", sep = "_")))
    return(outdt)
  }
  
  # Load data
  load(aseddsf)
  
  # Do for each strain separately
  out<-rbindlist(lapply(gstrains[,unique(strains)], function(strn){
    # Narrow to strain
    sdat<-colData(dds.ase)[colData(dds.ase)$Allele2==strn, ]
    gdat<-counts(dds.ase)[gstrains[strain==strn, gene_id], rownames(sdat)]
    
    # Do computations & return
    asedat<-onease(gdat, sdat, zeropad = zeropad)
    asedat[, strain:=strn]
    setcolorder(asedat, "strain")
    return(asedat)
  }))
  
  return(out)
}



#### Arguments and inputs ####
p<-arg_parser("Examine cis-trans estimates for artifactual correlation (see Fraser 2019; Zhang & Emerson 2019), experiment with corrections/ways to think about this
FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R and initial cis/trans classifications from ase_de_cistransclassifications.R", 
              name = "cistrans_estimates_explorations.R", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--resmatinput",
                help = "Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
                type = "character")
p<-add_argument(p, "--regpats",
                help = "Regulatory pattern classifications for all genes, strains for which classification could be made - columns strain, gene_id, regclass.
                (*_regpattern_per1plusunqalnggenestrain_<conf int>confdeaseoverlap.txt.gz output of ase_de_cistransclassifications.R)",
                type = "character")
p<-add_argument(p, "--aseddsf",
                help = "path to R saved object containing dds.ase, the DESeq2 ASE testing object made by ase_de_annotategenes_deseq2_fromemaseout.R",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files. Recommend including the informative threshold here if that is ever a thing that could change.",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally.",
                type = "character")
p<-add_argument(p, "--strains",
                help = "Strains to read in results for and process together. Either comma-separated (no spaces) list or path to no-header file with one line per strain.
                Must match how strains are named in input filenames.",
                type = "character")
p<-add_argument(p, "--refstrain",
                help = "Name of reference strain, matching name used in columns of input. Also used for plot labelling.",
                default = "N2")

# Threshold arguments
p<-add_argument(p, "--informthresh",
                help = "Gene must have this or more unique alignments in each sample to be considered informative for ASE/cis-trans analyses",
                default = 5)

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# Output directory
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

# --- Read in data
# combining strains and adding strain column
cat("....Reading in data....\n")
dat<-freadcombstrains(exampinput = p$resmatinput, strains)
setkey(dat, strain, gene_id)

# Add regulatory pattern
rp<-fread(p$regpats, header = T)
setkey(rp, strain, gene_id)
dat<-rp[dat]

# **Narrow to only informative genes for this analysis [new workflow]
dat<-dat[unqalnmts.min >= p$informthresh,]

# Make strains factors in order of interest for plotting
dat[,strain:=factor(strain, levels = strains)]

# Regulatory patterns** change if change something upstream **
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")

#### Examine initial estimates' resulting cis-trans estimates ####
# cis estimate is just ASE in log2 space, trans estimate is parental log2FC - ASE log2FC
cat("....Examining cis/trans correlations from provided log2FC estimates....\n")
# Get & save correlations
cistranscorrs<-rbindlist(lapply(strains, function(x){
  cistranscorr(onedat = dat[strain==x, ], ciscol = "log2FoldChange.ASE", 
               parentcol = paste0("log2FoldChange.ParentVsParent", p$refstrain),
               metadat = data.table(strain = x))
}))
write.table(cistranscorrs, file.path(p$outdir, paste0(p$baseoutname, "_cistranscorrelations_inputLFCs.txt")),
            sep = "\t", quote = F, row.names = F)

#--- Plot scatter plots of cis trans estimates
# get labels for facets
cistranscorrs[, label:=paste("R:", round(r, 2), "\nRho:", round(rho, 2), "\nTau:", round(tau, 2))]
dat[, transe:=get(paste0("log2FoldChange.ParentVsParent", p$refstrain)) - log2FoldChange.ASE] # makes easier for plotting

# Begin plot
pdf(file.path(p$outdir, paste0(p$baseoutname, "_cisVtrans_scatters_inputLFCs.pdf")), 12, 12)
## All genes etc
myscatter(dat, xname = "log2FoldChange.ASE", yname = "transe", xlabel = "Cis estimate (log2FC ASE)",
          ylabel = "Trans estimate (log2FC parental - log2FC ASE)",
          mytitle = paste0("All informative genes (", p$informthresh, "+ unique alignments)"),
          mysubt = "Axes span full range of data",
          eqaxes = T) +
  facet_wrap(~strain) + geom_label(data = cistranscorrs, mapping = aes(x = Inf, y = Inf, label = label), hjust = 1, vjust = 1)

## Restrict axes to -5, 5
myscatter(dat, xname = "log2FoldChange.ASE", yname = "transe", xlabel = "Cis estimate (log2FC ASE)",
          ylabel = "Trans estimate (log2FC parental - log2FC ASE)", 
          mytitle = paste0("All informative genes (", p$informthresh, "+ unique alignments)"),
          mysubt = "Axes restricted - genes with more extreme estimates excluded",
          eqaxes = F, giveaxes = c(-5, 5)) +
  facet_wrap(~strain) 

## Restrict axes to -2, 2
myscatter(dat, xname = "log2FoldChange.ASE", yname = "transe", xlabel = "Cis estimate (log2FC ASE)",
          ylabel = "Trans estimate (log2FC parental - log2FC ASE)", 
          mytitle = paste0("All informative genes (", p$informthresh, "+ unique alignments)"),
          mysubt = "Axes restricted - genes with more extreme estimates excluded",
          eqaxes = F, giveaxes = c(-2, 2)) +
  facet_wrap(~strain)

# End PDF
invisible(dev.off())

#### Cross-replicate stuff ####
cat("....Looking at cis/trans when the 'cis' estimate used to calculate 'trans' comes from a different replicate....\n")
# --- Generate relevant data: ASE from individual samples related
# Get ASE estimates from each sample individually (from each gene of interest)
persampdat<-asepersamp(aseddsf = p$aseddsf, gstrains = dat[,.(strain, gene_id)], zeropad = 0.1)
setkey(persampdat, strain, gene_id)
dat<-persampdat[dat] # combine all together for ease

# Get trans estimates from all sample pairs (cross and self) 
nsamps<-sum(grepl("samp", names(persampdat))) #       might as well make flexible for multiple samples....
asesampcols <- paste("samp", 1:nsamps, "aselfc", sep = "_") # change if change earlier function...
persamptr<-dat[, lapply(asesampcols, function(x) get(paste0("log2FoldChange.ParentVsParent", p$refstrain)) - get(x))]
names(persamptr)<-paste("samp", 1:nsamps, "transParentLFC", sep = "_")
dat<-data.table(persamptr, dat)
setcolorder(dat, c("strain", "gene_id"))

# SAVE the relevant/new columns
scols<-c("strain", "gene_id", "regclass", "log2FoldChange.ASE", paste0("log2FoldChange.ParentVsParent", p$refstrain), "transe",
         asesampcols, paste("samp", 1:nsamps, "transParentLFC", sep = "_"))
sdat<-dat[, scols, with = F]
setnames(sdat, "transe", "transestimate_fromorigLFCs")
write.table(sdat, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_cistransestimates_replicatepairsetc.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# --- Look at these: correlations & plots. CORRELATION is where the cross vs. self pairs come in
# All pairs to compute (& plot)
pairs2do<-data.table(rbind(matrix(rep(1:nsamps, each = 2), nrow = 3, byrow = T), # self pairs
                           t(combn(1:nsamps, 2)))) # cross pairs
setnames(pairs2do, c("samp_cis", "samp_fortrans"))
pairs2do[, pairclass:=ifelse(samp_cis==samp_fortrans, "Same", "Cross")]
pairs2do[, pairname:=paste0(pairclass, " replicate (", samp_cis, " for cis; ", samp_fortrans, " for trans)")]

# Get correlation estimates..... [long way format?? - possibly re-format long ways]
corrpairedsamps<-rbindlist(lapply(strains, function(strn){
  rbindlist(lapply(1:nrow(pairs2do), function(x){
    cistranscorr(onedat = dat[strain==strn, ], ciscol = paste("samp", pairs2do[x, samp_cis], "aselfc", sep = "_"),
                 parentcol = NA,
                 transcol = paste("samp", pairs2do[x, samp_fortrans], "transParentLFC", sep = "_"),
                 metadat = data.table(strain = strn, pairclass = pairs2do[x, pairclass], pairname = pairs2do[x, pairname]))
  }))
}))
write.table(corrpairedsamps, file.path(p$outdir, paste0(p$baseoutname, "_cistranscorrelations_allreplicatepairs.txt")),
            sep = "\t", quote = F, row.names = F)


# Plots
## Make data long format
pdat<-rbind(dat[, .(numsfrom = "All samples, DESeq2", pairname = "All samples, DESeq2", pairclass = "All samples, DESeq2",
                        strain, gene_id, cis = log2FoldChange.ASE, trans = transe)],
                rbindlist(lapply(1:nrow(pairs2do), function(x){
                  dat[, .(numsfrom = "Replicate pairs", pairname = pairs2do[x, pairname], pairclass = pairs2do[x, pairclass], 
                          strain, gene_id, cis = get(paste("samp", pairs2do[x, samp_cis], "aselfc", sep = "_")),
                          trans = get(paste("samp", pairs2do[x, samp_fortrans], "transParentLFC", sep = "_")))]
                })))

pdat[,strain:=factor(strain, levels = strains)]
## Correlations together too
cistranscorrs[, `:=`(numsfrom = "All samples, DESeq2", pairname = "All samples, DESeq2", pairclass = "All samples, DESeq2")]
corrpairedsamps[, `:=`(numsfrom = "Replicate pairs", label=paste("R:", round(r, 2), "\nRho:", round(rho, 2), "\nTau:", round(tau, 2)))]
setcolorder(corrpairedsamps, names(cistranscorrs))
allcorrs<-rbind(cistranscorrs, corrpairedsamps)

## Plot: all together (big ugly plot)
pdf(file.path(p$outdir, paste0(p$baseoutname, "_cisVtrans_scatters_replicatepairs.pdf")), 17, 16)
## Full axes
myscatter(pdat, xname = "cis", yname = "trans", xlabel = "Cis estimate (log2FC ASE)", 
          ylabel = "Trans estimate (log2FC parental - log2FC ASE)",
          mytitle = paste0("All informative genes (", p$informthresh, "+ unique alignments)"),
          mysubt = "Axes span full range of data",
          eqaxes = T) + facet_grid(strain~pairname) +
  geom_label(data = allcorrs, mapping = aes(x = Inf, y = Inf, label = label), hjust = 1, vjust = 1, size = 3)
## Restrict to -5, 5
myscatter(pdat, xname = "cis", yname = "trans", xlabel = "Cis estimate (log2FC ASE)", 
          ylabel = "Trans estimate (log2FC parental - log2FC ASE)",
          mytitle = paste0("All informative genes (", p$informthresh, "+ unique alignments)"),
          mysubt = "Axes restricted - genes with more extreme estimates excluded",
          eqaxes = F, giveaxes = c(-5, 5)) + facet_grid(strain~pairname)
invisible(dev.off())
  
## Plot how correlations stack up
### make long format (so can split different corr metrics)
longcorr<-rbindlist(list(
  allcorrs[, .(strain, numsfrom, pairname, pairclass, name = "Pearson's r", estimate = r)],
  allcorrs[, .(strain, numsfrom, pairname, pairclass, name = "Spearman's rho", estimate = rho)],
  allcorrs[, .(strain, numsfrom, pairname, pairclass, name = "Kendall's tau", estimate = tau)]
))
longcorr[grepl("Same replicate", pairname), numsfrom:="Same-replicate pair"]
longcorr[grepl("Cross replicate", pairname), numsfrom:="Cross-replicate pair"]

### make plot
pdf(file.path(p$outdir, paste0(p$baseoutname, "_cistranscorrdotplots_allmethods.pdf")), 8, 16)
ggplot(longcorr, aes(numsfrom, estimate)) + geom_dotplot(binaxis = "y", stackdir = "center") +
  xlab("Method for estimating cis/trans") + ylab("Correlation estimate") +
  theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 90), 
                     axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                     title = element_text(size = 17), legend.text = element_text(size = 13),
                     plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 12),
                     strip.text.y = element_text(size = 15)) +
  facet_grid(strain~name) 
invisible(dev.off())

#### Script completion message & session information ####
cat("....cistrans_estimates_explorations.R processing complete! Session information:....\n")
sessionInfo()