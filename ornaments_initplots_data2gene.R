#! /usr/bin/env/ Rscript
# Examine reference bias in ornaments output and get output per transcript and gene (in simplified way - looking at best covered variant)
# by Avery Davis Bell, begun 2024.01.25

# Get library location
if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
cat(paste("Library location:", mylibloc, "\n"))

# Load packages - trying to get to play nice on PACE
require(data.table, lib.loc = mylibloc) # trying to make this play nice on PACE...
require(argparser)
require(ggplot2, lib.loc =  mylibloc) # trying to make this play nice on PACE...

#### Functions ####
getlocusmetadata<-function(gff3file, chrexcl = c(NA), biotypes = c("protein_coding")){
  # Gets locus or coding sequence name, location, strand for all genes, narrowing to those on desired chromosomes 
  #     and with desired biotypes based on inputs
  # In: gff3file, path to *genes only* gff3 file containing info on all genes
  #     chrexcl, vector of chromosomes to exclude (or vector of NA) - matching however they're included in 
  # Out: data.table of gene information from gff. Columns:
  #   gene_id, gene ID, from Name field of gff. e.g. wormbase ID
  #   display_name, nicer-format display name to use. This is the 'locus' name if one exists, otherwise the sequence_name
  #   locus, locus name
  #   sequence_name, sequence name
  #   biotype, biotype of gene
  #   chr, chromosome
  #   start, start position
  #   end, end position
  #   strand, coding strand
  
  # Subfunctions
  formatgff<-function(gff, namesget = c("Name", "locus", "sequence_name", "biotype")){ # originally written in formatgff3_vits_20200803.R
    # Formats GFF3 to have useful columns only; only one column per gene
    # Input: data.table, no column names, of GFF information. *Needs to be pre-Tested for each entry to be GENE ONLY*
    #         namesget, vector of names of information to get from last column of GFF
    # Output: data.table with columns namesget, chr, start, end, strand. Keyed by gene_id.
    
    # subfunctions
    procinfo<-function(oneinfo, namesget){
      # Breaks down last column of GFF. In: one string of this information.
      # namesget is vector of names of information from info to get. e.g. ID, Name, locus, sequence_name, etc
      eachinfo<-strsplit(oneinfo,";")[[1]]
      eachinfo<-rbindlist(lapply(strsplit(eachinfo,"="), as.list))
      setnames(eachinfo,c("name", "value"))
      out<-data.table(matrix(eachinfo[match(namesget, name), value], nrow = 1))
      setnames(out, namesget)
      return(out)
    }
    
    setnames(gff, c("chr", "source", "type", "start", "end", "exclude1", "strand", "exclude2", "info"))
    gff<-data.table(gff, gff[, procinfo(oneinfo = info, namesget = namesget),by = 1:nrow(gff)])
    return(gff[,c(namesget, "chr", "start", "end", "strand"), with = F])
  }
  
  gffinfo<-formatgff(fread(gff3file, skip = 8, header=F))[!chr%in%chrexcl & biotype%in%biotypes, ]
  gffinfo[,display_name:=ifelse(is.na(locus), sequence_name, locus)]
  setnames(gffinfo, "Name", "gene_id")
  setcolorder(gffinfo, c("gene_id", "display_name"))
  
  setkey(gffinfo, gene_id)
  return(gffinfo)
}

freadornf<-function(alctf, tx2gene, geneinfo){
  # Reads one ornaments allele counts output file (allele_counts.txt<.gz>) and annotates with gene metadata.
  # Keeps only genes in geneinfo
  # *** if file doesn't exist, prints a warning message and makes a 0 row output 
  # In: alctf, path to ornaments allele_counts.txt<.gz output. Currently, this file has no header, 6 columns
  #     tx2gene, data.table mapping transcript IDs to gene IDs. columns transcript_id, gene_id **keyed by transcript id
  #     geneinfo, data.table with one row per gene to include for downstream analyses. Columns must include gene_id and display_name
  #         **keyed by gene_id
  # Out: data.table of allele counts for variants in genes included in geneinfo. Columns:
  #     transcript_id - transcript ID of transcript variant is in. From ornaments input.
  #   gene_id - Gene id of transcript this variant is in.
  #   variant_pos - transcript position (I think) of this variant - from ornaments input
  #   ref_allele - ref allele of this variant - from ornaments input
  #   alt_allele - alt allele of this variant - from ornaments input
  #   other_allele - column for possible other allele, almost (always?) NA - from ornaments input
  #   ref_count - counts of ref allele of this variant - from ornaments input
  #   alt_count - counts of alt allele of this variant - from ornaments input
  #   other_count - counts of possible other allele, almost (always?) NA - from ornaments input
  
  if(!file.exists(alctf)){
    # file doesn't exists, print warning and output dummy. Not sure this should happen especially with way script currently formatted
    msg<-paste("WARNING: allele counts file", alctf, "does not exist. Continuing without it, but you should stop this run if you expect this file to exist.\n")
    cat(msg)
    cts<-data.table(transcript_id = character(), variant_pos = integer(),
                    ref_allele = character(), alt_allele = character(), other_allele = character(),
                    ref_count = numeric(), alt_count = numeric(), other_count = numeric())
  }else{ # file exists
    # Get in counts       # XXXXX PROBABLY SHOULD MAKE THIS MORE FLEXIBLE IF FILE FORMAT CHANGES...
    cts<-fread(alctf, header = F)
    setnames(cts, c("transcript_id", "variant_pos", "ref_allele", "alt_allele", "other_allele", "ref_count", "alt_count", "other_count"))
    setkey(cts, transcript_id)
    
    # Add gene id; Narrow to only genes in gene_info, add display name
    cts<-tx2gene[cts][gene_id%in%geneinfo$gene_id]
  }

  # Return
  return(cts)
}

ctlist2longform<-function(cts.list, sampinfo){
  # Rbinds elements of cts.list after including relevant parts of sampinfo. ROUGHED OUT/SPECIFIC FOR CURRENT DATA
  # cts.list must be named for SampleID that's in sampinfo
  cts<-rbindlist(lapply(names(cts.list), function(sampid){
    dat<-cts.list[[sampid]]
    dat[, `:=`(total_count = alt_count + ref_count, SampleID = sampid, Generation = sampinfo[SampleID==sampid, Generation],
               Strain = ifelse(test = (sampinfo[SampleID==sampid, Generation]=="Parent"),
                               yes = sampinfo[SampleID==sampid, Strain],
                               no = sampinfo[SampleID==sampid, Allele2]),
               RNAExtractRep =  sampinfo[SampleID==sampid, RNAExtractRep])]
    setcolorder(dat, c("SampleID", "Generation", "Strain", "RNAExtractRep"))
    return(dat)
  }))
  return(cts)
}

npersampsumm<-function(cts, metacols = c("Generation", "Strain", "RNAExtractRep"), rthreshs, colcheck = "total_count", n_what = "variants"){
  # Gets number of rows in cts dt with value above each value in rthreshes in colcheck column
  # In: cts, data.table to examine. Must have columns SampleID; all in metacols; value of colcheck
  #     metacols, non-sampleID names of cts columns to include in output
  #     rthreshs, numerical threshold to check colcheck for being >= to
  #     n_what, value for n_what column of output
  # Out: data.table with one row per sample, columns SampleID, metacols, n_what, at_least (threshold this col is for), N (actual number), 
  
  setkey(cts, SampleID)
  # Numbers
  ns<-rbindlist(lapply(rthreshs, function(x) cts[, .(at_least = x,
                                                     n_what = n_what,
                                                     N = sum(get(colcheck)>=x)), 
                                                 by = SampleID]))
  # Metadata
  setkey(ns, SampleID)
  ns<-unique(cts[,c("SampleID", metacols), with = F])[ns]
  
  # Return
  return(ns)
}

refprop<-function(cts, metacols = c("Generation", "Strain", "RNAExtractRep"), rthresh, colcheck = "total_count", combined_across = "variants"){
  # For each sample, calculates TOTAL n reference alleles, n alt alleles, and reference allele proportion and CI
  #     from all the rows that have colcheck >= threshold in rthreshs
  # In: cts, data.table with one row per observation of interest (e.g., SNP, gene, transcript). Must have columns:
  #           SampleID, <any in metacols>, ref_count, alt_count, <value in colcheck>
  #     metacols,  non-sampleID names of cts columns to include in output
  #     rthresh, numerical threshold (ONE) colcheck must meet to be included 
  #     combined_across, value for combined_across column of output
  # Out: data.table with one row per sample per rthreshs value. Columns:
  #     SampleID, <any in metacols>, at_least (threshold value), combined_across (description from input),
  #     ref_count (summed across all input rows meeting criteria),
  #     alt_count (summed across all input rows meeting criteria), count_total (total allele counts), 
  #     prop_ref_alleles (from summed numbers), prop.low95 (lower 95th binomial CI bound on prop ref alleles), 
  #     prop.high95 (upper 95th binomial CI bound on prop ref alleles)
  #       ****NB PROPS, CIS FOR ROUNDED NUMBERS
  
  # Subfunctions
  propcidt<-function(x, n){
    # Runs binomial test and gets proportion and CI in data.table format
    # output: one row data.table with columns prop, lowci, highci
    res<-binom.test(x, n)
    out<-data.table(prop = res$estimate, lowci = res$conf.int[1], highci = res$conf.int[2])
    return(out)
  }
  
  # Narrow to appropriate data
  dat<-cts[get(colcheck) >= rthresh, ]
  setkey(dat, SampleID)
  
  # Get counts, prop & CI
  out<-dat[, .(at_least = rthresh, combined_across = combined_across, ref_count = sum(ref_count),
               alt_count = sum(alt_count), count_total = sum(ref_count + alt_count)),
           by = SampleID]
  props<-out[, propcidt(x = round(ref_count[1]), n = round(count_total[1])), by = SampleID]
  setnames(props, c("SampleID", "prop_ref_alleles", "prop.low95", "prop.high95"))
  
  # Re-annotate w/ metadata and Return
  outret<-unique(cts[,c("SampleID", metacols), with = F])[out][props]
  return(outret)
}

forestplot<-function(data, ptcol = "unnorm_pref", ptlowcol = "unnorm_pref.low95ci",
                     pthighcol = "unnorm_pref.high95ci", colby = "AltStrain", labcolby = "Strain",
                     shapeby = "Treatment", labshapeby = NA, vertline = 0.5, mytitle = "", myxlab = "",
                     myylab = ""){
  # Makes a forest plot (stacked lines with CIs), e.g. for a proportion
  # In: data, data.table with all the columns specified by other variables
  #     ptcol, column with main value to plot
  #     ptlowcol, column with low CI/error bar end to plot
  #     pthighcol, column with high CI/error bar end to plot
  #     colby, column to color points by
  #     labcolby, label for color legend
  #     shapeby, column to shape points by OR NA for no shaping!!
  #     labcolby, label for shape legend. if NA, shapeby column title will be used
  #     vertline, where to place vertical line if one is desired
  #     mytitle, title
  #     myxlab, xlabel
  #     myylab, ylabel
  # out: ggplot object
  
  # Data set up
  pdata<-copy(data)
  if(!is.na(shapeby)){
    pdata<-pdata[order(get(colby), get(shapeby), decreasing = T),]
  }else{
    pdata<-pdata[order(get(colby), decreasing = T),]
  }
  pdata[,ct:=1:nrow(pdata)]
  
  if(is.na(labshapeby) & !is.na(shapeby)){
    labshapeby<-shapeby
  }
  
  # Main plot
  if(!is.na(shapeby)){
    plt<-ggplot(pdata, aes(eval(as.name(ptcol)), ct, xmin = eval(as.name(ptlowcol)),
                           xmax = eval(as.name(pthighcol)))) + 
      geom_pointrange(aes(color = eval(as.name(colby)), shape = eval(as.name(shapeby)))) +
      ggtitle(mytitle) + xlab(myxlab) + ylab(myylab) +
      labs(color = labcolby, shape = labshapeby)
  }else{
    plt<-ggplot(pdata, aes(eval(as.name(ptcol)), ct, xmin = eval(as.name(ptlowcol)),
                           xmax = eval(as.name(pthighcol)))) + 
      geom_pointrange(aes(color = eval(as.name(colby)))) +
      ggtitle(mytitle) + xlab(myxlab) + ylab(myylab) +
      labs(color = labcolby)
  }
  plt<-plt + theme_bw() +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(), axis.title.x = element_text(size = 15),
                   axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), 
                   legend.text=element_text(size = 13), legend.title = element_text(size = 15),
                   title = element_text(size = 17), strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15)) 

  # Add vertical line if desired
  if(is.numeric(vertline)){
    plt<-plt + geom_vline(xintercept = 0.5, col = "darkgray", lty = "dashed")
  }
  
  return(plt)
}

#### Arguments & Input ####
# --- Command line arguments
p<-arg_parser("Examine reference bias in ornaments output and get output per transcript and gene (in simplified way - looking at best covered variant)", 
              name = "ornaments_initplots_data2gene.R ", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--sampinfo",
                help = "Path to sample information file containing information on parental and F1 samples. Columns should include:
                SampleID, Generation (Parental or F1), Allele1 (NA for parental, reference strain for F1), Allele2 (NA for parental, non-reference strain for F1 - used as Strain for F1 in generation/strain model),
                Strain (NA for F1, strain for parental. All strains/alleles should be the way you want them to show up in outputs)",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character")
p<-add_argument(p, "--allelecounts",
                help = "example filepath to ornaments allele_counts.txt (or allele_counts.txt.gz) RNA quantifiaction file for one sample. Transcripts in name-sorted order. Where each Sample ID goes, needs to have _SAMPID_ (e.g. path/to/file/_sampid_/allele_counts.txt.gz).
                For any other differences in filepath, include * for interpolation.",
                type = "character")
p<-add_argument(p, "--strains",
                help = "Strains in order you'd like them to be plotted - as in sampinfo. ***INCLUDE REF STRAIN***. Either comma separated or path to no-header file",
                type = "character")


# Arguments about gene information/inclusion
p<-add_argument(p, "--tx2genef",
                help = "Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.",
                type = "character")
p<-add_argument(p, "--genegff",
                help = "Path to *genes only* gff3 file containing info on all gene_ids present in input counts file; includes gene location, name, biotype and other information.",
                type = "character")
p<-add_argument(p, "--exclchrs",
                help = "Optional comma-separated, no space list of chromosomes to exclude entirely (named as in genegff). E.g. MtDNA cat be smart to exclude for ASE analysis.",
                default = NA,
                type = "character")
p<-add_argument(p, "--inclbiotype",
                help = "Comma-separated list of biotypes (as in genegff) to include in processing",
                default = "protein_coding",
                type = "character")

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# Output directory
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# other business
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

if(!is.na(p$exclchrs)){
  chrexcl<-strsplit(p$exclchrs, split = ",", fixed = T)[[1]]
}else{
  chrexcl<-NA
}

biotypes<-strsplit(p$inclbiotype, split = ",", fixed = T)[[1]]

# --- Read in gene information genome-wide for only genes to be retained downstream
cat("....Reading in gene information and narrowing to genes from user-specified chromosomes, biotypes....\n")
geneinfo<-getlocusmetadata(gff3file = p$genegff, chrexcl = chrexcl, biotypes = biotypes)
setkey(geneinfo, gene_id)
tx2gene<-fread(p$tx2genef, header = F)
setnames(tx2gene, c("transcript_id", "gene_id"))
setkey(tx2gene, transcript_id)

#### Get data in all required formats ####
# --- Read in allele counts
cat("....Reading in per-sample allele counts....\n")
sampinfo<-fread(p$sampinfo, header = T)
cts.list<-lapply(sampinfo$SampleID, function(sampid){
  myfile<-Sys.glob(gsub("_SAMPID_", sampid, p$allelecounts))
  if(length(myfile==1)){
    out<-freadornf(alctf = myfile, tx2gene, geneinfo)
  }else{
    msg<-paste("WARNING: no single allele counts file for sample", sampid, "!!! Processing will continue without it, interrupt if that seems wrong.\n")
    cat(msg)
    out<-NULL
  }
  return(out)
})
names(cts.list)<-sampinfo$SampleID

# --- Keep variant with highest number reads per transcript
cat("....Collapsing to one variant per transcript (that with highest coverage)....\n")
# Generate data - currently this is SLOW
tcts.list<-lapply(cts.list, function(ctsone){
  out<-rbindlist(lapply(ctsone[, unique(transcript_id)], function(tid){
    dat<-copy(ctsone[tid, ]) # make sure don't mess anything internal
    out<-dat[which.max(ref_count + alt_count), 
             .(transcript_id, gene_id, variant_pos, ref_allele, alt_allele, ref_count, alt_count)]
    return(out)
  }))
  return(out)
})

# Save out data
tdatdir<-file.path(p$outdir, "data/pertranscriptallelecounts")
if(!dir.exists(tdatdir)){dir.create(tdatdir, recursive = T)}
invisible(lapply(names(tcts.list), function(sampid){
  write.table(tcts.list[[sampid]],
              gzfile(file.path(tdatdir, paste0(sampid, "_pertranscript_allelecounts.txt.gz"))),
              sep = "\t", quote = F, row.names = F)
}))


# --- Keep transcript with highest number reads per gene
cat("....Collapsing to one variant (and transcript) per gene (that with highest coverage)....\n")
# Generate data
gcts.list<-lapply(tcts.list, function(tctsone){
  out<-rbindlist(lapply(tctsone[, unique(gene_id)], function(gid){
    dat<-copy(tctsone[gene_id==gid, ]) # make sure don't mess anything internal
    out<-dat[which.max(ref_count + alt_count), 
             .(transcript_id, gene_id, variant_pos, ref_allele, alt_allele, ref_count, alt_count)]
    return(out)
  }))
  return(out)
})

# Save out data
gdatdir<-file.path(p$outdir, "data/pergeneallelecounts")
if(!dir.exists(gdatdir)){dir.create(gdatdir, recursive = T)}
invisible(lapply(names(gcts.list), function(sampid){
  write.table(gcts.list[[sampid]],
              gzfile(file.path(gdatdir, paste0(sampid, "_pergene_allelecounts.txt.gz"))),
              sep = "\t", quote = F, row.names = F)
}))

# --- Summary numbers (& combining into long format)
cat("....Generating per-sample summary numbers (numbers of variants etc covered)....\n")
cts<-ctlist2longform(cts.list = cts.list, sampinfo = sampinfo)
tcts<-ctlist2longform(cts.list = tcts.list, sampinfo = sampinfo)
gcts<-ctlist2longform(cts.list = gcts.list, sampinfo = sampinfo)

# Get number of each thing covered per sample at multiple read thresholds
rthreshs<-c(0, 2, 5, 10, 20) ## add here if want to test different thresholds

nobspersamp<-rbindlist(list(npersampsumm(cts = cts, rthreshs = rthreshs, colcheck = "total_count", n_what = "variants"),
                            npersampsumm(cts = tcts, rthreshs = rthreshs, colcheck = "total_count", n_what = "transcripts"),
                            npersampsumm(cts = gcts, rthreshs = rthreshs, colcheck = "total_count", n_what = "genes")))
## Save
write.table(nobspersamp, file.path(p$outdir, paste0(p$baseoutname, "_numberobsdiffreadcounts.txt")),
            sep = "\t", quote = F, row.names = F)

## Plot summaries
nobspersamp[,Strain:=factor(Strain, levels = strains)]
### threshold vs number of genes
nplt<-ggplot(nobspersamp, aes(at_least, N)) + geom_point(aes(color = Strain)) + 
  xlab("Number of reads (threshold)") + ylab("N") +
  facet_grid(Generation~n_what, scales = "free_y") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
        title = element_text(size = 17), legend.text = element_text(size = 13),
        plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))
pdf(file.path(p$outdir, paste0(p$baseoutname, "_numberobsdiffreadcounts_threshVN.pdf")), 7, 5.5)
print(nplt)
invisible(dev.off())
### Split up by threshold
nplt2<-lapply(c("genes", "transcripts", "variants"), function(x){
  ggplot(nobspersamp[n_what==x], aes(Strain, N)) + geom_dotplot(aes(fill = Strain), binaxis = "y", stackdir = "center", dotsize = 2) +
    facet_grid(at_least~Generation) + ggtitle(paste("N", x)) + ylab(paste("N", x)) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 15), 
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))
})
pdf(file.path(p$outdir, paste0(p$baseoutname, "_numberobsdiffreadcounts_dotplot.pdf")),8, 8)
invisible(lapply(nplt2, print))
invisible(dev.off())

# ....maybe do # shared across samples in same strain too? And if so save that? [can come back later - that's for informative-ness for downstream analyses]

#### Examine reference bias ####
cat("....Analyzing and plotting reference bias metrics (preliminary)....\n")

# Reference allele proportion at multiple threshold across categories
refpersamp<-rbindlist(lapply(rthreshs, function(x){
  rbindlist(list(refprop(cts = cts, rthresh = x, colcheck = "total_count", combined_across = "variants"),
                 refprop(cts = tcts, rthresh = x, colcheck = "total_count", combined_across = "transcripts"),
                 refprop(cts = gcts, rthresh = x, colcheck = "total_count", combined_across = "genes")))
}))
## Save
write.table(refpersamp, file.path(p$outdir, paste0(p$baseoutname, "_proprefallelesdiffthresholds.txt")),
            sep = "\t", quote = F, row.names = F)

# Plots!
## data etc set up
refpersamp[,Strain:=factor(Strain, levels = strains)]

pasteminreads<-function(string){
  return(paste0(string, "+ reads"))
}
  
## Make plot
fplts<-lapply(c("genes", "transcripts", "variants"), function(x){
  forestplot(data = refpersamp[combined_across==x], ptcol = "prop_ref_alleles", ptlowcol = "prop.low95",
             pthighcol = "prop.high95", colby = "Strain", labcolby = "Strain",
             shapeby = NA, labshapeby = NA, vertline = 0.5, 
             mytitle = paste0(x, ": Allele skew across all reads summed"),
             myxlab = "Proportion reference alleles",
             myylab = "") + facet_grid(at_least~Generation, scales = "free_x", labeller = labeller(at_least = pasteminreads, Generation = label_value))
  # getting there. need to update theme in function (port to other script!!!) and do these more separately
})
## Save plots
pdf(file.path(p$outdir, paste0(p$baseoutname, "_proprefallelesdiffthresholds_forestplot.pdf")),8, 8)
invisible(lapply(fplts, print))
invisible(dev.off())


#### Script completion message & session information ####
cat("....ornaments_initplots_data2gene.R processing complete! Session information:....\n")
sessionInfo()
