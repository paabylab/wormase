#! /usr/bin/env/ Rscript
# Allele-specific expression analysis: reference bias, ASE genes across different gene categories
# FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R
# by Avery Davis Bell, begun 2022.03.22

#### System set up - especially required for playing nicely with PACE ####
# Get library location
if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
cat(paste("Library location:", mylibloc, "\n"))

require(data.table, lib.loc = mylibloc)
require(argparser)
require(ggplot2, lib.loc = mylibloc)
require(cowplot)

cat("Session info for early troubleshooting -- \n")
sessionInfo() #### tmp print for troubleshooting

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

# Proportion reference allele related plots
getprefgenesinderes<-function(res, mytitle = "", genesetdescrip = "All genes"){
  # this is a silly big function
  # Looks to quantify reference bias in another way - by looking at genes' log2FoldChanges/altVtotal allele proportions
  # also, looking at genes with log2FoldChanges/altVtotal allele proportions *above or equal to (not-ref-biased)* vs *less than (ref-biased)* 50%
  # Note: no longer information per sample, but rather per test/result given in (usually for a set of samples)
  # In: res, ASE results for one test for all genes of interest. Must be a data.table, must have columns:
  #              altVtotal - % alternate alleles as calculated from DESeq2 log2FC; <0.5 is ref-biased & >0.5 is alt-biased
  #              signifAtThresholds.ASE - T or F, is this gene significant ASE (parameters/what counts as significant separately set up)
  #     mytitle, title for plots, e.g., describing what result this is. ALSO used in output data.tables 'result' column
  #     genesetdescrip, description of genes contained in res - used for plots etc. E.g. "All genes" if all genes included; "Informative genes" if informative genes included
  # Out: list of:
  #   $statres, one-row data.table summarizing statistical test results AND number with ASE. Columns:
  #     result, mytitle input (describing result)
  #     geneset, genesetdescrip input
  #     ngenes, number of genes included
  #     nASE, number of genes with ASE (signifAtThresholds.ASE was T)
  #     propASE, proportion of genes with ASE (nASE/ngenes)
  #     nASE.ref, number of genes with ref-biased ASE (signifAtThresholds.ASE was T and altVtotal < 0.5)
  #     propASE.ref, proportion genes with ref-biased ASE
  #     nASE.alt, number of genes with alt-biased ASE (signifAtThresholds.ASE was T and altVtotal > 0.5)
  #     propgenesMoreRefAlleles, proportion genes in this set with more reference alleles per DESeq2 log2FC
  #     propASE.alt, proportion genes with ref-biased ASE
  #     propASE.low95ci, 95% binomial confidence interval on ASE proportion (lower bound; mostly for plots)
  #     propASE.high95ci, 95% binomial confidence interval on ASE proportion (upper bound; mostly for plots)
  #     propASE.ref.low95ci, 95% binomial confidence interval on ref-biased ASE proportion (lower bound; mostly for plots)
  #     propASE.ref.high95ci, 95% binomial confidence interval on ref-biased ASE proportion (upper bound; mostly for plots)
  #     propASE.alt.low95ci, 95% binomial confidence interval on alt-biased ASE proportion (lower bound; mostly for plots)
  #     propASE.alt.high95ci, 95% binomial confidence interval on alt-biased ASE proportion (upper bound; mostly for plots)
  #     propgenesMoreRefAlleles.low95ci, 95% confidence interval on above proportion (lower bound)
  #     propgenesMoreRefAlleles.high95ci, 95% confidence interval on above proportion (upper bound)
  #     binomtest.pval.2sided, binomial test p-value on this proportion vs. 0.5 (two-sided)
  #     binomtest.pval.refgreater, binomial test p-value on this proportion vs. 0.5 - for the hypothesis that more are reference biased (one-sided)
  #     wilcoxtest.diff50.pval.2sided, wilcox rank test p-value for differences from 0.5 (even alleles) for genes split by whether they have more reference alleles called by DESeq2 or not - about MAGNITUDE of skew. Two sided (testing hypothesis that there is no median difference)
  #     wilcoxtest.diff50.pval.refgreater, wilcox rank test p-value for differences from 0.5 (even alleles) for genes split by whether they have more reference alleles called by DESeq2 or not - about MAGNITUDE of skew. One sided (testing hypothesis that there is more REFERENCE skew)
  #     propRefAlleles.median, median proportion reference alleles across genes
  #     wilcoxtest.medianV50.pval.2sided, wilcox one sample test p-value for median prop ref alleles being different from 0.5
  #     wilcoxtest.medianV50.pval.refgreater, wilcox one sample test p-value for median prop ref alleles being greater than 0.5
  #  $areaplots, list of ggplot overlapped histogram plots showing distance from 0.5 for genes with more ref alleles vs. equal or more alt alleles:
  #     $raw, raw data
  #     $logy, log-scale Y axis
  #   $hist, ggplot basic histogram plots - prop ref alleles estimated per gene. 
  #     $raw, raw data
  #     $logy, log-scale Y axis
  
  # prop ASE & CI
  pase<-binom.test(res[,sum(signifAtThresholds.ASE)], nrow(res))
  prefase<-binom.test(res[,sum(signifAtThresholds.ASE & altVtotal < 0.5)], nrow(res))
  paltase<-binom.test(res[,sum(signifAtThresholds.ASE & altVtotal > 0.5)], nrow(res))
  
  # Sign test - # ref biased (<50% alt alleles) vs. (# ref + alt biased)
  # ONE TAILED to look for ref bias and TWO TAILED to look for any bias
  bres.all.2side<-res[,binom.test(sum(altVtotal < 0.5), sum(altVtotal >= 0.5) + sum(altVtotal < 0.5))]
  bres.all.refbias<-res[,binom.test(sum(altVtotal < 0.5), sum(altVtotal >= 0.5) + sum(altVtotal < 0.5), alternative = "greater")]
  
  # Test median vs. 50% ref alleles (added 2/2024)
  # ONE TAILED to look for ref bias and TWO TAILED to look for any bias
  w50.2side<-wilcox.test(res[, 1-altVtotal], mu = 0.5, alternative = "two.sided")
  w50.refbias<-wilcox.test(res[, 1-altVtotal], mu = 0.5, alternative = "greater")
  
  # Rank test -  take DIFF from 50% for each gene in above and below category & compare them??
  abs.diff50<-res[,abs(altVtotal - 0.5)]
  wres.all.2side<-wilcox.test(abs.diff50[res[,which(altVtotal < 0.5)]], abs.diff50[res[,which(altVtotal >= 0.5)]])
  wres.all.refgreater<-wilcox.test(abs.diff50[res[,which(altVtotal < 0.5)]], abs.diff50[res[,which(altVtotal >= 0.5)]],
                                   alternative = 'greater')
  
  # AREA PLOT difference from 0.5 on top of each other! 1 per result (set of samples). As histograms.
  pdata<-copy(res)
  pdata[,absdiff:=abs.diff50]
  ##  All genes
  diff50.all<-ggplot(pdata, aes(absdiff)) + geom_histogram(aes(fill = altVtotal<0.5), alpha = 0.4, stat = "bin", bins=50, position = 'identity') +
    xlab("Difference from 0.5 each allele (absolute value)") + ylab("Number of genes") + 
    ggtitle(mytitle, subtitle = paste0(genesetdescrip, " (", nrow(pdata), " genes)")) + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.text=element_text(size = 13), legend.title = element_text(size = 15),
          title = element_text(size = 17)) + labs(fill = ">50% ref.\nalleles")
  diff50.all.logy<-ggplot(pdata, aes(absdiff)) + geom_histogram(aes(fill = altVtotal<0.5), alpha = 0.4, stat = "bin", bins=50, position = 'identity') +
    xlab("Difference from 0.5 each allele (absolute value)") + ylab("log10(Number of genes)") + scale_y_log10() +
    ggtitle(mytitle, subtitle = paste0(genesetdescrip, " (", nrow(pdata), " genes)")) + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.text=element_text(size = 13), legend.title = element_text(size = 15),
          title = element_text(size = 17)) + labs(fill = ">50% ref.\nalleles")
  
  # **new - classic histogram too
  myhist<-ggplot(res, aes(1 - altVtotal)) + geom_histogram(stat = "bin", bins = 50, position = 'identity') +
    xlab("Proportion reference alleles called") + ylab("Number of genes") +
    ggtitle(mytitle, subtitle = paste0(genesetdescrip, " (", nrow(res), " genes)")) + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.text=element_text(size = 13), legend.title = element_text(size = 15),
          title = element_text(size = 17)) 
  myhist.logy<-ggplot(res, aes(1 - altVtotal)) + geom_histogram(stat = "bin", bins = 50, position = 'identity') +
    xlab("Proportion reference alleles called") + ylab("log10(Number of genes)") + scale_y_log10() +
    ggtitle(mytitle, subtitle = paste0(genesetdescrip, " (", nrow(res), " genes)")) + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.text=element_text(size = 13), legend.title = element_text(size = 15),
          title = element_text(size = 17)) 
  
  # Format test results as data.tables (one row per group of samples)
  statres<-data.table(result = mytitle, geneset = genesetdescrip, ngenes = nrow(res),
                      nASE = res[,sum(signifAtThresholds.ASE)],
                      propASE = pase$estimate,
                      nASE.ref = prefase$statistic,
                      propASE.ref = prefase$estimate,
                      nASE.alt = paltase$statistic,
                      propASE.alt = paltase$estimate,
                      propgenesMoreRefAlleles = bres.all.2side$estimate,
                      propASE.low95ci = pase$conf.int[1],
                      propASE.high95ci = pase$conf.int[2],
                      propASE.ref.low95ci = prefase$conf.int[1],
                      propASE.ref.high95ci = prefase$conf.int[2],
                      propASE.alt.low95ci = paltase$conf.int[1],
                      propASE.alt.high95ci = paltase$conf.int[2],
                      propgenesMoreRefAlleles.low95ci = bres.all.2side$conf.int[1],
                      propgenesMoreRefAlleles.high95ci = bres.all.2side$conf.int[2],
                      binomtest.pval.2sided = bres.all.2side$p.value,
                      binomtest.pval.refgreater = bres.all.refbias$p.value, 
                      wilcoxtest.diff50.pval.2sided = wres.all.2side$p.value,
                      wilcoxtest.diff50.pval.refgreater = wres.all.refgreater$p.value,
                      propRefAlleles.median = res[,median(1-altVtotal)],
                      wilcoxtest.medianV50.pval.2sided = w50.2side$p.value,
                      wilcoxtest.medianV50.pval.refgreater = w50.refbias$p.value)
  
  # Return
  return(list(statres = statres, areaplots = list(raw = diff50.all, logy = diff50.all.logy),
         hists = list(raw = myhist, logy = myhist.logy)))
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
                   title = element_text(size = 17)) 
  
  # Add vertical line if desired
  if(is.numeric(vertline)){
    plt<-plt + geom_vline(xintercept = 0.5, col = "darkgray", lty = "dashed")
  }
  
  return(plt)
}

violplot<-function(data, strainorder, mytitle = "", mysubt = "", colstrain = T){
  # makes violin plot of proportion ref alleles. All rows in input data are included
  # Inputs: data, data.table with columns strain and altVtotal
  #         strainorder, vector order of strains
  #       mytitle, mysubt - titling
  #       colstrain, fill in strains' violins
  # Outputs: plot
  
  pdat<-copy(data[,.(strain, altVtotal)])
  pdat[,strain:=factor(strain, levels = rev(strainorder))]
  
  plt<-ggplot(pdat, aes(1-altVtotal, strain)) + geom_violin() +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
    ggtitle(mytitle, subtitle = mysubt) + ylab("Strain") + xlab("Proportion reference alleles called at gene") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.position = "none", legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
  
  if(colstrain){
    cols<-scales::hue_pal()(length(strainorder)) # ensure same colors as are otherwise default assigned are used
    names(cols)<-strainorder
    plt<-plt + geom_violin(aes(fill = strain)) + scale_fill_manual(values = cols)
  }
  
  return(plt)
}

volcanoplot.ase<-function(res, mytitle, myalpha = 0.1, mythresh = 0.6, facetby = NA,
                          sizeby = NA, outlineby = NA, labsizeby = NA, laboutlineby = NA,truncaxis=1e-50){
  # Makes a volcano plot of proportion alt alleles vs. -log10(p value). Should provide only informative genes/those of interest for ASE.
  # Informative genes have different point shapes than non-informative.
  # In: res, data.table of ASE results for one test. Must have columns:
  #         padj.ASE, signifAtThresholds.ASE, altVtotal, <any specified by other arguments e.g. facetby>
  #     mytitle, title for plot
  #     myalpha, p-value threshold to draw on plot
  #     mythresh, magnitude threshold to draw on plot (0.6 will be drawn at 0.4, 0.6)
  #     facetby, optional column of res to facet plots by (NA to include all points on same plot)
  #     sizeby, optional column of res to size points by. ONLY TESTED FOR dnacoverage values!!
  #     outlineby, optional column of res containing TRUE and FALSE as only values to have TRUE values outlined (NA to not separate points by outline)
  #     labsizeby, What to label size legend. If NA and sizeby provided, sizeby used
  #     laboutlineby, What to label outline legend. If NA and outlineby provided, outlineby used
  #     truncaxis, p value below which to truncate axis
  # Out: volc, ggplot of the plot
  
  plt.data<-copy(res) # going to monkey with it - don't want to break anything
  myshapes<-rep(paste0('p>', truncaxis), nrow(plt.data))
  # Restrict p-values if needed; will make point type to show this
  if(plt.data[,min(padj.ASE, na.rm = T)<truncaxis]){
    myshapes[plt.data[,which(padj.ASE<truncaxis)]] <- paste0('p<', truncaxis)
    plt.data[padj.ASE<truncaxis ,padj.ASE:=truncaxis]
  }
  myshapes<-factor(myshapes)
  shape_scale<-c( 21, 24)
  names(shape_scale)<-c(paste0('p>', truncaxis), paste0('p<', truncaxis))
  
  # Plot
  if(is.na(sizeby)){
    volc<-ggplot(data = plt.data, aes(1 - altVtotal, -log10(padj.ASE))) + 
      geom_point(aes(shape = myshapes), stroke = 0.01,
                 fill = ifelse(plt.data$signifAtThresholds.ASE, "blue", "darkgray"), alpha = 0.5) +
      scale_shape_manual(values = shape_scale) + labs(shape = "") +
      geom_hline(yintercept = -log10(myalpha), lty = "dashed", col = "blue") + 
      geom_vline(xintercept = c(mythresh, 1 - mythresh), lty = "dashed", col = "darkgray") +
      ggtitle(mytitle) + 
      xlab("Proportion reference alleles") + ylab(bquote(~-Log[10]~italic(P))) + theme_bw() +
      theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
            axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
            title = element_text(size = 17), legend.position = "bottom", legend.text = element_text(size = 13),
            plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  }else{
    labsizeby<-ifelse(is.na(labsizeby), sizeby, labsizeby)
    
    volc<-ggplot(data = plt.data, aes(1 - altVtotal, -log10(padj.ASE))) + 
      geom_point(aes(shape = myshapes, size = eval(as.name(sizeby))), stroke = 0.01,
                 fill = ifelse(plt.data$signifAtThresholds.ASE, "blue", "darkgray"), alpha = 0.5) +
      scale_shape_manual(values = shape_scale) + labs(shape = "", size = labsizeby) +
      scale_size(range=c(0.1, 6)) +
      geom_hline(yintercept = -log10(myalpha), lty = "dashed", col = "blue") + 
      geom_vline(xintercept = c(mythresh, 1 - mythresh), lty = "dashed", col = "darkgray") +
      ggtitle(mytitle) + 
      xlab("Proportion reference alleles") + ylab(bquote(~-Log[10]~italic(P))) + theme_bw() +
      theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
            axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
            title = element_text(size = 17), legend.text = element_text(size = 13),
            plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  }
  
  # Other optional parts
  if(!is.na(facetby)){
    volc<- volc + facet_wrap(~eval(as.name(facetby)))
  }
  
  if(!is.na(outlineby)){
    laboutlineby<-ifelse(is.na(laboutlineby), outlineby, laboutlineby)
    volc<-volc + geom_point(aes(shape = myshapes, stroke = 0.7, color = eval(as.name(outlineby))), alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = 'red', "FALSE" = NA)) + labs(color = laboutlineby)
  }
  
  return(volc)
}

#### Arguments & inputs ####
# --- Command line arguments
p<-arg_parser("# Allele-specific expression analysis: reference bias, ASE genes across different gene categories.
              Data FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R", 
              name = "ase_refbias_withinstrain.R", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--input",
                help = "Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R or ase_de_annoategenes_deseq2_ornaments.R. (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
                type = "character")
p<-add_argument(p, "--aseformat",
                help = "'emase' or 'ornaments': how was ASE data in above matrix generated (with ase_de_annotategenes_deseq2_fromemaseout.R or ase_de_annoategenes_deseq2_ornaments.R?)",
                default = "emase")
p<-add_argument(p, "--strains",
                help = "Strains to read in results for and process together. Either comma-separated (no spaces) list or path to no-header file with one line per strain.
                Must match how strains are named in input filenames.",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character")

# Thresholds for plotting etc
p<-add_argument(p, "--alpha",
                help = "Alpha p-value threshold USED for FDR-like filtering in input signifAtThresholds.ASE column. (Used here for drawing appropriate plot thresholds)",
                type = "numeric",
                default = 0.05)
p<-add_argument(p, "--alleleskewthresh",
                help = "proportion alleles from one haplotype OR the other for results to be considered significant  USED in input signifAtThresholds.ASE column. (Used here for drawing appropriate plot thresholds).
                I.e., if 0.6, genes with >=60% of one allele (<=40% or >=60% alt. alleles) AND significant p-values were considered significant.",
                type = "numeric",
                default = 0.6)

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# data format
if(!p$aseformat%in%c("emase", "ornaments")){
  stop("--aseformat option must be either 'emase' or 'ornaments'")
}

# outdir
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# strains
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

# --- Read in data, combining strains and adding strain column
cat("....Reading in data....\n")
dat<-freadcombstrains(exampinput = p$input, strains)

#### Set up gene sets and get numerical summaries ####
cat("....Setting up gene sets and generating numerical summaries....\n")

# Gene sets. **ADD HERE IF WANT TO TRY A NEW ONE** 2/2024 updated for ornaments inputs then to be flexible depending on ornaments vs. emase input data
## Set up the categories to intersect
if(p$aseformat=="ornaments"){
  infos<-c("Informative, 1+ ornaments alignments", "Informative, 2+ ornaments alignments", 
           "Informative, 5+ ornaments alignments", "Informative, 10+ ornaments alignments") # 4 separate informative gene sets
  infos.short<-c("1plusOrnAlns", "2plusOrnAlns", "5plusOrnAlns", "10plusOrnAlns")
  infos.bool<-c("minorncount >= 1", "minorncount >= 2", "minorncount >= 5", "minorncount >= 10")
}else if(p$aseformat=="emase"){
  infos<-c("Informative, 2+ unique alignments", 
           "Informative, 5+ unique alignments", "Informative, 10+ unique alignments") # 4 separate informative gene sets
  infos.short<-c("2plusUnqAlns", "5plusUnqAlns", "10plusUnqAlns")
  infos.bool<-c("unqalnmts.min >= 2", "unqalnmts.min >= 5", "unqalnmts.min >= 10")
}


gens<-c("excl. hyperdivergent, coverage issues") # excluding hypdiv and no coverage (or high coverage)
gens.short<-c("exclhypdivbadcov")
gens.bool<-c("hypdiv == F & lowDNACov == F & highDNACov == F")

## Get all combinations in data.table
gsets<-data.table(descrip = c(infos, paste(infos, rep(gens, each = length(infos)), sep = "; ")),
                    # descrip is description for plot titles, understanding, etc
                  shortname = c(infos.short, paste(infos.short, rep(gens.short, each = length(infos.short)), sep = "_")),
                    # shortname is name for filenaming etc
                  exprtxt = c(infos.bool, paste(infos.bool, rep(gens.bool, each = length(infos.bool)), sep = " & ")))
                    # exprtxt is how to pull out these rows (with eval(parse(text=exprtxt)))

# Numerical summaries
nums.list<-lapply(1:nrow(gsets), function(x){
  perstrn<-lapply(strains, function(y){
    getprefgenesinderes(res = dat[eval(parse(text = gsets[x, exprtxt])) & strain == y, ],
                        mytitle = y, genesetdescrip = gsets[x, descrip])
  })
  names(perstrn)<-strains
  return(perstrn)
})
names(nums.list)<-gsets$shortname
## combine across strains, genesets & save
ngenes<-rbindlist(lapply(nums.list, function(x) rbindlist(lapply(x, function(y) y$statres))))
write.table(ngenes, file.path(p$outdir, paste0(p$baseoutname, "_numbers_ase_genewisealleleskew.txt")), quote = F, sep = "\t", row.names = F) 
## Add other columns in case helpful (for needing i.e. shorter plot names for facetting)
ngenes<-merge(gsets[,.(descrip, shortname)], ngenes, by.x = "descrip", by.y = "geneset", all.y = T)
setnames(ngenes, "descrip", "geneset")

#### Plots ####
cat("....Generating summary plots....\n")
# ---- set up output directories
vplotdir<-file.path(p$outdir, "volcanoplots")
alskewdir<-file.path(p$outdir, "globalalleleskew")
proprefdir<-file.path(p$outdir, "proprefpergene")
asedir<-file.path(p$outdir, "aseplots")
if(!dir.exists(vplotdir)){dir.create(vplotdir)}
if(!dir.exists(alskewdir)){dir.create(alskewdir)}
if(!dir.exists(proprefdir)){dir.create(proprefdir)}
if(!dir.exists(asedir)){dir.create(asedir)}

# ---- all gene sets separately
# Save 'folded over' histograms that're already created
invisible(
lapply(names(nums.list), function(x){
  pdf(file.path(alskewdir, paste0(p$baseoutname, "_geneskewfrom50plots", x, ".pdf")), 16, 12)
  print(plot_grid(plotlist = lapply(nums.list[[x]], function(y) y$areaplots$raw)))
  print(plot_grid(plotlist = lapply(nums.list[[x]], function(y) y$areaplots$logy)))
  invisible(dev.off())
  return(NULL)
})
)

# Forest plots
ngenes[,result:=factor(result, levels = strains)] ## plot in order they were given
## One per page
pdf(file.path(alskewdir, paste0(p$baseoutname, "_propgenesrefskew_forestplots.pdf")), 7, 5.5)
invisible(
lapply(1:nrow(gsets), function(x){

  print(
  forestplot(data = ngenes[geneset==gsets[x, descrip], ], ptcol = "propgenesMoreRefAlleles",
             ptlowcol = "propgenesMoreRefAlleles.low95ci", pthighcol = "propgenesMoreRefAlleles.high95ci",
             colby = "result", labcolby = "Strain", shapeby = NA, mytitle = gsets[x, descrip], 
             myxlab = "Proportion genes with nominally more reference alleles called", myylab = "")
  )
  return(NULL)
})
)
invisible(dev.off())
## All together - faceted by gene set

pdf(file.path(alskewdir, paste0(p$baseoutname, "_propgenesrefskew_forestplots_onepage.pdf")), 12, 8)
forestplot(data = ngenes, ptcol = "propgenesMoreRefAlleles",
           ptlowcol = "propgenesMoreRefAlleles.low95ci", pthighcol = "propgenesMoreRefAlleles.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", 
           myxlab = "Proportion genes with nominally more reference alleles called", myylab = "") + 
  facet_wrap(~factor(shortname, levels = gsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# Histograms of allelic proportion (do w/r/t ref, probably)
invisible(
  lapply(names(nums.list), function(x){
    pdf(file.path(proprefdir, paste0(p$baseoutname, "_proprefalleleshists", x, ".pdf")), 16, 12)
    print(plot_grid(plotlist = lapply(nums.list[[x]], function(y) y$hists$raw)))
    print(plot_grid(plotlist = lapply(nums.list[[x]], function(y) y$hists$logy)))
    invisible(dev.off())
    return(NULL)
  })
)


# Violin plots of allelic proportion per gene - this is new, goes back to actual data
pdf(file.path(proprefdir, paste0(p$baseoutname, "_proprefalleles_violplots.pdf")), 7, 5.5)
invisible(
  lapply(1:nrow(gsets), function(x){
    print(
      violplot(data = dat[eval(parse(text = gsets[x, exprtxt])), ], strainorder = strains,
                 mytitle = gsets[x, descrip], colstrain = T)
    )
    return(NULL)
  })
)
invisible(dev.off())

#           LABEL WITH MW P VALUE? make labeled and unlabeled versions??
#           pass through gene set?



# ---- some gene sets on same plot (some categories differentiated with shape, color etc)
# Volcano plots

## All informative thresholds, with points outlined by hyperdivergence & sized by DNA coverage [only some of these are super relevant]
pdf(file.path(vplotdir, paste0(p$baseoutname, "_volcanoplots_allgenesets_covsizehypdivoutline.pdf")), 10, 7)
invisible(lapply(1:nrow(gsets), function(x){
  print(volcanoplot.ase(res = dat[eval(parse(text = gsets[x, exprtxt])), ], mytitle = gsets[x, descrip], 
                        myalpha = p$alpha, mythresh = p$alleleskewthresh, facetby = "strain", 
                        sizeby = "dnacoverage", outlineby = "hypdiv", labsizeby = "DNA cov. (med. norm.)", 
                        laboutlineby = "Hyperdivergent"), truncaxis = 1e-20)
}))
invisible(dev.off())

## T/F DNA coverage...
pdf(file.path(vplotdir, paste0(p$baseoutname, "_volcanoplots_allgenesets_lowcovoutline.pdf")), 8, 7)
invisible(lapply(1:nrow(gsets), function(x){
  print(volcanoplot.ase(res = dat[eval(parse(text = gsets[x, exprtxt])), ], mytitle = gsets[x, descrip], 
                        myalpha = p$alpha, mythresh = p$alleleskewthresh, facetby = "strain", 
                        outlineby = "lowDNACov", laboutlineby = "Low DNA coverage"), truncaxis = 1e-20)
}))
invisible(dev.off())

# maybe - mean expression [from parentals or...?] vs. skew or the like? probably don't go to extra trouble

# ---- ASE plots

# proportion ASE plots. 
#     Plot % called ASE vs gene set on same plots - lines/colors connecting strain
#     w/95%CI (binomial) on ASE...but note the CI is BINOMIAL assuming the cutoff between ASE and not is right/binary/etc
#   If among strain relationship holds, that tells us something (not scaling with divergence, etc) [to notebook]
#   also do same sort of plot for ref skew!

# Forest plot of proportion genes with ASE - all gene sets as facets
pdf(file.path(asedir, paste0(p$baseoutname, "_propgenesASE_forestplots_onepage.pdf")), 12, 8)
forestplot(data = ngenes, ptcol = "propASE",
           ptlowcol = "propASE.low95ci", pthighcol = "propASE.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", vertline = F,
           myxlab = "Proportion genes with ASE called", myylab = "") + 
  facet_wrap(~factor(shortname, levels = gsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# Forest plot: proportion of ASE genes that are ref-biased
## Get data formatted this way
ngenes<-data.table(ngenes, rbindlist(lapply(1:nrow(ngenes), function(x){
  res<-ngenes[x, binom.test(nASE.ref, nASE.ref + nASE.alt)]
  return(data.table(pRefOfASE = res$estimate,
                    pRefOfASE.low95ci = res$conf.int[1],
                    pRefOfASE.high95ci = res$conf.int[2]))
})))
## Make plot
pdf(file.path(asedir, paste0(p$baseoutname, "_propgenesrefbiasedOfASE_forestplots_onepage.pdf")), 12, 8)
forestplot(data = ngenes, ptcol = "pRefOfASE",
           ptlowcol = "pRefOfASE.low95ci", pthighcol = "pRefOfASE.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", 
           myxlab = "Proportion ASE genes that are reference-biased", myylab = "") + 
  facet_wrap(~factor(shortname, levels = gsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# Forest plot showing proportion of genes in set that are ASE; ref-biased ASE; alt-biased ASE all together
## Reformat data
bydirns<-copy(ngenes)
bydirns<-rbindlist(list(bydirns, bydirns, bydirns))
bydirns[, dir:=rep(c("Total", "Ref-biased", "Alt-biased"), each = nrow(bydirns)/3)] # ngenes evaluates inside, breaks everything
bydirns[dir=="Total", `:=`(prop = propASE, lowci = propASE.low95ci, highci = propASE.high95ci)]
bydirns[dir=="Ref-biased", `:=`(prop = propASE.ref, lowci = propASE.ref.low95ci, highci = propASE.ref.high95ci)]
bydirns[dir=="Alt-biased", `:=`(prop = propASE.alt, lowci = propASE.alt.low95ci, highci = propASE.alt.high95ci)]

## Make plot
pdf(file.path(asedir, paste0(p$baseoutname, "_propgenesDiffDirsASE_forestplots_onepage.pdf")), 12, 8)
### Total, ref-biased, alt-biased
forestplot(data = bydirns, ptcol = "prop", ptlowcol = "lowci", pthighcol = "highci",
           colby = "result", labcolby = "Strain", shapeby = "dir", labshapeby = "ASE direction",
           mytitle = "", myxlab = "Proportion genes with given ASE", myylab = "",
           vertline = F) +
  facet_wrap(~factor(shortname, levels = gsets$shortname)) + theme(strip.text.x = element_text(size = 15))
### Just ref & alt biased (total is naturally 2x larger than these)
forestplot(data = bydirns[dir%in%c("Ref-biased", "Alt-biased")],
           ptcol = "prop", ptlowcol = "lowci", pthighcol = "highci",
           colby = "result", labcolby = "Strain", shapeby = "dir", labshapeby = "ASE direction",
           mytitle = "", myxlab = "Proportion genes with given ASE", myylab = "",
           vertline = F) +
  facet_wrap(~factor(shortname, levels = gsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# XXX SANDBOX - PLOT VS. DIVERGENCE FROM N2??
#     but might want more thinking; this doesn't exclude # variants in hypdiv, low coverage genes, for example - might want to do that
# nSNPsIndels = c(176595, 127926, 79823, 366034) # VERY one off; from worm21 poster
# nsnpsindels<-data.table(n = nSNPsIndels, strain = c("CB4856", "EG4348", "JU1088", "QX1211"))
# tmp<-copy(ngenes)
# tmp<-merge(nsnpsindels, tmp, by.x = "strain", by.y = "result", all.y = T)
# ggplot(tmp, aes(n, propASE)) + geom_pointrange(aes(ymin = propASE.low95ci, ymax = propASE.high95ci, color =strain)) + facet_wrap(~shortname) + xlab("Number SNPs + INDELs between genome and N2") + ylab("Proportion genes called ASE")
# 

#### Look within hyperdivergent genes - a bit of a tack on/one off (don't want plots to be in facets with others) ####
cat("....Looking within genes in hyperdivergent haplotypes only....\n")
# directory
hypdir<-file.path(p$outdir, "hypdiv")
if(!dir.exists(hypdir)){dir.create(hypdir)}

# Set up gene sets
if(p$aseformat=="ornaments"){
  nunqs<-c(1, 2, 5, 10)
  hgsets<-data.table(descrip = paste0("Hyperdivergent, ", nunqs, "+ ornaments alignments"),
                     shortname = paste0("hypdiv", nunqs, "plusOrnAlns"),
                     exprtxt = paste("hypdiv == T &", "minorncount >=", nunqs))
}else if(p$aseformat=="emase"){
  nunqs<-c(2, 5, 10)
  hgsets<-data.table(descrip = paste0("Hyperdivergent, ", nunqs, "+ unique alignments"),
                     shortname = paste0("hypdiv", nunqs, "plusUnqAlns"),
                     exprtxt = paste("hypdiv == T &", "unqalnmts.min >=", nunqs))
}

h.nums.list<-lapply(1:nrow(hgsets), function(x){
  perstrn<-lapply(strains, function(y){
    getprefgenesinderes(res = dat[eval(parse(text = hgsets[x, exprtxt])) & strain == y, ],
                        mytitle = y, genesetdescrip = hgsets[x, descrip])
  })
  names(perstrn)<-strains
  return(perstrn)
})
names(h.nums.list)<-hgsets$shortname
## combine across strains, genesets & save
hyp.ngenes<-rbindlist(lapply(h.nums.list, function(x) rbindlist(lapply(x, function(y) y$statres))))
write.table(hyp.ngenes, file.path(hypdir, paste0(p$baseoutname, "_hypdivgenes_numbers_ase_genewisealleleskew.txt")), quote = F, sep = "\t", row.names = F) 
## Add other columns in case helpful (for needing i.e. shorter plot names for facetting)
hyp.ngenes<-merge(hgsets[,.(descrip, shortname)], hyp.ngenes, by.x = "descrip", by.y = "geneset", all.y = T)
setnames(hyp.ngenes, "descrip", "geneset")

# --- Plots - unrelated to ASE calls
# Save 'folded over' histograms that're already created
invisible(
  lapply(names(h.nums.list), function(x){
    pdf(file.path(hypdir, paste0(p$baseoutname, "_geneskewfrom50plots_", x, ".pdf")), 12, 8)
    print(plot_grid(plotlist = lapply(h.nums.list[[x]], function(y) y$areaplots$raw)))
    print(plot_grid(plotlist = lapply(h.nums.list[[x]], function(y) y$areaplots$logy)))
    invisible(dev.off())
    return(NULL)
  })
)

# Forest plot - total # genes ref, alt biased
pdf(file.path(hypdir, paste0(p$baseoutname, "_hypdivgenes_propgenesrefskew_forestplots_onepage.pdf")), 12, 5)
forestplot(data = hyp.ngenes, ptcol = "propgenesMoreRefAlleles",
           ptlowcol = "propgenesMoreRefAlleles.low95ci", pthighcol = "propgenesMoreRefAlleles.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", 
           myxlab = "Proportion genes with nominally more reference alleles called", myylab = "") + 
  facet_wrap(~factor(shortname, levels = hgsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# --- ASE call related ---
# Forest plot of proportion genes with ASE - all gene sets as facets
pdf(file.path(hypdir, paste0(p$baseoutname, "_hypdivgenes_propgenesASE_forestplots_onepage.pdf")), 12, 5)
forestplot(data = hyp.ngenes, ptcol = "propASE",
           ptlowcol = "propASE.low95ci", pthighcol = "propASE.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", vertline = F,
           myxlab = "Proportion genes with ASE called", myylab = "") + 
  facet_wrap(~factor(shortname, levels = hgsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# Forest plot: proportion of ASE genes that are ref-biased
## Get data formatted this way
hyp.ngenes<-data.table(hyp.ngenes, rbindlist(lapply(1:nrow(hyp.ngenes), function(x){
  res<-hyp.ngenes[x, binom.test(nASE.ref, nASE.ref + nASE.alt)]
  return(data.table(pRefOfASE = res$estimate,
                    pRefOfASE.low95ci = res$conf.int[1],
                    pRefOfASE.high95ci = res$conf.int[2]))
})))
## Make plot
pdf(file.path(hypdir, paste0(p$baseoutname, "_hypdivgenes_propgenesrefbiasedOfASE_forestplots_onepage.pdf")), 12, 5)
forestplot(data = hyp.ngenes, ptcol = "pRefOfASE",
           ptlowcol = "pRefOfASE.low95ci", pthighcol = "pRefOfASE.high95ci",
           colby = "result", labcolby = "Strain", shapeby = NA, mytitle = "", 
           myxlab = "Proportion ASE genes that are reference-biased", myylab = "") + 
  facet_wrap(~factor(shortname, levels = hgsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())

# Forest plot showing proportion of genes in set that are ASE; ref-biased ASE; alt-biased ASE all together
## Reformat data
h.bydirns<-copy(hyp.ngenes)
h.bydirns<-rbindlist(list(h.bydirns, h.bydirns, h.bydirns))
h.bydirns[, dir:=rep(c("Total", "Ref-biased", "Alt-biased"), each = nrow(h.bydirns)/3)] # hyp.ngenes evaluates inside, breaks everything
h.bydirns[dir=="Total", `:=`(prop = propASE, lowci = propASE.low95ci, highci = propASE.high95ci)]
h.bydirns[dir=="Ref-biased", `:=`(prop = propASE.ref, lowci = propASE.ref.low95ci, highci = propASE.ref.high95ci)]
h.bydirns[dir=="Alt-biased", `:=`(prop = propASE.alt, lowci = propASE.alt.low95ci, highci = propASE.alt.high95ci)]

## Make plot
pdf(file.path(hypdir, paste0(p$baseoutname, "_hypdivgenes_propgenesDiffDirsASE_forestplots_onepage.pdf")), 12, 5)
### Total, ref-biased, alt-biased
forestplot(data = h.bydirns, ptcol = "prop", ptlowcol = "lowci", pthighcol = "highci",
           colby = "result", labcolby = "Strain", shapeby = "dir", labshapeby = "ASE direction",
           mytitle = "", myxlab = "Proportion genes with given ASE", myylab = "",
           vertline = F) +
  facet_wrap(~factor(shortname, levels = hgsets$shortname)) + theme(strip.text.x = element_text(size = 15))
### Just ref & alt biased (total is naturally 2x larger than these)
forestplot(data = h.bydirns[dir%in%c("Ref-biased", "Alt-biased")],
           ptcol = "prop", ptlowcol = "lowci", pthighcol = "highci",
           colby = "result", labcolby = "Strain", shapeby = "dir", labshapeby = "ASE direction",
           mytitle = "", myxlab = "Proportion genes with given ASE", myylab = "",
           vertline = F) +
  facet_wrap(~factor(shortname, levels = hgsets$shortname)) + theme(strip.text.x = element_text(size = 15))
invisible(dev.off())


#### Script completion message & session information ####
cat("....ase_refbias_withinstrain.R processing complete! Session information:....\n")
sessionInfo()
