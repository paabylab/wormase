#! /usr/bin/env/ Rscript
# Compare F1 and parental expression within-strain
# FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R
# by Avery Davis Bell, begun 2022.03.30 - but not come back to until 2022.06.07
require(data.table, quietly = T)
require(argparser, quietly = T)
require(ggplot2, quietly = T)
require(cowplot, quietly = T)
require(RColorBrewer, quietly = T)

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

inhmode<-function(dat, refstrain = "N2", alpha = 0.1, lfcthresh = 0.5849625){
  # Initial classification of inheritance mode using F1 vs. parental gene expression results
  #   USES provided 'signifAtThresholds' columns - doesn't currently look at p-values etc
  # In: dat, data.table with one row for every gene for which to classify inheritance patterns. Can be multiple rows (e.g. strains) per genes - this isn't checked.
  #       Necessary columns pertain to F1 vs. parental DE *and* inter-parental DE:
  #           each of the following with suffixes .F1VsParent<refstrain>, .F1VsParentNonRef, and .ParentVsParent<N2>
  #         log2FoldChange.<suff above>, padj.<suff above>
  #     refstrain, name of reference strain as in column names
  #     alpha, p-value threshold for classifying significance - used alone sometimes, in concert with lfctrhesh sometimes (pvalue.ASE and pvalue.ParentVsParent<ref> columns)
  #     lfcthresh, log2FoldChange threshold for classifying significance 
  # Out: none - modifies dat in place to add 'inhmode' column with the following values (see notebook for full classification info):
  #   no_change: no DE parents, F1 and either parent
  #   additive: F1 has opposite-dir DE called between both parents (+/- or -/+), with parents having DE too.
  #         ONLY significance threshold is used to determine F1 vs. parental in this case.
  #   <refstrain>_dominant: F1 has same expression as N2, while N2 and F1 are significantly lower expressed than other parent OR
  #       F1 has same expression as N2, while N2 and F1 are significantly higher expressed than other parent
  #   alt_dominant: F1 has 'same expression as alt, while alt and F1 are significantly lower expressed than N2 parent OR
  #       F1 has same' expression as alt, while alt and F1 are significantly higher expressed than N2 parent
  #   overdominant: F1 more expressed than either parent (both signifAtThresholds; log2FC positive)
  #   underdominant: F1 less expressed than either parent (both signifAtThresholds; log2FC negative)
  #   ambiguous: any log2FC patterns that don't fit in above categories, e.g. parents aren't different but F1 is called different from only one of them
  
  # get column names in easy format (for those that aren't standard)
  parenlfc<-paste0("log2FoldChange.ParentVsParent", refstrain)
  parenp<-paste0("padj.ParentVsParent", refstrain)
  f1reflfc<-paste0("log2FoldChange.F1VsParent", refstrain)
  f1refp<-paste0("padj.F1VsParent", refstrain)
  
  # Seed all with ambiguous designation
  dat[, inhmode:="ambiguous"]
  
  # One rule categories
  ## No change
  dat[(get(parenp) >= alpha | abs(get(parenlfc)) <= lfcthresh) & (get(f1refp) >= alpha | abs(get(f1reflfc)) <= lfcthresh) &
        (padj.F1VsParentNonRef >= alpha | abs(log2FoldChange.F1VsParentNonRef) <= lfcthresh), inhmode:="no_change"]
  ## Overdominant
  dat[(get(f1refp) < alpha & get(f1reflfc) > lfcthresh) & 
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef > lfcthresh), inhmode:="overdominant"]
  ## Underdominant
  dat[(get(f1refp) < alpha & get(f1reflfc) < -1*lfcthresh) & 
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef < -1*lfcthresh), inhmode:="underdominant"]
  
  # Two rule categories (need to check direction)
  ## ref dominant
  dat[(get(parenp) < alpha & get(parenlfc) > lfcthresh) & get(f1refp) > alpha &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef < 0), inhmode:=paste0(refstrain, "_dominant")]
  dat[(get(parenp) < alpha & get(parenlfc) < -1*lfcthresh) & get(f1refp) > alpha &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef > 0), inhmode:=paste0(refstrain, "_dominant")]
  dat[(get(parenp) < alpha & get(parenlfc) > 0) & get(f1refp) > alpha &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef < -1*lfcthresh), inhmode:=paste0(refstrain, "_dominant")]
  dat[(get(parenp) < alpha & get(parenlfc) < 0) & get(f1refp) > alpha &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef > lfcthresh), inhmode:=paste0(refstrain, "_dominant")]
  ## alt dominant
  dat[(get(parenp) < alpha & get(parenlfc) < -1*lfcthresh) & (get(f1refp) < alpha & get(f1reflfc) < 0) & 
        padj.F1VsParentNonRef > alpha, inhmode:="alt_dominant"]
  dat[(get(parenp) < alpha & get(parenlfc) > lfcthresh) & (get(f1refp) < alpha & get(f1reflfc) > 0) & 
        padj.F1VsParentNonRef > alpha, inhmode:="alt_dominant"]
  dat[(get(parenp) < alpha & get(parenlfc) < 0) & (get(f1refp) < alpha & get(f1reflfc) < -1*lfcthresh) & 
        padj.F1VsParentNonRef > alpha, inhmode:="alt_dominant"]
  dat[(get(parenp) < alpha & get(parenlfc) > 0) & (get(f1refp) < alpha & get(f1reflfc) > lfcthresh) & 
        padj.F1VsParentNonRef > alpha, inhmode:="alt_dominant"]
  ## additive
  dat[(get(parenp) < alpha & get(parenlfc) < -1*lfcthresh) & (get(f1refp) < alpha & get(f1reflfc) < 0) &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef > 0), inhmode:="additive"]
  dat[(get(parenp) < alpha & get(parenlfc) > lfcthresh) & (get(f1refp) < alpha & get(f1reflfc) > 0) &
        (padj.F1VsParentNonRef < alpha & log2FoldChange.F1VsParentNonRef < 0), inhmode:="additive"]
}

nsummf1parentde<-function(onedat, refstrain = "N2", strainres, genesetin,
                          inhmodes = c("no_change", "additive", paste0(refstrain, "_dominant"), "alt_dominant", "overdominant", "underdominant")){
  # Gets number summaries of # genes DE between F1 and parents (ref and non-ref parent; overlapping; etc)
  # In: onedat, data.table containing all results to summarize at once (i.e. one strain, one geneset) - one row per gene,
  #        of format generated with ase_de_annotategenes_deseq2_fromemaseout.R AND run through inhmode() function
  #         Must have columns: gene_id, signifAtThresholds.F1VsParent<refstrain>, signifAtThresholds.F1VsParentNonRef, log2FoldChange.F1VsParent<refstrain>, log2FoldChange.F1VsParentNonRef,
  #           inhmode (as added with inhmode() function - values in inhmodes)
  #     refstrain, name of reference strain as in column names of onedat
  #     strainres, name of strain these results came from, used for output
  #     genesetin, name or description of this geneset, used for output
  #     inhmodes, vector of all inheritance modes to count up (from inhmode input column). This is also ORDER they'll be output in.
  # Out: list of data.tables, $ns, $ps,  $ns contains counts of each category, $ps contains proportion of genes in the gene set that fall into this category,
  #     and also non-no_change-inheritance mode genes' breakdown of non-no_change genes 
  # - in LONG format [see details].
  #   Columns of $ns:
  #     strain, strain this result is for
  #     geneset, name of geneset passed from input
  #     ngenes, # genes in this gene set
  #     ngenes_changed, # genes with non-no_change inheritance pattern [so can get breakdown of this]
  #     DEv<refstrain>, # genes with F1 vs. refstrain DE called
  #     DEvNonRef, # genes with F1 vs. alt strain DE called
  #     DEvEither, # genes with F1 vs. either parent DE called
  #     DEvBoth, # genes with F1 vs. both parent DE called
  #     <inheritance mode>, # classified as each inheritance mode
  #   Columns of $ps:
  #     strain, strain this result is for
  #     geneset, name of geneset passed from input
  #     ngenes, # genes in this gene set
  #     ngenes_changed, # genes with non-no_change inheritance pattern [so can get breakdown of this]
  #     propName, which N is numerator of proportion (see descriptions above). 
  #     propDenom, ngenes or ngenes_changed: what is proportion OF? ***CAREFUL with ngenes_changed proportions - squishes out any real difference in total # not no_change genes, which is relevant to many things
  #         only non-no_change inheritance modes have the latter proportion
  #     prop, value of proportion
  #     lowci, binomial 95% CI lower bound
  #     highci, binomial 95% CI upper bound

  sigref<-paste0("signifAtThresholds.F1VsParent", refstrain)
  log2ref<-paste0("log2FoldChange.F1VsParent", refstrain)
  
  # Numbers
  ns<-onedat[,.(strain = strainres,
                geneset = genesetin,
                ngenes = length(inhmode),
                ngenes_changed = sum(inhmode!="no_change"),
                # direction-less counts of DE genes
                DEvRef = sum(get(sigref)), # Differentially expressed vs. reference-strain parent
                DEvNonRef = sum(signifAtThresholds.F1VsParentNonRef),
                DEvEither = sum(get(sigref) | signifAtThresholds.F1VsParentNonRef), # DE vs. ref, non ref, or both
                DEvBoth = sum(get(sigref) & signifAtThresholds.F1VsParentNonRef), # DE vs. both parents
                # Including direction counts of DE genes: inheritance mode
                as.data.table(t(as.matrix(table(inhmode))))
                )]
  ## make sure all inheritance mode columns are included (table() will of course not count 0s) 
  sapply(inhmodes, function(x){
    if(!x%in%names(ns)){
      ns[,newcol:=0]
      setnames(ns, "newcol", x)
    }
  })
  setcolorder(ns, inhmodes, after = 8) # re-order inheritance classifications to be easily interpretable
  setnames(ns, "DEvRef", paste0("DEv", refstrain))
  
  # Proportions (long-ways for plotting)
  ## Of total N genes
  ps<-rbindlist(lapply(names(ns)[!names(ns)%in%c("strain", "geneset", "ngenes", "ngenes_changed")], function(x){
    cis<-data.table(t(ns[, binom.test(get(x), ngenes)$conf.int]))
    setnames(cis, c("lowci", "highci"))
    out<-data.table(ns[, .(strain, geneset, ngenes)],
                    propName = x, propDenom = "ngenes",
                    prop = ns[,get(x)/ngenes],
                    cis)
    return(out)
  }))
  ## Non-no_change genes' proportion
  ps_notcons<-rbindlist(lapply(inhmodes[!inhmodes=="no_change"], function(x){
    cis<-data.table(t(ns[, binom.test(get(x), ngenes_changed)$conf.int]))
    setnames(cis, c("lowci", "highci"))
    out<-data.table(ns[, .(strain, geneset, ngenes)],
                    propName = x, propDenom = "ngenes_changed",
                    prop = ns[,get(x)/ngenes_changed],
                    cis)
    return(out)
  }))
  
  ## combine
  ps<-rbind(ps, ps_notcons)
  
  # Return
  return(list(ns = ns, ps = ps))
}

inhmodescat<-function(pdat, inhmodes, xcol = "log2FoldChange.F1VsParentN2", ycol = "log2FoldChange.F1VsParentNonRef",
                      mytitle = "", mysubt = "", xlabel = bquote(log[2]~"(Fold Change: F1 vs. ref strain)"), ylabel = bquote(log[2]~"(Fold Change: F1 vs. non-ref strain)"),
                      labcolby = "Inheritance mode", outlineby = NA, laboutlineby = NA){
  # Makes a scatter plot of F1s log2FCs vs. parents colored by inheritance mode
  # In: pdat, data with one row per gene/point to plot. Columns xcol value, ycol value, inhmode [containing inhmodes values]
  #     inhmodes, vector of inheritance modes in inhmode column of pdat in order to color them in. **Currently assumes contains 'no_change' and 'ambiguous'**
  #         used to level factor
  #     xcol, character name of column for x axis values, typically F1 vs. reference parent log2FCs
  #     ycol, character name of column for y axis values, typically F1 vs. non-reference parent log2FCs
  #     mytitle, title for plot
  #     mysubt, subtitle for plot
  #     xlabel, x axis label for plot
  #     ylabel, y axis label for plot
  #     labcolby, title of color legend label
  #     outlineby, optional column of res containing TRUE and FALSE as only values to have TRUE values outlined in black (NA to not separate points by outline)
  #     laboutlineby, What to label outline legend. If NA and outlineby provided, outlineby used
  
  # Format data
  pdat<-copy(pdat) # don't want to mess with original data
  ## Make inhmode factor of desired order
  pdat[, inhmode:=factor(inhmode, levels = inhmodes)]
  ## colors
  mycols<-rep(NA, length(inhmodes))
  mycols[inhmodes=="no_change"]<-"lightgray"
  mycols[inhmodes=="ambiguous"]<-"darkgray"
  mycols[!inhmodes%in%c("no_change", "ambiguous")]<-brewer.pal(sum(!inhmodes%in%c("no_change", "ambiguous")), "Set2")
  names(mycols)<-inhmodes
  
  
  # Make plot
  plt<-ggplot(pdat, aes(eval(as.name(xcol)), eval(as.name(ycol)))) + geom_point(aes(fill = inhmode), pch = 21, alpha = 0.6, stroke = 0.001) +
    scale_fill_manual(values = mycols) + labs(fill = labcolby) +
    ggtitle(mytitle, subtitle = mysubt) + xlab(xlabel) + ylab(ylabel) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  if(!is.na(outlineby)){
    laboutlineby<-ifelse(is.na(laboutlineby), outlineby, laboutlineby)
    plt<-plt + geom_point(aes(stroke = 0.7, color = eval(as.name(outlineby))), pch = 21, alpha = 0.6, fill = NA) +
      scale_color_manual(values = c("TRUE" = 'black', "FALSE" = NA)) + labs(color = laboutlineby) +
      guides(color = guide_legend(override.aes = list(size = 5)))
  }
  
  return(plt)
}

stackedbar<-function(pinhclass, xcol = "strain", stackcol = "propName", stacknumcol = "prop", stackorder = NA,
                     barcols,
                     legendlabel = "", myxlab = "Strain", myylab = "Proportion genes", mytitle = "", mysubt = ""){
  # Makes stacked bar chart, one stacked bar per category in x col
  # In: pinhclass, data.table with all data to plot. Must have columns named with values of xcol, stackcol, stacknumcol
  #     xcol, column containing category for X axis of bar chart
  #     stackcol, column containing category to split bars into
  #     stacknumcol, column containing numbers to plot as stacked bar
  #     stackorder, OPTIONAL order in which to plot stackcol values top to bottom
  #     barcols, COLORS for each in stack order
  #     legendlabel
  #     myxlab, myylab, mytitle, mysubt: plot labels as intuitively named
  # Out: ggplot2 barplot object. Can facet this externally!   
  
  pdata<-copy(pinhclass)
  if(!is.na(stackorder[1])){
    # Re-level factor 
    pdata[,mystackcol:=factor(pdata[,get(stackcol)], levels = stackorder)]
  }
  
  plt<-ggplot(pdata, aes(eval(as.name(xcol)), eval(as.name(stacknumcol)))) +
    geom_bar(aes(fill = mystackcol), stat = "identity", position = "stack") + labs(fill = legendlabel) +
    xlab(myxlab) + ylab(myylab) + ggtitle(mytitle, subtitle = mysubt) +
    scale_fill_manual(values = barcols) + # better colors
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  return(plt)
}

divplot<-function(onedat, xcol = "nvars.thou", xlabel = "SNVs/INDELs vs. N2 (thousands)", 
                  ylabel = "", mytitle = "", mysize = 1){
  # One off just to make this plot
  plt<-ggplot(onedat, aes(eval(as.name(xcol)), prop)) +
    geom_pointrange(aes(ymin = lowci, ymax = highci, color = strain), size = mysize) +
    xlab(xlabel) + ylab(ylabel) + ggtitle(mytitle) + labs(color = "Strain") + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 9))
  
  return(plt)
}

propvprop<-function(pdata, prop1name = "N2_dominant", prop1denom = "ngenes",  prop2name = "alt_dominant", prop2denom = "ngenes",
                    colby = "strain", namecolby = "F1", myxlab = "Proportion N2 dominant", myylab = "Proportion alt. dominant",
                    mytitle = "", mysubt = ""){
  # Plots two proportions in pdata vs. each other [converts from longformat]
  # In:
  #   pdata: data.table with columns propName, propDenom, colby value, prop (data), lowci, highci (confidence intervals for prop)
  #   prop1name: propName value of proportion to plot on X axis
  #   prop1denom: propDenom value of proportion to plot on X axis
  #   prop2name: propName value of proportion to plot on Y axis
  #   prop2denom: propDenom value of proportion to plot on Y axis
  #   colby, column of pdata to color points by
  #   namecolby, how to label color legend
  #   myxlab, myylab, mytitle, mysubt: axis titles, main title, subtitle for plot
  # Out: ggplot
  
  # Format data
  dat1<-pdata[propName==prop1name & propDenom==prop1denom, .(prop, lowci, highci)]
  setnames(dat1, c("prop", "lowci", "highci"), paste(c("prop", "lowci", "highci"), "1", sep = "."))
  dat2<-pdata[propName==prop2name & propDenom==prop2denom, .(prop, lowci, highci)]
  setnames(dat2, c("prop", "lowci", "highci"), paste(c("prop", "lowci", "highci"), "2", sep = "."))
  pdat<-cbind(pdata[propName==prop1name & propDenom==prop1denom, .SD, .SDcols = -c("prop", "lowci", "highci")], dat1, dat2)
  
  myaxlim<-c(pdat[,min(c(lowci.1, lowci.2), na.rm = T)], pdat[,max(c(highci.1, highci.2), na.rm = T)])
  
  # Make plot
  plt<-ggplot(pdat, aes(prop.1, prop.2)) + geom_pointrange(aes(ymin = lowci.2, ymax = highci.2, color = eval(as.name(colby)))) +
    geom_segment(aes(x = lowci.1, y = prop.2, xend = highci.1, yend = prop.2, color = eval(as.name(colby)))) + # adding 95% CIs on x axis
    xlim(myaxlim) + ylim(myaxlim) + geom_abline(slope = 1, intercept = 0, lty = "dashed", color = "gray") +
    labs(color = namecolby) + xlab(myxlab) + ylab(myylab) + ggtitle(mytitle, subtitle = mysubt)  + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 9))
  
  return(plt)
}

#### Arguments & inputs ####
# --- Command line arguments
p<-arg_parser("Compare F1 and parental expression within-strain.
              Data FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R", 
              name = "f1_parental_comparisons_withinstrain.R", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--input",
                help = "Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
                type = "character")
p<-add_argument(p, "--aseformat",
                help = "'emase' or 'ornaments': how was ASE data in above matrix generated (with ase_de_annotategenes_deseq2_fromemaseout.R or ase_de_annoategenes_deseq2_ornaments.R?)",
                default = "emase")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
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

# Gene sets to plot for. Automatically does all, excluding hypdiv etc; informative [using threshold below] and excluding hypdiv etc
p<-add_argument(p, "--informthresh",
              help = "Gene must have this or more unique alignments (EMASE; or this or more Ornaments counts) in each sample to be considered informative for ASE/cis-trans analyses",
              default = 5)

# Thresholding for inheritance mode
p<-add_argument(p, "--alpha",
                help = "P-value threshold for considering gene DE for *adjusted p-value*. Combined with magnitude threshold (--lfcthresh) for many inheritance mode classifications, used on its own in some cases where other DE used as a prior.", 
                default = 0.1)
p<-add_argument(p, "--lfcthresh",
                help = "log2FoldChange threshold for considering gene DE when combined with --alpha.", 
                default = 0.5849625)

# Related to plotting vs. divergence
p<-add_argument(p, "--genomebp",
                help = "Length of reference genome in bp. Used to get variants vs. reference per kb.
                Default is genome length from NCBI https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6/ - 
                for ws235 but unlikely to have changed much (and not enough to be meaningful for our purposes)",
                default = 100286401)
p<-add_argument(p, "--varsvsref",
                help = "Path to file containing columns strain (one row for each strain in main input), nvars (number variants vs. reference genome)",
                type = "character")

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# data format
if(!p$aseformat%in%c("emase", "ornaments")){
  stop("--aseformat option must be either 'emase' or 'ornaments'")
}

# Output directory
f1vsparentdir<-file.path(p$outdir)
if(!dir.exists(f1vsparentdir)){dir.create(f1vsparentdir, recursive = T)}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

# --- Read in data, combining strains and adding strain column
cat("....Reading in data....\n")
dat<-freadcombstrains(exampinput = p$input, strains)

#### F1 vs. parents - numerical summaries of DE ####
cat("....Analyzing F1 vs. parental strains DE - inheritance mode classification and numerical summaries....\n")

# Annotate each gene with inheritance mode classification
dat<-inhmode(dat, refstrain = p$refstrain,  alpha = p$alpha, lfcthresh = p$lfcthresh)
inhmodes<-c("no_change", "ambiguous", "additive", paste0(p$refstrain, "_dominant"), "alt_dominant", "overdominant", "underdominant")
    # update if change inhmode function
  ## COLORS
inhcols<-rep(NA, length(inhmodes))
names(inhcols)<-inhmodes
inhcols["no_change"]<-"lightgray"
inhcols["ambiguous"]<-"darkgray"
inhcols[c("additive", paste0(p$refstrain, "_dominant"), "alt_dominant", "overdominant", "underdominant")]<-rev(brewer.pal(5, "Paired"))
## Save out. NB thought about not saving, saving with all data, etc - decided to save just this so can merge back in later and change this script if needed
write.table(dat[,.(strain, gene_id, inhmode)], gzfile(file.path(p$outdir, paste0(p$baseoutname, "_inheritancemode_pergenestrain.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# Gene sets. **ADD HERE IF WANT TO TRY A NEW ONE**
if(p$aseformat=="emase"){
  thrcol<-"unqalnmts.min"
  infdescrip<-"unique alignments"
}else if(p$aseformat=="ornaments"){
  thrcol<-"minorncount"
  infdescrip<-"ornaments alignments"
}

f1gsets<-data.table(descrip = c("All", "Excl. hyperdivergent, coverage issues", paste0("ASE informative (", p$informthresh, "+", infdescrip, "excl. hypdiv and coverage issues)")),
                  shortname = c("all", "exclhypdivbadcov", paste0("ase", p$informthresh, "plusinfoexcls")),
                  exprtxt = c("hypdiv == T | hypdiv == F", "hypdiv == F & lowDNACov == F & highDNACov == F", 
                              paste0(thrcol, ">=", p$informthresh," & hypdiv == F & lowDNACov == F & highDNACov == F")))
# exprtxt is how to pull out these rows (with eval(parse(text=exprtxt))) - first is just one of many ways to make all rows come out

# Get & save numerical summaries per strain, gene set
nspslist<-lapply(1:nrow(f1gsets), function(x){
  lapply(strains, function(y){
    nsummf1parentde(onedat = dat[eval(parse(text = f1gsets[x, exprtxt])) & strain == y, ],
                   refstrain = p$refstrain, strainres = y, genesetin = f1gsets[x, shortname],
                   inhmodes = inhmodes)
  })
})
## Combine
ninhmode<-rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$ns))))
pinhmode<-rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$ps))))
## Save
write.table(ninhmode, file.path(p$outdir, paste0(p$baseoutname, "_numbers_f1parentsdeinhmode.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(pinhmode, file.path(p$outdir, paste0(p$baseoutname, "_proportion_f1parentsdeinhmode.txt")),
            sep = "\t", quote = F, row.names = F)

## Factor updating for plots
pinhmode[,geneset:=factor(geneset, levels = f1gsets$shortname)] # so will plot in order specified

#### bar plots of breakdown of inh mode [one bar per strain, stacked colors of inh mode??] ####
# Proportion of all genes
## All gene sets on one page faceted [useful for seeing all at once, notebook]
pdf(file.path(p$outdir, paste0(p$baseoutname, "_inhmode_propofallgenes_allsets.pdf")), 12, 6)
plt<-stackedbar(pinhclass = pinhmode[propName%in%inhmodes & propDenom=="ngenes", ],
                xcol = "strain", stackcol = "propName", stacknumcol = "prop", stackorder = inhmodes,
                barcols = inhcols,
                legendlabel = "Inheritance mode", myxlab = "Strain", myylab = "Proportion of all genes in set")
print(plt + facet_wrap(~geneset) + theme(axis.text.x = element_text(size = 14, angle = 90)))
invisible(dev.off())
## Gene sets individually, one plot per page [useful for presentations, zooming in]
pdf(file.path(p$outdir, paste0(p$baseoutname, "_inhmode_propofallgenes_indivplots.pdf")), 8, 5.5)
invisible(lapply(1:nrow(f1gsets), function(x){
  plt<-stackedbar(pinhclass = pinhmode[propName%in%inhmodes & propDenom=="ngenes" & geneset ==f1gsets[x, shortname], ],
                  xcol = "strain", stackcol = "propName", stacknumcol = "prop", stackorder = inhmodes,
                  barcols = inhcols,
                  legendlabel = "Inheritance mode", myxlab = "Strain", myylab = "Proportion of all genes in set",
                  mytitle = f1gsets[x, descrip])
  print(plt)
  return(NULL)
}))
invisible(dev.off())

# Changed expression as proportion of changed expressiion
#       NOTE But this is normalized to # expression changed so be careful as that may skew interpretation
## All gene sets on one page faceted [useful for seeing all at once, notebook]
pdf(file.path(p$outdir, paste0(p$baseoutname, "_inhmode_propofchangedgenes_allsets.pdf")), 12, 4)
plt<-stackedbar(pinhclass = pinhmode[propName%in%inhmodes & propDenom=="ngenes_changed", ],
                xcol = "strain", stackcol = "propName", stacknumcol = "prop", stackorder = inhmodes,
                barcols = inhcols,
                legendlabel = "Inheritance mode", myxlab = "Strain", myylab = "Proportion of changed expr. genes in set")
print(plt + facet_wrap(~geneset) + theme(axis.text.x = element_text(size = 14, angle = 90)))
invisible(dev.off())
## Gene sets individually, one plot per page [useful for presentations, zooming in]
pdf(file.path(p$outdir, paste0(p$baseoutname, "_inhmode_propofchangedgenes_indivplots.pdf")), 8, 5.5)
invisible(lapply(1:nrow(f1gsets), function(x){
  plt<-stackedbar(pinhclass = pinhmode[propName%in%inhmodes & propDenom=="ngenes_changed" & geneset ==f1gsets[x, shortname], ],
                  xcol = "strain", stackcol = "propName", stacknumcol = "prop", stackorder = inhmodes,
                  barcols = inhcols,
                  legendlabel = "Inheritance mode", myxlab = "Strain", myylab = "Proportion of changed expr. genes in set",
                  mytitle = f1gsets[x, descrip])
  print(plt)
  return(NULL)
}))
invisible(dev.off())

#### Scatter plots ####
cat("....Making F1 vs. parental strains DE scatter plots....\n")
lfcpltdir<-file.path(p$outdir, "lfcVlfcplots")
if(!dir.exists(lfcpltdir)){dir.create(lfcpltdir)}

# Make basic plots each gene set
pdf(file.path(lfcpltdir, paste0(p$baseoutname, "_lfcVlfcplots_allgenesets_plain.pdf")), 18, 12) # size is for 4 strains - not currently updated for more!
invisible(lapply(1:nrow(f1gsets), function(x){
  plts<-lapply(strains, function(str){
    inhmodescat(pdat = dat[strain==str & eval(parse(text = f1gsets[x, exprtxt]))],
                inhmodes, xcol = paste0("log2FoldChange.F1VsParent", p$refstrain), ycol = "log2FoldChange.F1VsParentNonRef",
                mytitle = str, mysubt = paste0("Genes: ", f1gsets[x, descrip]))
  })
  print(plot_grid(plotlist = plts))
}))
invisible(dev.off())

# Hyperdivergent region genes outlined [irrelevant for hypdiv-excluded geneset]
pdf(file.path(lfcpltdir, paste0(p$baseoutname, "_lfcVlfcplots_allgenesets_hypdivoutline.pdf")), 18, 12) # size is for 4 strains - not currently updated for more!
invisible(lapply(1:nrow(f1gsets), function(x){
  plts<-lapply(strains, function(str){
    inhmodescat(pdat = dat[strain==str & eval(parse(text = f1gsets[x, exprtxt]))],
                inhmodes, xcol = paste0("log2FoldChange.F1VsParent", p$refstrain), ycol = "log2FoldChange.F1VsParentNonRef",
                mytitle = str, mysubt = paste0("Genes: ", f1gsets[x, descrip]), outlineby = "hypdiv", laboutlineby = "Hyperdivervent")
  })
  print(plot_grid(plotlist = plts))
}))
invisible(dev.off())

# ASE genes outlined [irrelevant for ASE genes set]
pdf(file.path(lfcpltdir, paste0(p$baseoutname, "_lfcVlfcplots_allgenesets_aseinformoutline.pdf")), 18, 12) # size is for 4 strains - not currently updated for more!
invisible(lapply(1:nrow(f1gsets), function(x){
  plts<-lapply(strains, function(str){
    pdat = copy(dat[strain==str & eval(parse(text = f1gsets[x, exprtxt]))]) ## HARD CODING gene set info used...
    pdat[,asegene:=factor(ifelse(eval(parse(text = f1gsets[shortname=="ase2plusinfoexcls", exprtxt])), TRUE, FALSE), levels = c(TRUE, FALSE))] 
    inhmodescat(pdat = pdat,
                inhmodes, xcol = paste0("log2FoldChange.F1VsParent", p$refstrain), ycol = "log2FoldChange.F1VsParentNonRef",
                mytitle = str, mysubt = paste0("Genes: ", f1gsets[x, descrip]), outlineby = "asegene", laboutlineby = "ASE inform.")
  })
  print(plot_grid(plotlist = plts))
}))
invisible(dev.off())

#### Vs divergence plots - to be able to easily see proportions across strains, even if divergence isn't key here ####
cat("....Plotting vs. divergence....\n")
divpltdir<-file.path(p$outdir, "propvdivergenceplots")
if(!dir.exists(divpltdir)){dir.create(divpltdir)}

# Read in divergence info
ndiv<-fread(p$varsvsref)
ndiv[,`:=`(nvars.thou = nvars/1e03, nvars.perkb = nvars/(p$genomebp/1e03))] 
setkey(ndiv, strain)
setkey(pinhmode, strain)
pinhmode<-ndiv[pinhmode]

# plotting info
toplot<-nspslist[[1]][[1]]$ps[,.(propName, propDenom)] # in right order within one strain, gene set
toplot[propDenom=="ngenes", niceden:="total n genes"]
toplot[propDenom=="ngenes_changed", niceden:="n expression changed genes"]

# plots faceted by gene set (useful for full notebooking, seeing everything at once)
pdf(file.path(divpltdir, paste0(p$baseoutname, "_propvdivergence_allinhplots_allsetsfacets.pdf")), 12, 4)
invisible(
  lapply(1:nrow(toplot), function(x){
    plt<-divplot(pinhmode[propName==toplot[x, propName] & propDenom==toplot[x, propDenom],], 
                 ylabel = toplot[x, paste(propName, "of", niceden)], mytitle = toplot[x, paste(propName, "of", niceden)])
    print(plt + facet_wrap(~geneset))
    return(NULL)
  })
)
invisible(dev.off())

# Individual plots per gene set (useful for presentations, streamlining, etc) 
invisible(lapply(1:nrow(f1gsets), function(y){
  pdf(file.path(divpltdir, paste0(p$baseoutname, "_propvdivergence_allinhplots_", f1gsets[y, shortname], ".pdf")), 7, 5.5)
  invisible(
    lapply(1:nrow(toplot), function(x){
      plt<-divplot(pinhmode[propName==toplot[x, propName] & propDenom==toplot[x, propDenom] & geneset==f1gsets[y, shortname],], 
                   ylabel = toplot[x, paste(propName, "of", niceden)], mytitle = toplot[x, paste(propName, "of", niceden)])
      print(plt)
      return(NULL)
    })
  )
  invisible(dev.off())
}))

#### More specific analyses ####
cat("....Final specific analyses....\n")
# --- Reference-dominant vs. alt-dominant proportion plots ---
subandir<-file.path(p$outdir, "subanalyses")
if(!dir.exists(subandir)){dir.create(subandir)}

# plots faceted by gene set (useful for full notebooking, seeing everything at once)
pdf(file.path(subandir, paste0(p$baseoutname, "_refdominantValtdominant_allsetsfacets.pdf")), 12, 4)
print(propvprop(pinhmode) + facet_wrap(~geneset))
invisible(dev.off())

# Individual plots per gene set (useful for presentations, streamlining, etc) 
pdf(file.path(subandir, paste0(p$baseoutname, "_refdominantValtdominant_indivplots.pdf")), 7, 5.5)
invisible(lapply(1:nrow(f1gsets), function(y){
  print(propvprop(pinhmode[geneset == f1gsets[y, shortname]], mytitle = f1gsets[y, descrip]))
}))
invisible(dev.off())


#### Script completion message & session information ####
cat("....f1_parental_comparisons_withinstrain.R processing complete! Session information:....\n")
sessionInfo()