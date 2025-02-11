#! /usr/bin/env/ Rscript
# Looks at mean expression vs. ASE classifications in drosophila - data from Glaser-Schmitt et al 2024 (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011257)
# by Avery Davis Bell, begun 2025.01.14
require(data.table)
require(argparser)
require(DESeq2)
require(cqn) # for length normalization; bioconductor package
require(vsn) # for mean sd plots; bioconductor package
require(ggplot2)
require(formattable)
require(ggforce)
require(RColorBrewer)

# plotting theme (newly added 2025.01)
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11),
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14),
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

#### Functions ####
allcors<-function(x, y, info = ""){
  p<-cor.test(x, y, method = "pearson")
  s<-cor.test(x, y, method = "spearman")
  k<-cor.test(x, y, method = "kendall")
  return(data.table(info = info,
                    pearson.r = p$estimate,
                    pearson.r2 = p$estimate^2,
                    pearson.pval = p$p.value,
                    spearman.rho = s$estimate,
                    spearman.pval = s$p.value,
                    kendall.tau = k$estimate,
                    kendall.pval = k$p.value))
}

getrank<-function(val, vec, decile = T, suff = ""){
  # Get rank, quantile of val in terms of vector. ties are averaged
  # suffix added to output
  myranks<-rank(vec)
  valrank<-myranks[which(vec==val)][1]
  out<-data.table(nrank = valrank, proprank = valrank/length(vec))
  setnames(out, paste0(c("nrank", "proprank"), suff))

  if(decile==T){ # add decile
    out[, decile:=(floor(get(paste0("proprank", suff))*10) + 1)]
    out[decile>10, decile:=10]
    setnames(out, "decile", paste0("decile", suff))
  }

  return(out)
}

ginfofromtinfo<-function(tinfof){
  # Gets a gene-level estimate of transcript length, GC content from bioawk processing of FlyBase transcript FASTA
  # In: tinfof, path to input file. No header. Format (columns):
  #       transcript ID, comment - all the other info from FASTA header (includes gene ID as 'parent'), transcript length, transcript GC%
  # Out: data.table with one row per gene in tinfof. Columns:
  #       FBgn, gene ID
  #       ntranscripts, number transcripts from FlyBase for this gene
  #       length.median, median length of this gene's transcripts
  #       GC.median, median GC % (well, proportion) for this gene's transcripts

  # Read in, organize, name, etc
  tinfo<-fread(tinfof, header = F, sep = "\t")
  setnames(tinfo, c("FBtr", "Comment", "length", "GC"))
  tinfo[, FBgn:=tinfo[, tstrsplit(Comment, ";")][, tstrsplit(V8, "=")]$V2]
  setkey(tinfo, FBgn)
  tinfo<-tinfo[!is.na(FBgn)] # remove any where transcript doesn't map to a gene

  # Get per gene
  #       here is where figure out how to combine across transcripts. Just median probably OK.
  glgc<-tinfo[, .(ntranscripts = .N, length.median = median(length), GC.median = median(GC)), by = FBgn]

  return(glgc)
}

sampinfofromnames<-function(namevec){
  # Creates sample information data.table from vector of counts' names
  # In: namevec, names of samples from GEO counts matrices. Format of each is <tissue signifier>_<strain/crossed strain>_<biological replicate N>
  # Out: data.table with one row per name in initial namevec. Columns:
  #       SampleID, Sample ID as in namevec
  #       Tissue, tissue from input - first part of name
  #       Replicate, replicate from input - last part of name
  #       Generation, generation (F1 or Parent) inferred from name
  #       Strain, strain or strainxstrain from input
  #       Allele1, If an F1, first strain parent
  #       Allele2, If an F1, second strain parent

  # Split existing name nicely
  sdt<-data.table(SampleID = namevec)
  sdt<-data.table(sdt, sdt[,tstrsplit(SampleID, "_")])
  setnames(sdt, c("SampleID", "Tissue", "Strain", "Replicate"))

  # Add other info
  sdt[grep("x", Strain), Generation:="F1"] # danger, but works here with these exact data
  sdt[is.na(Generation), Generation:="Parent"]
  setkey(sdt, SampleID)
  ## Alleles for F1s, just in case
  als<-sdt[Generation=="F1", tstrsplit(Strain, "x"), by = SampleID]
  setnames(als, c("SampleID", "Allele1", "Allele2"))
  sdt<-als[sdt]
  setcolorder(sdt, c("SampleID", "Tissue", "Replicate", "Generation", "Strain"))

  # Return (in same order as input if have changed!!)
  sdt<-sdt[namevec, ]
  return(sdt)
}

cqn4deseq<-function(sinfo, ginfo, cts){
  # Uses cqn to get length (and library size) normalizations that can be used by DESeq2
  # In: sinfo, sample information in same order as samples are in expression data. Required columns:
  #               SampleID (as in expression data), nreads (total n reads - library size)
  #     ginfo, gene information data.table with one row per gene in expression data (in that order). Required columns:
  #             FBgn, gene ID (just nice to have)
  #             length.median, length to use for normalization
  #             GC.median, GC content to use for normalization
  #     counts, count matrix! rows are genes (same order as ginfo), columns are samples (same order as sinfo)
  # Out: matrix (genex x samples) of normalization factors generated by cqn, converted to normal scale for use by DESeq2

  # Get data set up as CQN likes it where needed
  sizeFs<-sinfo$nreads
  names(sizeFs)<-sinfo$SampleID

  # Run cqn
  cqnfit<-cqn(cts, lengths = ginfo$length.median, x = ginfo$GC.median,
              sizeFactors = sizeFs)

  # Get for use with DESeq2 & return
  cqnNormFactors <- exp(cqnfit$glm.offset)
  return(cqnNormFactors)
}

procasefiles<-function(exampfile, scombs, lfcthresh.asede = 0, alpha.asede = 0.05, lfcthresh.owncol
                       ){ # lfcthresh.regclass =  0.5849625, alpha.regclass = 0.05, rpci = 0.999
  # Reads in ASE data for all strain crosses, sanitizes names of columns & keeps only relevant data;
  # **based on exactly how the S3 and S4 data were in the xslx from PLoS Gen
  # In: exampfile, EXAMPLE Path to midgut ASE (subsampled) output, from Data S3. Where specific strain cross is, should have STRAIN
  #     scombs, vector of strain combinations - these are in file names and used for output
  #     lfcthresh.asede, threshold for log2FC for gene to be called significantly ASE/DE
  #     alpha.asede, threshold for adjusted p value to be called significantly ASE/DE
  #     lfcthresh.owncol, tthreshold for log2FC for gene to be called signif ASE/DE in additional column (signifAtThr...ASE.thresh)
  # Out: data.table with one row per gene-strain combination from input. ** different numbers of genes per strain **. Columns:
  #   FBgn: flybase gene ID
  #   strain: strain combination this result is for - <strain>x<strain>
  #   log2FoldChange.ParentVsParent: LFC for DE (from input)
  #   padj.ParentVsParent: padj  for DE (from input)
  #   log2FoldChange.ASE: LFC for ASE (from input)
  #   padj.ASE: padj for ASE (from input)
  #   padj.CHM: padj for test of ASE vs DE (from input)
  #   signifAtThresholds.ASE: is ASE significant at LFC and alpha thresholds given (added here)
  #   signifAtThresholds.ParentVsParent: is DE significant at LFC and alpha thresholds given (added here)
  #   signifAtThresholds.ASE.thresh, is ASE significant at alfa and lfctrhesh.owncol
  #   signifAtThresholds.ParentVsParent.thresh, is DE significant at alpha and lfcthresh.owncol
  #   regclass.orig: called regulatory class (from input)
  #   inform: obligately T - just fill in T if gene is included here (not sure about this)
  #   regclass.update: called regulatory class - uses same data as input, but uses my terms, calls any without ASE and DE as conserved, fixes a mixup with cis+trans/cisxtrans not taking into account directionality
  #   ciscompensated, for genes with cis, T or F for if compensated. F: cis and enhancing classes; T: compensatory and overcompensating classesd ('compensating' means DE lower magnitude than ASE, excluding this from this column entirely)
  #   regclass.update.comb: regclass.update with cis-trans opposing combined together: compensatory, overcompensating, and compensating are termed 'cis-trans opposing'

  # Initial read in
  fs<-sapply(scombs, gsub, pattern = "STRAIN", x = exampfile)

  dat<-rbindlist(lapply(1:length(fs), function(i){
    strn<-scombs[i]
    out<-fread(fs[[i]], header = T)[, .SD, .SDcols = -c(2:13)][!FBgn==""] # some extra space weirdly appended for some reason
    setnames(out, c("FBgn", "log2FoldChange.ParentVsParent", "padj.ParentVsParent", "log2FoldChange.ASE", "padj.ASE", "padj.CHM", "regclass.orig"))

    out[, strain:=strn]
    setcolorder(out, c("FBgn", "strain"))


    return(out)
  }))

  # Add other relevant columns: inform; signifAtThresholds...
  dat[, `:=`(inform = T, # just saying it's informative if it's in here...not exactly sure how to handle...I think this is right based on reading the paper
             signifAtThresholds.ASE = (padj.ASE < alpha.asede & abs(log2FoldChange.ASE)>=lfcthresh.asede),
             signifAtThresholds.ParentVsParent = (padj.ParentVsParent < alpha.asede & abs(log2FoldChange.ParentVsParent)>=lfcthresh.asede))]

  # Re-classify regclass using their p-values, LFCs, but fixing some issues!!
  # (I still think they're probably not calling enough as cis? )
  dat[regclass.orig=="all cis", regclass.update:="cis"]
  dat[regclass.orig=="all trans", regclass.update:="trans"]
  dat[regclass.orig=="ambiguous", regclass.update:="ambiguous"]
  dat[regclass.orig=="compensatory", regclass.update:="compensatory"]
  dat[regclass.orig=="conserved", regclass.update:="conserved"]
  dat[is.na(regclass.update) & padj.CHM< 0.05 & abs(log2FoldChange.ParentVsParent) > abs(log2FoldChange.ASE) &
        ((log2FoldChange.ParentVsParent < 0  & log2FoldChange.ASE < 0) | (log2FoldChange.ParentVsParent > 0  & log2FoldChange.ASE > 0)),
      regclass.update:="enhancing"] # enhancing looks like DE is more than ASE and in same direction
  dat[is.na(regclass.update) & padj.CHM< 0.05 & abs(log2FoldChange.ParentVsParent) < abs(log2FoldChange.ASE) &
        ((log2FoldChange.ParentVsParent < 0  & log2FoldChange.ASE < 0) | (log2FoldChange.ParentVsParent > 0  & log2FoldChange.ASE > 0)),
      regclass.update:="compensating"] # less DE but ASE in same direction === compensatory - trans compensating for but not knocking out cis
  dat[is.na(regclass.update) & padj.CHM < 0.05 & ((log2FoldChange.ParentVsParent < 0  & log2FoldChange.ASE > 0) | (log2FoldChange.ParentVsParent > 0  & log2FoldChange.ASE < 0)),
      regclass.update:="overcompensating"] # DE & ASE are in opposite directions
  ## Ambiguous called in cases where ASE & DE not called - if both ASE & DE not called, I am going to classify as conserved!!
  dat[regclass.orig=="ambiguous" & signifAtThresholds.ASE==F & signifAtThresholds.ParentVsParent==F, regclass.update:="conserved"]

  # Add cis compensated info from these new categories (of ones with clear cis effects)
  dat[regclass.update%in%c("cis", "enhancing"), ciscompensated:=F]
  dat[regclass.update%in%c("compensatory", "overcompensating"), ciscompensated:=T]

  # Add reg class with umbrella cis-trans opposing
  dat[, regclass.update.comb:=regclass.update]
  dat[regclass.update%in%c("compensatory", "overcompensating", "compensating"), regclass.update.comb:="cis-trans opposing"]

  # ASE, DE with log2FC thresholds. NA out ones that are signif but below threshold
  dat[, signifAtThresholds.ASE.thresh:=signifAtThresholds.ASE]
  dat[abs(log2FoldChange.ASE)<lfcthresh.owncol & signifAtThresholds.ASE==T, signifAtThresholds.ASE.thresh:=NA]
  dat[, signifAtThresholds.ParentVsParent.thresh:=signifAtThresholds.ParentVsParent]
  dat[abs(log2FoldChange.ParentVsParent)<lfcthresh.owncol & signifAtThresholds.ParentVsParent==T, signifAtThresholds.ParentVsParent.thresh:=NA]

  ### other thoughts....
  # ALSO probably try with LFC threshold!!! like, a column that puts ambiguous with lower FCs as 1.5 or something....
  #     (for reg class; maybe expand conserve)
  #
  # # ***ASE vs not where ASE = p < 0.05 & FC > 1.5, no ASE = p > 0.05 [& FC lower? not sure] ***. Maybe DE, too.
  #     JUST ADD AS MANY COLUMNS AS I MIGHT NEED HERE, INCLUDING COMP VS NOT ETC, AND DOCUMENT THEM HERE

  setcolorder(dat, c("signifAtThresholds.ASE", "signifAtThresholds.ParentVsParent", "signifAtThresholds.ASE.thresh" , "signifAtThresholds.ParentVsParent.thresh"),
              before = "regclass.orig")
  return(dat)
}

# --- stats, plotting etc from my earlier scripts
anovatuk<-function(dat, mystrain, testinforow, testagainstcol = "testDat", colorvec){
  # Runs an ANOVA of category vs. test data (all in all date) for one strain, one category
  # In: dat - data with columns strain, any specified in multiwaytests (columnname, mynarrow), testagainstcol
  #     mystrain - strain to narrow to for this test
  #     testinforow: one row data.table with columns myname (long format name), shortname (short format name),
  #             columnname (column CATEGORY data is in alldat), mynarrow (text logical expression of how to further narrow data for this test, e.g. 'inform==T'),
  #     testagisntcol - column name in alldat that's the y axis/continuous data to test against
  #     colorvec - named in order vector of all category values; used for ordering results
  # Out: list of data.tables -
  #       anout - ANOVA results. One row data.table with columns
  #               Df_resid, DF residuals from ANOVA out
  #               SumSq_resid, sum squares residuals from ANOVA out
  #               MeanSq_resid, mean squares residuals from ANOVA out
  #               Df_category,  DF for test category from ANOVA out
  #               SumSq_category, sum squares for test category from ANOVA out
  #               MeanSq_category, mean squares for test category from ANOVA out
  #               Fvalue, ANOVA F value
  #               pvalue, ANOVA p value
  #       tukout - Tukey's HSD results. One row per inter-category comparison. Columns:
  #               comparison, cat1-cat2 detail of comparison : default way tukey output displays
  #               cat1, category this test is in (vs cat 2)
  #               cat2, category this test is in (vs cat1)
  #               diff, Tukey output
  #               lwr, Tukey output
  #               upr,  Tukey output
  #               tukey.padj, Tukey output (adjusted p)
  #       tuklabs - data.table with columns <column name in testinfo row>, label - < or > if this category is p < 0.05 different from FIRST category in this datat
  #           (for plotting)
  #       ns - data.table with columns categorylabel, n: number of observations in each category (that are non-NA!!)

  # Narrow data to that of interest
  thisdat<-dat[strain==mystrain & eval(parse(text = testinforow[, mynarrow])), ]
  yvals<-thisdat[, get(testagainstcol)] # naming for prettier stats calls, returns
  catvals<-factor(thisdat[, get(testinforow[, columnname])], # naming for prettier stats calls, returns
                  levels = names(colorvec)) # leveling so comparisons ordered as sensibly as automatically possible
  # ns
  ns<-data.table(as.matrix(data.table(catvals, yvals)[!is.na(yvals), table(catvals)]), keep.rownames = T)
  setnames(ns, c("categorylabel", "n"))

  # ANOVA
  myan<-anova(lm(yvals ~ catvals))
  anout<-data.table( Df_resid = myan$Df[2], SumSq_resid = myan$`Sum Sq`[2], MeanSq_resid = myan$`Mean Sq`[2],
                     Df_category = myan$Df[1], SumSq_category = myan$`Sum Sq`[1],
                     MeanSq_category = myan$`Mean Sq`[1], Fvalue = myan$`F value`[1], pvalue = myan$`Pr(>F)`[1])
  # one way so just saving in one row

  # Tukey's HSD
  mytuk<-TukeyHSD(aov(yvals ~ catvals)) # NOTE my design likely not balanced
  ## get categories separate: UPDATED way that doesn't rely on a given character not being in names. **Confirmed that given that levels = names(colorvec), this is true
  ##   also updated to work if not all categories are actually in data...
  myusedlevs<-sort(factor(names(table(catvals)[table(catvals)!=0]), levels = names(colorvec)))
  mycomps<-unlist(
    lapply(1:(length(myusedlevs)-1), function(i){
      lapply((i+1):length(myusedlevs), function(j){
        c(myusedlevs[j], myusedlevs[i])
      })
    }),
    recursive = F)

  tukout<-data.table(comparison = rownames(mytuk[[1]]), # maybeeee add column with rownames from tuk so if something's weird, can confirm
                     cat1 = sapply(mycomps, function(x) x[1]),
                     cat2 = sapply(mycomps, function(x) x[2]),
                     mytuk[[1]])
  setnames(tukout, "p adj", "tukey.padj")

  # LABELS based on ANOVA and TUKEY results for comparison vs. first (reference level) - not for saving, for plotting
  ## asterisk if p < 0.05; next line is '>' if mean is larger than reference category, '<' if less than reference category
  mycats<-names(colorvec)
  if(anout$pvalue<0.05){ # only make non-blank labels/look at comparisons if ANOVA is nominally significant
    tuklabs<-tukout[cat2==mycats[1], .(cat1, tukey.padj, diff)]
    tuklabs[tukey.padj >= 0.05, label:= ""]
    tuklabs[tukey.padj < 0.05 & diff < 0, label := "<"]
    tuklabs[tukey.padj < 0.05 & diff > 0, label := ">"]
    # tuklabs[,label := ifelse(tukey.padj<0.05, "*", "")]
    tuklabs<-tuklabs[,.(cat1, label)]
    ## include any not observed categories (including first level itself!) and order appropriately
    tuklabs<-rbind(data.table(cat1 = mycats[!mycats%in%tuklabs$cat1], label = ""),
                   tuklabs)
    # order in order of input vector
    tuklabs<-tuklabs[order(match(cat1, mycats))]
  }else{ # all are blank/no asterisk
    tuklabs<-data.table(cat1 = names(colorvec),
                        label = "")
  }
  setnames(tuklabs, "cat1",testinforow[, columnname]) # Name so it's same column name as in overall data

  # Return
  return(list(anout = anout, tukout = tukout, tuklabs = tuklabs, ns = ns))
}

sinawmed<-function(datin, xcol = "signifAtThresholds.ASE", ycol = "pSegSites", facrow = "strain", faccol = "sites",
                   colorcol, colorvec, boxcol = rgb(0, 0, 1, 0.6), faclabels = F,
                   myxlab = "ASE called", myylab = "Proportion segregating sites", mytitle = "", mysubt = "", myscales = "fixed",
                   facvec = ""){
  # Modified from many previous scripts...
  # Makes a faceted sina plot with boxplot (IQR only!) overlaid
  # In: datin, data, with columns titled by values of xcol, ycol, facrow, and faccol
  #     xcol, character name of column for x axis/to split data by
  #     ycol, character name of column for y axis/of values
  #     facrow, charactter name of column whose values will be row-wise facets; if same as faccol, will be facet wrapped
  #     faccol, charactter name of column whose values will be col-wise facets; if same as facrow, will be facet rapped
  #     colorcol, name of column containing values to color by (factors)
  #     colorvec, named vector specifying colors for each value of colorcol in order they should be on plot. Names are values, values are colors
  #     boxcol, color for boxplot outline. Recommend including transparency IN this color
  #     faclabels, optional (put F to exclude) data.table specifying labels for each facet, i.e. p-values or the like.
  #         must have columns facrow value, faccol value, label (label has actual text to include)
  #     myxlab, x axis label
  #     myylab, y axis label
  #     mytitle, plot title
  #     mysubt, subtitle
  #     myscales, passed to facetting scales=
  #     facvec, vector of values facet take to ORDER FACCROW BY (i.e., makes this a factor leveled this way)

  # Make sure stuff is in right order
  pdat<-copy(datin)
  ## try if factors are here - they are!
  pdat[, mycol:=factor(get(colorcol), levels = names(colorvec))]
  if(colorcol == xcol){ # factor/order xcol as above
    pdat[, myxcol:=factor(get(xcol), levels = names(colorvec))]
  }else{
    pdat[, myxcol:=get(xcol)]
  }

  # Plot
  plt<-ggplot(pdat, aes(myxcol, eval(as.name(ycol)))) +
    ggforce::geom_sina(aes(color = mycol), alpha = 0.3, size = 0.2) +
    geom_boxplot(outlier.color = NA, alpha = 0, col = boxcol, coef = 0, lwd = 0.2) + # this boxplot shows just median and IQR
    scale_color_manual(values = colorvec) +
    xlab(myxlab) + ylab(myylab) + ggtitle(mytitle, subtitle = mysubt) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14),
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.position = "none",
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
                       strip.text.y = element_text(size = 14))
  # facets
  if(facrow==faccol){
    if(facvec[1]!=""){ # need to add leveling in here for it to propagate through
      plt<-plt + facet_wrap(~factor(eval(as.name(facrow)), levels = facvec), scales = myscales)
    }else{
      plt<-plt + facet_wrap(~eval(as.name(facrow)), scales = myscales)
    }
  }else{
    plt<-plt + facet_grid(eval(as.name(facrow))~eval(as.name(faccol)), scales = myscales)
  }

  # Add labels if provided
  if(is.data.table(faclabels)){ # thanks to https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
    plt<-plt + geom_text(
      size    = 3,
      data    = faclabels,
      mapping = aes(x = Inf, y = Inf, label = label),
      hjust   = 1.05,
      vjust   = 1.5
    )
  }

  # Other code
  # stat_summary(fun = "median", geom = "point", col = "red") +# this auto-computes and adds the medians! but couldn't figure out making it horizontal bar easily
  # geom_point(stat = "summary", fun = median, color = "red") + # this auto-computes and adds the medians! but couldn't figure out making it horizontal bar easily


  # Deal with if want scales to be free or not, or to do both
  # May want to re-order strains, sites - provide something with it factor-ized if so
  return(plt)
}

#### Arguments & Inputs ####
p<-arg_parser("Looks at mean expression vs. ASE classifications in drosophila - data from Glaser-Schmitt et al 2024 (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011257)",
              name = "glaserschmitt_drosase_meanexprvsaseetc.R", hide.opts = TRUE)

# Inputs - gene information
p<-add_argument(p, "--transcriptinfo",
                help = "Path to file with transcript name, comment, length, GC% - from using bioawk on transcript fasta",
                type = "character",
                default = "FlyBase_transcript_LengthGC_20220115.txt.gz")

# Inputs - sample and expression-level related (from GEO)
p<-add_argument(p, "--gcountsMG",
                help = "Path to midgut gene counts for each sample (GSE263264_MG_gcounts_allgenotypes_allreads.tsv.gz from GEO)",
                type = "character",
                default = "GSE263264_MG_gcounts_allgenotypes_allreads.tsv.gz")
p<-add_argument(p, "--gcountsHG",
                help = "Path to hindgut gene counts for each sample (GSE263264_HG_gcounts_allgenotypes_allreads.tsv.gz from GEO)",
                type = "character",
                default = "GSE263264_HG_gcounts_allgenotypes_allreads.tsv.gz")

# Inputs - ASE/analysis-related (from PLoS Genetics)
p<-add_argument(p, "--aseMG",
                help = "EXAMPLE Path to midgut ASE (subsampled) output, from Data S3. Where specific strain cross is, should have STRAIN",
                type = "character",
                default = "s3data_ASESTRAIN.txt")
p<-add_argument(p, "--aseHG",
                help = "EXAMPLE Path to hindgut ASE (subsampled) output, from Data S4. Where specific strain cross is, should have STRAIN",
                type = "character",
                default = "s4data_ASESTRAIN.txt")


# Outputs
p<-add_argument(p, "--outdir",
                help = "Output directory",
                type = "character",
                default = "out")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "drosase_glaserschmitt2024")

p<-parse_args(p)

#### Get in and format gene data ####
cat(".....Reading in and format relevant gene data (lengths, GC content)....\n")
glgc<-ginfofromtinfo(p$transcriptinfo)

#### Get in and format expression data ####
cat(".....Reading in and formatting expression level data....\n")
# --- Read in counts: COMBINING THEM (have checked they're in same order)
mg.cts<-fread(p$gcountsMG, header = T)
hg.cts<-fread(p$gcountsHG, header = T)

print("Gene info is identical between MG and HG gene counts files?")
all(mg.cts[, .(FBgn, CGnum, Chrom, Length)]==hg.cts[, .(FBgn, CGnum, Chrom, Length)]) # TRUE, means can just retain gene info on one

## Extract gene info
ginfo.cts<-mg.cts[, .(FBgn, CGnum, Chrom, Length)]
mg.cts<-mg.cts[, .SD, .SDcols = -c(1:4)]
hg.cts<-hg.cts[, .SD, .SDcols = -c(1:4)]

## Get sample information from counts headers
print("Sample info is identical between MG and HG gene counts files?")
all(names(mg.cts)==names(hg.cts)) # FALSE - because they have the data collection site in names

sampinfo.mg<-sampinfofromnames(names(mg.cts)) # make sure put sampleID as rows too
sampinfo.hg<-sampinfofromnames(names(hg.cts)) # make sure put sampleID as rows too - LATER when use

## Add in library size TOTALS as sample info (before narrowing to higher-expressed genes)
sampinfo.mg[, nreads:=unlist(mg.cts[, lapply(.SD, sum)])]
sampinfo.hg[, nreads:=unlist(hg.cts[, lapply(.SD, sum)])]

## Parent-only sample info, likely what I'll use downstream
sinfo.mg.paren<-sampinfo.mg[Generation=="Parent", ]
sinfo.hg.paren<-sampinfo.hg[Generation=="Parent", ]

# --- Get normalized data [for PARENTS]: step one, restrict to relevant genes
# First, restrict to genes that they included in any analysis: those with >= 15 reads in EACH sample
#       AND, for my purposes, only care about parents
## Find genes to keep: >= 15 reads AND with length, GC content info
gkeep.mg<-(rowMins(as.matrix(mg.cts))>=15 & ginfo.cts[, FBgn%in%glgc$FBgn] )
ginfo.mg<-data.table(ginfo.cts[gkeep.mg ,], glgc[ginfo.cts[gkeep.mg , FBgn], .(ntranscripts, length.median, GC.median)]) # add in length & GC content here, too

gkeep.hg<-(rowMins(as.matrix(hg.cts))>=15 & ginfo.cts[, FBgn%in%glgc$FBgn])
ginfo.hg<-data.table(ginfo.cts[gkeep.hg ,], glgc[ginfo.cts[gkeep.hg , FBgn], .(ntranscripts, length.median, GC.median)]) # add in length & GC content here, too

## Keep only these genes for parents; convert counts to matrices
mg.cts<-as.matrix(mg.cts[gkeep.mg, sinfo.mg.paren$SampleID, with = F]) # now lower count ones are DROPPED
rownames(mg.cts)<-ginfo.mg$FBgn
hg.cts<-as.matrix(hg.cts[gkeep.hg, sinfo.hg.paren$SampleID, with = F]) # now lower count ones are DROPPED
rownames(hg.cts)<-ginfo.hg$FBgn

## Print out number genes retained, dropped, etc
cat(paste("\t\t\t MIDGUT:", nrow(mg.cts), "genes retained for analysis of", nrow(ginfo.cts), "original genes. (those with 15+ reads in each sample and FlyBase transcript info).\n"))
cat(paste("\t\t\t HINDGUT:", nrow(hg.cts), "genes retained for analysis of", nrow(ginfo.cts), "original genes. (those with 15+ reads in each sample and FlyBase transcript info).\n"))

# --- Get normalized parent data
# Get normalization factors taking into account length, gc (using cqn)
norm.mg<-cqn4deseq(sinfo = sinfo.mg.paren, ginfo = ginfo.mg, cts = mg.cts)
norm.hg<-cqn4deseq(sinfo = sinfo.hg.paren, ginfo = ginfo.hg, cts = hg.cts)

# Sample info in DESeq2 format
colinfo.mg<-data.frame(sinfo.mg.paren)
rownames(colinfo.mg)<-sinfo.mg.paren$SampleID
colinfo.hg<-data.frame(sinfo.hg.paren)
rownames(colinfo.hg)<-sinfo.hg.paren$SampleID

# DESeq2 set up
dds.mg.parents<-DESeqDataSetFromMatrix(
  countData = mg.cts,
  colData = colinfo.mg,
  design = ~ Replicate + Strain # doesn't super matter
)
dds.hg.parents<-DESeqDataSetFromMatrix(
  countData = hg.cts,
  colData = colinfo.hg,
  design = ~ Replicate + Strain # doesn't super matter
)
## Add in normalization factors
normalizationFactors(dds.mg.parents)<-norm.mg
normalizationFactors(dds.hg.parents)<-norm.hg

# DESeq2 vst normalization
vsd.mg<-vst(dds.mg.parents, blind = T)
vsd.hg<-vst(dds.hg.parents, blind = T)

# DESeq2 shifted log transform (log2(cts norm + 1))
ntd.mg <- normTransform(dds.mg.parents)
ntd.hg<-normTransform(dds.hg.parents)

#### make sure norm factors get used here! make sure don't need to run something first: THEY DO

# --- Get mean, variance, SD per gene OVERALL for quick tests (& in case want to use)
all.mean.mg<-data.table(ginfo.mg,
                        mean.vst = rowMeans(assay(vsd.mg)),
                        median.vst = rowMedians(assay(vsd.mg)),
                        var.vst = rowVars(assay(vsd.mg)),
                        sd.vst = rowSds(assay(vsd.mg)),
                        mean.shlog = rowMeans(assay(ntd.mg)),
                        median.shlog = rowMedians(assay(ntd.mg)),
                        var.shlog = rowVars(assay(ntd.mg)),
                        sd.shlog = rowSds(assay(ntd.mg)))
setkey(all.mean.mg, FBgn)
all.mean.mg<-all.mean.mg[all.mean.mg[, getrank(sd.vst, all.mean.mg$sd.vst, suff = ".sd.vst"), by = FBgn]][all.mean.mg[, getrank(sd.shlog, all.mean.mg$sd.shlog, suff = ".sd.shlog"), by = FBgn]]
    # can add more here if want more ranks - not sure that's needed

all.mean.hg<-data.table(ginfo.hg,
                        mean.vst = rowMeans(assay(vsd.hg)),
                        median.vst = rowMedians(assay(vsd.hg)),
                        var.vst = rowVars(assay(vsd.hg)),
                        sd.vst = rowSds(assay(vsd.hg)),
                        mean.shlog = rowMeans(assay(ntd.hg)),
                        median.shlog = rowMedians(assay(ntd.hg)),
                        var.shlog = rowVars(assay(ntd.hg)),
                        sd.shlog = rowSds(assay(ntd.hg)))
setkey(all.mean.hg, FBgn)
all.mean.hg<-all.mean.hg[all.mean.hg[, getrank(sd.vst, all.mean.hg$sd.vst, suff = ".sd.vst"), by = FBgn]][all.mean.hg[, getrank(sd.shlog, all.mean.hg$sd.shlog, suff = ".sd.shlog"), by = FBgn]]

# DIRECTORY for any quality check plots I may save
qcdir<-file.path(p$outdir, "qualitychecks")
if(!dir.exists(qcdir)){dir.create(qcdir, recursive = T)}

# look at mean-var relationship?  why not
## BEFORE vst (just log'ing)

ma.ntd.mg<-vsn::meanSdPlot(assay(ntd.mg))$gg + ggtitle("Midgut parental data shifted log") + myggtheme # *** come back here - why is this looking better??
ma.ntd.hg<-vsn::meanSdPlot(assay(ntd.hg))$gg + ggtitle("Hindgut parental data shifted log") + myggtheme # *** come back here - why is this looking better??

## AFTER vst
ma.mg<-vsn::meanSdPlot(assay(vsd.mg))$gg + ggtitle("Midgut parental data post-VST") + myggtheme
ma.hg<-vsn::meanSdPlot(assay(vsd.hg))$gg + ggtitle("Hindgut parental data post-VST") + myggtheme
## Save all
pdf(file.path(qcdir, paste0(p$baseoutname, "_parentalprepostVST_MAPlots.pdf")), 7, 6)
print(ma.ntd.mg)
print(ma.ntd.hg)
print(ma.mg)
print(ma.hg)
invisible(dev.off())

# see if can see length not an issue??
mvl.mg<-ggplot(all.mean.mg, aes(mean.vst, length.median)) + geom_hex() + geom_smooth(color = "red") +
  xlab("Mean normalized expression (vst)") + ylab("Median transcript length") +
  ggtitle("Midgut parental data post-VST") + myggtheme
mvl.hg<-ggplot(all.mean.hg, aes(mean.vst, length.median)) + geom_hex() + geom_smooth(color = "red") +
  xlab("Mean normalized expression (vst)") + ylab("Median transcript length") +
  ggtitle("Hindgut parental data post-VST") + myggtheme

pdf(file.path(qcdir, paste0(p$baseoutname, "_parentalVST_MeanvLengthPlots.pdf")), 7, 6)
print(mvl.mg)
print(mvl.hg)
invisible(dev.off())

# --- Get per-strain-pair averages to use downstream
    # in allele1 and Allele2 of the F1s
alcombos.mg<-unique(sampinfo.mg[Generation=="F1", .(Allele1, Allele2)])
expr2parens.mg<-rbindlist(lapply(1:nrow(alcombos.mg), function(strnind){
  samps<-rownames(colData(vsd.mg))[colData(vsd.mg)$Strain%in%c(alcombos.mg[strnind, Allele1], alcombos.mg[strnind, Allele2])]
  out<-data.table(strain = alcombos.mg[strnind, paste(Allele1, Allele2, sep = "x")],
                  ginfo.mg,
                  meanexp.vst = rowMeans(assay(vsd.mg)[, samps]),
                  meanexp.shlog = rowMeans(assay(ntd.mg)[, samps]))
  return(out)
}))

alcombos.hg<-unique(sampinfo.hg[Generation=="F1", .(Allele1, Allele2)])
expr2parens.hg<-rbindlist(lapply(1:nrow(alcombos.hg), function(strnind){
  samps<-rownames(colData(vsd.hg))[colData(vsd.hg)$Strain%in%c(alcombos.hg[strnind, Allele1], alcombos.hg[strnind, Allele2])]
  out<-data.table(strain = alcombos.hg[strnind, paste(Allele1, Allele2, sep = "x")],
                  ginfo.hg,
                  meanexp.vst = rowMeans(assay(vsd.hg)[, samps]),
                  meanexp.shlog = rowMeans(assay(ntd.hg)[, samps]))
  return(out)
}))

#### New: look at gross expression value vs variance patterns ####
cat(".....Aside: looking at gross expression valule vs. variance patterns....\n")
# DIRECTORY for this
mvdir<-file.path(p$outdir, "parentmeanvar")
if(!dir.exists(mvdir)){dir.create(mvdir, recursive = T)}

# --- Get data into long format (tissue, expression normalization)
allexp<-rbind(data.table(tissue = "Midgut", all.mean.mg),
              data.table(tissue = "Hindgut", all.mean.hg))
allexp.vst<-allexp[, .(tissue, FBgn, CGnum, Chrom, Length, ntranscripts, length.median, GC.median,
                       mean.vst, median.vst, var.vst, sd.vst, nrank.sd.vst, proprank.sd.vst, decile.sd.vst)]
allexp.shlog<-allexp[, .(tissue, FBgn, CGnum, Chrom, Length, ntranscripts, length.median, GC.median,
                       mean.shlog, median.shlog, var.shlog, sd.shlog, nrank.sd.shlog, proprank.sd.shlog, decile.sd.shlog)]
setnames(allexp.vst, c("mean.vst", "median.vst", "var.vst", "sd.vst", "nrank.sd.vst", "proprank.sd.vst", "decile.sd.vst"),
         c("mean.exp", "median.exp", "var.exp", "sd.exp", "nrank.sd.exp", "proprank.sd.exp", "decile.sd.exp"))
setnames(allexp.shlog, c("mean.shlog", "median.shlog", "var.shlog", "sd.shlog", "nrank.sd.shlog", "proprank.sd.shlog", "decile.sd.shlog"),
         c("mean.exp", "median.exp", "var.exp", "sd.exp", "nrank.sd.exp", "proprank.sd.exp", "decile.sd.exp"))

allexp.long<-rbind(data.table(expnormmethod = "vst", allexp.vst),
                   data.table(expnormmethod = "shlog", allexp.shlog))
rm(allexp.vst, allexp.shlog) # clean up

# Save this
write.table(allexp.long, gzfile(file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# --- compare & contrast the two metrics?? just for kicks
normcomp<-ggplot(allexp, aes(mean.vst, mean.shlog)) + geom_hex() + myggtheme +
  ggtitle("Comparing mean with VST vs. log2 normalizations")
normcomp.rank<-ggplot(allexp, aes(rank(mean.vst), rank(mean.shlog))) + geom_hex() + myggtheme +
  ggtitle("Comparing mean with VST vs. log2 normalizations - ranks")
normcomp.sd<-ggplot(allexp, aes(sd.vst, sd.shlog)) + geom_hex() + myggtheme +
  ggtitle("Comparing SD expression with VST vs. log2 normalizations")
normcomp.sd.rank<-ggplot(allexp, aes(nrank.sd.vst, nrank.sd.shlog)) + geom_hex() + myggtheme +
  ggtitle("Comparing SD expression with VST vs. log2 normalizations - ranks")

pdf(file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_plots.pdf")),
              7, 6)
print(normcomp)
print(normcomp.rank)
print(normcomp.sd)
print(normcomp.sd.rank)
invisible(dev.off())

# --- General correlations (plots, tests - abs & ranks)
setkey(allexp.long, expnormmethod, tissue)
allexp.corrs<-allexp.long[, allcors(x = mean.exp, y = sd.exp, info = "Mean vs SD"), by = .(expnormmethod, tissue)]
write.table(allexp.corrs, file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_corrs.txt")),
            sep = "\t", quote = F, row.names = F)

# --- Decile explorations
deccol<-colorRampPalette(c("lightblue", "darkblue"))(10)
names(deccol)<-1:10

dectestinfo<-data.table(myname = c("Mean expression"),
                        shortname = c("exp"),
                        columnname = c("decile.sd.exp"),
                        mynarrow = rep("!is.na(FBgn)", 1))

tissues<-c("Midgut", "Hindgut")
expnorms<-c("vst", "shlog")
anres.dec<-lapply(tissues, function(tis){
  lapply(expnorms, function(myexp){
    thisdat<-allexp.long[tissue==tis & expnormmethod==myexp, ]
    thisdat[, strain:=""] # just to not have to change function

    out<-anovatuk(dat = thisdat, mystrain = "", testinforow = dectestinfo, testagainstcol = "mean.exp",
                  colorvec = deccol)

    out$anout[, `:=`(tissue = tis, expnormmethod = myexp)]
    out$tukout[, `:=`(tissue = tis, expnormmethod = myexp)]
    out$tuklabs[, `:=`(tissue = tis, expnormmethod = myexp)]
    out$ns[, `:=`(tissue = tis, expnormmethod = myexp)]

    return(out)
  })
})

# COMBINE & save results
andec<-rbindlist(lapply(anres.dec, function(x) rbindlist(lapply(x, function(y) y$anout))))
setcolorder(andec, c("tissue", "expnormmethod"))
tukdec<-rbindlist(lapply(anres.dec, function(x) rbindlist(lapply(x, function(y) y$tukout))))
setcolorder(tukdec, c("tissue", "expnormmethod"))
tuklabdec<-rbindlist(lapply(anres.dec, function(x) rbindlist(lapply(x, function(y) y$tuklabs))))
setcolorder(tuklabdec, c("tissue", "expnormmethod"))
ndec<-rbindlist(lapply(anres.dec, function(x) rbindlist(lapply(x, function(y) y$ns))))
setcolorder(ndec, c("tissue", "expnormmethod"))
## save
write.table(andec, file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_SDDeciles_ANOVA.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(tukdec, file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_SDDeciles_Tukey.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(ndec, file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_SDDeciles_Ns.txt")),
            sep = "\t", quote = F, row.names = F)

#  Sina decile plots
andec[, label:=paste("ANOVA p =", ifelse(pvalue < 0.001, paste(formattable::scientific(pvalue, digits = 1)),
                                         round(pvalue, digits = 3)))]
andec.sina<-sinawmed(datin = allexp.long,
         xcol = "decile.sd.exp", ycol = "mean.exp", facrow = "expnormmethod",
         faccol = "tissue", colorcol = "decile.sd.exp",
         colorvec = deccol,
         boxcol = rgb(0, 0, 0, 0.6), faclabels = andec, myxlab = "SD of expression",
         myylab = "Normalized mean expression",  mytitle = "Parent strains' expression across samples",
         myscales = "free", facvec = c("vst", "shlog")) +

pdf(file.path(mvdir, paste0(p$baseoutname, "_overallparentalmeanvarexp_multiplenormalizations_SDDecilesSina.pdf")),
      7, 6)
print(andec.sina)
invisible(dev.off())

#### Get in & format ASE-related data ####
cat(".....Reading in and formatting ASE/analysis data....\n")
# --- Read in ASE data
scombs<-alcombos.hg[, paste(Allele1, Allele2, sep = "x")]
ase.mg<-procasefiles(exampfile = p$aseMG, scombs = scombs, lfcthresh.asede = 0, alpha.asede = 0.05, lfcthresh.owncol = log2(1.5)) # currently these threshs are from paper
ase.hg<-procasefiles(exampfile = p$aseHG, scombs = scombs, lfcthresh.asede = 0, alpha.asede = 0.05, lfcthresh.owncol = log2(1.5))

ase.mg[, tissue:="Midgut"]
setcolorder(ase.mg, c("strain", "tissue", "FBgn"))
ase.hg[, tissue:="Hindgut"]
setcolorder(ase.hg, c("strain", "tissue", "FBgn"))

ase<-rbind(ase.mg, ase.hg)
setkey(ase, strain, tissue, FBgn)
rm(ase.mg, ase.hg) # clean up

# Quick PLOT to help me figure out their definitions
devase<-ggplot(ase[!is.na(regclass.orig)], aes(log2FoldChange.ParentVsParent, log2FoldChange.ASE)) +
  geom_point(aes(color = regclass.orig, fill = regclass.orig), alpha = 0.4) +
  xlab("DE (log2FC)") + ylab("ASE (log2FC)") + ggtitle("Regulatory classes from paper") + myggtheme +
  facet_grid(strain~tissue)
pdf(file.path(qcdir, paste0(p$baseoutname, "_regclass_scatters.pdf")), 8, 10)
print(devase)
invisible(dev.off())

## same plot, but for my updates of their definitions
devase.new<-ggplot(ase[!is.na(regclass.update)], aes(log2FoldChange.ParentVsParent, log2FoldChange.ASE)) +
  geom_point(aes(color = regclass.update, fill = regclass.update), alpha = 0.4) +
  xlab("DE (log2FC)") + ylab("ASE (log2FC)") + ggtitle("Updated regulatory classes") + myggtheme +
  facet_grid(strain~tissue)
pdf(file.path(qcdir, paste0(p$baseoutname, "_regclassNEW_scatters.pdf")), 8, 10)
print(devase.new)
invisible(dev.off())

# --- Add in expression level value! (key by strain and gene); combine tissues
# Format expr level
expr2parens.mg[, tissue:="Midgut"]
setcolorder(expr2parens.mg,  c("strain", "tissue", "FBgn"))
expr2parens.hg[, tissue:="Hindgut"]
setcolorder(expr2parens.hg,  c("strain", "tissue", "FBgn"))
expr2parens<-rbind(expr2parens.mg, expr2parens.hg)
setkey(expr2parens, strain, tissue, FBgn)
rm(expr2parens.mg, expr2parens.hg) # clean up

# Combine and clean up
dat<-ase[expr2parens]
rm(expr2parens, ase) # clean up

# Any not inform = T are inform = F
dat[is.na(inform), inform:=F]

# --- Add 'all' combined strain (within-tissue currently): genes with results in all strains (so each gene represented multiple times)
# Get N strains for each gene
npres<-dat[, sum(!is.na(log2FoldChange.ASE)), by = .(FBgn, tissue)]
setnames(npres, "V1", "NTestedCrosses")
npres[, testedInAll:=ifelse(NTestedCrosses==length(scombs), T, F)]

setkey(dat, FBgn, tissue)
setkey(npres, FBgn, tissue)

dat<-npres[dat]
setcolorder(dat, c("NTestedCrosses", "testedInAll"), after = "meanexp.vst")
setcolorder(dat, "strain")
# get 'all' subset
allc.dat<-dat[testedInAll==T,]
allc.dat[, strain:="all"]
# add in
dat.wall<-rbind(dat, allc.dat)
# Strain list in desired order, with 'all'
combswall<-c("all", scombs)


# ---- SAVE long-format data that has everyone, expression data (means) & ASE data
write.table(dat.wall, gzfile(file.path(p$outdir, paste0(p$baseoutname, "data_inclcombostrain_ExprAndASE.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# Inform - not in input doesn't mean not informative, necessarily; could have been excluded at any point...
#       closer read of paper to see how they defined, see if can find other data....

#### Run tests of expression vs. ASE etc ####
cat(".....Performing analyses....\n")
# --- colors, categories, etc
## Reg pattern, COLORS
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")
regcols<-rep(NA, length(regpats))
names(regcols)<-regpats
regcols["conserved"]<-"lightgray"
regcols["ambiguous"]<-"darkgray"
regcols[c("compensating", "compensatory", "overcompensating")]<-brewer.pal(9, "Blues")[c(3, 6, 9)] # wider dynamic range of color than doing brewer.pal(3, "Blues")
regcols[is.na(regcols)]<-brewer.pal(4, "Set2")[c(1, 2, 4)]
## simplified
cistransop<-c("compensating", "compensatory", "overcompensating")
combinedrps<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "cis-trans opposing")
combinedrpcols<-c(na.omit(regcols[combinedrps]), "cis-trans opposing" = "#4292C6")
## Original Glaser-Schmitt
gsregs<-c("conserved", "ambiguous", "all cis", "cis+trans", "all trans", "cisxtrans", "compensatory")
gsregcols<-rep(NA, length(gsregs))
names(gsregcols)<-gsregs
gsregcols["conserved"]<-"lightgray"
gsregcols["ambiguous"]<-"darkgray"
gsregcols[("compensatory")]<-"#4292C6"
gsregcols[is.na(gsregcols)]<-brewer.pal(4, "Set2")[c(1, 2, 4, 3)]
## changed only categories (all and simplified)
regcols.ch<-regcols[!names(regcols)%in%c("conserved", "ambiguous")]
combinedrpcols.ch<-combinedrpcols[!names(combinedrpcols)%in%c("conserved", "ambiguous")]

# T/F colors
asecols<-c("gray", "red") # ASE yes/no
names(asecols)<-c(FALSE, TRUE)

# --- set up tests to do [NEED TO FIX FOR THESE DATA]
# Two way comparison - ASE vs. not - style comparisons
twowaytests<-data.table(myname = c("Used for ASE/DE testing vs not", "ASE vs not", "DE vs not", "ASE vs not - FC > 1.5", "DE vs not - FC > 1.5",
                                   "Cis-regulatory change: compensated in trans or not"),
                        shortname = c("inform", "ASE", "DE", "ASE_FC", "DE_FC", "ciscompensated"), # for file etc out
                        myxname = c("Used for ASE/DE testing", "ASE called", "DE called", "ASE called (FC > 1.5)", "DE called (FC > 1.5)",
                                    "Cis change is compensated"), # for plotting
                        columnname = c("inform", "signifAtThresholds.ASE", "signifAtThresholds.ParentVsParent","signifAtThresholds.ASE.thresh", "signifAtThresholds.ParentVsParent.thresh",
                                       "ciscompensated"), # data test: column with T/Fs for this test
                        mynarrow = c("inform==T | inform==F", "inform==T & !is.na(signifAtThresholds.ASE)", "inform==T & !is.na(signifAtThresholds.ParentVsParent)",
                                     "inform==T&!is.na(signifAtThresholds.ASE.thresh)", "inform==T&!is.na(signifAtThresholds.ParentVsParent.thresh)",
                                     "inform==T & !is.na(ciscompensated)"),  # narrow data further than it already is for this test?
                        colorvec = c("asecols", "asecols", "asecols", "asecols", "asecols", "asecols") # what colors to use: so far always red/black for T/F
)

# Multi-way comparison - regulatory pattern-style comparisons
multiwaytests<-data.table(myname = c("Regulatory pattern (original, issues)", "Regulatory pattern (updated, all)",
                                     "Regulatory pattern (updated, simplified)", "Regulatory pattern (updated, changed only)",
                                     "Regulatory pattern (updated, simplified, changed only)"),
                          shortname = c("regpat_orig", "regpat_update", "regpat_update_simple", "regpat_update_changed", "regpat_update_simplechanged"),
                          myxname = rep("", 5), # filler - so matches twowaytests (sometimes that's helpful)
                          columnname = c("regclass.orig", "regclass.update", "regclass.update.comb", "regclass.update", "regclass.update.comb"),
                          mynarrow = c("!is.na(regclass.orig)", "!is.na(regclass.update)", "!is.na(regclass.update.comb)",
                                       "!is.na(regclass.update) & !regclass.update%in%c('ambiguous', 'conserved')",
                                       "!is.na(regclass.update) & !regclass.update%in%c('ambiguous', 'conserved')"),
                          colorvec = c("gsregcols", "regcols", "combinedrpcols", "regcols.ch", "combinedrpcols.ch"))

# expression info, descriptions (NEW)
expinfo<-data.table(colname = c("meanexp.vst", "meanexp.shlog"),
                    ylabel = c("VST-normalized gene expression\n(parental pair)",
                               "Log2-normalized gene expression\n(parental pair)"),
                    subt = c("VST-normalized gene expression data", "Log2-normalized gene expression data"),
                    fname = c("vstnorm", "shlognorm"))

# --- Do 2-way tests, plots
cat("--....Doing two-way Mann-Whitney type analyses of quantitative test data vs. various gene categorizations....\n")
tissues<-c("Hindgut", "Midgut")

invisible(lapply(1:nrow(expinfo), function(ei){
  mycol<-expinfo[ei, colname]
  myy<-expinfo[ei, ylabel]
  mys<-expinfo[ei, subt]
  myf<-expinfo[ei, fname]

  mwout<-lapply(1:nrow(twowaytests), function(ind){
    # Test
    mymw<-rbindlist(lapply(tissues, function(tis){
      rbindlist(lapply(combswall, function(strn){

        thisdat<-dat.wall[tissue==tis & strain==strn & eval(parse(text = twowaytests[ind, mynarrow])), ]

        if(length(thisdat[get(twowaytests[ind, columnname])==T, get(mycol)])>0 &
           length(thisdat[get(twowaytests[ind, columnname])==F, get(mycol) ])>0){
          res<-wilcox.test(thisdat[get(twowaytests[ind, columnname])==T, get(mycol)],
                           thisdat[get(twowaytests[ind, columnname])==F, get(mycol) ])

          out<-data.table(tissue = tis,
                          strain = strn,
                          expmeasure = mycol,
                          test = twowaytests[ind, myname],
                          med.T = thisdat[get(twowaytests[ind, columnname])==T, median(get(mycol), na.rm = T)],
                          med.F = thisdat[get(twowaytests[ind, columnname])==F, median(get(mycol), na.rm = T)],
                          T.greater = thisdat[get(twowaytests[ind, columnname])==T, median(get(mycol), na.rm = T)] >
                            thisdat[get(twowaytests[ind, columnname])==F, median(get(mycol), na.rm = T)],
                          ngenes.T = thisdat[, sum(get(twowaytests[ind, columnname])==T, na.rm = T)],
                          ngenes.F = thisdat[, sum(get(twowaytests[ind, columnname])==F, na.rm = T)],
                          W = res$statistic,
                          wilcox.p.value = res$p.value)
        }else{ # not enough obs - skip this one
          out<-NULL
        }
        return(out)
      }))
    }))

    # Plot
    ## Format labels
    mymw[, label:= paste(ifelse(T.greater, paste(twowaytests[ind, shortname], "> non; p ="),
                                paste(twowaytests[ind, shortname], "< non; p =")),
                         ifelse(wilcox.p.value < 0.001, paste(formattable::scientific(wilcox.p.value, digits = 1)),
                                round(wilcox.p.value, digits = 3)))]
    ## Make plot
    plt<-sinawmed(datin = dat.wall[eval(parse(text = twowaytests[ind, mynarrow])),],
                  xcol = twowaytests[ind, columnname], ycol = mycol, facrow = "strain", faccol = "tissue",
                  colorcol = twowaytests[ind, columnname], colorvec = eval(as.name(twowaytests[ind, colorvec])),
                  boxcol = rgb(0, 0, 1, 0.6), faclabels = mymw, myxlab= twowaytests[ind, myxname],
                  myylab = myy, mytitle = twowaytests[ind, myname], mysubt = mys, myscales = "fixed", facvec = combswall)

    # Return
    return(list(mw = mymw[, .SD, .SDcols = -c(ncol(mymw))], plt = plt))
  })

  # Save test results
  write.table(rbindlist(lapply(mwout, function(x) x$mw)),
              file.path(p$outdir, paste0(p$baseoutname, "_", myf, "_twowaytests_mannwhitneyresults.txt")),
              sep = "\t", quote = F, row.names = F)

  # Save plots (in one PDF)
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_", myf, "_twowaytests_sinaplots.pdf")),
      8, 8)
  invisible(lapply(mwout, function(x) print(x$plt)))
  invisible(dev.off())

  # --- Do multi-way tests, plots
  cat(paste("--....Doing multi-way ANOVA type analyses of quantitative test data vs. various gene categorizations (", mys, ")....\n"))
  # Perform tests & make plots

  multwayout<-lapply(1:nrow(multiwaytests), function(ind){
    myaov<-lapply(tissues, function(tis){
      oneav<-lapply(combswall, function(strn){

        allouts<-anovatuk(dat = dat.wall[tissue==tis, ], mystrain = strn, testinforow = multiwaytests[ind, ],
                          testagainstcol = mycol, colorvec = eval(as.name(multiwaytests[ind, colorvec])))

        # Format for return
        anout<-data.table(strain = strn, tissue = tis,  expmeasure = mycol, category = multiwaytests[ind, myname],
                          allouts$anout)
        tukout<-data.table(strain = strn, tissue = tis, expmeasure = mycol, category = multiwaytests[ind, myname],
                           allouts$tukout)
        tuklabs<-data.table(strain = strn,  tissue = tis, expmeasure = mycol, category = multiwaytests[ind, myname],
                            allouts$tuklabs)
        ns<-data.table(strain = strn, tissue = tis,  expmeasure = mycol, category = multiwaytests[ind, myname],
                       allouts$ns)

        return(list(antab = anout, tukres = tukout, tuklabs = tuklabs, ns = ns))
      })
      antab<-rbindlist(lapply(oneav, function(x) x$antab))
      tukres<-rbindlist(lapply(oneav, function(x) x$tukres))
      tuklabs<-rbindlist(lapply(oneav, function(x) x$tuklabs))
      ns<-rbindlist(lapply(oneav, function(x) x$ns))

      return(list(antab = antab, tukres = tukres, tuklabs = tuklabs, ns = ns))
    }) # end lapply over tissues
    # Format test results together (again)
    antab<-rbindlist(lapply(myaov, function(x) x$antab))
    tukres<-rbindlist(lapply(myaov, function(x) x$tukres))
    tuklabs<-rbindlist(lapply(myaov, function(x) x$tuklabs))
    ns<-rbindlist(lapply(myaov, function(x) x$ns))

    # --- PLOTS
    # Make plots
    ## ANOVA p-value labels
    antab[, label:=paste("ANOVA p =", ifelse(pvalue < 0.001, paste(formattable::scientific(pvalue, digits = 1)),
                                             round(pvalue, digits = 3)))]

    ### Absolute y axes
    anplt<-sinawmed(datin = dat.wall[eval(parse(text = multiwaytests[ind, mynarrow])) , ],
                    xcol = multiwaytests[ind, columnname], ycol = mycol, facrow = "strain",
                    faccol = "tissue", colorcol = multiwaytests[ind, columnname],
                    colorvec = eval(as.name(multiwaytests[ind, colorvec])),
                    boxcol = rgb(0, 0, 0, 0.6), faclabels = antab, myxlab = multiwaytests[ind, myname],
                    myylab = myy, mysubt = mys, mytitle =multiwaytests[ind, myname],
                    myscales = "fixed", facvec = combswall) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

    # Return
    return(list(antab = antab, tukres = tukres, ns = ns, plt = anplt))
  })

  # Save test results
  names(multwayout)<-multiwaytests$shortname
  ## Combine
  ansave<-rbindlist(lapply(multwayout, function(x) x$antab))
  tuksave<-rbindlist(lapply(multwayout, function(x) x$tukres))
  nsave<-rbindlist(lapply(multwayout, function(x) x$ns))
  ## Write out
  write.table(ansave, file.path(p$outdir, paste0(p$baseoutname, "_", myf, "_multiwaytests_ANOVAresults.txt")),
              sep = "\t", quote = F, row.names = F)
  write.table(tuksave, file.path(p$outdir, paste0(p$baseoutname,"_", myf, "_multiwaytests_TukeyHSDresults.txt")),
              sep = "\t", quote = F, row.names = F)
  write.table(nsave, file.path(p$outdir, paste0(p$baseoutname,"_", myf, "_multiwaytests_nspercategory.txt")),
              sep = "\t", quote = F, row.names = F)

  # Save plots
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_", myf, "_multiwaytests_sinaplots.pdf")),
      10, 10)
  invisible(lapply(multwayout, function(x) print(x$plt)))
  invisible(dev.off())


}))# BIG NEW LAPPLY: diff expression normalizations


#### Script completion message & session information ####
cat(".....glaserschmitt_drosase_meanexprvsaseetc.R processing complete! Session information:....\n")
sessionInfo()
