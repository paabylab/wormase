#! /usr/bin/env/ Rscript
# Looks at mean expression vs. variation data from Wolf et al 2023: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010833
# by Avery Davis Bell, begun 2024.05.20
require(data.table)
require(argparser)
require(ggplot2)

#### Functions [from prev scripts] ####
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

getrank<-function(val, vec){
  # Get rank, quantile of val in terms of vector. ties are averaged
  myranks<-rank(vec)
  valrank<-myranks[which(vec==val)][1]
  return(data.table(nrank = valrank, proprank = valrank/length(vec)))
}

anovatuk<-function(dat, columnname, testagainstcol = "testDat", colorvec){
  # Runs an ANOVA of category vs. test data (all in all date) for one strain, one category
  # In: dat - data with columns any specified in multiwaytests (columnname, mynarrow), testagainstcol
  #     mystrain - strain to narrow to for this test
  #             columnname (column CATEGORY data is in alldat),
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
  thisdat<-copy(dat)
  yvals<-thisdat[, get(testagainstcol)] # naming for prettier stats calls, returns
  catvals<-factor(thisdat[, get(columnname)], # naming for prettier stats calls, returns
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
  setnames(tuklabs, "cat1", columnname) # Name so it's same column name as in overall data

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
  if(!is.na(facrow) & !is.na(faccol)){
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
  }

  # Other code
  # stat_summary(fun = "median", geom = "point", col = "red") +# this auto-computes and adds the medians! but couldn't figure out making it horizontal bar easily
  # geom_point(stat = "summary", fun = median, color = "red") + # this auto-computes and adds the medians! but couldn't figure out making it horizontal bar easily


  # Deal with if want scales to be free or not, or to do both
  # May want to re-order strains, sites - provide something with it factor-ized if so
  return(plt)
}

#### Arguments & Inputs ####
p<-arg_parser("Looks at mean expression vs. variation data from Wolf et al 2023: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010833",
              name = "wolf2023humexpanalyses.R", hide.opts = TRUE)

# Inputs
p<-add_argument(p, "--wolfsupp4",
                help = "Filepath to S4 Data from Wolf et al 2023 (https://doi.org/10.1371/journal.pgen.1010833.s013), across study rank xlsx data tab as tab-delimited file",
                type = "character")

# Outputs
p<-add_argument(p, "--outdir",
                help = "Output directory",
                type = "character",
                default = "out")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "humangeneexp_wolf2023")

# Parse
p<-parse_args(p)

#### Do the work ####
expdat<-fread(p$wolfsupp4, header = T)
# --- correlation of SD rank & mean rank
# All correlations
meansdcor<-allcors(x = expdat$Mean.Rank, y = expdat$SD.Rank, info = "Mean rank vs SD rank")

# Plot
sctplot<-ggplot(expdat, aes(SD.Rank, Mean.Rank)) + geom_point(alpha = 0.2) +
  ggtitle("Human gene expression mean, variation rankings", subtitle = "from Wolf et al 2023") +
  theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14),
                     axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                     title = element_text(size = 16), legend.position = "none",
                     plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
                     strip.text.y = element_text(size = 14))

# Save these results/plots
write.table(meansdcor,
            file.path(p$outdir, paste0(p$baseoutname, "_cortestresults_meanrankvssdrank.txt")),
            sep = "\t", quote = F, row.names = F)
pdf(file.path(p$outdir, paste0(p$baseoutname, "_scatterplot_meanrankvssdrank.pdf")), 6, 6)
print(sctplot)
invisible(dev.off())

# --- take top X% variation as DE, bottom X% as not-DE, see if there are differences in these groups
# Set up data
## Add percentile of SD; create bins from this
setkey(expdat, Gene)
expdat[, SD.proprank:=getrank(SD.Rank, expdat$SD.Rank)$proprank , by = Gene]
expdat[, decile:=(floor(SD.proprank*10) + 1)] # works but sometimes last one will be nudged up
expdat[decile>10, decile:=10]
## Other
deccol<-colorRampPalette(c("lightblue", "darkblue"))(10)
names(deccol)<-1:10

## ANOVAs - not best stat but quickest here
anres<-anovatuk(dat = expdat, columnname ="decile", testagainstcol = "Mean.Rank", colorvec = deccol)
write.table(anres$anout,
            file.path(p$outdir, paste0(p$baseoutname, "_ANOVAresults_sddecileVsmeanrank.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(anres$tukout,
            file.path(p$outdir, paste0(p$baseoutname, "_ANOVAresults_tukey_sddecileVsmeanrank.txt")),
            sep = "\t", quote = F, row.names = F)

## Sina plot
pdf(file.path(p$outdir, paste0(p$baseoutname, "_sinaplot_sddecileVsmeanrank.pdf")), 7.5, 6)
sinawmed(datin = expdat, xcol = "decile", ycol = "Mean.Rank", facrow = NA, faccol = NA, colorcol = "decile",
         colorvec = deccol, myxlab = "Decile of gene expression variability", myylab = "Gene expression (rank, mean)",
         mytitle = "Human gene expression mean, variation rankings", mysubt = "from Wolf et al 2023",
         boxcol = "red4")
invisible(dev.off())


# maybe do a couple different %s [function that'll do this]
#     MW test AND sina plot
# If I add rank to each row this'll be super easy...

#### Script completion message & session information ####
cat("....wolf2023humexpanalyses.R processing complete! Session information:....\n")
sessionInfo()
