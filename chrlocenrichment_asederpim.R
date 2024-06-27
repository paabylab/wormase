#! /usr/bin/env/ Rscript
# For examining chromosomal location/enrichment of ASE genes and other ASE-related gene classifications. Using:
# annotated ASE results from ase_de_annotategenes_deseq2_fromemaseout.R; 
# regulatory patterns from ase_de_cistransclassifications.R; 
# inheritance modes from f1_parental_comparisons_withinstrain.R

# by Avery Davis Bell, begun 2022.11.01. Updated 2024.03
require(argparser)
require(data.table)
require(ggplot2)
require(RColorBrewer)

#### Functions ####
freadcombstrainsimrpdom<-function(exampinput, strains, inhmode, regpats, chrreg){
  # Reads in exampinput file for each strain in strains; combines the data.tables with added (first) column specifying strain
  #  Adds inheritance mode, regulatory pattern info.
  # In: exampinput, example filepath to file to read in as data.table, with strain identified in filename. Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz
  #     strains, vector of strain IDs to substitute into exampinput and read in
  #     inhmode, path to file containing Inheritance mode classifications for all genes and strains. Columns strain, gene_id, inhmode
  #     regpats, path to file containing Regulatory pattern classifications for all genes, strains for which classification could be made - columns strain, gene_id, regclass
  #     chrreg, data.table containing domain info (for assigning each gene to its domain). Columns:
  #       chr, domain [tip, arm, center], subdomain [left or right where relevant], start [bp position, included], end [bp position, included]
  # Out: data.table with columns strain, all columns in exampinput, regclass, inhmode, domain, subdomain [chr region input-related]. Keyed by strain, gene_id
  
  out<-rbindlist(lapply(strains, function(x){
    cat(paste("Reading in file for strain", x, "; recommend checking that this in fact is for the right strain! Path is: \n",
              file.path(gsub("STRAIN", replacement = x, exampinput, fixed = T), "\n")))
    out1<-fread(file.path(gsub("STRAIN", replacement = x, exampinput, fixed = T)))
    out1[,strain:=x]
  }))
  setcolorder(out, "strain")
  setkey(out, strain, gene_id)
  
  # Add inheritance mode
  im<-fread(inhmode, header = T)
  setkey(im, strain, gene_id)
  out<-im[out]
  
  # Add regulatory pattern
  rp<-fread(regpats, header = T)
  setkey(rp, strain, gene_id)
  out<-rp[out]
  
  # Assign each gene to domain
  out[,`:=`(genepos = start + (end-start)/2, genechr = chr)] # genepos is midpoint of gene; genechr is to differentiate when combine with other DT
  dom<-out[, chrreg[chr==genechr & start<=genepos & end>=genepos, .(domain, subdomain)], by = .(gene_id, strain)]
  setkey(dom, strain, gene_id)
  out<-dom[out[, .SD, .SDcols = -"genechr"]]

  # Reorder and return
  setcolorder(out, c("strain", "gene_id", "display_name", "chr", "start", "end", "genepos", "domain", "subdomain"))
  return(out)
}

posdensplot<-function(datin, chrs, poscol = "genepos", colorcol = "signifAtThresholds.ASE", colorvec,
                      stackfacet = "strain", labcolor = "", myylab = "Density (gene midpoints)", mytitle = "", 
                      mysubt = "", myscales = "free", chrreg = F){
  # Makes density plots of gene location on chromosome - split by some categorical information (e.g. ASE gene vs not, provided in fillcol)
  # In: datin, data.table of data to plot. Must have columns chr; named with values of poscol, fillcol, stackfacet
  #     chrs, chromosome information. Chromosomes will be plotted in this order!! Columns chr, start, end. (used for ordering and to make sure plots extend to beginnings, ends of chrs)
  #     poscol, name of column with basepair genomic position to plot on X/get density for
  #     colorcol, name of column with categorical info - each category will have different line on density plot
  #     colorvec, named vector specifying colors for each value of colorcol in order they should be in legend. Names are values, values are colors
  #           This order is also used to plot lines - data re-ordered to be in this same order; last values plotted on top
  #     labcolor, what to label the color legend
  #     myylab, y axis label
  #     mytitle, plot title
  #     mysubt, plot subtitle
  #     myscales, passed to scales= of faceting. **MUST be "free" or "free_x" - x axis needs to vary
  #     chrreg, F or data.table with tip, arm, center chromosomal region start and end positions (columns chr, domain, start, end)
  # Out: the ggplot
  
  # Set up data (extend to ends of chromosomes, etc)
  pdat<-copy(datin[,.(chr, get(poscol), get(colorcol), get(stackfacet))]) # don't want to modify original dt
  setnames(pdat, c("chr", "pltpos", "splitcol", "pltstack"))
  
  cdat<-rbindlist(lapply(1:pdat[, length(unique(pltstack))], function(x){
    rbind(chrs[, .(chr, pltpos = start, splitcol = pdat[,unique(splitcol)], pltstack = pdat[,unique(pltstack)][x])],
          chrs[, .(chr, pltpos = end, splitcol = pdat[,unique(splitcol)],  pltstack = pdat[,unique(pltstack)][x])])
  }))
  pdat<-rbind(pdat, cdat)
  pdat[,chr:=factor(chr, levels = chrs$chr)]
  pdat<-pdat[order(match(splitcol, names(colorvec)))]
  
  # Make plot
  plt<-ggplot(pdat, aes(pltpos/1e06)) + geom_density(aes(color = splitcol)) +
    facet_grid(pltstack~chr, scales = myscales) +
    scale_color_manual(values = colorvec) + labs(color = labcolor) +
    ggtitle(mytitle, subtitle = mysubt) + xlab("Genomic position (Mb)") + ylab(myylab) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13), legend.position = "bottom",
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
                       strip.text.y = element_text(size = 14), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
  
  # Do:
  #   overlay with center vs arms vs etc - maybe draw those first? and do line color instead of fill for density?
  if(is.data.table(chrreg)){
    plt<-plt + geom_rect(data = chrreg, inherit.aes = F, aes(xmin = start/1e06, xmax = end/1e06, 
                                                             ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
      scale_fill_manual(values = c("arm" = "white", "center" = "gray", "tip" = "gray")) 
  }
  
  return(plt)
  
  # Dev notes
  # how does it work with facet? can you compare across chr or not with density
  # density area under curve is 1. doesn't show diffs in absolute number across chr.
  # If do a plain histogram, wildly different numbers across i.e. gene categories
}

chisqfull<-function(datin, catrowcol = "signifAtThresholds.ASE", rowcats = c(T, F),
                    catcol = "domain", colcats = c("arm", "center", "tip"), annotdt = NULL){
  # Runs chi-sq test and reports these results as well as underlying counts, proportions in rowwise-focused way
  # In: datin, data where rows are observations. Must have columns named with values of catrowcol, catcol
  #     catrowcol, name of column in datin that contains categories to use as rows of contingency table, output
  #     rowcats, values of catrowcol to tabulate 
  #     catcol, name of column in datin that contains categories to use as columns of contingency table, output
  #     colcats, values of catcol to tabulate
  #     annotdt, optional one-row data.table with further information to include in output
  # Out: data.table with as many rows as rowcats; columns:
  #     <any in annotdt>
  #     category, rowcats value for this row
  #     <One for each of the column categories (colcats)> - number of observations with this row, column value
  #     <One for each of the column categories (colcats)>_pOfTotal - proportion of TOTAL observations with this column value - same for each row
  #     <One for each of the column categories (colcats)>_pOfThisRowCategory - proportion of this row category with each column value
  #     ChiSq - test statistic for this test. Same for each row!
  #     Df - Chisq test Degrees of freedom for this test. Same for each row!
  #     pvalue - ChiSq test p-value for this test. Same for each row!
  
  # Create contingency table & run chi-sq
  cttab<-sapply(colcats, function(x) datin[, sapply(rowcats, function(y) sum(get(catcol)==x & get(catrowcol)==y))])
  res<-chisq.test(cttab)
  
  # Format output counts, proportions
  propAll<-data.table(t(as.matrix(colSums(cttab)/sum(cttab)))) # proportion of TOTAL genes that each column category is (will be same across output rows)
  setnames(propAll, paste0(names(propAll), "_pOfTotal"))
  propInCat<-data.table(cttab/rowSums(cttab))
  setnames(propInCat, paste0(names(propInCat), "_pOfThisRowCategory"))
  
  # Annotate with other output info & return
  out<-data.table(annotdt, category = rowcats, cttab, propAll, propInCat,
                  ChiSq = res$statistic, Df = res$parameter, pvalue = res$p.value)
  return(out)
}

longformnsps<-function(datin, catrowcol = "signifAtThresholds.ASE", rowcats = c(T, F),
                         catcol = "domain", colcats = c("arm", "center", "tip"), annotdt = NULL){
  # Gets categories x chromosome domain numbers, proportions. In long format oriented way - one row per 
  # In: datin, data where rows are observations. Must have columns named with values of catrowcol, catcol
  #     catrowcol, name of column in datin that contains categories to use as rows of contingency table, output
  #     rowcats, values of catrowcol to tabulate 
  #     catcol, name of column in datin that contains categories to use as columns of contingency table, output
  #     colcats, values of catcol to tabulate
  #     annotdt, optional one-row data.table with further information to include in output
  # Out: data.table with one row per combination of rowcat and column cat. (e.g., TRUE ASE and center domain)
  #        columns (last few not in order for ease of documentation)::
  #       <any in annotdt>
  #       categ1name, catrowcol input - what is name of the category 1 column
  #       category1, Actual category of categ1 that this row describes (e.g. TRUE for ASE)
  #       categ2name, catcol input - what is name of the category 2 column
  #         category2, Actual category of categ2 that this row describes (e.g. center for domain)
  #       total.n, total # genes - input number
  #       total.thiscateg1, total # genes with categ1name (e.g., TRUE for ASE)
  #       total.thiscateg2, total # genes with categ2name at category2 (e.g., center for domain)
  #       n.thiscombo, number of genes that have categ1name at category 1 and categ2name at category 2 (e.g., TRUE for ASE and in center domain)
  #       prop.<total, thiscateg1, thiscateg2>, proportion this number comprises: of all genes, of genes with category 1 having categ1 name (e.g., of ASE = TRUE genes),
  #             of genes with category 2 having categ2 name (e.g., of domain = center genes)
  #       low95.<total, thiscateg1, thiscateg2>, lower binomial 95% confidence interval bound for the proportion specified
  #       high95.<total, thiscateg1, thiscateg2>, upper binomial 95% confidence interval bound for the proportion specified
  
  # subfunctions
  pci<-function(x, n, outnamesuff = ""){
    # Gets x/n proportion and 95% binomial CI. VECTORWISE. As data.table
    # In: x, numeric (or int) vector of numerators
    #     n, numeric (or int) vector of denominators
    #     outnamesuff, optional suffix for data.table output names
    # Out: data.table with one row for each entry of x/n
    if(length(x)!=length(n)){
      stop("x and n must be the same length!")
    }
    res<-lapply(1:length(x), function(i){
      if(n[i] > 0){
        binom.test(x[i], n[i])
      }else{list(estimate = NaN, conf.int = c(NaN, NaN))}
    })
    out<-data.table(prop = sapply(res, function(x) x$estimate),
                    low95ci = sapply(res, function(x) x$conf.int[1]),
                    high95ci = sapply(res, function(x) x$conf.int[2]))
    setnames(out, paste0(names(out), outnamesuff))
    return(out)
  }
  
  # Get numbers
  ndt<-rbindlist(lapply(colcats, function(x){
    rbindlist(lapply(rowcats, function(y){
      data.table(categ1name = catrowcol,
                 category1 = y,
                 categ2name = catcol,
                 category2 = x,
                 total.n = nrow(datin),
                 total.thiscateg1 = datin[, sum(get(catrowcol)==y)],
                 total.thiscateg2 = datin[, sum(get(catcol)==x)],
                 n.thiscombo = datin[, sum(get(catcol)==x & get(catrowcol)==y)])
    }))
  })) 
  
  # Get proportions
  ndt<-data.table(ndt,
                  ndt[, pci(x = n.thiscombo, n = total.n, outnamesuff = ".total")],
                  ndt[, pci(x = n.thiscombo, n = total.thiscateg1, outnamesuff = ".thiscateg1")], # proportion of those that have ASE (EG)
                  ndt[, pci(x = n.thiscombo, n = total.thiscateg2, outnamesuff = ".thiscateg2")]) # proportion of those that are in the arm domain (EG)
  
  # Add annotating columns & return
  out<-data.table(annotdt, ndt)
  return(out)
}

stackedbar<-function(pinhclass, mycolors, xcol = "strain", stackcol = "propName", stacknumcol = "prop", xorder = NA,
                     legendlabel = "", myxlab = "Strain", myylab = "Proportion genes", mytitle = "", mysubt = ""){
  # Makes stacked bar chart, one stacked bar per category in x col. Copied from other scripts.
  # In: pinhclass, data.table with all data to plot. Must have columns named with values of xcol, stackcol, stacknumcol
  #     mycolors, colors for bars: vector of length of categories in data; also used as ORDER for bars!! Names are categories, values are colors
  #     xcol, column containing category for X axis of bar chart
  #     stackcol, column containing category to split bars into
  #     stacknumcol, column containing numbers to plot as stacked bar
  #     xorder, OPTIONAL order in which to arrange categories on x axis (in xcol)
  #     legendlabel
  #     myxlab, myylab, mytitle, mysubt: plot labels as intuitively named
  # Out: ggplot2 barplot object. Can facet this externally!   
  
  # Format data
  pdata<-copy(pinhclass)
  pdata[,mystackcol:=factor(pdata[,get(stackcol)], levels = names(mycolors))] # relevel factor
  if(!is.na(xorder[1])){
    # Re-level factor 
    pdata[,xcol:=factor(pdata[,get(xcol)], levels = xorder)]
  }else{pdata[,xcol:=get(xcol)]}
  
  plt<-ggplot(pdata, aes(xcol, eval(as.name(stacknumcol)))) +
    geom_bar(aes(fill = mystackcol), stat = "identity", position = "stack") + labs(fill = legendlabel) +
    xlab(myxlab) + ylab(myylab) + ggtitle(mytitle, subtitle = mysubt) +
    scale_fill_manual(values = mycolors) + # provided colors
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  return(plt)
}

divplot<-function(onedat, xcol = "nvars.thou", xlabel = "SNVs/INDELs vs. N2 (thousands)", 
                  ylabel = "", mytitle = "", mysubt = "", mysize = 1, ymincol = "lowci", ymaxcol = "highci"){
  # One off just to make this plot. Copied from ase_de_cistransclassification.R
  plt<-ggplot(onedat, aes(eval(as.name(xcol)), prop)) +
    geom_pointrange(aes(ymin = eval(as.name(ymincol)), ymax = eval(as.name(ymaxcol)), color = strain), size = mysize) +
    xlim(c(0, onedat[,max(get(xcol))] + 0.1*onedat[,max(get(xcol))])) +
    xlab(xlabel) + ylab(ylabel) + ggtitle(mytitle, subtitle = mysubt) + labs(color = "Strain") + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 9))
  
  return(plt)
}


#### Arguments & inputs ####
p<-arg_parser("For examining chromosomal location/enrichment of ASE genes and other ASE-related gene classifications. Using:
      annotated ASE results from ase_de_annotategenes_deseq2_fromemaseout.R; 
      regulatory patterns from ase_de_cistransclassifications.R; 
      inheritance modes from f1_parental_comparisons_withinstrain.R",
              name = "chrlocenrichment_asederpim.R", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--resmatinput",
                help = "Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
                type = "character")
p<-add_argument(p, "--regpats",
                help = "Regulatory pattern classifications for all genes, strains for which classification could be made - columns strain, gene_id, regclass.
                (*_regpattern_per1plusunqalnggenestrain_<conf int>confdeaseoverlap.txt.gz output of ase_de_cistransclassifications.R)",
                type = "character")
p<-add_argument(p, "--inhmode",
                help = "Inheritance mode classifications for all genes and strains. Columns strain, gene_id, inhmode
                (*_inheritancemode_pergenestrain.txt.gz output from f1_parental_comparisons_withinstrain.R)",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally.",
                type = "character")
p<-add_argument(p, "--strains",
                help = "Strains to read in results for and process together. Either comma-separated (no spaces) list or path to no-header file with one line per strain.
                Must match how strains are named in input filenames. Strains will be plotted/leveled in this order.",
                default = "JU1088,EG4348,CB4856,QX1211")
p<-add_argument(p, "--refstrain",
                help = "Name of reference strain, matching name used in columns of resmatinput. Also used for plot labelling.",
                default = "N2")
p<-add_argument(p, "--chrfile",
                help = "File with columns chr, start, end. Used for plotting coordinates",
                type = "character")
p<-add_argument(p, "--chrregions",
                help = "File containing chromosome tip, arm, center start/end positions; columns chr, domain [tip, arm, center], 
                subdomain [left or right where relevant], start [bp position, included], end [bp position, included]",
                type = "character")

# Threshold arguments
p<-add_argument(p, "--informthresh",
                help = "Gene must have this or more unique alignments in each sample to be considered informative for ASE/cis-trans analyses",
                default = 2) #NB setting default to what pilot run was, not suggested/best default value!

# Related to plotting vs. divergence. **NEW 2024.03
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

# Output directory
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

# chromosomes, regions
chrs<-fread(p$chrfile, header = T)
chrreg<-fread(p$chrregions, header = T)

# --- Read in expression-based data, add domain info
cat("....Reading in data....\n")
dat<-freadcombstrainsimrpdom(exampinput = p$resmatinput, strains, p$inhmode, p$regpats, chrreg)
## Save out genes annotated with their location/domain for other analyses! New 2024.03
write.table(dat[strain==strains[1], .(gene_id, display_name, chr, start, end, genepos, domain, subdomain)],
            gzfile(file.path(p$outdir, paste0(p$baseoutname, "_geneswithdomain.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

## add informative column
dat[,inform:=ifelse(unqalnmts.min>=p$informthresh, T, F)]

# Set up colors
asecols<-c("black", "red") # ASE yes/no
names(asecols)<-c(FALSE, TRUE)
# Fancier
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")
inhmodes<-c("no_change", "ambiguous", "additive", paste0(p$refstrain, "_dominant"), "alt_dominant", "overdominant", "underdominant")
## Reg pattern, inheritance mode COLORS
regcols<-rep(NA, length(regpats))
names(regcols)<-regpats
regcols["conserved"]<-"lightgray"
regcols["ambiguous"]<-"darkgray"
regcols[c("compensating", "compensatory", "overcompensating")]<-brewer.pal(9, "Blues")[c(3, 6, 9)] # wider dynamic range of color than doing brewer.pal(3, "Blues")
regcols[is.na(regcols)]<-brewer.pal(4, "Set2")[c(1, 2, 4)]
## colors for inheritance modes
inhcols<-rep(NA, length(inhmodes))
names(inhcols)<-inhmodes
inhcols["no_change"]<-"lightgray"
inhcols["ambiguous"]<-"darkgray"
inhcols[c("additive", "N2_dominant", "alt_dominant", "overdominant", "underdominant")]<-rev(brewer.pal(5, "Paired"))

#### Make density plots ####
# Set up data information: what to plot, include, exclude, etc
## Genes: hyperdivergent vs. not
hypdivinfo<-data.table(myname = c("inclhypdiv", "nonhypdiv"), # ***hypdiv & low coverage excluded!!
                       descrip = c("Including hyperdivergent genes", "Excluding hyperdivergent & low coverage genes"),
                       vals = c("(hypdiv == T | hypdiv == F)", "(hypdiv == FALSE & lowDNACov == FALSE)"))
## Categories: ASE, DE, reg pattern, inh mode
##        include inform == T vals; descrips & short labels
coldata<-data.table(shortname = c("aseVnot", "deVnot", "deVnot", "regpattern", "inhmode", "inhmode"),
                    legendname = c("ASE", "DE", "DE", "Reg.\npattern", "Inh.\nmode", "Inh.\nmode"),
                    descrip = c("ASE (informative genes)", "DE (all genes)", "DE (ASE informative genes)",
                                "Regulatory pattern (ASE informative genes)", "Inheritance mode (all genes)",
                                "Inheritance mode (ASE informative genes)"),
                    colname = c("signifAtThresholds.ASE", paste0("signifAtThresholds.ParentVsParent", p$refstrain),
                                paste0("signifAtThresholds.ParentVsParent", p$refstrain), "regclass", "inhmode", "inhmode"),
                    colorvec = c("asecols", "asecols", "asecols", "regcols", "inhcols", "inhcols"),
                    datsub =c("inform == T", "(inform == T | inform == F)", "inform==T", "inform == T",
                              "(inform == T | inform == F)", "inform== T"))

# --- Make plots
perdatplts<-lapply(1:nrow(coldata), function(datind){
  pergset<-(lapply(1:nrow(hypdivinfo), function(hypind){
    plt<-posdensplot(datin = dat[eval(parse(text = hypdivinfo[hypind, vals])) & eval(parse(text = coldata[datind, datsub])),],
                     chrs = chrs, colorcol = coldata[datind, colname], colorvec = eval(as.name(coldata[datind, colorvec])),
                     labcolor = coldata[datind, legendname], mytitle = coldata[datind, descrip], mysubt = hypdivinfo[hypind, descrip], 
                     chrreg = chrreg, myscales = "free")
    return(plt)
  }))
  names(pergset)<-hypdivinfo$myname
  return(pergset)
})

# --- Save plots
invisible(lapply(unique(coldata$shortname), function(sname){
  myplts<-perdatplts[which(coldata$shortname==sname)]
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_", sname, "_genedensityplot.pdf")), 
      10, max(4, 1.8*length(strains)))
  invisible(lapply(myplts, function(x){
    print(x$inclhypdiv)
    print(x$nonhypdiv)
  }))
  invisible(dev.off())
  return(NULL)
}))

# Can do stats first & then annotate them on density plots if want to once get to stats...

#### ChiSq tests ####
chisqres<-rbindlist(lapply(1:nrow(coldata), function(datind){
  rbindlist(lapply(1:nrow(hypdivinfo), function(hypind){
    dat[eval(parse(text = hypdivinfo[hypind, vals])) & eval(parse(text = coldata[datind, datsub])),
        chisqfull(datin = .SD, catrowcol = coldata[datind, colname], rowcats = names(eval(as.name(coldata[datind, colorvec]))),
                  catcol = "domain", colcats = c("arm", "center", "tip"),
                  annotdt = data.table(testedCategory = coldata[datind, descrip], 
                                       hypdivgenes = hypdivinfo[hypind, myname])),
        by = strain]
  }))
}))
write.table(chisqres, file.path(p$outdir, paste0(p$baseoutname, "_chrdomainVmulticategorychisqresults.txt")),
            sep = "\t", quote = F, row.names = F)


#### Make bar plots - new 2024.03 ####
cat("....Making barplot summaries of proportions in different regions....\n")
barpltdir<-file.path(p$outdir, "barplots")
if(!dir.exists(barpltdir)){dir.create(barpltdir)}

# --- Get and save numbers
allprops<-rbindlist(lapply(1:nrow(coldata), function(datind){
  rbindlist(lapply(1:nrow(hypdivinfo), function(hypind){
    dat[eval(parse(text = hypdivinfo[hypind, vals])) & eval(parse(text = coldata[datind, datsub])),
        longformnsps(datin = .SD, catrowcol = coldata[datind, colname], rowcats = names(eval(as.name(coldata[datind, colorvec]))),
                     catcol = "domain", colcats = c("arm", "center", "tip"),
                     annotdt = data.table(testedCategory = coldata[datind, descrip], 
                                          hypdivgenes = hypdivinfo[hypind, myname])),
        by = strain]
  }))
}))

## Save
write.table(allprops, file.path(barpltdir, paste0(p$baseoutname, "_longformproportions_chrdomainXmulticategories.txt")),
            sep = "\t", quote = F, row.names = F)

# --- barplots of all of these
# Update colors to look OK with bars
asecols<-c("gray", "black") # ASE yes/no
names(asecols)<-c(FALSE, TRUE)
domaincolors<-c("gray90", "gray40", "gray10")
names(domaincolors)<-c("arm", "center", "tip")

# Make a zillion plots
invisible(lapply(coldata$shortname, function(sname){
  # Bar plots with domain on x axis, segmented into other category
  pdf(file.path(barpltdir, paste0(p$baseoutname, "_domainVs", sname, "_barplots_domainx.pdf")),
      max(8, 1.75*length(strains)), max(4, 1.6*length(strains)))
  invisible(lapply(coldata[shortname==sname, descrip], function(mydescrip){
    invisible(lapply(1:nrow(hypdivinfo), function(hypind){
      pdat<-allprops[testedCategory==mydescrip & hypdivgenes==hypdivinfo[hypind, myname], ]
      pdat[,strain:=factor(strain, levels = strains)]
      ## NUMBERS
      print(
      stackedbar(pinhclass = pdat,
                 mycolors = eval(as.name(coldata[shortname==sname, unique(colorvec)])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "n.thiscombo", xorder = c("arm", "center", "tip"),
                 legendlabel = "", myxlab = "Chromosomal domain", myylab = "Number of genes", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
      
      ## Global proportions
      print(
      stackedbar(pinhclass = allprops[testedCategory==mydescrip & hypdivgenes==hypdivinfo[hypind, myname], ],
                 mycolors = eval(as.name(coldata[shortname==sname, unique(colorvec)])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "prop.total", xorder = c("arm", "center", "tip"),
                 legendlabel = "", myxlab = "Chromosomal domain", myylab = "Proportion of all genes", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
      ## ** proportion within domain
      print(
      stackedbar(pinhclass = allprops[testedCategory==mydescrip & hypdivgenes==hypdivinfo[hypind, myname], ],
                 mycolors = eval(as.name(coldata[shortname==sname, unique(colorvec)])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "prop.thiscateg2", xorder = c("arm", "center", "tip"),
                 legendlabel = "", myxlab = "Chromosomal domain", myylab = "Proportion of genes in domain", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
    }))
  }))
  invisible(dev.off())
  
  # Bar plots with other category on x axis, segmented into 
  pdf(file.path(barpltdir, paste0(p$baseoutname, "_domainVs", sname, "_barplots.pdf")),
      max(8, 1.75*length(strains)), max(4, 1.6*length(strains)))
  invisible(lapply(coldata[shortname==sname, descrip], function(mydescrip){
    invisible(lapply(1:nrow(hypdivinfo), function(hypind){
      pdat<-allprops[testedCategory==mydescrip & hypdivgenes==hypdivinfo[hypind, myname], ]
      pdat[,strain:=factor(strain, levels = strains)]
      ## NUMBERS
      print(
      stackedbar(pinhclass = pdat,
                 mycolors = domaincolors,
                 xcol = "category1", stackcol = "category2", stacknumcol = "n.thiscombo",
                 xorder = names(eval(as.name(coldata[shortname==sname, unique(colorvec)]))),
                 legendlabel = "", myxlab = mydescrip, myylab = "Number of genes", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
      
      ## Global proportions
      print(
      stackedbar(pinhclass = pdat,
                 mycolors = domaincolors,
                 xcol = "category1", stackcol = "category2", stacknumcol = "prop.total",
                 xorder = names(eval(as.name(coldata[shortname==sname, unique(colorvec)]))),
                 legendlabel = "", myxlab = mydescrip, myylab = "Proportion of all genes", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
      
      ## ** proportion within category
      print(
      stackedbar(pinhclass = pdat,
                 mycolors = domaincolors,
                 xcol = "category1", stackcol = "category2", stacknumcol = "prop.thiscateg1",
                 xorder = names(eval(as.name(coldata[shortname==sname, unique(colorvec)]))),
                 legendlabel = "", myxlab = mydescrip, myylab = "Proportion of genes in this bar category", 
                 mytitle = mydescrip, mysubt = hypdivinfo[hypind, descrip]) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      )
    }))
  }))
  invisible(dev.off())
})) # end lapply through data/categories

#### Proportion in different domains vs. divergence - new 2024.03 ####
cat("....Plotting proportion of genes in chromosomal domain that are various categories vs. divergence from N2....\n")
divdir<-file.path(p$outdir, "divergence")
if(!dir.exists(divdir)){dir.create(divdir)}

# Get in divergence info (so can order strains this way even before plotting divergence)
ndiv<-fread(p$varsvsref)
ndiv[,`:=`(nvars.thou = nvars/1e03, nvars.perkb = nvars/(p$genomebp/1e03))] 
setkey(ndiv, strain)

# Make plots
invisible(lapply(coldata$shortname, function(sname){
  pdf(file.path(divdir, paste0(p$baseoutname, "_propindomainvsdivergencefromN2_", sname, ".pdf")),
      7.5, max(6, length(eval(as.name(coldata[shortname==sname, colorvec]))) * 1.5))
  invisible(lapply(coldata[shortname==sname, descrip], function(mydescrip){
    invisible(lapply(1:nrow(hypdivinfo), function(hypind){
      # Data set up
      pdat<-allprops[testedCategory==mydescrip & hypdivgenes==hypdivinfo[hypind, myname],
                     .(strain, testedCategory, hypdivgenes, categ1name, category1, categ2name, 
                       category2, n.thiscombo, prop.thiscateg2, low95ci.thiscateg2, high95ci.thiscateg2)]
      setkey(pdat, strain)
      pdat<-pdat[ndiv]
      pdat[,strain:=factor(strain, levels = ndiv[order(nvars.thou), strain])]
      setnames(pdat, "prop.thiscateg2", "prop")
      pdat[, category2:=factor(category2, levels = c("arm", "center", "tip"))]
      pdat[, category1:=factor(category1, levels=names(eval(as.name(coldata[descrip==mydescrip, colorvec]))))]
      
      # Make plot
      print(
      divplot(onedat = pdat, ylabel = "Proportion of genes in given domain/column facet\nthat are the row facet's category", mytitle = paste(mydescrip, "- proportion of genes in domain that are each category"),
              mysubt = hypdivinfo[hypind, descrip], ymincol = "low95ci.thiscateg2", ymaxcol = "high95ci.thiscateg2") + 
        facet_grid(category1~category2, scales = "free_y")
      )
      }))
    }))
  invisible(dev.off())
}))

# ideally would have umbrealla categories here (transgressive, compensatory)...but don't

#### Save session info ####
cat("....chrlocenrichment_asederpim.R complete! Session information:....\n")
sessionInfo()
