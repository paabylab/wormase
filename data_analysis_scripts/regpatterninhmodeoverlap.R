#! /usr/bin/env/ Rscript
# Investigate overlap between regulatory patterns and inheritance modes
# FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R; 
#     regulatory patterns from ase_de_cistransclassifications.R; 
#     inheritance modes from f1_parental_comparisons_withinstrain.R
# by Avery Davis Bell, begun 2022.07.28
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

regxinhnums<-function(onedat, regpats, inhmodes, dtadd = NULL){
  # Creates regulatory pattern x inheritance mode counts tables, as well as % of all genes each regulatory pattern-inh. mode classification is
  #     And % of each regulatory pattern each inheritance mode is and vice versa
  # In: onedat, data.table with one row per gene, columns regclass (containing regpats values) and inhmode (containing inhmodes values) 
  #     regpats, character vector of regulatory patterns to count up/co-classify (in desired order)
  #     inhmodes, character vector of inheritance modes to count up/co-classify (in desired order)
  #     dtadd, any information to add to data.table outputs in data.table format - i.e. one row containing strain, gene set information
  # Out: list of lists. $matrices - list of MATRICES:
  #     $ns, regulatory pattern x inheritance modes matrix with counts of genes in each reg pattern-inh mode pair
  #     $pAll, ns/total number of input genes: proportion of genes in each reg pattern-inh mode pair
  #     $pWitninReg, regulatory pattern x inheritance modes matrix containing proportion of each regulatory pattern each inheritance mode is (rows sum to 1)
  #     $pWithinInh, *inheritance modes x regulatory patterns (transposed vs others!) matrix containing proportion of each inheritance mode each regulatory pattern is (rows sum to 1)
  # $dts - list of datatables, same names as matrices, but now first columns are any provided in dtadd,
  #   next column contains regulatory patterns (or inh modes for pWithinInh); following are the numbers
  # $longform - contains all aforementioned information plus some in long format (good for plotting etc). One row per regpat/inhmode combination.
  #   Columns (last few not in order for ease of documentation):
  #   <any in dtadd>
  #     regpat, regulatory pattern for this row (from regpats)
  #     inhmode, inheritance mode for this row (from inhmodes)
  #     n, # genes with this regulatory pattern - inheritance mode combination
  #     totalngenes, total # genes analyzed (i.e. sum of all regpat-inhmode combinations) - same for every row
  #     totalnregpat, total # genes with the given regulatory pattern - same for each row with this reg. pattern
  #     totalninhmode, total # genes with the given inheritance mode - same for each row with this inheritance mode.
  #     prop<allgenes, nregpat, ninhmode> - proportion of each category this row (regpat/inhmode combination) is. .all is of all genes;
  #         .nregpat is of totalnregpat; .ninhmode is of totalninhmode
  #     low95ci<allgenes, nregpat, ninhmode> - lower bound binomial 95% confidence interval on proportion
  #     high95ci<allgenes, nregpat, ninhmode> - upper bound binomial 95% confidence interval on proportion
  
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
  
  # Numbers
  ns<-sapply(inhmodes, function(ih){ sapply(regpats, function(rp) {onedat[,sum(regclass==rp & inhmode==ih)]})})
  imsums<-colSums(ns)
  rpsums<-rowSums(ns)
  # ns<-cbind(ns, rowSums(ns))
  # rownames(ns)<-c(regpats, "total")
  # colnames(ns)<-c(inhmodes, "total")
  
  # Proportions
  pAll<-ns/sum(ns) # Of all
  pWithinReg<-ns/rpsums # Within regulatory pattern (e.g. % of cis that are each inheritance mode)
  pWithinInh<-t(ns)/imsums # Within inheritancde mode (e.g. % of additive that are each regulatory pattern)
  
  # Long format numbers and proportions (with CI?). Very ugly but useful for plotting
  longform<-rbindlist(lapply(inhmodes, function(x){
    data.table(regpat = rownames(ns), inhmode = x, n = ns[,x], totalngenes = sum(ns),
               totalnregpat = rpsums[rownames(ns)], totalninhmode = imsums[x])
  }))
  longform<-data.table(longform, longform[, .(pci(x = n, n = totalngenes, outnamesuff = ".allgenes"),
                                              pci(x = n, n = totalnregpat, outnamesuff = ".nregpat"),
                                              pci(x = n, n = totalninhmode, outnamesuff = ".ninhmode"))])
  longform<-data.table(dtadd, longform) # add dtadd to this too
  
  # Return
  return(list(matrices = list(ns = ns, pAll = pAll, pWithinReg = pWithinReg, pWithinInh = pWithinInh),
              dts = list(ns = data.table(dtadd, as.data.table(ns, keep.rownames = "regpat")),
                         pAll = data.table(dtadd, as.data.table(pAll, keep.rownames = "regpat")), 
                         pWithinReg = data.table(dtadd, as.data.table(pWithinReg, keep.rownames = "regpat")),
                         pWithinInh = data.table(dtadd, as.data.table(pWithinInh, keep.rownames = "inhmode"))),
              longform = longform))
}

stackedbar<-function(pinhclass, mycolors, xcol = "strain", stackcol = "propName", stacknumcol = "prop", xorder = NA,
                     legendlabel = "", myxlab = "Strain", myylab = "Proportion genes", mytitle = "", mysubt = ""){
  # Makes stacked bar chart, one stacked bar per category in x col
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
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
                       strip.text.y = element_text(size = 14))
  
  return(plt)
}

heatplot<-function(pdat, xcat = "inhmode", ycat = "regpat", numcol = "n", 
                   yaxlab = "", xaxlab = "", leglab = "", mytitle = "", mysubt = ""){
  # Makes a heat plot where tiles represent number in ncol, with tiles' x-axis categories from xcat and y axis categories from ycat
  # In: pdat, data to plot - must contain numbers falling into 2 overlapping categories, with columns for each category and the number
  #     xcat, character name of column containing categorical variable for x axis
  #     ycat, character name of column containing categorical variable for y axis
  #     numcol, character name of column containing number by which to color the tiles (so a number matching the specified xcat and ycat categories)
  #     yaxlab, y axis title (optional)
  #     xaxlab, x axis title (optional)
  #     leglab, legend label (optional)
  # Out: Plot

  plt<-ggplot(pdat, aes(eval(as.name(xcat)), eval(as.name(ycat)))) + geom_tile(aes(fill = eval(as.name(numcol)))) + 
    xlab(xaxlab) + ylab(yaxlab) + scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
    ggtitle(mytitle, subtitle = mysubt) + labs(fill = leglab) + theme_bw() + 
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
          strip.text.y = element_text(size = 14))
  
  return(plt)
}

#### Arguments & inputs ####
# --- Command line arguments
p<-arg_parser("Investigate overlap between regulatory patterns and inheritance modes. FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R; 
    regulatory patterns from ase_de_cistransclassifications.R; 
    inheritance modes from f1_parental_comparisons_withinstrain.R", 
              name = "regpatterninhmodeoverlap.R", hide.opts = TRUE)

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
p<-add_argument(p, "--informthresh",
                help = "Gene must have this or more unique alignments (EMASE; or this or more Ornaments counts) in each sample to be considered informative for ASE/cis-trans analyses",
                default = 2) # 2 is default because with pilot data I ran 2 (before I gave this as an option....)
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

# Add inheritance mode
im<-fread(p$inhmode, header = T)
setkey(im, strain, gene_id)
dat<-im[dat]

# Add regulatory pattern
rp<-fread(p$regpats, header = T)
setkey(rp, strain, gene_id)
dat<-rp[dat]

# Regulatory patterns, inheritance mode **in order of interest** [maybe set up colors here too?] ** change if change something upstream **
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")
inhmodes<-c("no_change", "ambiguous", "additive", paste0(p$refstrain, "_dominant"), "alt_dominant", "overdominant", "underdominant")

# Add umbrella categories ??
#         XXXXXXX

#### inheritance mode/regulatory pattern overlap: counts and associated summaries ####
# --- Gene sets: NARROWING from all, currently just doing 2+ unique alignments incl/excly hypdiv/lowcov
## Set up the categories to intersect
infos<-c(paste0("Informative, " , p$informthresh,"+ unique alignments"))
infos.short<-c(paste0(p$informthresh, "plusUnqAlns"))
infos.bool<-c(paste0("unqalnmts.min >= ", p$informthresh))

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

# --- Get counts (& proportions)
nspslist<-lapply(1:nrow(gsets), function(x){
  out<-lapply(strains, function(str){
    regxinhnums(onedat = dat[strain == str & eval(parse(text = gsets[x, exprtxt]))],
                regpats = regpats, inhmodes = inhmodes, dtadd = data.table(strain = str, geneset = gsets[x,shortname]))
  })
  names(out)<-strains
  return(out)
})
names(nspslist)<-gsets$shortname

# Save [combining together though it's kinda gnarly..]
write.table(rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$dts$ns)))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_counts.txt")), sep = "\t", row.names = F, quote = F)
write.table(rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$dts$pAll)))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_props_ofall.txt")), sep = "\t", row.names = F, quote = F)
write.table(rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$dts$pWithinReg)))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_props_ofeachregpat.txt")), sep = "\t", row.names = F, quote = F)
write.table(rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$dts$pWithinInh)))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_props_ofeachinhmode.txt")), sep = "\t", row.names = F, quote = F)

#### --- Bar plots ####
# COLOR SET UP
## pick colors used for each reg pattern for all plots. All the ones with cis-trans opposing in diff shades of same color
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

# MAKE PLOTS
bpdir<-file.path(p$outdir, "barplots")
if(!dir.exists(bpdir)){dir.create(bpdir)}
## Set up data
longfprops<-rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) y$longform))))
longfprops[, strain:=factor(strain, levels = strains)]
longfprops[, inhmode:=factor(inhmode, levels = inhmodes)]
longfprops[, regpat:=factor(regpat, levels = regpats)]

#--- All on one page (strains x gene sets)
pdf(file.path(bpdir, paste0(p$baseoutname, "_strainsxgenesets_globalnsps_barplots.pdf")), 8, 10) # PDF with global Ns, global ps
# Ns
## Inheritance mode bars, reg pattern colors
stackedbar(longfprops, mycolors = regcols, xcol = "inhmode", stackcol = "regpat", 
           stacknumcol = "n", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
           myylab = "Number genes", mytitle = "Number genes", mysubt = "all inheritance mode categories") + 
  facet_grid(strain~geneset)
  # exclude no change!!
stackedbar(longfprops[!(inhmode=="no_change" & regpat=="conserved"),], mycolors = regcols, xcol = "inhmode", stackcol = "regpat", 
           stacknumcol = "n", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
           myylab = "Number genes", mytitle = "Number genes", mysubt = "excluding no_change-conserved genes!") + 
  facet_grid(strain~geneset)
## Reg pattern bars, inheritance mode colors
stackedbar(longfprops, mycolors = inhcols, xcol = "regpat", stackcol = "inhmode", 
           stacknumcol = "n", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
           myylab = "Number genes", mytitle = "Number genes", mysubt = "all regulatory pattern categories") + 
  facet_grid(strain~geneset)
  # exclude conserved
stackedbar(longfprops[!(inhmode=="no_change" & regpat=="conserved"), ], mycolors = inhcols, xcol = "regpat", stackcol = "inhmode", 
           stacknumcol = "n", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
           myylab = "Number genes", mytitle = "Number genes", mysubt = "excluding  no_change-conserved genes!") + 
  facet_grid(strain~geneset)

# Proportion of all genes within strain/gene set
## Inheritance mode bars, reg pattern colors
stackedbar(longfprops, mycolors = regcols, xcol = "inhmode", stackcol = "regpat", 
           stacknumcol = "prop.allgenes", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
           myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set", mysubt = "all inheritance mode categories") + 
  facet_grid(strain~geneset)
# exclude no change!!
stackedbar(longfprops[!(inhmode=="no_change" & regpat=="conserved"),], mycolors = regcols, xcol = "inhmode", stackcol = "regpat", 
           stacknumcol = "prop.allgenes", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
           myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set", mysubt = "excluding no_change-conserved genes!") + 
  facet_grid(strain~geneset)
## Reg pattern bars, inheritance mode colors
stackedbar(longfprops, mycolors = inhcols, xcol = "regpat", stackcol = "inhmode", 
           stacknumcol = "prop.allgenes", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
           myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set", mysubt = "all regulatory pattern categories") + 
  facet_grid(strain~geneset)
# exclude conserved
stackedbar(longfprops[!(inhmode=="no_change" & regpat=="conserved"), ], mycolors = inhcols, xcol = "regpat", stackcol = "inhmode", 
           stacknumcol = "prop.allgenes", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
           myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set", mysubt = "excluding no_change-conserved genes!") + 
  facet_grid(strain~geneset)

invisible(dev.off()) # end PDF with global ns/ps

pdf(file.path(bpdir, paste0(p$baseoutname, "_strainsxgenesets_propswithincategs_barplots.pdf")), 8, 10) # all proportions within categories pdf
# Proportion genes WITHIN x axis categories
## Inheritance mode bars, reg pattern colors
stackedbar(longfprops, mycolors = regcols, xcol = "inhmode", stackcol = "regpat", 
           stacknumcol = "prop.ninhmode", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
           myylab = "Proportion (within inheritance mode)", mytitle = "Proportion within inheritance mode", mysubt = "all inheritance mode categories") + 
  facet_grid(strain~geneset)
## Reg pattern bars, inheritance mode colors
stackedbar(longfprops, mycolors = inhcols, xcol = "regpat", stackcol = "inhmode", 
           stacknumcol = "prop.nregpat", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
           myylab = "Proportion (within regulatory pattern)", mytitle = "Proportion within regulatory pattern", mysubt = "all regulatory pattern categories") + 
  facet_grid(strain~geneset)
invisible(dev.off()) # end PDF with proportions within categories

# --- strains are x categories, colors are reg pattern, one facet per inh mode? This may be useful!!
lapply(1:nrow(gsets), function(x){
  pdf(file.path(bpdir, paste0(p$baseoutname, "_strainbars_globalnsps_barplots_", gsets[x, shortname], ".pdf")), 9, 8) # global ns/ps pdf - within gene set
  # Ns (free y!). This is great!!
  ## Inheritance mode bars, reg pattern colors
  print(stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = regcols, xcol = "strain", stackcol = "regpat", 
             stacknumcol = "n", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
             myylab = "Number genes", mytitle = "Number genes, inheritance mode facets", 
             mysubt = gsets[x, descrip]) + 
    facet_wrap(~inhmode, scales = "free_y"))
  ## Reg pattern bars, inheritance mode colors
  print(stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = inhcols, xcol = "strain", stackcol = "inhmode", 
             stacknumcol = "n", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
             myylab = "Number genes", mytitle = "Number genes, regulatory pattern facets",
             mysubt = gsets[x, descrip]) + 
    facet_wrap(~regpat, scales = "free_y"))
  
  
  # Proportion all genes (free y!)
  ## Inheritance mode bars, reg pattern colors
  print(stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = regcols, xcol = "strain", stackcol = "regpat", 
             stacknumcol = "prop.allgenes", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
             myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set", 
             mysubt = gsets[x, descrip]) + 
    facet_wrap(~inhmode, scales = "free_y"))
  ## Reg pattern bars, inheritance mode colors
  print(stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = inhcols, xcol = "strain", stackcol = "inhmode", 
             stacknumcol = "prop.allgenes", legendlabel = "inh. mode",  myxlab = "Regulatory pattern",
             myylab = "Proportion (all in set)", mytitle = "Proportion of all in strain, gene set",
             mysubt = gsets[x, descrip]) + 
    facet_wrap(~regpat, scales = "free_y"))
  invisible(dev.off()) # close global ns/ps pdf - within gene set
  
  pdf(file.path(bpdir, paste0(p$baseoutname, "_strainbars_propswithincategs_barplots_", gsets[x, shortname], ".pdf")), 9, 8) # all proportions within categories pdf - within gene set
  # Proportion genes WITHIN x axis categories
  ## Inh mode facets, reg pattern colors
  print(
    stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = regcols, xcol = "strain", stackcol = "regpat", 
               stacknumcol = "prop.ninhmode", legendlabel = "reg. pattern", myxlab = "Inheritance mode", 
               myylab = "Proportion (within inheritance mode)", mytitle = "Proportion within inheritance mode", 
               mysubt = gsets[x, descrip]) + 
      facet_wrap(~inhmode)
  )
  ## Reg pattern bars, inheritance mode colors
  print(
    stackedbar(longfprops[geneset==gsets[x, shortname],], mycolors = inhcols, xcol = "strain", stackcol = "inhmode", 
               stacknumcol = "prop.nregpat", legendlabel = "inh. mode", myxlab = "Regulatory pattern", 
               myylab = "Proportion (within regulatory pattern)", mytitle = "Proportion within regulatory pattern", 
               mysubt = gsets[x, descrip]) + 
      facet_wrap(~regpat)
  )
  invisible(dev.off()) # end PDF with proportions within categories within gene set
  
  return(NULL)
})

# for above - consider rotate y ax labels further? (so would be perpendicular and fit?)

#--- Gene sets individually, one page per type of plot - NAH, not doing for now

# Label with Ns somehow? that'd be great [though busy]

#     Umbrella categories too? transgressive = under/overdominant (? others); etc (cis/trans opposing at all?)
#         Could simply add columns that have bigger categories and do the same procedure 2x, that likely makes the same sense

#### Heat plots ####
hpdir<-file.path(p$outdir, "heatplots")
if(!dir.exists(hpdir)){dir.create(hpdir)}

#--- Total numbers
pdf(file.path(hpdir, paste0(p$baseoutname, "_strainsxgenesets_globalns_heatplots.pdf")), 12, 6)
heatplot(pdat = longfprops, numcol = "n", leglab = "N genes", mytitle = "All categorized genes") +
  facet_grid(geneset~strain)
## **excluding only conserved, no change combination!!
heatplot(pdat = longfprops[!(inhmode=="no_change" & regpat=="conserved"), ], numcol = "n", leglab = "N genes", 
         mytitle = "Excluding conserved-no_change combination genes") +
  facet_grid(geneset~strain)
## Excluding conserved, no_change genes [EITHER]
heatplot(pdat = longfprops[inhmode!="no_change" & regpat!="conserved", ], numcol = "n", leglab = "N genes", 
         mytitle = "Excluding conserved, no_change genes") +
  facet_grid(geneset~strain)
## Also exclude ambiguous
heatplot(pdat = longfprops[!inhmode%in%c("no_change", "ambiguous") & !regpat%in%c("conserved", "ambiguous"), ], numcol = "n", leglab = "N genes", 
         mytitle = "Excluding conserved, no_change, ambiguous genes") +
  facet_grid(geneset~strain)
invisible(dev.off())

# --- Total props [same plots basically, different legend #] - useful because now normalized within strain to have same span
pdf(file.path(hpdir, paste0(p$baseoutname, "_strainsxgenesets_globalps_heatplots.pdf")), 12, 6)
heatplot(pdat = longfprops, numcol = "prop.allgenes", 
         leglab = "P genes", mytitle = "All categorized genes") +
  facet_grid(geneset~strain)
## **excluding only conserved, no change combination!!
heatplot(pdat = longfprops[!(inhmode=="no_change" & regpat=="conserved"), ], numcol = "prop.allgenes", leglab = "P genes", mytitle = "Excluding conserved-no_change combination genes") +
  facet_grid(geneset~strain)
## Excluding conserved, no_change genes [EITHER]
heatplot(pdat = longfprops[inhmode!="no_change" & regpat!="conserved", ], numcol = "prop.allgenes", leglab = "P genes", 
         mytitle = "Excluding conserved, no_change genes") +
  facet_grid(geneset~strain)
## Also exclude ambiguous
heatplot(pdat = longfprops[!inhmode%in%c("no_change", "ambiguous") & !regpat%in%c("conserved", "ambiguous"), ], numcol = "prop.allgenes", leglab = "P genes", 
         mytitle = "Excluding conserved, no_change, ambiguous genes") +
  facet_grid(geneset~strain)
invisible(dev.off())

# Props w/in categories? NAH

# --- By chromosome
# Get #s broken down by chromosome
cechrs<-c("I", "II", "III", "IV", "V", "X")
nspslistchr<-lapply(1:nrow(gsets), function(x){
  out<-lapply(strains, function(str){
    out1<-lapply(cechrs, function(mychr){
      regxinhnums(onedat = dat[strain == str & eval(parse(text = gsets[x, exprtxt])) & chr==mychr],
                  regpats = regpats, inhmodes = inhmodes, dtadd = data.table(strain = str, geneset = gsets[x,shortname], chr = mychr))
    })
    names(out1)<-cechrs
    return(out1)
  })
  names(out)<-strains
  return(out)
})
names(nspslistchr)<-gsets$shortname

## Save [combining together though it's kinda gnarly..]
write.table(rbindlist(lapply(nspslistchr, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$dts$ns)))))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_counts_bychr.txt")), sep = "\t", row.names = F, quote = F)
write.table(rbindlist(lapply(nspslistchr, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$dts$pAll)))))),
            file.path(p$outdir, paste0(p$baseoutname, "_regxinh_props_ofall_bychr.txt")), sep = "\t", row.names = F, quote = F)


## Longform for plots
longfpropschr<-rbindlist(lapply(nspslistchr, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$longform))))))
longfpropschr[, strain:=factor(strain, levels = strains)]
longfpropschr[, inhmode:=factor(inhmode, levels = inhmodes)]
longfpropschr[, regpat:=factor(regpat, levels = regpats)]
longfpropschr[, chr:=factor(chr, levels = cechrs)]

# Make plots 
## Ns
pdf(file.path(hpdir, paste0(p$baseoutname, "_strainsxgenesets_globalns_bychr_heatplots.pdf")), 12, 10)
invisible(lapply(1:nrow(gsets), function(x){
  print(heatplot(pdat = longfpropschr[geneset==gsets[x, shortname], ], numcol = "n", leglab = "N genes", mytitle = "All categorized genes",
           mysubt = gsets[x, descrip]) +
    facet_grid(strain~chr))
  ## **excluding only conserved, no change combination!!
  print(heatplot(pdat = longfpropschr[geneset==gsets[x, shortname] & !(inhmode=="no_change" & regpat=="conserved"), ], 
           numcol = "n", leglab = "N genes", mytitle = "Excluding conserved-no_change combination genes", mysubt = gsets[x, descrip]) +
    facet_grid(strain~chr))
  return(NULL)
}))
invisible(dev.off())

## Overall Ps
pdf(file.path(hpdir, paste0(p$baseoutname, "_strainsxgenesets_globalps_bychr_heatplots.pdf")), 12, 10)
invisible(lapply(1:nrow(gsets), function(x){
  print(heatplot(pdat = longfpropschr[geneset==gsets[x, shortname], ], numcol = "prop.allgenes", leglab = "P genes", mytitle = "All categorized genes",
                 mysubt = gsets[x, descrip]) +
          facet_grid(strain~chr))
  ## **excluding only conserved, no change combination!!
  print(heatplot(pdat = longfpropschr[geneset==gsets[x, shortname] & !(inhmode=="no_change" & regpat=="conserved"), ], 
                 numcol = "prop.allgenes", leglab = "P genes", mytitle = "Excluding conserved-no_change combination genes", mysubt = gsets[x, descrip]) +
          facet_grid(strain~chr))
  return(NULL)
}))
invisible(dev.off())

# try using ggsave for kicks. WHEN ONLY ONE PLOT THOUGH

#### --- Statistical testing... ####
# Fishers or chi? Fishers. BUT need to strip out no change, conserved stuff? Hmm...like, there will be obligate enrichments due to definitions of underlying stuff?
# May need to test more specific hypothesis/es - enrichment of over/underdominant in compensatory stuff? Maybe look back at Cutter preprint

#### Script completion message & session information ####
cat("....regpatterninhmodeoverlap.R processing complete! Session information:....\n")
sessionInfo() 