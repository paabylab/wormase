#! /usr/bin/env/ Rscript
# Classify genes as cis, trans, compensatory, etc, and initial summaries. [After working/thinking on this in scripts like ase_de_comparisons_withinstrain.R]
# FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R
# by Avery Davis Bell, begun 2022.07.05
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

regpattern<-function(onedat, alpha = 0.1, lfcthresh = 0.5849625, refstrain = "N2",
                     ci = 0.99){
  # Classifies each gene in onedat as its regulatory pattern [currently see notebook for details?]. *Genes without unique alignments can't be informative and are classified NA
  #     This would be prettier with one complicated logic flow looping through each row...
  # In: onedat, data.table with one row for every gene for which to classify regulatory pattern[s]. Columns: 
  #         unqalnmts.min
  #         log2FoldChange.<suff>, lfcSE.<suff>, padj.<suff>, with <suff> being ASE or ParentVsParent<refstrain>
  #     alpha, p-value threshold for classifying significance - used alone sometimes, in concert with lfctrhesh sometimes (pvalue.ASE and pvalue.ParentVsParent<ref> columns)
  #     lfcthresh, log2FoldChange threshold for classifying significance - log2FoldChange.ASE and log2FoldChange.ParentVsParent<ref> columns
  #     refstrain, name of reference strain in onedat column headers
  #     ci, confidence interval for seeing if ASE, DE overlap - 0.99 constructs 99% confidence intervals
  # Out: input data.table COPY [not efficient but haven't figured out how to get around with how I'm calling this] with new columns:
  #       ase.lowci, ase.highci - lower and upper confidence interval boundaries on ASE log2FoldChange. Confidence level is provided as input.
  #       deparents.lowci, deparents.highci - lower and upper confidence interval boundaries on parent vs. parent DE log2FoldChange. Confidence level is provided as input.
  #       asedeoverlap, logical - do ASE and DE CIs overlap?
  #       regclass, KEY output: what regulatory mode was assigned to this gene? Character one word description of regulatory class. See notebook/documentation for full details. Possible values:
  #         NA - if unqalnmts.min<1, get NA (doesn't classify if no chance this is informative)
  #         conserved - no DE or ASE
  #         ambiguous -  ambiguous category, currently when ASE and DE have overlapping CIs but opposite directions of effect
  #         cis - ASE and DE of similar magnitudes (overlapping CIs, same direction)
  #         trans - DE but no ASE
  #         compensatory - ASE but no DE
  #         compensating - ASE and DE; trans in opposite direction of cis (but cis winning out/ASE 'bigger')
  #         overcompensating - ASE and DE; trans in opposite direction of cis so much so that DE is opposite dir from ASE
  #         enhancing - cis and trans effects acting in same dir (DE and ASE unequal, DE/ASE > 1)
  
  # Subfunction for one row [prettier, clearer way to do this I think...]
  regone<-function(onerow, alpha, lfcthresh, log2ref, seref, pref){
    # Classifies regulatory pattern for one row (gene).
    # In: onerow, one row data.table with ASE, DE information to use to classify regulatory pattern. Columns:
    #       unqalnmts.min, 
    #       columns named log2ref [log2FoldChange for parent vs parent DE], seref [lfcSE for parent vs parent DE], pref [adjusted p-value for parent vs parent]
    #       log2FoldChange.ASE, lfcSE.ASE, padj.ASE, ASE results;
    #       asedeoverlap - have ASE and DE log2FCs been deemed to overlap?
    #     alpha, p-value threshold for classifying significance - used alone sometimes, in concert with lfctrhesh sometimes (pvalue.ASE and pvalue.ParentVsParent<ref> columns)
    #     lfcthresh, log2FoldChange threshold for classifying significance - log2FoldChange.ASE and log2FoldChange.ParentVsParent<ref> columns
    # Out: character one word description of regulatory class. See notebook/documentation for full details. Options:
    #         NA - if unqalnmts.min<1, get NA (doesn't classify if no chance this is informative)
    #         conserved - no DE or ASE
    #         ambiguous - ambiguous category, currently when ASE and DE have overlapping CIs but opposite directions of effect
    #         cis - ASE and DE of similar magnitudes (same direction, overlapping CIs)
    #         trans - DE but no ASE
    #         compensatory - ASE but no DE
    #         compensating - ASE and DE; trans in opposite direction of cis (but cis winning out/ASE 'bigger')
    #         overcompensating - ASE and DE; trans in opposite direction of cis so much so that DE is opposite dir from ASE
    #         enhancing - cis and trans effects acting in same dir (DE and ASE unequal, DE/ASE > 1)

    # NA p-values are a pain. Rather than write to deal with NA p-values at every step, going to separate p-values; put NA p-values as 1
    pase<-ifelse(onerow[,is.na(padj.ASE)], 1, onerow[,padj.ASE])
    pde<-ifelse(onerow[,is.na(get(pref))], 1, onerow[,get(pref)])
    
    if(onerow[,unqalnmts.min]<1){ 
      reg<-as.character(NA) # Ones without any informative reads are obligately NA
    }else{ # at least 1 informative read per sample - go ahead and classify
      if(pase < alpha & onerow[, (abs(log2FoldChange.ASE) > lfcthresh)]){ # ASE at sig & mag thresholds
        if(pde >= alpha){ # definitively no DE
          reg<-"compensatory"
        }else{ # DE - classify if overlaps with ASE
          if(onerow[,asedeoverlap]){
            if(onerow[,get(log2ref)/log2FoldChange.ASE > 0]){ # overlaps and direction is same
              reg<-"cis"
            }else if(onerow[,get(log2ref)/log2FoldChange.ASE <= 0]){ # overlaps but direction is different - ambiguous
              reg<-"ambiguous"
            } # end if/else DE/ASE >0 inside of specification that ASE and DE CIs overlap
          }else{ # ase and de don't overlap
            if(onerow[,get(log2ref)/log2FoldChange.ASE > 1]){
              reg<-"enhancing"
            }else{ # compensating - split further
              if(onerow[,get(log2ref)/log2FoldChange.ASE > 0]){
                reg<-"compensating"
              }else if(onerow[,get(log2ref)/log2FoldChange.ASE < 0]){
                reg<-"overcompensating"
              } # end if/else DE/ASE > 0
            } # end if/else DE/ASE > 1
          } # end if/else asedeoverlap
        } # end if/else DE p >= alpha
        
      }else if(pase >= alpha | onerow[, abs(log2FoldChange.ASE) <= lfcthresh]){ # no ASE at sig & mag thresholds
        if(pde >= alpha | onerow[, abs(get(log2ref)) <= lfcthresh]){ # no ASE, no DE
          reg<-"conserved"
        }else if(pde < alpha & onerow[, abs(get(log2ref)) > lfcthresh]){ # definitive DE [inside no ASE]
          if(pase < alpha){ # ASE passing p but not magnitude cutoff - want to classify as ASE given that there's DE
            if(onerow[,asedeoverlap]){ # ASE and DE overlap: cis
              reg<-"cis"
            }else{ # no asedeoverlap: process what kind of cis/trans this is
              if(onerow[,get(log2ref)/log2FoldChange.ASE > 1]){
                reg<-"enhancing"
              }else{ # compensating - split further
                if(onerow[,get(log2ref)/log2FoldChange.ASE > 0]){
                  reg<-"compensating"
                }else if(onerow[,get(log2ref)/log2FoldChange.ASE < 0]){
                  reg<-"overcompensating"
                } # end if/else DE/ASE > 0
              } # end if/else DE/ASE > 1
            } # end if/else asedeoverlap
          }else{reg<-"trans"} # truly DE, truly no ASE
        } # end if else DE/no DE inside no ASE
      } # end if else ASE/no ASE
    } # end else [1+ unq alignments]
    
    
    return(reg)
  }
  
  # copy
  onedat<-copy(onedat)
  
  # column names
  log2ref<-paste0("log2FoldChange.ParentVsParent", refstrain)
  seref<-paste0("lfcSE.ParentVsParent", refstrain)
  pref<-paste0("padj.ParentVsParent", refstrain)
  
  # Construct CIs, determine overlap
  effp<-1-ci
  qtiles<-c(effp/2, 1-effp/2)
  onedat[unqalnmts.min >=1, `:=`(ase.lowci = log2FoldChange.ASE + qnorm(qtiles[1])*lfcSE.ASE, ase.highci = log2FoldChange.ASE + qnorm(qtiles[2])*lfcSE.ASE,
                                 deparents.lowci = get(log2ref) + qnorm(qtiles[1])*get(seref), deparents.highci = get(log2ref) + qnorm(qtiles[2])*get(seref))]
  onedat[unqalnmts.min >=1, asedeoverlap:=(ase.lowci>=deparents.lowci & ase.lowci<=deparents.highci)|
           (ase.highci>=deparents.lowci & ase.highci<=deparents.highci)|
           (deparents.lowci>=ase.lowci & deparents.lowci<=ase.highci)|
           (deparents.highci>=ase.lowci & deparents.highci<=ase.highci)]
  
  # Classify. Definitely SLOWER to do this way than to update the data.table itself with own logic flows, but that was very ugly and confusing code...
  onedat[, regclass:=regone(.SD, alpha, lfcthresh, log2ref, seref, pref), by = seq_len(nrow(onedat))]
  
 # return [modifying in place not a winner with lapply'ing through several CIs?]
  return(onedat)
}

nsummreg<-function(onedat, descripcols, regpats = c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating"),
                   umbrellacats = NA, refstrain = "N2"){
  # Gets number summary of # genes assigned to each of the regulatory patterns (by regpattern())
  # In: onedat, data.table containing results to summarize into one row (i.e. one strain, one gene set). One row per gene.
  #       Only columns used are: regclass, which must contain values present in regpats (or subset thereof)
  #                               signifAtThresholds.ASE, T or F for gene considered ASE
  #                               signifAtThresholds.ParentVsParent<refstrain>, T or F for gene considered DE
  #     descripcols, data.table that will be appended to beginning of outputs - include any description of note (strain, gene set, etc)
  #     regpats, vector of all regulatory patterns to count up (from regclass input column). This is also ORDER they'll be output in.
  #           must include cis, trans, compensating, compensatory, overcompensating
  #     umbrellacats, optional categories to combine to get counts of larger 'umbrella' categories; NA for none. If provide, format is a list - names of elements are names of umbrella category, value is vector of the regpats that fall into this category
  #         e.g. not conserved category
  #     refstrain, name of reference strain.
  # Out: list of data.tables, $ns, $ps,  $ns contains counts of each category + any umbrella categories, $ps contains proportion of genes in the gene set that fall into this category,
  #         ***AND PS NEWLY CONTAINS proportions of ASE and DE genes that each relevant category is
  #   Columns of one-row data.table ns:
  #     all columns of descripcols
  #     ngenes, number of genes in this gene set
  #     one column per regulatory pattern from regpats. Contains number of genes with this regulatory pattern.
  #     [if umbrellacats provided - column for every umbrella category, sum of regulatory patterns that fall into this]
  #   Columns of ps [ps is in long format for plotting]
  #     all columns of descripcols
  #     ngenes,  number of genes in this gene set
  #     propName, what proportion is in this row? All values are regulatory pattern classifiers (and any umbrella categories if provided)
  #     denom, what is denominator for this row: 'totalngenes' means total genes in input, 
  #             'ASEgenes' means genes with signifAtThresholds.ASE==T,
  #             'DEgenes' means genes with signifAtThresholds.ParentVsParent<refstrain>==T
  #     prop, proportion of genes in this regulatory class (of denominator)
  #     lowci, binomial 95% confidence interval lower bound on proportion
  #     highci, binomial upper 95% confidence interval lower bound on proportion
  
  # Numbers
  ns<-onedat[,.(ngenes = length(regclass), as.data.table(t(as.matrix(table(regclass)))))]
  ## make sure all inheritance mode columns are included (table() will of course not count 0s) 
  invisible(sapply(regpats, function(x){
    if(!x%in%names(ns)){
      ns[,newcol:=0]
      setnames(ns, "newcol", x)
    }
  }))
  setcolorder(ns, regpats, after = 1)
  ## Umbrella categories, if any
  if(!is.na(umbrellacats[1])){
    invisible(lapply(names(umbrellacats), function(x){
      ns[,newcol:=onedat[,sum(regclass%in%umbrellacats[[x]])]]
      setnames(ns, "newcol", x)
    }))
  }
  
  # Proportions of total (long-ways for plotting)
  ps<-rbindlist(lapply(names(ns)[names(ns)!="ngenes"], function(x){
    cis<-data.table(t(ns[, binom.test(get(x), ngenes)$conf.int]))
    setnames(cis, c("lowci", "highci"))
    out<-data.table(ngenes = ns[, ngenes], 
                    propName = x,
                    denom = "totalngenes",
                    prop = ns[,get(x)/ngenes],
                    cis)
    return(out)
  }))
  
  # proportions of ASE, DE genes for relevant categories: added later!. DO make sure ASE/DE called with category for this overlap (not always the case)
  sumase<-onedat[, sum(signifAtThresholds.ASE)]
  testase<-c("cis", "compensatory")
  if("cis-trans opposing"%in%names(umbrellacats)){ ## If have cis-trans opposing umbrella cat, add that [it's maybe most relevant one]
    testase<-c("cis", "compensatory", "UMBRELLA")
  }
  
  ps.ase<-rbindlist(lapply(testase, function(x){
    if(x=="UMBRELLA"){
      regcheck<-umbrellacats[["cis-trans opposing"]]
      regname<-"cis-trans opposing"
    }else{
      regcheck<-c(x)
      regname<-x
    }
    
    myn<-onedat[, sum(signifAtThresholds.ASE==T & regclass%in%regcheck, na.rm = T)]
    cis<-binom.test(myn, sumase)$conf.int
    out<-data.table(ngenes = myn,
                    propName = regname,
                    denom = "ASEgenes",
                    prop = myn/sumase,
                    lowci = cis[1],
                    highci = cis[2])
    return(out)
  }))
  
  sumde<-onedat[,sum(get(paste0("signifAtThresholds.ParentVsParent", refstrain)))]
  testde<-c("cis", "trans")
  ps.de<-rbindlist(lapply(testde, function(x){
    myn<-onedat[,sum(get(paste0("signifAtThresholds.ParentVsParent", refstrain))==T & regclass==x, na.rm = T)]
    cis<-binom.test(myn, sumde)$conf.int
    out<-data.table(ngenes = myn,
                    propName = x, 
                    denom = "DEgenes",
                    prop = myn/sumde,
                    lowci = cis[1],
                    highci = cis[2])
  }))
  ps<-rbind(ps, ps.ase, ps.de)
  
  # Final format & return
  ns<-data.table(descripcols, ns)
  ps<-data.table(descripcols, ps)
  return(list(ns = ns, ps = ps))
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
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  return(plt)
}

regpatscat<-function(pdat, mycolors, colcol = "regclass", xcol = "log2FoldChange.ParentVsParentN2", ycol = "log2FoldChange.ASE",
                      mytitle = "", mysubt = "", xlabel = bquote(log[2]~"(Fold Change: alt strain vs. ref strain) (DE)"),
                     ylabel = bquote(log[2]~"(Fold Change: alt allele vs. ref allele) (ASE)"),
                      labcolby = "", outlineby = NA, laboutlineby = NA, axmin = NA, axmax = NA,
                     annotlines = T){
  # Makes a scatter plot of DE vs ASE colored by regulatory pattern (more flexible than this in fact)
  # In: pdat, data with one row per gene/point to plot. Columns xcol value, ycol value, colcol value 
  #     mycolors, colors for bars: vector of length of categories in data; also used as ORDER for bars!! Names are categories, values are colors
  #         *also used to level factor
  #         Must have all values that are in colcol column of pdat
  #     colcol, name of column with data to color by
  #     xcol, character name of column for x axis values, typically DE
  #     ycol, character name of column for y axis values, typically ASE
  #     mytitle, title for plot
  #     mysubt, subtitle for plot
  #     xlabel, x axis label for plot
  #     ylabel, y axis label for plot
  #     labcolby, title of color legend label
  #     outlineby, optional column of res containing TRUE and FALSE as only values to have TRUE values outlined in black (NA to not separate points by outline)
  #     laboutlineby, What to label outline legend. If NA and outlineby provided, outlineby used
  #     axmin, minimum value for both X and y axes. If NA, this will be minimum value of data in either xcol or ycol.
  #     axmax, maximum value for both x and y axes, if NA, this will be maximum value of data in either xcol or ycol
  #     annotlines, T or F: draw x = 0, y = 0 lines on plot?
  # Out: ggplot
  
  # Format data
  pdat<-copy(pdat) # don't want to mess with original data
  ## Make color col factor of desired order
  pdat[, colbycol:=factor(get(colcol), levels = names(mycolors))]
  ## Axis values
  if(is.na(axmin)){
    axmin<-min(c(pdat[, min(get(xcol), na.rm = T)], pdat[, min(get(ycol), na.rm = T)]))
  }
  if(is.na(axmax)){
    axmax<-max(c(pdat[, max(get(xcol), na.rm = T)], pdat[, max(get(ycol), na.rm = T)]))
  }
  # nexcl<-plt.data[,sum((get(xcol) < axmin | get(ycol) < axmin | 
  #                         get(xcol) > axmax | get(ycol) > axmax), na.rm = T)]
  
  # Make plot
  plt<-ggplot(pdat, aes(eval(as.name(xcol)), eval(as.name(ycol)))) + geom_point(aes(fill = colbycol), pch = 21, alpha = 0.7, stroke = 0.001) +
    scale_fill_manual(values = mycolors) + labs(fill = labcolby) +
    xlim(c(axmin, axmax)) + ylim(c(axmin, axmax)) +
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
      xlim(c(axmin, axmax)) + ylim(c(axmin, axmax)) +
      guides(color = guide_legend(override.aes = list(size = 5)))
  }
  
  if(annotlines){
    plt<-plt + geom_vline(xintercept = 0, col = "darkgray", lty = "dashed") + geom_hline(yintercept = 0, col = "darkgray", lty = "dashed")
  }
  
  return(plt)
}

divplot<-function(onedat, xcol = "nvars.thou", xlabel = "SNVs/INDELs vs. N2 (thousands)", 
                  ylabel = "", mytitle = "", mysubt = "", mysize = 1){
  # One off just to make this plot
  plt<-ggplot(onedat, aes(eval(as.name(xcol)), prop)) +
    geom_pointrange(aes(ymin = lowci, ymax = highci, color = strain), size = mysize) +
    xlab(xlabel) + ylab(ylabel) + ggtitle(mytitle, subtitle = mysubt) + labs(color = "Strain") + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 9))
  
  return(plt)
}

#### Arguments & inputs ####
# --- Command line arguments
p<-arg_parser("Classify genes as cis, trans, compensatory, etc, and initial summaries. [After working/thinking on this in scripts like ase_de_comparisons_withinstrain.R]
              FROM annotated results generated with ase_de_annotategenes_deseq2_fromemaseout.R", 
              name = "ase_de_cistransclassifications.R", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--input",
                help = "Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
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
                Must match how strains are named in input filenames.",
                type = "character")
p<-add_argument(p, "--refstrain",
                help = "Name of reference strain, matching name used in columns of input. Also used for plot labelling.",
                default = "N2")

# Thresholding for cis/trans etc
p<-add_argument(p, "--alpha",
                help = "P-value threshold for considering gene ASE or DE for *adjusted p-value*. Combined with magnitude threshold (--lfcthresh) for first identification, used on its own if other expression difference already established.", 
                default = 0.1)
p<-add_argument(p, "--lfcthresh",
                help = "log2FoldChange threshold for considering gene ASE or DE. Combined with significance threshold (--alpha) to establish first expression diff.", 
                default = 0.5849625)
p<-add_argument(p, "--asedecis",
                help = "Path to file specifying confidence intervals to use to determine if ASE and DE magnitudes are different. Must have one entry for each strain for each desired threshold; thresholds can be consistent or different across strains. Columns:
                strain (name of strain), CI (confidence interval - 0.95 means 95%, etc), descrip (description of CI for output files etc. Each strain should have a row for each description. Order of descriptions in this file is considered factor level order!)",
                type = "character")

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

# Output directory
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}

# --- Read in data, combining strains and adding strain column
cat("....Reading in data....\n")
dat<-freadcombstrains(exampinput = p$input, strains)

#### Do classifications ####
cat("....Classifying regulatory patterns; generating numerical summaries ....\n")
# Possible regulatory patterns, in order for factor etc. ***CHANGE IF CHANGE THIS.****
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")
umbrellacats<-list("cis-trans opposing" = c("compensating", "compensatory", "overcompensating"),
                   "not conserved" = c("cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")) # Any that want to see COMBINED!
#     NB compensatory, compensating, overcompensating are all cis-trans opposing categories. 

# --- Confidence intervals
confs<-fread(p$asedecis, header = T)
confs.ord<-as.character(confs[,unique(descrip)])

# --- Get each gene's reg pattern several ways. Currently NOT annotating in data.table with CI etc...
regdat<-lapply(confs.ord, function(x){
  out<-lapply(strains, function(str){ # have to break down by strain because CI numbers may differ across strains
    regpattern(onedat = dat[strain==str, ], alpha = p$alpha, lfcthresh = p$lfcthresh, 
               refstrain = p$refstrain, ci = confs[strain==str & descrip==x, CI])
  })
  names(out)<-strains
  return(out)
})
names(regdat)<-confs.ord

# --- Save out. NB saving for *only genes with non-NA classification for space - NAs auto added when merge back in.* Also saving only NEW columns.
invisible(lapply(confs.ord, function(x){
  write.table(rbindlist(lapply(regdat[[x]], function(y) y[unqalnmts.min>=1, .(strain, gene_id, regclass)])), 
              gzfile(file.path(p$outdir, paste0(p$baseoutname, "_regpattern_per1plusunqalnggenestrain_", x, "confdeaseoverlap.txt.gz"))),
              sep = "\t", quote = F, row.names = F)
}))

# --- Numerical summaries
# Gene sets. **ADD HERE IF WANT TO TRY A NEW ONE** [will be this x strain x CI. these are the informative gene sets used previously]
## Set up the categories to intersect
infos<-c("Informative, 1+ unique alignments", "Informative, 2+ unique alignments", "Informative, 5+ unique alignments") # 3 separate informative gene sets
infos.short<-c("1plusUnqAlns", "2plusUnqAlns", "5plusUnqAlns")
infos.bool<-c("unqalnmts.min >= 1", "unqalnmts.min >= 2", "unqalnmts.min >= 5")

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

# Numbers per set separately
nspslist<-lapply(confs.ord, function(x){
  lapply(1:nrow(gsets), function(y){
    lapply(strains, function(str){
      nsummreg(onedat = regdat[[x]][[str]][eval(parse(text = gsets[y, exprtxt])),], 
               descripcols = data.table(strain = str, CI = x, geneset = gsets[y, shortname]),
               regpats = regpats, umbrellacats = umbrellacats)
    })
  })
})

# Combine numbers per set
nregpat<-rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$ns))))))
pregpat<-rbindlist(lapply(nspslist, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$ps))))))
  
# Save
write.table(nregpat, file.path(p$outdir, paste0(p$baseoutname, "_numbers_regpatterns.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(pregpat, file.path(p$outdir, paste0(p$baseoutname, "_proportion_regpatterns.txt")),
            sep = "\t", quote = F, row.names = F)

#### Plot set up related ####
# pick colors used for each reg pattern for all plots. All the ones with cis-trans opposing in diff shades of same color
regcols<-rep(NA, length(regpats))
names(regcols)<-regpats
regcols["conserved"]<-"lightgray"
regcols["ambiguous"]<-"darkgray"
regcols[c("compensating", "compensatory", "overcompensating")]<-brewer.pal(9, "Blues")[c(3, 6, 9)] # wider dynamic range of color than doing brewer.pal(3, "Blues")
regcols[is.na(regcols)]<-brewer.pal(4, "Set2")[c(1, 2, 4)]

# Get in divergence info (so can order strains this way even before plotting divergence)
ndiv<-fread(p$varsvsref)
ndiv[,`:=`(nvars.thou = nvars/1e03, nvars.perkb = nvars/(p$genomebp/1e03))] 
setkey(ndiv, strain)

# Make geneset a factor in the order I want it
pregpat[,geneset:=factor(geneset, levels = gsets$shortname)]

#### Stacked bar plots showing breakdown ####
cat("....Making barplot summaries of classification proportions....\n")
barpltdir<-file.path(p$outdir, "stackedbarplots")
if(!dir.exists(barpltdir)){dir.create(barpltdir)}

# All together on one page [useful for seeing all at once, notebook]
pdf(file.path(barpltdir, paste0(p$baseoutname, "_regpatterns_proportions_allsets.pdf")), 16, 16)
stackedbar(pregpat[propName%in%regpats & denom=="totalngenes", ], mycolors = regcols, xorder = ndiv[order(nvars), strain]) + facet_grid(geneset~CI)
invisible(dev.off())

# Each CI has its own PDF, each gene set has its own page [useful for presentations, zooming in]
invisible(lapply(confs.ord, function(x){
  pdf(file.path(barpltdir, paste0(p$baseoutname, "_regpatterns_proportions_", x, "confdeaseoverlap.pdf")), 9, 6)
  lapply(1:nrow(gsets), function(y){
    plt<-stackedbar(pregpat[CI==x & geneset==gsets[y, shortname] & propName%in%regpats & denom=="totalngenes", ], mycolors = regcols, xorder = ndiv[order(nvars), strain],
                    mysubt = gsets[y, descrip], mytitle = paste("Confidence interval:", x))
    print(plt)
    return(NULL)
  })
  invisible(dev.off())
}))

#### Scatter plots showing all genes ####
cat("....Plotting DE vs ASE scatter plots....\n")
scatpltdir<-file.path(p$outdir, "devsaseplots")
if(!dir.exists(scatpltdir)){dir.create(scatpltdir)}

# Faceted by strain; each CI has its own PDF and each gene set has its own page
## No outline
invisible(lapply(confs.ord, function(x){
  pdf(file.path(scatpltdir, paste0(p$baseoutname, "_asevsdeplot_regpatterns_", x, "confdeaseoverlap.pdf")), 16, 15) # **PDF dimensions only tested with 4 strains
  lapply(1:nrow(gsets), function(y){
   # set up data with strain facets ordered the way I want
    pdat <- rbindlist(lapply(regdat[[x]], function(z) z[eval(parse(text = gsets[y, exprtxt])), ]))
    pdat[, strain:=factor(strain, levels = ndiv[order(nvars), strain])]
    # make plot
    plt<-regpatscat(pdat = pdat, mycolors = regcols, mytitle = paste("Confidence interval:", x), mysubt = gsets[y, descrip]) +
      facet_wrap(~strain)
    print(plt)
    return(NULL)
  })
  invisible(dev.off())
  return(NULL)
}))
## Hypdiv outline [irrelevant for hypdiv-excluded genesets] 
invisible(lapply(confs.ord, function(x){
  pdf(file.path(scatpltdir, paste0(p$baseoutname, "_asevsdeplot_regpatterns_", x, "confdeaseoverlap_hypdivoutline.pdf")), 16, 15) # **PDF dimensions only tested with 4 strains
  lapply(1:nrow(gsets), function(y){
    # set up data with strain facets ordered the way I want
    pdat <- rbindlist(lapply(regdat[[x]], function(z) z[eval(parse(text = gsets[y, exprtxt])), ]))
    pdat[, strain:=factor(strain, levels = ndiv[order(nvars), strain])]
    # make plot
    plt<-regpatscat(pdat = pdat, mycolors = regcols, mytitle = paste("Confidence interval:", x), mysubt = gsets[y, descrip],
                    outline = "hypdiv", laboutlineby = "Hyperdivergent") +
      facet_wrap(~strain)
    print(plt)
    return(NULL)
  })
  invisible(dev.off())
  return(NULL)
}))

# Each strain has its own page [helpful for zooming in; variable axes]
#     NOT currently doing this - a LOT of nesting etc. Can make individual ones as one off if/when actually desired.

#### Plot proportions vs. divergence ####
cat("....Plotting vs. divergence....\n")
divpltdir<-file.path(p$outdir, "propvdivergenceplots")
if(!dir.exists(divpltdir)){dir.create(divpltdir)}

setkey(pregpat, strain)
pregpat<-ndiv[pregpat]
# Many on same page
## One page per proportion, grid of CIs vs. gene sets
pdf(file.path(divpltdir, paste0(p$baseoutname, "_regpatternsvsdivergence_allsets.pdf")), 16, 16)
invisible(lapply(c(regpats, names(umbrellacats)), function(x){
  print(
  divplot(pregpat[propName==x & denom=="totalngenes",], xcol = "nvars.thou", xlabel = paste("SNVs/INDELs vs.", p$refstrain, "(thousands)"), 
          ylabel = paste("Proportion genes classified as", x), mytitle = x, mysize = 0.7) + 
    facet_grid(geneset~CI)
  )
  return(NULL)
}))
invisible(dev.off())
    # Use colorblind friendly colors!!! [? not sure if want diff from existing plots though]

## any combo with props on same page...? might be nice to see?
##  # 9 regulatory patterns/umbrella cats; can be broken into 7 + 2 #

# separate CIs (PDFs) and gene sets (pages); 9 regulatory patterns/umbrella cats on same page
invisible(lapply(confs.ord, function(x){
  pdf(file.path(divpltdir, paste0(p$baseoutname, "_regpatternsvsdivergence_", x, "confdeaseoverlap.pdf")), 16, 10)
  invisible(lapply(1:nrow(gsets), function(y){
    # Order proportion names for facets
    pdat<-pregpat[CI==x & geneset==gsets[y, shortname] & denom=="totalngenes",]
    pdat[,propName:=factor(propName, levels = c(regpats, names(umbrellacats)))]
    pdat[, strain:=factor(strain, levels = ndiv[order(nvars), strain])]
    # Plot
    plt<-divplot(pdat, xcol = "nvars.thou", xlabel = paste("SNVs/INDELs vs.", p$refstrain, "(thousands)"), 
                 mytitle = paste("Confidence interval for ASE/DE overlap:", x), mysubt = gsets[y, descrip], mysize = 0.7) +
      facet_wrap(~propName, scales = "free_y")
    print(plt)
    return(NULL)
  }))
  invisible(dev.off())
}))

# --- NEW: plot proportion of ASE/DE genes that are various categories
## ASE
invisible(lapply(confs.ord, function(x){
  pdf(file.path(divpltdir, paste0(p$baseoutname, "_regpatternsvsdivergence_", x, "confdeaseoverlap_PropOfASEGenes.pdf")), 10, 4)
  invisible(lapply(1:nrow(gsets), function(y){
    # Order proportion names for facets
    pdat<-pregpat[CI==x & geneset==gsets[y, shortname] & denom=="ASEgenes",]
    pdat[,propName:=factor(propName, levels = c(regpats, names(umbrellacats)))]
    pdat[, strain:=factor(strain, levels = ndiv[order(nvars), strain])]
    # Plot
    plt<-divplot(pdat, xcol = "nvars.thou", xlabel = paste("SNVs/INDELs vs.", p$refstrain, "(thousands)"), 
                 ylabel = "Proportion of ASE genes in this category",
                 mytitle = paste("Confidence interval for ASE/DE overlap:", x), mysubt = gsets[y, descrip], mysize = 0.7) +
      facet_wrap(~propName, scales = "free_y")
    print(plt)
    return(NULL)
  }))
  invisible(dev.off())
}))
## DE
invisible(lapply(confs.ord, function(x){
  pdf(file.path(divpltdir, paste0(p$baseoutname, "_regpatternsvsdivergence_", x, "confdeaseoverlap_PropOfDEGenes.pdf")), 8, 4)
  invisible(lapply(1:nrow(gsets), function(y){
    # Order proportion names for facets
    pdat<-pregpat[CI==x & geneset==gsets[y, shortname] & denom=="DEgenes",]
    pdat[,propName:=factor(propName, levels = c(regpats, names(umbrellacats)))]
    pdat[, strain:=factor(strain, levels = ndiv[order(nvars), strain])]
    # Plot
    plt<-divplot(pdat, xcol = "nvars.thou", xlabel = paste("SNVs/INDELs vs.", p$refstrain, "(thousands)"), 
                 ylabel = "Proportion of DE genes in this category",
                 mytitle = paste("Confidence interval for ASE/DE overlap:", x), mysubt = gsets[y, descrip], mysize = 0.7) +
      facet_wrap(~propName, scales = "free_y")
    print(plt)
    return(NULL)
  }))
  invisible(dev.off())
}))

#### Script completion message & session information ####
cat("....ase_de_cistransclassifications.R processing complete! Session information:....\n")
sessionInfo()
