#! /usr/bin/env/ Rscript
# Summarize & plot counts of unique eq class alignments - summarize/plot *_eqclasses_genes.txt.gz salmonalleleeqclasses.py across samples
# by Avery Davis Bell, begun 2021.05.04
require(data.table, quietly = T)
require(ggplot2, quietly = T)
require(argparser, quietly = T)

#### Functions ####
freadcollapsegneq<-function(eqgenef, nalnsincl = 0, nthreshs = c(1:5, 10, 20), htitle = ""){
  # Reads one _eqclasses_genes.txt.gz input (output of salmonalleleeqclasses.py); summarizes naln_hapgenespec column
  #   Makes histogram
  # In: eqgenef, path to *_eqclasses_genes.txt.gz input (output of salmonalleleeqclasses.py)
  #     nalnsincl, integer - genes with >= this number of alignments will be summarized in outputs (others excluded early!). Should never be higher than 1!
  #     nthreshs, vector of integers. For each, a column saving the count of genes with THIS MANY OR MORE haplotype & gene-specific alignments will be output, in addition to other columns
  #     htitle, title for histogram that's generated - e.g. including sample information & any relevant description
  # Out: $dt, one-row data.table. Columns - ALL are summaries of the naln_hapgenespec column of input:
  #         n_genes - total # genes included - before thresholding on nalnsincl
  #         n_genes_any<nalnsincl>plusalignments - number genes left after thresholding on nalnsincl
  #         median_uniqalns - median unique (haplotype & gene specific) alignments for genes remaining after thresholding on nalnsincl
  #         mean_uniqalns - mean unique (haplotype & gene specific) alignments for genes remaining after thresholding on nalnsincl
  #         max_uniqalns - maxiumum unique (haplotype & gene specific) alignments for genes remaining after thresholding on nalnsincl
  #         n_genes_0uniqalns - # genes with 0 unique (haplotype & gene specific) alignments of genes remaining after thresholding on nalnsincl
  #         n_genes_<nthreshs>plusuniqalns - # genes with <each of nthreshes values> OR MORE unique (haplotype & gene specific) alignments of genes remaining after thresholding on nalnsincl
  #         genelists - LIST column. List of length one, INNER list is of length nthreshs - one vector of gene_ids matching that make up the n_genes_<nthreshs>plusuniqalns.
  #            names(genelists[[1]]) are genes_<nthreshs>plusuniqalns>
  #            To keep track of which genes met which thresholds!
  #      $hists, ggplot2 histograms for this sample:
  #       $all, all genes included. List of 2 hists. Number of unique alignments plotted. First element of internal list is just raw - not that useful; second has # unique alignments (X) on LOG SCALE - EXCLUDES 0s
  #       $alnmin, only genes with >= nalnsincl alignments included. List of 2 hists. Number of unique alignments plotted. First element is just raw - not that useful; second has # unique alignments (X) on LOG SCALE - EXCLUDES 0s
  
  # Read in
  gndt<-fread(eqgenef, select = c("gene_id", "naln_all", "naln_hapgenespec"))

  # Summary data.table
  summrow<-data.table(gndt[,.(length(gene_id), sum(naln_all>=nalnsincl),
                   median(naln_hapgenespec[naln_all>=nalnsincl]), mean(naln_hapgenespec[naln_all>=nalnsincl]), max(naln_hapgenespec[naln_all>=nalnsincl]),
                   sum(naln_hapgenespec[naln_all>=nalnsincl]==0))],
                   gndt[,lapply(nthreshs, function(x) sum(naln_hapgenespec[naln_all>=nalnsincl]>=x))],
                   list(lapply(nthreshs, function(x) gndt[naln_hapgenespec>=x & naln_all>=nalnsincl, gene_id]))) # internal list of vectors of gene names
  setnames(summrow, c("n_genes", paste0("n_genes_any", nalnsincl,"plusalignments"),
                      "median_uniqalns", "mean_uniqlans", "max_uniqalns",
                      "n_genes_0uniqalns", paste0("n_genes_", nthreshs, "plusuniqalns"),
                      "genelists"))
  names(summrow$genelists[[1]])<-paste0("genes_", nthreshs, "plusuniqalns") # internal list of vectors of gene names! Name each of them
  
  # Histograms 
  h.all<-ggplot(gndt, aes(naln_hapgenespec)) + geom_histogram() +
    xlab("Number of gene & haplotype-specific alignments") + ylab("Number of genes") +
    ggtitle(htitle) +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17))
  h.all.log10<-h.all + scale_x_log10() + xlab("Number of gene & haplotype-specific alignments (log10 scale, 0s excluded!)")
  
  h.alnmin<-ggplot(gndt[naln_all>nalnsincl,], aes(naln_hapgenespec)) + geom_histogram() +
    xlab("Number of gene & haplotype-specific alignments") + ylab("Number of genes") +
    ggtitle(htitle, subtitle = paste("Only genes with at least", nalnsincl, "alignments (non-unique)")) +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17))
  h.alnmin.log10<-h.alnmin + scale_x_log10() + xlab("Number of gene & haplotype-specific alignments (log10 scale, 0s excluded!)")
  
  # Return
  return(list(dt = summrow, hists = list(all = list(h.all, h.all.log10),
                                         alnmin = list(h.alnmin, h.alnmin.log10))))
}

freadsummeq<-function(sampinfo, exampinput, descripcols, outstem, nthreshs = c(1:5, 10, 20)){
  # Reads in & summarizes eq class information from _eqclasses_genes.txt.gz for all samples, focusing on haplotype & gene-specific alignments (naln_hapgenespec column)
  # Also plots histograms for each sample, outputs together into one big old PDF
  # In: sampinfo, data.table with columns SampleID, all those named in descripcols
  #     exampinput, Example filepath to *_eqclasses_genes.txt.gz input (output of salmonalleleeqclasses.py). Where sample name goes, should have _samp_ to be substituted!  If anything else differs across samples (e.g. strain), use a glob (*) character in its place.
  #     descripcols, names of columns of sampinfo to use for histogram titles & to include in output data.table - e.g. Strain, Treatment, or the like
  #     outstem, output filestem (for histgram PDFs - descriptor and suffix will be added)
  #     nthreshs, vector of integers. For each, a column saving the count of genes with THIS MANY OR MORE haplotype & gene-specific alignments will be output, in addition to other columns
  # Out: data.table with one row per sample. Columns SampleID, all those named in descripcols (from sampinfo),
  #       n_genes - total # genes included - before thresholding on 1 alignment
  #         n_genes_any1plusalignments - number genes left after thresholding on 1 alignment
  #         median_uniqalns - median unique (haplotype & gene specific) alignments for genes remaining after thresholding on1 alignment
  #         mean_uniqalns - mean unique (haplotype & gene specific) alignments for genes remaining after thresholding on 1 alignment
  #         max_uniqalns - maxiumum unique (haplotype & gene specific) alignments for genes remaining after thresholding on 1 alignment
  #         n_genes_0uniqalns - # genes with 0 unique (haplotype & gene specific) alignments of genes remaining after thresholding on 1 alignment
  #         n_genes_<1:5, 10, 20>plusuniqalns - # genes with <1-5, 10, 20> OR MORE unique (haplotype & gene specific) alignments of genes remaining after thresholding on nalnsincl 
  #         genelists - LIST column. List of length one, INNER list is of length nthreshs - one vector of gene_ids matching that make up the n_genes_<nthreshs>plusuniqalns.
  #            names(genelists[[1]]) are genes_<nthreshs>plusuniqalns>
  #            To keep track of which genes met which thresholds!
  # Side effects: 4 PDFs containing histograms of number of gene & haplotype specific alignments (x axis) vs. # genes (y axis). Subset and axis as noted in file names:
  #   <outstem>_histograms_genesxnuniqalignments_all.pdf, 
  #   <outstem>_histograms_genesxnuniqalignments_all_logx.pdf
  #   <outstem>_histograms_genesxnuniqalignments_atleast1alignment.pdf
  #   <outstem>_histograms_genesxnuniqalignments_atleast1alignment_logx.pdf
  
  # Get data & hists
  dtshists.l<-lapply(1:nrow(sampinfo), function(x){
    out<-freadcollapsegneq(eqgenef = Sys.glob(gsub("_samp_", sampinfo[x, SampleID], exampinput, fixed = T)),
                           nalnsincl = 1, nthreshs = nthreshs,
                           htitle = paste(unlist(sampinfo[x, c("SampleID", descripcols), with = F]), collapse = " - "))
    out$dt<-data.table(sampinfo[x, c("SampleID", descripcols), with = F], out$dt)
    return(out)
  })
  
  # Format output data
  eqgenects<-rbindlist(lapply(dtshists.l, function(x) x$dt))
  
  # Plot histograms
  ## All genes
  pdf(paste0(outstem, "_histograms_genesxnuniqalignments_all.pdf"), 7.75, 5.5)
  invisible(lapply(dtshists.l, function(x) print(x$hists$all[[1]])))
  invisible(dev.off())
  ## All genes, log scale X excluding those with 0 unique alignments
  pdf(paste0(outstem, "_histograms_genesxnuniqalignments_all_logx.pdf"), 7.75, 5.5)
  invisible(lapply(dtshists.l, function(x) print(x$hists$all[[2]])))
  invisible(dev.off())
  ## Genes with 1 or more alignments
  pdf(paste0(outstem, "_histograms_genesxnuniqalignments_atleast1alignment.pdf"), 7.75, 5.5)
  invisible(lapply(dtshists.l, function(x) print(x$hists$alnmin[[1]])))
  invisible(dev.off())
  ## Genes with 1 or more alignments, log scale X excluding those with 0 unique alignments
  pdf(paste0(outstem, "_histograms_genesxnuniqalignments_atleast1alignment_logx.pdf"), 7.75, 5.5)
  invisible(lapply(dtshists.l, function(x) print(x$hists$alnmin[[2]])))
  invisible(dev.off())
  
  # Return
  return(eqgenects)
}

makedotplot<-function(dt, groupby, plotcol, facetby = NA, mytitle = "", myy = ""){
  # Makes a dotplot of plotcol in dt, grouped and colored by column groupby, optionally faceted by facetby
  # In: dt, data.table containing columns groupby, plotcol, and facetby if facetby isn't NA
  #     groupby, name of column to group the points by in a dotplot (character)
  #     plotcol, name of column with data to actually plot on y axis
  #     facetby, optional (if don't want, leave NA) name of column in dt to facet (separate) the plots by
  #     mytitle, plot title
  #     myy, y axis label
  # Out: ggplot object
  
  plt<-ggplot(dt, aes(eval(as.name(groupby)), eval(as.name(plotcol)))) +
    geom_dotplot(binaxis = "y", stackdir = "center", aes(fill = eval(as.name(groupby)))) +
    ggtitle(mytitle) + xlab(groupby) + ylab(myy) + labs(fill = "") + 
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text=element_text(size = 13),
          strip.text.x = element_text(size = 15)) #l egend.title = element_text(size = 15),
  if(!is.na(facetby)){
    plt<-plt + facet_grid(.~eval(as.name(facetby)))
  }
  return(plt)
}

countovids<-function(ids){
  # Counts the total number of observations in ids; breaks down as to whether present in 1 element of ids, 2 elements...all n elements
  # In: ids, LIST of vectors to overlap
  # Out: one row data.table. Columns:
  #       n_any - # of ids in any one of input lists
  #       n_<1- number of elements of ids> - # of ids in exactly <#> elements of list
  #       n_<2 - one fewer than number of elements of ids>plus - # of ids in <#> or more elements of list
  #       p_<1- number of elements of ids> - proportion of all ids seen in exactly  <#> elements of list
  #       p_<2 - one fewer than number of elements of ids>plus - proportion of all ids in <#> or more elements of list
  ids_ct<-table(unlist(ids))
  out<-setDT(as.list(c(length(ids_ct), sapply(1:length(ids), function(x) sum(ids_ct==x)),
                       sapply(2:(length(ids)-1), function(x) sum(ids_ct>=x)),
                       sapply(1:length(ids), function(x) sum(ids_ct==x)/length(ids_ct)),
                       sapply(2:(length(ids)-1), function(x) sum(ids_ct>=x)/length(ids_ct)))))
  setnames(out, c("n_any", paste("n", 1:length(ids), sep = "_"), 
                  paste0("n_", 2:(length(ids)-1), "plus"),
                  paste("p", 1:length(ids), sep = "_"),
                  paste0("p_", 2:(length(ids)-1), "plus")))
  return(out)
}

getovids<-function(ids, nthresh){
  # Counts the total number of observations in ids; returns those present in nthresh or more
  # In: ids, LIST of vectors to overlap
  #     nthresh, if ID is present this OR MORE times, it is returned as part of 
  # Out: vector of IDs that were present >=nthresh time in input
  ids_ct<-table(unlist(ids))
  return(names(ids_ct)[ids_ct>=nthresh])
}

#### Parse arguments ####
p<-arg_parser("Summarize & plot counts of unique eq class alignments - summarize/plot *_eqclasses_genes.txt.gz salmonalleleeqclasses.py across samples", 
              name = "eqclassalnmtsummary_multsamples.R", hide.opts = TRUE)
p<-add_argument(p, "--sampinfo",
                help = "REQUIRED. Path to sample information file for all samples to process. Must include columns SampleID, any values passed as --splitby, --facetby, --groupby",
                type = "character")
p<-add_argument(p, "--exampinput",
                help = "REQUIRED. Example filepath to *_eqclasses_genes.txt.gz input (output of salmonalleleeqclasses.py).
                Where sample name goes, should have _samp_ to be substituted! If anything else differs across samples (e.g. strain), use a glob (*) character in its place.",
                type = "character")
p<-add_argument(p, "--outstem",
                help = "Output filestem (no directory info)",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Output directory. Defaults to current directory. (If --genelists provided, internal directory is created.)",
                default = NA, type = "character")
p<-add_argument(p, "--groupby",
                help = "REQUIRED. Column name (in --sampinfo file) to group samples by in output plots - samples will be split by this on the same plot.",
                default = "Strain")
p<-add_argument(p, "--splitby",
                help = "Optional. Column name (in --sampinfo file) to split samples by, keeping all plots separate between these groups.",
                default = NA, type = "character")
p<-add_argument(p, "--facetby", 
                help = "Optional. Column name (in --sampinfo file) to facet samples by, keeping them on separate plots but in the same PDF in these groups.",
                default = NA, type = "character")
p<-add_argument(p, "--genelists",
                help = "Optional. If included, informative gene lists PER GROUPBY/SPLITBY/FACETBY CATEGORY (smallest possible given those inputs) are written.
                Format: each gene list to output is :-delimited (number of samples or more of groupby/splitby/facetby category that must have unique alignments:number of unique alignments or more required),
                if more than one desired, these are ,-delimited.  E.g., <--genelists 3:1,1:1> will output 2 gene lists per groupby/splitby/facetby category:
                first is genes where 3 or more samples have 1 or more unique alignments (naln_hapgenespec column of --exampinput >=1), second is where 1 or more samples have 1 or more unique alignments.
                Will be written to <outdir>/genelists/<outstem>_genelist*",
                default = NA, type = "character")
p<-parse_args(p)

# Output directory/file naming
if(is.na(p$outdir)){
  outdir<-getwd()
}else{
  outdir<-p$outdir
  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
  }
}
baseout<-paste0(outdir, "/", p$outstem)

sampinfo<-fread(p$sampinfo)
if(!p$groupby%in%names(sampinfo)){
  stop("The argument --groupby needs to specify a column in --sampinfo")
}
if(!is.na(p$splitby) & p$splitby=="NA"){facetby<-NA}
if(!is.na(p$splitby) & !p$splitby%in%names(sampinfo)){
  stop("The argument --splitby needs to specify a column in --sampinfo if it is included")
}
if(!is.na(p$facetby) & p$facetby=="NA"){p$facetby<-NA}
if(!is.na(p$facetby) & !p$facetby%in%names(sampinfo)){
  stop("The argument --facetby needs to specify a column in --sampinfo if it is included")
}

if(!is.na(p$genelists) & p$genelists!="NA"){
  # Parse early so make sure get n alignments desired propagated through all summaries etc
  genelistspecs<-as.data.table(do.call(rbind, lapply(strsplit(p$genelists, ",")[[1]],
                                                     function(x) t(as.matrix(as.integer(strsplit(x, ":")[[1]]))))))
  setnames(genelistspecs, c("nsamples", "nuniqalns"))
  nthreshs.forout<-sort(unique(c(c(1:5, 10, 20), genelistspecs$nuniqalns)))
}else{
  genelistspecs<-NULL
  nthreshs.forout<-c(1:5, 10, 20)
}

#### Read in & set up data ####
# Sample information
if(!is.na(p$splitby)){ # split sampinfo and lapply everything else by LIST of sample infos
  cats<-unique(sampinfo[,get(p$splitby)])
  sampinfo.l<-lapply(cats, function(x) sampinfo[get(p$splitby)==x])
  names(sampinfo.l)<-cats # names for plots etc
  outnames<-paste0("_", cats) # include UNDERSCORE for using this in output file naming
}else{
  sampinfo.l<-list(sampinfo)
  names(sampinfo.l)<-""
  outnames<-"" # no category split, so no extra info in output filenames
}

# Sample-level information about genes with unique alignments
## Histograms are written out during this!
eqgenects.l<-lapply(1:length(sampinfo.l), function(x){
  freadsummeq(sampinfo = sampinfo.l[[x]], exampinput = p$exampinput, descripcols = na.omit(c(p$groupby, p$facetby)),
              outstem = paste0(baseout, outnames[x]), nthreshs = nthreshs.forout)
})

## Save summary numbers
invisible(lapply(1:length(sampinfo.l), function(x){
  write.table(eqgenects.l[[x]][,.SD, .SDcols = -c(ncol(eqgenects.l[[x]]))],  # LAST column contains gene lists in not nice writeable format - exclude
              paste0(baseout, outnames[x], "_geneuniquealnsummaries.txt"),
              quote = F, row.names = F, sep = "\t")
}))

#### By-sample plots ####
# Dot plots of # genes at various thresholds, split & grouped & faceted as described by input options
nthreshs <- nthreshs.forout # Thresholds used in eqgenects titles
nicetitles<-paste("Genes with at least", nthreshs, "gene/haplotype-unique alignments")
invisible(lapply(1:length(sampinfo.l), function(x){
  pdf(paste0(baseout, outnames[x], "_ngeneuniquealnplots.pdf"), 9, 6)
  sapply(1:length(nthreshs), function(y){
    print(makedotplot(eqgenects.l[[x]], groupby = p$groupby, plotcol = paste0("n_genes_", nthreshs[y], "plusuniqalns"), 
                facetby = p$facetby, mytitle = nicetitles[y], myy = paste("N", nicetitles[y])))
  })
  invisible(dev.off())
}))

# WIDTH BASED on # facets? not currently, probably a good idea longterm


#### Overlap within groupby, optionally facetby: number genes at a threshold SHARED within a GROUP of samples ####
# Get & save tables: number of genes with different thresholds in different N samples
ovtabs.l<-lapply(1:length(sampinfo.l), function(x){
  ## Within facetby, if there is facetby
  if(!is.na(p$facetby)){
    fac<-rbindlist(lapply(sampinfo.l[[x]][,unique(get(p$groupby))], function(y){
        rbindlist(lapply(sampinfo.l[[x]][,unique(get(p$facetby))], function(z){
          rbindlist(lapply(nthreshs, function(a){
            out<-data.table(y, z, a, countovids(lapply(eqgenects.l[[x]][get(p$groupby)==y & get(p$facetby)==z, genelists],
                                                    function(b) b[paste0("genes_", a, "plusuniqalns")])))
            setnames(out, c("y", "z", "a"), c(p$groupby, p$facetby, "unique_aln_threshold"))
          }))
        }))
    }))
    write.table(fac, paste0(baseout, outnames[x], "_uniqealngenesampleoverlap_by", p$facetby, ".txt"),
                quote = F, row.names = F, sep = "\t")
  }else{fac<-NA}
  ## All together within groupby
  allov<-rbindlist(lapply(sampinfo.l[[x]][,unique(get(p$groupby))], function(y){
      rbindlist(lapply(nthreshs, function(a){
        out<-data.table(y, a, countovids(lapply(eqgenects.l[[x]][get(p$groupby)==y, genelists],
                                                   function(b) b[paste0("genes_", a, "plusuniqalns")])))
        setnames(out, c("y", "a"), c(p$groupby, "unique_aln_threshold"))
      }))
  }))
  write.table(allov, paste0(baseout, outnames[x], "_uniqealngenesampleoverlap_allwithin", p$groupby, ".txt"),
              quote = F, row.names = F, sep = "\t")
  return(list(faceted = fac, all = allov))
})

# Plot these? Meh.

##### Save gene lists, if desired ####
if(!is.null(genelistspecs)){
  geneldir<-paste0(outdir, "/genelists")
  if(!dir.exists(geneldir)){dir.create(geneldir, recursive = T)}
  genelout<-paste0(geneldir, "/", p$outstem)
  
  invisible(lapply(1:length(sampinfo.l), function(x){
    ## Within faceting if --facetby provided
    if(!is.na(p$facetby)){
      lapply(sampinfo.l[[x]][,unique(get(p$groupby))], function(y){ # Specific within groupby
        lapply(sampinfo.l[[x]][,unique(get(p$facetby))], function(z){ # Specific within facetby
          lapply(1:nrow(genelistspecs), function(a){ # Specific within genelistspecs
            glist<-getovids(ids = lapply(eqgenects.l[[x]][get(p$groupby)==y & get(p$facetby)==z, genelists],
                                   function(b) b[paste0("genes_", genelistspecs[a, nuniqalns], "plusuniqalns")]),
                            nthresh = genelistspecs[a, nsamples])
            write.table(as.matrix(glist),
                        gzfile(paste0(genelout, outnames[x], "_genelist_", y, "_", z, "_",
                                      genelistspecs[a, nsamples], "plussamples_", 
                                      genelistspecs[a, nuniqalns], "plusuniqlns", ".txt.gz")),
                        quote = F, row.names = F, col.names = F)
          }) 
        })
      })
    }else{ ## not including faceting if --facetbyprovided
      lapply(sampinfo.l[[x]][,unique(get(p$groupby))], function(y){ # Specific within groupby
          lapply(1:nrow(genelistspecs), function(a){ # Specific within genelistspecs
            glist<-getovids(ids = lapply(eqgenects.l[[x]][get(p$groupby)==y, genelists],
                                         function(b) b[paste0("genes_", genelistspecs[a, nuniqalns], "plusuniqalns")]),
                            nthresh = genelistspecs[a, nsamples])
            write.table(as.matrix(glist),
                        gzfile(paste0(genelout, outnames[x], "_genelist_", y, "_",
                                      genelistspecs[a, nsamples], "plussamples_", 
                                      genelistspecs[a, nuniqalns], "plusuniqlns", ".txt.gz")),
                        quote = F, row.names = F, col.names = F)
          }) 
      })
    }
  }))
}

cat("....eqclassalnmtsummary_multsamples.R processing complete! Session information:....\n")
sessionInfo()
