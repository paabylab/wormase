#! /usr/bin/env/ Rscript
# Get nucleotide diversity in C ele genome via PopGenome for
#  annotate with gene information from GFF

# by Avery Davis Bell, begun 2022.09.06
require(argparser)
require(data.table)
require(PopGenome)
require(ggplot2)

#### Functions ####
geneidfromfeature<-function(gene_ids_plus, spltstring = "gene-->gene_id "){
  # Gets the gene_id from a long string returned by get.feature.names for wormbase GTF
  # In: gene_ids_plus, vector of results from calling get.feature.names for genes from wormbase GTF. example entry: 
  #           "gene-->gene_id WBGene00022278; gene_source WormBase; gene_biotype protein_coding;;transcript-->gene_id WBGene00022278; transcript_id Y74C9A.4b.2; gene_source WormBase; gene_biotype protein_coding; transcript_source WormBase; transcript_biotype protein_coding;"
  # Out: character vector of same length of gene_ids_plus, with just the gene_id directly following the gene-->gene_id in input
  return(
    sapply(strsplit(gene_ids_plus, spltstring, fixed = T), function(x){
      strsplit(x[2], ";")[[1]][1]
    })
  )
}

posvsvalplot<-function(dat, chrs, poscol = "pos.mid", ycol = "pSegSites_ofTotalN", stackfacet = "population",
                       myylab = "", mytitle = "", mysubt = "", myscales = "free", addsmooth = T){
  # Makes your basic plot with genomic position on X axis, value of something of interest on Y. 
  # In: dat, data.table of data to plot. Must have columns chr; named with values of poscol, ycol, stackfacet
  #     chrs, chromosome information. Chromosomes will be plotted in this order!! Columns chr, start, end. (used for ordering and to make sure plots extend to beginnings, ends of chrs)
  #     poscol, name of column with basepair genomic position to plot
  #     ycol, name of column with values to plot on y axis
  #     stackfacet, name of column for vertical facets/rows (plot columns will be chromosome)
  #     myylab, y axis label
  #     mytitle, plot title
  #     mysubt, plot subtitle
  #     myscales, passed to scales= of faceting. **MUST be "free" or "free_x" - x axis needs to vary
  #     addsmooth, whether to add loess-smoothed line on top
  
  # Set up data (chromosome stuff)
  pdat<-copy(dat[,.(chr, get(poscol), get(ycol), get(stackfacet))]) # don't want to modify original dt
  setnames(pdat, c("chr", "pltpos", "plty", "pltstack"))
  
  cdat<-rbindlist(lapply(1:pdat[, length(unique(pltstack))], function(x){
    rbind(chrs[, .(chr, pltpos = start, plty = as.numeric(NA), pltstack = pdat[,unique(pltstack)][x])],
          chrs[, .(chr, pltpos = end, plty = as.numeric(NA),  pltstack = pdat[,unique(pltstack)][x])])
    }))
  pdat<-rbind(pdat, cdat)
  
  pdat[,chr:=factor(chr, levels = chrs$chr)]
  
  
  # Make plot
  plt<-ggplot(pdat, aes(pltpos/1e06, plty)) + geom_point(alpha = 0.4, color = "gray") +
    facet_grid(pltstack~chr, scales = myscales) +
    ggtitle(mytitle, subtitle = mysubt) + xlab("Genomic position (Mb)") + ylab(myylab) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 13),
                       strip.text.y = element_text(size = 14))
  
  # add line smoothing if desired
  if(addsmooth){
    plt<-plt + geom_smooth(method = "loess")
  }
  
  return(plt)
}

#### Arguments and input ####
p<-arg_parser("Get nucleotide diversity per gene in C ele genome via PopGenome; annotate with gene information from GFF
              and hyperdivergent haplotype status across strains of interest", 
              name = "allgenes_nucdiv_cendr.R ", hide.opts = TRUE)

p<-add_argument(p, "--vcf",
                help = "Path to VCF to process. Must have matching tabix index in same directory.
                Note: ***popgenome/called packages for VCF processing do NOT do well with the ~ shortcut in filenames - need to write it all out!",
                type = "character")
p<-add_argument(p, "--gtf",
                help = "Path to *unzipped* GTF matching input VCF.",
                type = "character")
p<-add_argument(p, "--outdir",
                help = "Output directory. Currently tmp files written here too.",
                default = "out")
p<-add_argument(p, "--out",
                help = "Output file prefix. Data, plots will be written to this _<descriptive suffixes>",
                default = "out")
p<-add_argument(p, "--pops",
                help = "File containing sub-populations for which to compute statistics. **ALL samples will be added to a population called 'all' automatically.
                Format: 2 columns titled PopName and ID. PopName contains the name for the population (same for all IDs for that population), ID the sample ID matching VCF. As many rows per population as samples in that population.",
                type = "character")
p<-add_argument(p, "--chrfile",
                help = "File with columns chr, start, end. One line per chromosome to process; start and end should match GTF/GFF for these data",
                type = "character")
p<-add_argument(p, "--reffa",
                help = "Example filepath to reference genome fasta for one chromosome matching GTF, VCF build! Chromosome in name must be as provided in --chrfile.
                Where chromosome ID is in filename, replace with -CHR-, e.g. c_elegans.PRJNA13758.WS276.genomic.-CHR-.fa.
                NOTE: need to be unzipped for PopGenome!",
                type = "character")
p<-add_argument(p, "--windowbp",
                help = "Size of window to bin genome into for sliding window analyses. Non-overlapping windows of this size will be considered.",
                default = 10e03)

cat("....Parsing arguments....\n")
p<-parse_args(p)

#### PopGenome processing and data formatting ####
# Housekeeping & reading in files used for all chrs
setwd(p$outdir) # Needed for tmp stuff
chrs<-fread(p$chrfile, header = T)
# GTF - used for CDS-level processing
cdsgtf<-fread(p$gtf, header = F, skip = "#")[V3=="CDS", ]
cdsgtf_geneids<-unlist(strsplit(geneidfromfeature(cdsgtf$V9, spltstring = "gene_id \""), 
                                split = "\"", fixed = T)) # need to strip escaped quotes before & after

## non-all population set up
popsdt<-fread(p$pops, header = T)
setkey(popsdt, PopName)
pops<-lapply(popsdt[, PopName, by = PopName][, PopName], function(x) popsdt[x, ID])
names(pops)<-popsdt[, PopName, by = PopName][, PopName]

# Set up chromosome
popgenout<-lapply(1:nrow(chrs), function(chrind){
  mychr<-chrs[chrind, chr]
  mystart<-chrs[chrind, start]
  myend<-chrs[chrind, end]
  cat(paste0("....PopGenome formatting chromosome ", mychr, "....\n"))

# ---- Read in for all samples, set appropriate populations and variant classifications ----
#       ***popgenome/called packages do NOT do well with the ~ shortcut in filenames - need to write it all out!
  GENOME.class<-readVCF(p$vcf, numcols = 1e03, tid = mychr, frompos = mystart, topos = myend,
                        gffpath = p$gtf, include.unknown = T, approx = F)
  
  # Set populations
  myinds<-get.individuals(GENOME.class)[[1]][seq(1, length(get.individuals(GENOME.class)[[1]]), by = 2)] # each second entry is just '.2' of ID - second haplotype
  GENOME.class<-set.populations(GENOME.class, append(list(all = myinds), pops), diploid = T)

  # Get syn/nonsyn
  chrfasta<-gsub("-CHR-", mychr, p$reffa, fixed = T)
  GENOME.class<-set.synnonsyn(GENOME.class, ref.chr = chrfasta) # GENOME.class@region.data@synonymous is 0 for NON-synonymous, 1 for synonymous
  
# ---- Split into various genome intervals (genes, windows, ??coding seqs[but maybe that's WITHIN the other 2 instead]) ----
  cat("    ....splitting into genes....\n")
  # Genes - used only for getting # genic SNPs
  genes <- splitting.data(GENOME.class, subsites="gene")
  ## Get gene_ids, positions for annotating output
  gene_ids_plus<-get.feature.names(genes, p$gtf, chr = mychr) # this pulls whole string for any matching regions.
  gene_info<-data.table(gene_id = geneidfromfeature(gene_ids_plus),
                        rbindlist(lapply(strsplit(genes@region.names, " - "), function(x){
                          pos<-as.numeric(x)
                          return(data.table(pos.start = pos[1], pos.end = pos[2], pos.mid = mean(pos))) 
                        })),
                        genelength = genes@n.sites)
  
  # CDS - used for doing coding, syn, nonsyn; have to do within gene because any CDS overlapping a gene gets assigned to that gene
  gene_ids<-geneidfromfeature(gene_ids_plus) # this is within-chromosome
  tmpdir<-file.path(p$outdir, "tmpcdsgtfs")
  dir.create(tmpdir)
  
  cat("    ....splitting into CDS within-gene: this can take some time!....\n")
  cds.pergenesplit<-lapply(gene_ids, function(x){ # This takes quite a while for each chromosome (an hour-ish for chrII).
    
    if(x%in%cdsgtf_geneids){ # Only do if this gene has CDS - otherwise popgenome can't handle
      # Make tmp GTFs. Lots of I/O, but saves loading, and deals with letting you know if it's protein coding
      mycdsdat<-cdsgtf[cdsgtf_geneids==x,]
      ## duplicate last CDS - last exon (by position) is getting omitted for unknown internal reasons
      mycdsdat<-mycdsdat[order(V5),]
      mycdsdat<-rbind(mycdsdat, mycdsdat[nrow(mycdsdat), ])
      ## save
      mygtf<-file.path(tmpdir, paste0(x, "_cds.gtf"))
      write.table(mycdsdat, mygtf, sep = "\t", col.names = F, row.names = F, quote = F)
      
      # Pull gene with popgenome 
      capture.output(onegene<-split_data_into_GFF_attributes(GENOME.class, gff.file = mygtf, 
                                                             chr = mychr, attribute = x))
      
      # clean up
      file.remove(mygtf)
    }else{ # this gene has no CDS
      onegene<-NULL
    }
    return(onegene)
  })
  
  gene_info[, cdsdata:=!sapply(cds.pergenesplit, is.null)] # save if has CDS info or not
  gene_info[, cdslength:=sapply(cds.pergenesplit, function(x){
    if(is.null(x)){
      return(as.numeric(NA))
    }else{
      return(x@n.sites)
    }
  })] # save CDS length 
  # Would ideally Combine into one object for more easily doing stats using concatenate.classes, but it seems this hasn't been well implemented
  #   tmp<-concatenate.classes(cds.pergenesplit[gene_info$cdsdata]) 
  unlink(tmpdir, recursive = T)
  
  # Windows [update: only doing OVERALL sites here, not subsites]
  cat("    ....splitting into genomic windows....\n")
  gwindows<-sliding.window.transform(GENOME.class, width = p$windowbp, jump = p$windowbp, type = 2)

   ## Get windows as numeric positions for annotating output
  gwindows_pos<-rbindlist(lapply(strsplit(gwindows@region.names, " "), function(x){
    pos<-as.numeric(x[c(1, 3)])
    return(data.table(pos.start = pos[1], pos.end = pos[2], pos.mid = mean(pos))) # mid likely most useful for plotting
  }))

# ---- Compute statistics of interest ---- 
  cat(paste0("...Computing diversity and neutrality stats for all populations, site categories for chromosome ", mychr, "....\n"))
  
  ## Genes -> only doing all sites within a gene (CDS was wonky)
  genes<-diversity.stats(genes, pi = T)
  genes<-neutrality.stats(genes, FAST = T)
  
  ## Coding regions per gene -> coding, syn, nonsyn
  sitesvals<-c("coding", "syn", "nonsyn")
  ### Narrow to non-null first: these are ones with gene_info[cdsdata==T]
  cds.pergenesplit<-cds.pergenesplit[!sapply(cds.pergenesplit, is.null)]
  cds.pergenesplit<-lapply(sitesvals, function(x){
    lapply(cds.pergenesplit, function(y){
      capture.output(myret<-diversity.stats(y, subsites = x, pi = T))
      capture.output(myret<-neutrality.stats(myret, subsites = x, FAST = T))
      return(myret)
    })
  })
  names(cds.pergenesplit)<-sitesvals
  
  ## Windows -> ONLY DOING ALL SITES
  gwindows<-diversity.stats(gwindows, pi = T)
  gwindows<-neutrality.stats(gwindows, FAST = T)
  
# ---- Format sensibly for combining across chromosomes ---- 
cat(paste0("...Formatting diversity and neutrality stats for all populations, site categories for chromosome ", mychr, "....\n"))
  # Genes. LOTS of rows!
  pergenepop<-rbindlist(lapply(1:length(genes@populations), function(x){ # x is population index
    rbind(
      # Do genes/overall
      data.table(chr = mychr, gene_info,
                 population = names(genes@populations)[x],
                 sites = "all",
                 nSegregatingSites = genes@n.segregating.sites[, x], 
                 pSegSites = genes@n.segregating.sites[, x]/gene_info[, genelength], 
                 pi_raw = genes@Pi[, x],
                 pi_persite = genes@Pi[, x]/gene_info[, genelength], 
                 tajimasD = genes@Tajima.D[, x]
                 # ***ANYTHING ELSE??***
      ),
      # Subsites within coding regions
      rbindlist(lapply(sitesvals, function(y){ # y is which site subset (char)
        rbindlist(
          lapply(1:length(cds.pergenesplit[[y]]), function(z){ # z is index of gene data
            mygenedat<-gene_info[cdsdata==T, ][z, ]
            data.table(chr = mychr, mygenedat,
                       population = names(genes@populations)[x],
                       sites = y,
                       nSegregatingSites = cds.pergenesplit[[y]][[z]]@n.segregating.sites[, x],
                       pSegSites = cds.pergenesplit[[y]][[z]]@n.segregating.sites[, x]/mygenedat[, cdslength],
                       pi_raw = cds.pergenesplit[[y]][[z]]@Pi[, x],
                       pi_persite = cds.pergenesplit[[y]][[z]]@Pi[, x]/mygenedat[, cdslength], # all sites divided by total gene length; others by CDS length
                       tajimasD = cds.pergenesplit[[y]][[z]]@Tajima.D[, x]
                       # ***ANYTHING ELSE??***
                       )
          })
        )
      })
    )
    )
  }))
  
  # Windows
  perwindpop<-rbindlist(lapply(1:length(gwindows@populations), function(x){ # x is population index
    data.table(chr = mychr, gwindows_pos,
               population = names(gwindows@populations)[x],
               nSegregatingSites = gwindows@n.segregating.sites[, x],
               pSegSites = gwindows@n.segregating.sites[, x]/gwindows@n.sites,
               pi_raw = gwindows@Pi[, x],
               pi_persite = gwindows@Pi[, x]/gwindows@n.sites,
               tajimasD = gwindows@Tajima.D[, x]
               # ***ANYTHING ELSE??***
               )
  }))
  
  # Return!
  return(list(genes = pergenepop, windows = perwindpop))
})

#### Format across chromosomes & save #### 
cat("....PopGenome processing complete; finalizing formatting and saving out summary statistics....\n")
pergene<-rbindlist(lapply(popgenout, function(x) x$genes))
perwindow<-rbindlist(lapply(popgenout, function(x) x$windows))

write.table(pergene, gzfile(file.path(p$outdir, paste0(p$out, "_nucdivpergene.txt.gz"))), sep = "\t", quote = F, row.names = F)
write.table(perwindow, gzfile(file.path(p$outdir, paste0(p$out, "_nucdivper", p$windowbp, "window.txt.gz"))), sep = "\t", quote = F, row.names = F)

#### Make a few automatic basic plots for fun ####
cat("....Generating location vs. nucleotide diversity plots....\n")

# Naming, formatting, etc
statsinfo<-data.table(colname = c("pSegSites", "pi_persite", "tajimasD", "nSegregatingSites", "pi_raw"),
                      descrip = c("Proportion segregating sites", "Pi (per site)",
                                  "Tajima's D", "N segregating sites", "Raw pi"))
subsitesinfo<-data.table(colname = c("all", "coding", "syn", "nonsyn"),
                         descrip = c("All segregating sites", "Coding region (CDS) segregating sites", 
                                     "Synonymous segregating sites", "Nonsynonymous segregating sites"))
datainfo<-data.table(datname = c("pergene", "perwindow"),
                     descrip = c("Genes", paste0("Genomic windows (", p$windowbp/1e03, " kb)")),
                     filedescrip = c("pergene", paste0("per", p$windowbp, "window")))

# Populations are stacked; one statistic/sites subset per page
#     PDF separated just by what regions are (genes or windows)
npops<-length(pergene[,unique(population)])
invisible(
lapply(1:nrow(datainfo), function(di){
  pdf(file.path(p$outdir, paste0(p$out, "_nucdivplots_", datainfo[di, filedescrip], "_stackedpops.pdf")), 10, max(4, 1.6*npops))
  lapply(1:nrow(statsinfo), function(x){
    if(datainfo[di, datname]=="pergene"){ # do subsites of pergene
      lapply(1:nrow(subsitesinfo), function(y){
        print(
          posvsvalplot(dat = eval(as.name(datainfo[di, datname]))[sites==subsitesinfo[y, colname], ], 
                       chrs = chrs, ycol = statsinfo[x, colname], stackfacet = "population",
                       myylab = statsinfo[x, descrip], mytitle = datainfo[di, descrip],
                       mysubt = paste(statsinfo[x, descrip], subsitesinfo[y, descrip], sep = " | "))
        )
        return(NULL)
      })
      }else{ # make only plot for perwindow
        print(
        posvsvalplot(dat = eval(as.name(datainfo[di, datname])), 
                     chrs = chrs, ycol = statsinfo[x, colname], stackfacet = "population",
                     myylab = statsinfo[x, descrip], mytitle = datainfo[di, descrip],
                     mysubt = statsinfo[x, descrip])
        )
    }
    return(NULL)
  })
  invisible(dev.off())
  return(NULL)
})
)

# Sites subsets are stacked; one statistic/population per page
di<-1
pdf(file.path(p$outdir, paste0(p$out, "_nucdivplots_", datainfo[di, filedescrip], "_stackedsubsites.pdf")), 10, 8)
invisible(
  lapply(1:nrow(statsinfo), function(x){
    lapply(eval(as.name(datainfo[di, datname]))[,unique(population)], function(y){
      print(
        posvsvalplot(dat = eval(as.name(datainfo[di, datname]))[population==y, ],
                     chrs = chrs, ycol = statsinfo[x, colname], stackfacet = "sites",
                     myylab = statsinfo[x, descrip], mytitle = paste("Population:", y))
      )
    })
    return(NULL)
  })
)
invisible(dev.off())
    
#### Save session info ####
cat("....nucdivcendr_geneswindows_allandasestrains.R complete! Session information:....\n")
sessionInfo()
