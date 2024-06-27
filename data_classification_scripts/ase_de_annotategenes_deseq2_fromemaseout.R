#! /usr/bin/env/ Rscript
# Allele-specific expression analysis + comparison to parental samples, annotated with various other experimental information about the genes
#   including number of unique alignments in underlying Salmon; whether gene is in CeNDR called hyperdivergent haplotype; etc
#    using DESeq2 to generate log2 fold changes, p-values, etc from EMASE output 
# Strains are included in same model/input data, but results are within-strain and included separately
# by Avery Davis Bell, begun 2021.02.28

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
require(DESeq2, lib.loc = mylibloc ) # trying to make this play nice on PACE...
require(ashr, lib.loc =  mylibloc ) # trying to make this play nice on PACE...
require(ggVennDiagram, lib.loc = mylibloc)
require(ggplot2, lib.loc =  mylibloc) # trying to make this play nice on PACE...

cat("Session info for early troubleshooting -- \n")
sessionInfo() #### tmp print for troubleshooting

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

fread.emase<-function(emfile){
  # Reads & properly column names emase output files where the first line is header but starts with comment character and the second line starts with comment character
  # In: emfile, path to emase-zero output file to read
  # Out: data.table of emase-zero output with column names as in that file
  
  # Get names
  cnames<-strsplit(readLines(emfile, n = 1), "\t")[[1]]
  if(startsWith(cnames[1], "#")){
    cnames[1] <- substr(cnames[1], 2, nchar(cnames[1]))
  }
  
  # Get data
  dt<-fread(emfile, skip = 2, header = F, col.names = cnames)
  return(dt)
}

totctmatrix_emasedir<-function(sampinfo, emasedir){
  # Reads in emase count matrix 'total' column - one column per sample.
  # In: sampinfo, data.table with sample information. Only required column is SampleID *in same format as in exampfile*.
  #     emasedir, Path to directory with one EMASE output .gene.counts (or .gene.counts.gz) per sample. Filenames must start with sample ID (from sampinfo)
  # Out: matrix of genes (rows, row names) x samples. Columns are named SampleID. Values are from 'total' column of emase input  
  
  # Read in
  emfiles<-list.files(emasedir)
  samp.emase<-lapply(sampinfo$SampleID, function(x){
    samp.files<-emfiles[grepl(x, emfiles)]
    emfile<-file.path(emasedir, samp.files[grepl("gene.counts", samp.files)])
    em<-fread.emase(emfile)
    setnames(em, names(em)[2:3], paste(x, names(em)[2:3], sep = "_"))
    return(em[,.(target_id, total)])
  })
  print("read all")
  
  # Combine
  ctmat<-as.matrix(do.call(cbind, lapply(samp.emase, function(x) x[,total])))
  row.names(ctmat)<-samp.emase[[1]][,target_id]
  colnames(ctmat)<-sampinfo$SampleID
  
  # Return
  return(ctmat)
}

alctmatrix_emasedir<-function(sampinfo, emasedir){
  # Combines emase outputs for all samples into count matrix with one column per allele. Generates sample info with one row per allele
  # In: sampinfo, data.table with sample information. Only required column is SampleID *in same format as in exampfile*. Order of samples here will be their order in output!
  #         Other optional columns will be included in sample info output!
  #     emasedir, Path to directory with one EMASE output .gene.counts (or .gene.counts.gz) per sample. Filenames must start with sample ID (from sampinfo)
  # Out: list. 
  #     $cts, matrix of genes (rows, row names) x sample/allele pairs. Two columns per sample (allele 1, allele 2). Named SampleID_<name of allele from EMASE>
  #     $samps, data.table of samples. Two rows for each input sample, in same order as $cts. Row names are same as column names of $cts
  #       SampleID and all other information from sampinfo are included
  #       critically, an Allele column is added: which allele (from emase) does row describe. For use in model to get ASE estimates.
  
  # Read in all samples' emase counts; re-naming columns to reflect sample ID
  emfiles<-list.files(emasedir)
  samp.emase<-lapply(sampinfo$SampleID, function(x){
    samp.files<-emfiles[grepl(x, emfiles)]
    emfile<-paste0(emasedir, "/", samp.files[grepl("gene.counts", samp.files)])
    em<-fread.emase(emfile)
    setnames(em, names(em)[2:3], paste(x, names(em)[2:3], sep = "_"))
    return(em[,.SD,.SDcols = -4]) # don't keep total column!
  })
  
  # Combine
  ctmat<-as.matrix(do.call(cbind, lapply(samp.emase, function(x) x[,.SD,.SDcols=-1])))
  row.names(ctmat)<-samp.emase[[1]][,target_id]
  
  # Sample info: 2 entries per sample (1 per allele)
  samps<-copy(sampinfo)
  samps<-rbind(samps, samps)
  setkey(samps, SampleID)
  samps<-samps[sampinfo$SampleID, ] # group sets of 2 sample together - in order of ctmat
  row.names(samps)<-colnames(ctmat)
  ## Add Allele column: this is Strain of the specific allele counted in the given column
  myals<-sapply(strsplit(row.names(samps), "_"), function(x) x[length(x)]) ## need to take end - if other underscores in sample name, hard coding to take x[2] doesn't work##
  samps[, Allele:=myals] 
  
  # Return
  return(list(cts = ctmat, samps = samps))
}

getaseres<-function(dds, resn, alpha = 0.1, skewthresh = 0.6, desuff = ""){
  # Extracts result of interest from dds, gets in terms of allele proportion. Assumes result of interest is allele effect vs. ref (optionally in a given condition)
  # log2 fold changes shrunken with ashr!
  # In: dds, DESeq2 object that includes results to extract from. 
  #     resn, name of result to get
  #     alpha, p-value threshold
  #     skewthresh, proportion alleles from one haplotype OR the other for results to be considered significant (signifAtThresholds output column). I.e., if 0.6, genes with <=0.4 or >=0.6 alt alleles are considered significant
  #         # ** NEW: also used for calculating lfcThreshold in ASHR call
  #     desuff, suffix for columns that come straight from DESeq2 output (useful if these will be combined in same DT with other results later)
  # Out: data.table with columns gene_id, (KEYED by gene_id)
  #       baseMean, log2FoldChange, lfcSE, pvalue, padj (from DESeq2 ASHR-shrunken results) - will end in desuff if provided
  #       altVref - proportion of ref allele that alt allele is (calculated from log2FC)
  #       altVtotal - proportion of total that alt allele is (calculated from log2FC)
  #       signifAtThresholds.<desuff>,  T or F: is this gene significant at the alpha and allele skew thresholds provided
  
  res<-results(dds, name = resn, alpha = alpha)
  res<-as.data.table(lfcShrink(dds, res = res, type = 'ashr', 
                               lfcThreshold = log2(-1*skewthresh/(skewthresh-1))), keep.rownames = T)
  res[, altVref:=2^log2FoldChange] # prop of ref that alt is
  res[, altVtotal:=altVref/(1 + altVref)] # prop of total alleles that alt is
  
  # Annotate with overlap with thresholds
  skthrs<-sort(c(skewthresh, 1 - skewthresh))
  res[, signifAtThresholds:=ifelse(!is.na(padj) & padj<alpha & (altVtotal<=skthrs[1] | altVtotal>=skthrs[2]), T, F)]
  
  # Name columns appropriately
  setnames(res, "rn", "gene_id")
  setnames(res, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "signifAtThresholds"),
           paste0(c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "signifAtThresholds"), desuff))
  
  setkey(res, "gene_id")
  return(res)
}

plotPCA_givePCs_colshape<-function (object, colgroup = "Strain", shapegroup = "Treatment", ntop = 500,
                                    returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # also expects to plot color as one variable, shape as another
  # In: see ?plotPCA except for:
  #   colgroup: what to color by
  #   shapegroup: what to shape by
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  rv <- rowVars(assay(object), useNames = T)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(c(colgroup, shapegroup) %in% names(colData(object)))) {
    stop("the argument colgroup and shapegroup should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, c(colgroup, shapegroup), 
                                               drop = FALSE])
  
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  plt<- ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = colgroup, shape = shapegroup)) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
  
  
  return(plt)
}

plotPCA_givePCs_colshapesz<-function (object, colgroup = "Strain", shapegroup = "Treatment",sizeby = "DropletBleach", sizeptval = "Yes",
                                           ntop = 500, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # also expects to plot color as one variable, shape as another, size as another. CLUNKY
  # In: see ?plotPCA except for:
  #   colgroup: what to color by
  #   shapegroup: what to shape by
  #   sizeby: column name - what to size by.CATEGORICAL column.
  #   sizeptval: value of sizeby to set as reference (have come first on legend, will be smallest point)
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  rv <- rowVars(assay(object), useNames = T)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(c(colgroup, shapegroup) %in% names(colData(object)))) {
    stop("the argument colgroup and shapegroup should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, c(colgroup, shapegroup, sizeby), 
                                               drop = FALSE])
  
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], 
                  intgroup.df, name = colnames(object))
  
  # x'ing point related
  d[,sizeby]<-relevel(d[,sizeby], ref = sizeptval)
    
  plt<-ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = colgroup, shape = shapegroup, size = sizeby)) + 
    geom_point(alpha = 0.6) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    labs(size = sizeby) +
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          title = element_text(size = 17), legend.text = element_text(size = 13),
          plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15))
}

getgenstrainres<-function(dds, strn, comp = "N2", alpha = 0.1, lfcthresh = 0.5849625){
  # Gets 3 sets of results involving strain provided from DESeq object where Generation-Strain combination was used as explanatory; outputs in same data.table
  #   1. Given strain (parent) vs. reference strain (comp)
  #   2. Given strain F1 hybrid with reference vs. reference parental strain (comp)
  #   3. Given strain F1 hybrid with reference vs. non-reference parental strain (strn)
  # In: dds, DESeq2 object with results for GenerationStrain column (of coldata) named Parent.<Strain - given strn and comp> and F1.<Strain - given strn and comp>
  #         i.e., must be able to call results(dds, contrast = c("GenerationStrain", "Parent.<strn>", "Parent.<comp>")) and other queries
  #     strn, Test strain to be in numerator of comparisons (to retrieve results for)
  #     comp, Comparison strain to compare to - should be reference strain, usually
  #     alpha, adjusted p-value threshold for calling genes significantly DE
  #     lfcthresh, log2 fold change threshold for calling genes significantly DE **included in call to ashr too**
  # Out: data.table of results, keyed by gene_id. Columns
  #       gene_id, gene ID
  #       baseMean.GenerationStrain, baseMean from this dds. (same for all results)
  #     *Each of the following sets of columns x3 suffixes:
  #                                                 first is <name>.ParentVsParent<comp> - results for Given strain (parent) vs. reference strain (comp)
  #                                                 second is <name>.F1VsParent<comp> - results for  Given strain F1 hybrid with reference vs. reference parental strain (comp)
  #                                                 third is <name>.F1VsParentNonRef - results for Given strain F1 hybrid with reference vs. non-reference parental strain (strn)
  #       Columns with these suffixes: log2FoldChange, lfcSE, pvalue, padj (all from DESeq2); signifAtThresholds (T or F - padj passes alpha threshold and log2FoldChange passes lfcthresh)                                              
  
  # Subfunctions
  formatgenstrnres<-function(dds, res, desuff = "", alpha = 0.1, lfcthresh = 0.5849625, 
                             exclcols = NA){
    # ASHR-shrinks results; Names columns of res as what they are plus desuff; adds column signifAtThresholds (T if alpha, lfcthresh passed)
    # In: dds, DESeq object these results came from
    #     res, DESeq results() output. Columns baseMean, log2FoldChange, lfcSE, pvalue, padj
    #     desuff, what to end column names with (delimiter not added - include it in input if you want it)
    #     alpha,  adjusted p-value threshold for calling genes significantly DE
    #     lfcthresh, log2 fold change threshold for calling genes significantly DE
    #     exclcols, name of column to exclude from output, or NA to keep all
    # Out: data.table of results, but with new column signifAtThresholds<desuff> and with all other columns having desuff appended
    #     log2FCs are ASHR-shrunken; currently not keeping gene_ids
    
    # Format results
    res<-as.data.table(lfcShrink(dds, res = res, type = 'ashr', lfcThreshold = lfcthresh))
    if(!(is.na(exclcols))){
      res<-res[,(exclcols):=NULL]
    }
    res[,signifAtThresholds:=ifelse(!is.na(padj) & padj<alpha & 
                                      (log2FoldChange > abs(lfcthresh) | log2FoldChange < -1*abs(lfcthresh)), T, F)]
    
    # Update names
    setnames(res, names(res), paste0(names(res), desuff))
    
    # Return
    return(res)
  }
  
  # Get all results & shared metadata
  out<-data.table(gene_id = rownames(dds), # same for all results
                  baseMean.GenerationStrain = mcols(dds)$baseMean, # same for all results
                  do.call(cbind, list(
                    formatgenstrnres(dds = dds,
                                     res = results(dds, contrast = c("GenerationStrain", paste0("Parent.", strn), paste0("Parent.", comp)), alpha = alpha, lfcThreshold = lfcthresh),
                                     desuff = paste0(".ParentVsParent", comp),
                                     alpha = alpha, 
                                     lfcthresh = lfcthresh, 
                                     exclcols = "baseMean"),
                    formatgenstrnres(dds = dds,
                                     res = results(dds, contrast = c("GenerationStrain", paste0("F1.", strn), paste0("Parent.", comp)), alpha = alpha,  lfcThreshold = lfcthresh),
                                     desuff = paste0(".F1VsParent", comp),
                                     alpha = alpha, 
                                     lfcthresh = lfcthresh, 
                                     exclcols = "baseMean"),
                    formatgenstrnres(dds = dds,
                                     res = results(dds, contrast = c("GenerationStrain", paste0("F1.", strn), paste0("Parent.", strn)), alpha = alpha,  lfcThreshold = lfcthresh),
                                     desuff = ".F1VsParentNonRef",
                                     alpha = alpha, 
                                     lfcthresh = lfcthresh, 
                                     exclcols = "baseMean")
                  )))
  
  # Key by gene_id & return
  setkey(out, gene_id)
  return(out)
}

freadeqctshapgene<-function(sampids, geneids, exampuniqeqcts){
  # Reads in the number of haplotype and gene specific alignments for the given samples, genes 
  # In: sampids, sample IDs to read files for. These will be used as column names in the output and *must* be in exampuniqeqcts file paths
  #     geneids, gene IDs to retain records for.
  #     exampuniqeqcts, Example filepath to per-sample counts of eq classes, haplotype-specific eq classes, gene-specific eq classes for each _gene_ (*eqclasses_genes.txt.gz output of salmonalleleeqclasses.py).
  #                   Where sample ID (as in --sampleinfo) goes in filename, should have string SAMP.
  #                   For any other differences among filenames (e.g. genome aligned to), include the glob '*' asterisk character
  # Out: data.table with one row per gene in geneids. Columns gene_id, one for each in sampids. Keyed by gene_id.
  
  out.l<-lapply(sampids, function(x){
    eqgenef<-Sys.glob(gsub("SAMP", x, exampuniqeqcts))
    onesamp<-fread(eqgenef, select = c("gene_id", "naln_hapgenespec"))[gene_id%in%geneids]
    setnames(onesamp, "naln_hapgenespec", x)
    setkey(onesamp, gene_id)
    return(onesamp)
  })
  out<-Reduce(function(x, y) merge(x, y, all = T), out.l)
  setkey(out, gene_id) # redundant
  return(out)
}

summndt<-function(ndt, prefix){
  # Summarizes a data.table with columns of numbers (see output description) - written for # unique alignments, so makes sense currently for use with positive integers
  # In: ndt, data.table with multiple columns; all must contain numbers
  #     prefix, what to start output column names with. No delimiter added, this used directly.
  # Out: data.tale with one row for each row in input. Columns:
  #   <prefix>nOver0 - # of columns with values > 0
  #   <prefix>allOver0 - T or F, are all input columns over 0
  #   <prefix>min - minimum value of that row
  #   <refix>med - median value of that row

  out<-ndt[, .(sum(unlist(.SD) > 0), sum(unlist(.SD) > 0)==ncol(ndt),
               min(unlist(.SD)), median(unlist(.SD))),
           by = .I]
  out[,I:=NULL]
  setnames(out, paste0(prefix, c("nOver0", "allOver0", "min", "med")))
  return(out)
}

prochypdivbed<-function(divregbed, geneinfo, isotypes){
  # For genes of interest, determines whether they overlap divergent regions in each of the isotypes of interest
  # In: divregbed,Filepath to hyperdivergent haplotypes BED file as downloaded from CeNDR 20210121 release (published as Lee et al 2020/2021 preprint/pub)
  #       BED format: no column names; columns are chr, start, end, isotype
  #     geneinfo, data.table with information on genes of interest. Needed columns are gene_id, chr, start, end (gene coordinates)
  #     isotypes, vector of isotype names (as in divregbed) for which to determine hyperdivergent status
  # Out: data.table with columns:
  #      gene_id, gene_id from geneinfo. KEYED BY THIS.
  #     one column per isotype, named <isotype from isotypes vector>.div containing T or F for if that gene overlaps hyperdivergent region
  
  # Read in divergent regions, keep only isotypes of interest
  divs<-fread(divregbed, header = F)
  setnames(divs, c("chr", "start", "end", "iso"))
  divs<-divs[iso%in%isotypes, ]
  setkey(divs, chr, start, end)
  
  # Determine divergent yes/no for each gene/strain pair
  ginfo<-copy(geneinfo) # don't want to modify input DT 
  setkey(ginfo, chr, start, end)
  gdiv<-data.table(ginfo[, gene_id],
                   (do.call(cbind, lapply(isotypes, function(x){
                     div1<-divs[iso==x, ]
                     setkey(div1, chr, start, end)
                     return(!is.na(foverlaps(ginfo, div1, type = "any", mult = "first", which = T)))
                   })))
  )
  setnames(gdiv, c("gene_id", paste(isotypes, "div", sep = ".")))
  
  # Return
  setkey(gdiv, gene_id)
  return(gdiv)
}

#### Arguments & Input ####
# --- Command line arguments
p<-arg_parser("Allele-specific expression analysis + comparison to parental samples, annotated with various other experimental information about the genes
        including number of unique alignments in underlying Salmon; whether gene is in CeNDR called hyperdivergent haplotype; etc
        using DESeq2 to generate log2 fold changes, p-values, etc from EMASE output 
        Strains are included in same model/input data, but results are within-strain and included separately", 
              name = "ase_de_annotategenes_deseq2_fromemaseout.R ", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--sampleinfo",
                help = "Path to sample information file containing information on parental and F1 samples. Columns should include:
                SampleID, Generation (Parental or F1), Allele1 (NA for parental, reference strain for F1), Allele2 (NA for parental, non-reference strain for F1 - used as Strain for F1 in generation/strain model),
                Strain (NA for F1, strain for parental. All strains/alleles should be the way you want them to show up in outputs),
                terms in --aseotherterms and --genstrainmodel should be column names if there are other terms provided not already accounted for.",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character")
p<-add_argument(p, "--emasedir",
                help = "Path to directory with one EMASE output .gene.counts (or .gene.counts.gz) per sample. Filenames must start with sample ID (from sampinfof1, sampinfoparent)",
                type = "character")
p<-add_argument(p, "--allelerename",
                help = "Optional - allele mapping/renaming between files. Path to file with columns emase, sample: initial should be alleles as in EMASE outputs; final should be alleles as in sample information!", 
                default = NA)

# Statistical testing-related arguments
p<-add_argument(p, "--aseotherterms",
                help = "(optional) Term(s) (columns of sampinfo) to include in ASE model before allele. ONLY TESTED WITH ONE PROVIDED. Comma-separated if multiple. If one provided,
                  design will be ~<this> + <this>:Sample[grouped as needed] + Allele",
                default = NA,
                type = "character")
p<-add_argument(p, "--refcategoryinfo", 
                help = "Path to matrix describing the reference level for ALL factors in the ASE model, even if you don't care (Allele, other terms). Columns 'colname', 'reflevel'. reflevel for allele should be as in sample info, not emase.
                **Important: ref level for Allele is assumed to be reference strain to be used as denominator for comparisons here and in generation/strain model",
                type = "character")
p<-add_argument(p, "--alpha",
                help = "Alpha p-value threshold for FDR-like filtering (used for ASE and Generation/Strain group comparisons)",
                type = "numeric",
                default = 0.1)
p<-add_argument(p, "--alleleskewthresh",
                help = "proportion alleles from one haplotype OR the other for results to be considered significant. (used to threshold calls that can be significant within DESeq2).
                I.e., if 0.6, genes with >=60% of one allele (<=40% or >=60% alt. alleles) AND significant p-values are considered significant.
                Leave at default 0.5 to just threshold on p-value.",
                type = "numeric",
                default = 0.5)
p<-add_argument(p, "--genstrainmodel",
                help = "Model design for model including all samples (parental and F1). Last term should be GenerationStrain - 
                this will be made by combining generation and strain (or allele2 for F1s) to be able to compare among parental strains and between F1s and their parents.
                Include any batch or other covariates to regress out! Quote-wrap if any spaces.",
                default = "~GenerationStrain")
p<-add_argument(p, "--genstrainlfc",
                help = "Log2FoldChange threshold for genes to be called significant in GenerationStrain model (so, for parents vs. each other and F1 vs. each parent).
                Used to threshold calls that can be significant within DESeq2 testing. Consider matching this to alleleskewthresh (i.e. an alleleskewthresh of 0.6 is equivalent to a fold-change threshold of 1.5, LFC threshold of 0.5849625",
                default = 0.5849625)

# Arguments about gene inclusion/information 
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

# Arguments about informative-ness/unique alignments
p<-add_argument(p, "--exampuniqeqcts",
                help = "Example filepath to per-sample counts of eq classes, haplotype-specific eq classes, gene-specific eq classes for each _gene_ (*eqclasses_genes.txt.gz output of salmonalleleeqclasses.py).
                Where sample ID (as in --sampleinfo) goes in filename, should have string SAMP. One file should exist for each F1 sample.
                For any other differences among filenames (e.g. genome aligned to), include the glob '*' asterisk character",
                type = "character")

# Arguments characterizing genes in other ways (hyperdivergence, gene coverage)
p<-add_argument(p, "--strain2iso",
                help = "Path to strain-to-isotype mapping file for all non-reference strains included here (even if strain and isotype are identical!). Columns strain, isotype; one row per included strain.",
                type = "character")
p<-add_argument(p, "--leehypdivbed",
                help = "Filepath to hyperdivergent haplotypes BED file as downloaded from CeNDR 20210121 release (published as Lee et al 2020/2021 preprint/pub)",
                type = "character")
p<-add_argument(p, "--genednacov",
                help = "Path to file containing information on DNA coverage per gene (e.g. *_genecoveragefromexons_mednorm.txt.gz output of exploregenecoverage_fromexons.R/mosdepthmergedexons.nf workflow).
                Columns must include gene_id, <named for all isotypes matching strains in RNA-seq data>. Each isotype-named column
                contains DNA coverage value to be added (e.g. gene median-normalized coverage) - no downstream processing/normalization done on these.",
                type = "character")
p<-add_argument(p, "--dnalowcovthresh",
                help = "Threshold below which DNA coverages (proved in --genednacov) will be flagged as low coverage in column of output (lowDNACov)",
                default = 0.25)
p<-add_argument(p, "--dnahighcovthresh",
                help = "Threshold above which DNA coverages (proved in --genednacov) will be flagged as high coverage in column of output (highDNACov)",
                default = 2)

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# Output directory
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

if(!is.na(p$exclchrs)){
  chrexcl<-strsplit(p$exclchrs, split = ",", fixed = T)[[1]]
}else{
  chrexcl<-NA
}

biotypes<-strsplit(p$inclbiotype, split = ",", fixed = T)[[1]]

# --- Read in gene information genome-wide for only genes to be retained downstream
cat("....Reading in gene information and narrowing to genes from user-specified chromosomes, biotypes....\n")
geneinfo<-getlocusmetadata(gff3file = p$genegff, chrexcl = chrexcl, biotypes = biotypes)

# --- Read in counts
# Sample information
samps<-fread(p$sampleinfo, header = T)

# Other information
if(!is.na(p$allelerename)){allelefix<-fread(p$allelerename, header = T)}
reflevels<-fread(p$refcategoryinfo, header = T)

    # --- Get in count data per sample (no allele split) --- #
cat("....Reading in and formatting gene expression data - total per sample....\n")
cts.tot<-totctmatrix_emasedir(samps, p$emasedir)

    # --- Get in count data per allele for F1 ASE --- #
cat("....Reading in and formatting gene expression data - F1s per allele....\n")
cts.ase<-alctmatrix_emasedir(samps[Generation=="F1",], p$emasedir)
if(!is.na(p$allelerename)){ # Fix alleles inherited from EMASE inputs
  for(i in 1:nrow(allelefix)){
    cts.ase$samps[Allele==allelefix[i,emase], Allele:=allelefix[i, sample]]
  }
}

# --- Narrow to genes to include for all downstream steps
# Exclude those excluded from gene information
ngs.cts<-nrow(cts.tot) # keep record so can say how many originally were in data
cts.tot<-cts.tot[geneinfo$gene_id, ] 
cts.ase$cts<-cts.ase$cts[geneinfo$gene_id, ] 

# Exclude based on low coverage (total reads < 10 across all samples)
ngs.prelowcov<-nrow(cts.tot)
keep<-which(rowSums(cts.tot)>=10)
cts.tot<-cts.tot[keep, ]
cts.ase$cts<-cts.ase$cts[keep, ]

# Print to log genes retained/excluded (Ns)
cat(paste("       ", length(keep), "genes retained for analysis:", ngs.cts - ngs.prelowcov, "excluded due to input chromosome and/or biotype restrictions;", 
          ngs.prelowcov - length(keep), "excuded due to low expression (sum of samples' read counts < 10)\n"))

#### Perform ASE testing: F1s only ####
cat("....Formatting data for DESeq2 ASE analysis....\n")
# Column metadata
coldata.ase<-data.frame(copy(cts.ase$samps))
row.names(coldata.ase)<-row.names(cts.ase$samps)
## Alleles
coldata.ase$Allele<-factor(coldata.ase$Allele)
## Sample groups
if(!is.na(p$aseotherterms)){
  # Figure out what to group within
  grp.within<-strsplit(p$aseotherterms, ",")[[1]]
  
  # Get unique combinations of these (assign groups within these combinations)
  vals.grp.within<-expand.grid(lapply(grp.within, function(x) unique(coldata.ase[,x])), stringsAsFactors = F)
  colnames(vals.grp.within)<-grp.within
  sampgrp<-rep(NA, nrow(coldata.ase))
  for(i in 1:nrow(vals.grp.within)){
    if(ncol(vals.grp.within) > 1){ # make sure all columns match
      samps.these<-apply(coldata.ase[,grp.within]==c(vals.grp.within[i, grp.within]), 1, all)
    }else{ # only one column needs match
      samps.these<-(coldata.ase[,grp.within]==c(vals.grp.within[i, grp.within]))
    }
    sampgrp[samps.these] <- rep(1:(sum(samps.these)/2), each = 2)
  }
  coldata.ase$SampleGrp<-factor(sampgrp)
}else{
  # make a SampleGrp column from samples for consistency of using same column name in design
  # Convert to simple number only in case SampleID column has non-safe values
  coldata.ase$SampleGrp<-factor(coldata.ase$SampleID, labels = 1:length(unique(coldata.ase$SampleID)))
}
## Reference levels
for(x in 1:nrow(reflevels)){
  coldata.ase[,reflevels[x, colname]]<-factor(coldata.ase[,reflevels[x, colname]])
                                             # levels = unique(c(reflevels[x, reflevel], unique(coldata.ase[,reflevels[x, colname]]))))
  coldata.ase[,reflevels[x, colname]]<-relevel(coldata.ase[,reflevels[x, colname]], ref = reflevels[x, reflevel])
}

# ASE model design
if(is.na(p$aseotherterms)){
  modeldesign<- "~ SampleGrp + Allele" # this is for NO resultswithin OR other terms
}else{
  # add in other term & blocking factor using this term
  if(length(grp.within==1)){ # currently implemented - 1 grp, 1 interaction
    modeldesign<-paste0("~ ", paste(grp.within, collapse = " + "), # main effects of other terms
                        " + ", paste0(grp.within[1], ":SampleGrp"), # sample blocking
                        " + Allele")
  }else if(length(grp.within)>1){ # UNTESTED: multiple grps, interactions
    modeldesign<-paste0("~ ", paste(grp.within, collapse = " + "), # main effects of other terms
                        " + ", paste(unlist(lapply(1:(length(grp.within) - 1), function(x){
                          sapply((x + 1):length(grp.within), function(y) paste(grp.within[x], grp.within[y], sep = ":"))
                        })), collapse = " + "), # pairwise interactions of other terms
                        # more than pairwise? urghhhh.  NOT IMPLEMENTED ****************************************************
                        " + ", paste(paste(grp.within, collapse = ":"), "SampleGrp", sep=":"), # sample blocking
                        " + Allele")
  }
} # end else meaning there are otherterms
cat(paste("    ASE model design is:", modeldesign, "\n")) 

# DESeq2 format.
## Make sure all combinations of factors exist
m1<-model.matrix(as.formula(modeldesign), data = coldata.ase)
all.zero <- apply(m1, 2, function(x) all(x==0))
if(sum(all.zero)==0){ # No all-zero columns; proceed with normal DESeq2 procedure
  dds.ase<-DESeqDataSetFromMatrix(countData = round(cts.ase$cts),
                                  colData = coldata.ase,
                                  design = as.formula(modeldesign))
  sizeFactors(dds.ase)<-rep(1, nrow(coldata.ase)) ## set size factors equal (since 2 observations come per individual and individual effects already blocked!)
  cat("....Running F1 DESeq2 ASE analysis ....\n")
  dds.ase<-DESeq(dds.ase)
  
}else{ # Not all combinations of factors exist
  cat("     NOTE: Levels or combinations of levels without any samples have resulted
      in column(s) of zeros in the model matrix, so model matrix with these columns manually removed
      will be used.\n")
  mUse <- m1[,-which(all.zero)]
  dds.ase<-DESeqDataSetFromMatrix(countData = round(cts.ase$cts),
                                  colData = coldata.ase,
                                  design = ~ Allele) # PLACEHOLDER design
  sizeFactors(dds.ase)<-rep(1, nrow(coldata.ase)) ## set size factors equal (since 2 observations come per individual and individual effects already blocked!)
  
  # Run Wald test providing own model matrix
  dds.ase<-DESeq(dds.ase, test = "Wald", full = mUse, betaPrior = F)
}

save(dds.ase, file = file.path(p$outdir, paste0(p$baseoutname, "_dds_ase.RData")), ascii = F)

# --- Get ASE results in long term format (+ combined with gene metadata) - these are the results that'll be carried forward
cat("....Pulling ASE results to start output results matrices....\n")
# Set up results naming info etc & pull results
resnames.ase<-resultsNames(dds.ase)[grepl("Allele", resultsNames(dds.ase))]
# resnaming.ase<-data.table(resname = resnames.ase,  # OLD; when have different models, result names can be slightly different for same result
#                           allele = sapply(resnames.ase, function(x) strsplit(x, "Allele")[[1]][-1]))
resnaming.ase<-data.table(resname = resnames.ase,
                         allele = levels(colData(dds.ase)$Allele)[2:nlevels(colData(dds.ase)$Allele)])
setkey(resnaming.ase, resname)

# Results pull
reslist<-lapply(resnaming.ase$resname, function(x){
  getaseres(dds = dds.ase, resn = x, alpha = p$alpha, skewthresh = p$alleleskewthresh, desuff = ".ASE")
})
names(reslist)<-resnaming.ase$allele # allele/strain is name that will carry through
  
# Beginning gene info annotation (from GFF)
reslist<-lapply(reslist, function(x){
  return(geneinfo[x])
})

#### Model with F1s and parentals ####
cat("....Formatting data for DESeq2 DE analyses (parentals & F1s in same model)....\n")
# Column data
coldata.all<-copy(samps)
coldata.all[Generation=="F1" & is.na(Strain), Strain:=Allele2] # if Strain is NA in F1s, use allele 2
coldata.all[, GenerationStrain:=factor(paste(Generation, Strain, sep = "."))] # . isn't great for filenames but works with factors
coldata.all<-data.frame(coldata.all, stringsAsFactors = T)
rownames(coldata.all)<-coldata.all$SampleID

# Full DESeq2
dds.all<-DESeqDataSetFromMatrix(countData = round(cts.tot),
                                colData = coldata.all,
                                design = as.formula(p$genstrainmodel))
# RUN
cat("....Running DESeq2 DE analyses (parentals & F1s in same model)....\n")
dds.all<-DESeq(dds.all)
# Save DESeq2 object (use file.path notation)
save(dds.all, file = file.path(p$outdir, paste0(p$baseoutname, "_dds_genstrain.RData")), ascii = F)

# --- PCA plot with all samples [later addition] ---
cat("..quick PCA plots from all samples together...\n")
vsd<-vst(dds.all, blind = T)
pcadir<-file.path(p$outdir, "pcaplots")
if(!dir.exists(pcadir)){dir.create(pcadir, recursive = T)}
# IF do make these, make bigger axis etc with theme
# Generation & Strain only
pdf(file.path(pcadir, paste0(p$baseoutname, "_dds_genstrain_vst_pcaplots.pdf")), 7, 5.5)
plotPCA_givePCs_colshape(vsd, colgroup = "Strain", shapegroup = "Generation", ntop = 500, 
                         xpc = 1, ypc = 2, mytitle = "vst-transformed 500 most variable genes")
plotPCA_givePCs_colshape(vsd, colgroup = "Strain", shapegroup = "Generation", ntop = 500, 
                         xpc = 2, ypc = 3, mytitle = "vst-transformed 500 most variable genes")
plotPCA_givePCs_colshape(vsd, colgroup = "Strain", shapegroup = "Generation", ntop = 500, 
                         xpc = 3, ypc = 4, mytitle = "vst-transformed 500 most variable genes")
invisible(dev.off())
# Also with one other variable if included in model (only implemented for one!)
if(!is.na(p$aseotherterms)){
  pdf(file.path(pcadir, paste0(p$baseoutname, "_dds_genstrain_vst_pcaplots_with", grp.within[1], ".pdf")), 7, 5.5)
  # this fn returns plot so it needs to be printed
  print(plotPCA_givePCs_colshapesz(vsd, colgroup = "Strain", shapegroup = "Generation", sizeby = grp.within[1], 
                                   sizeptval = reflevels[colname==grp.within[1], reflevel], ntop = 500, 
                                   xpc = 1, ypc = 2, mytitle = "vst-transformed 500 most variable genes"))
  print(plotPCA_givePCs_colshapesz(vsd, colgroup = "Strain", shapegroup = "Generation", sizeby = grp.within[1], 
                                   sizeptval = reflevels[colname==grp.within[1], reflevel], ntop = 500, 
                                   xpc = 2, ypc = 3, mytitle = "vst-transformed 500 most variable genes"))
  print(plotPCA_givePCs_colshapesz(vsd, colgroup = "Strain", shapegroup = "Generation", sizeby = grp.within[1], 
                                   sizeptval = reflevels[colname==grp.within[1], reflevel], ntop = 500, 
                                   xpc = 3, ypc = 4, mytitle = "vst-transformed 500 most variable genes"))
  invisible(dev.off())
}


# --- Get contrast results (DE among-parental & F1 vs. parental) & combine in ---
reslist<-lapply(names(reslist), function(x){
  return(reslist[[x]][getgenstrainres(dds = dds.all, strn = x, comp = reflevels[colname=="Allele", reflevel],
                                      alpha = p$alpha, lfcthresh = p$genstrainlfc)])
})
names(reslist)<-resnaming.ase$allele # allele/strain is name that will carry through

#### Determine informative-ness/unique alignments & annotate results with this info #### 
cat("....Parsing and adding information about unique alignments per gene....\n")
# Get in counts per gene per sample
unqcts<-freadeqctshapgene(sampids = samps[Generation=="F1", SampleID], geneids = reslist[[1]]$gene_id,
                          exampuniqeqcts = p$exampuniqeqcts)

# Summarize per strain & add in to results
reslist<-lapply(names(reslist), function(x){
  mysamps<-samps[Generation=="F1" & Allele2==x, SampleID]
  unctinfo<-summndt(ndt = unqcts[ , mysamps, with = F], prefix = "unqalnmts.")
  unctinfo[,gene_id:=unqcts$gene_id]
  setkey(unctinfo, gene_id)
  out<-unctinfo[reslist[[x]]]
  setcolorder(out, c(names(geneinfo), names(unctinfo)[1:(ncol(unctinfo)-1)]))
  return(out)
})
names(reslist)<-resnaming.ase$allele # allele/strain is name that will carry through

#### Annotate with hyperdiverged haplotype information #####
cat("....Parsing and adding information about hyperdivergent haplotypes....\n")
# Get hyperdivergent haplotype overlap
strain2iso<-fread(p$strain2iso, header = T)
gdiv<-prochypdivbed(divregbed = p$leehypdivbed, geneinfo = reslist[[1]][,.(gene_id, chr, start, end)],
                    isotypes = strain2iso$isotype)    
# Add in to results
reslist<-lapply(names(reslist), function(x){
  div.x<-gdiv[,.(gene_id, get(paste0(strain2iso[strain==x, isotype], ".div")))]
  setnames(div.x, c("gene_id", "hypdiv"))
  out<-div.x[reslist[[x]]]
  setcolorder(out, c(names(geneinfo), "hypdiv"))
  return(out)
})  
names(reslist)<-resnaming.ase$allele # allele/strain is name that will carry through

#### Annotate with gene DNA coverage information #####
cat("...Parsing and adding information about DNA coverage at genes....\n")
dnacov<-fread(p$genednacov, header = T)
setkey(dnacov, gene_id)
reslist<-lapply(names(reslist), function(x){
  cov.x<-dnacov[,.(gene_id, get(strain2iso[strain==x, isotype]))]
  setnames(cov.x, c("gene_id", "dnacoverage"))
  cov.x[, lowDNACov:=dnacoverage < p$dnalowcovthresh]
  cov.x[, highDNACov:=dnacoverage > p$dnahighcovthresh]
  out<-cov.x[reslist[[x]]]
  setcolorder(out, c(names(geneinfo), "hypdiv", "dnacoverage", "lowDNACov", "highDNACov"))
  return(out)
})
names(reslist)<-resnaming.ase$allele # allele/strain is name that will carry through

#### Save ####
cat("...Writing output files....\n")
invisible( # So NULLs not printed
lapply(names(reslist), function(x){
  write.table(reslist[[x]], file = gzfile(file.path(p$outdir, paste0(p$baseoutname, "_", x, "_annotatedASEDEresults.txt.gz"))),
              quote = F, row.names = F, sep = "\t")
})
)

#### Script completion message & session information ####
cat("....ase_de_annotategenes_deseq2_fromemaseout.R processing complete! Session information:....\n")
sessionInfo()