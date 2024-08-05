#! /usr/bin/env/ Rscript
# Multipurpose script to take any per-gene 'phenotype' (categorical or quantitative) and plot and test it against ASE-derived categories (ASE/not, DE/not, regulatory pattern, inheritance mode)
# by Avery Davis Bell, begun 2024.04.24

#### System set up - especially required for playing nicely with PACE ####
# Get library location
if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
cat(paste("Library location:", mylibloc, "\n"))

require(data.table, lib.loc = mylibloc)
require(argparser, lib.loc = mylibloc)
require(ggplot2, lib.loc = mylibloc)
require(RColorBrewer, lib.loc = mylibloc)
require(ggforce, lib.loc = mylibloc) # used to make geom_sina
require(formattable, lib.loc = mylibloc) # used to do scientific notation formatting
require(scales, lib.loc = mylibloc) # for pseudolog transform

#### Functions ####
# --- Data reading
genesfilter<-function(gfile){
  # Reads list of genes to retain; determines if should be per strain or global [or all strains retained]
  # In: OPTIONAL Filepath to list of genes to narrow to if desired [if omitted, all genes in --resmatinput will be used].
  #       Can be no-header list of gene_ids (same genes retained for all strains)
  #       OR two-column file with columns strain, gene_id (genes will be retained in strain-specific manner)
  #       OR empty string; if so 'all' tokens will be returned
  # Out: list of $gflag, 'all', 'keepsame', OR 'bystrain': keep all genes; the same genes for each strain; or keep genes by strain
  #               $glist, 'all', OR <vector of gene_ids to keep for all strain>, OR <two-column data table with columns strain and gene_id>

  if(gfile==""){ # Flag to keep all if no file provided
    gflag<-"all"
    glist<-"all"
  }else{ # if file name provided, make sure it exists then proceed based on what's in there
    if(!file.exists(gfile)){
      stop(paste("Gene list file", gfile, "does not exist; fix and re-try"))
    }
    # Check header
    header<-strsplit(readLines(gfile, n = 1), "\t")[[1]]
    if("gene_id"%in%header & "strain"%in%header){ # Do by strain
      gflag<-"bystrain"
      glist<-fread(gfile, header = T)
    }else{ # keep genes globally
      gflag<-"keepsame"
      glist<-fread(gfile, header = F)$V1
    } # end if/else gene_id and strain in header
  }

  return(list(gflag = gflag, glist = glist))
}

freadcombdat<-function(strains, exampaseinput, inhmode, regpats, gflag, glist, totestf, testcolumn, informthresh){
  # Reads in all the data used by this script into one data.table; narrows to only genes of interest if inputs specify to do so
  # In: strains, vector of strain IDs to substitute into exampaseinput and read in
  #     exampaseinput, example filepath to ASE data-related file to read in as data.table, with strain identified in filename. Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz
  #     inhmode, path to file containing Inheritance mode classifications for all genes and strains. Columns strain, gene_id, inhmode
  #      regpats, path to file containing Regulatory pattern classifications for all genes, strains for which classification could be made - columns strain, gene_id, regclass
  #     gflag, 'all', 'keepsame', OR 'bystrain': keep all genes; the same genes for each strain; or keep genes by strain
  #     glist, 'all' if gflag is 'all', OR if gflag is 'keepsame' <vector of gene_ids to keep for all strain>, OR if gflag is 'bystrain' <two-column data table with columns strain and gene_id>
  #     totestf, Path to file containing the data to test against ASE etc! Must have columns gene_id, testcolumn value;
  #                   **if it has column strain, data will be matched in strain-specific manner; if not, data will be assigned multiply to each strain
  #     testcolumn, character name of column in totestf that has the data to keep for output
  #     informthresh, Gene must have this or more unique alignments in each sample to be considered informative for ASE/cis-trans analyses - in column unqalnmts.min in exampaseinput
  # Out: data.table with columns strain, inform [informative for ASE - T or F], regclass, inhmode, testDat [testcolumn data from totestf], <all columns in exampaseinput>

  # --- Check all files of interest exist
  ## ASE files
  asefs<-sapply(strains, function(x) file.path(gsub("STRAIN", replacement = x, exampaseinput, fixed = T)))
  cat(paste("ASE-related results files (check that this is what you intend!:\n", paste(asefs, collapse = "\n"), "\n"))
  if(!all(file.exists(c(asefs, inhmode, regpats)))){
    stop(paste("The following user-proved ASE-related results files DO NOT EXIST:\n",
               paste(c(asefs, inhmode, regpats)[!file.exists(c(asefs, inhmode, regpats))], collapse = "\n")))
  }
  ## test against ASE data related stuff
  if(!file.exists(totestf)){
    stop(paste("Test-related results file DOES NOT EXIST: ", totestf))
  }else{
    # check it has gene_id and testcolumn here
    testdat<-fread(totestf, nrows = 2, header = T)
    if(!"gene_id"%in%names(testdat) | !testcolumn%in%names(testdat)){
      stop(paste("Test data file (path: ", totestf, ") does not contain required columns gene_id and", testcolumn))
    }
  }

  # --- Get all genes to keep if gflag is 'all'
  if(gflag=="all"){
    glist<-fread(asefs[1], select = "gene_id")$gene_id
  }

  # --- Read in ASE-associated data [for genes of interest]
  dat<-rbindlist(
    lapply(names(asefs), function(strn){
    # Read in all genes
    out1<-fread(asefs[[strn]], header = T)
    out1[, strain:=strn]
    # Narrow to genes of interest
    if(gflag%in%c("all", "keepsame")){ # keep same genes for all strains
      out1<-out1[gene_id%in%glist, ]
    }else if(gflag=="bystrain"){ # keep different genes per strain
      out1<-out1[gene_id%in%glist[strain==strn, gene_id], ]
    }
  })
  )
  setkey(dat, strain, gene_id)
  ## Inform threshold
  dat[, inform:=unqalnmts.min>=informthresh]
  ## Reg pattern, inheritance mode
  im<-fread(inhmode, header = T)
  setkey(im, strain, gene_id)
  dat<-im[dat]
  rp<-fread(regpats, header = T)
  setkey(rp, strain, gene_id)
  dat<-rp[dat]
  rm(im, rp) # just a lil clean up

  # --- Add in 'test' data [for genes of interest] - keys & joins by strain if it's there
  tdat<-fread(totestf, header = T)
  if("strain"%in%names(tdat)){
    tdat<-tdat[, .(gene_id, strain, get(testcolumn))]
    setnames(tdat, c("gene_id", "strain", "testDat"))
    setkey(tdat, strain, gene_id)
    dat<-tdat[dat] # this retains all genes in dat. Doing internal to if/else for correct joining
  }else{
    tdat<-tdat[, .(gene_id, get(testcolumn))]
    setnames(tdat, c("gene_id", "testDat"))
    setkey(tdat, gene_id)
    setkey(dat, gene_id)
    dat<-tdat[dat] # this retains all genes in dat. Doing internal to if/else for correct joining
  }

  # --- Return [with columns ordered nicely]
  setkey(dat, strain, gene_id)
  setcolorder(dat, c("strain", "gene_id", "inform", "regclass", "inhmode", "testDat"))
  return(dat)
}

# --- Stats/analysis
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

chisqfull<-function(datin, catrowcol = "signifAtThresholds.ASE", rowcats = c(T, F),
                    catcol = "domain", colcats = c("arm", "center", "tip"), annotdt = NULL,
                    colcats.names = colcats){
  # Origninally From chrlocenrichment_asederpim.R
  # Runs chi-sq test and reports these results as well as underlying counts, proportions in rowwise-focused way.
  #       **removes any category that are all 0s***
  # In: datin, data where rows are observations. Must have columns named with values of catrowcol, catcol
  #     catrowcol, name of column in datin that contains categories to use as rows of contingency table, output
  #     rowcats, values of catrowcol to tabulate
  #     catcol, name of column in datin that contains categories to use as columns of contingency table, output
  #     colcats, values of catcol to tabulate
  #     annotdt, optional one-row data.table with further information to include in output
  #     colcats.names, names of column categories (matching colcats) for output (may be prettier/have extra info vs. colcats alone)
  # Out: data.table with as many rows as rowcats; columns:
  #     <any in annotdt>
  #     category, rowcats value for this row
  #     <One for each of the column categories (colcats.names)> - number of observations with this row, column value
  #     <One for each of the column categories (colcats.names)>_pOfTotal - proportion of TOTAL observations with this column value - same for each row
  #     <One for each of the column categories (colcats.names)>_pOfThisRowCategory - proportion of this row category with each column value
  #     ChiSq - test statistic for this test. Same for each row!
  #     Df - Chisq test Degrees of freedom for this test. Same for each row!
  #     pvalue - ChiSq test p-value for this test. Same for each row!

  # Create contingency table & run chi-sq
  cttab<-sapply(colcats, function(x) datin[, sapply(rowcats, function(y) sum(get(catcol)==x & get(catrowcol)==y, na.rm = T))])
  colnames(cttab)<-colcats.names
  ## Remove any zeros
  rowcats.out<-rowcats[rowSums(cttab)!=0]
  cttab<-cttab[rowSums(cttab)!=0,]
  res<-chisq.test(cttab)

  # Format output counts, proportions
  propAll<-data.table(t(as.matrix(colSums(cttab)/sum(cttab)))) # proportion of TOTAL genes that each column category is (will be same across output rows)
  setnames(propAll, paste0(names(propAll), "_pOfTotal"))
  propInCat<-data.table(cttab/rowSums(cttab))
  setnames(propInCat, paste0(names(propInCat), "_pOfThisRowCategory"))

  # Annotate with other output info & return
  out<-data.table(annotdt, category = rowcats.out, cttab, propAll, propInCat,
                  ChiSq = res$statistic, Df = res$parameter, pvalue = res$p.value)
  return(out)
}

longformnsps<-function(datin, catrowcol = "signifAtThresholds.ASE", rowcats = c(T, F),
                       catcol = "domain", colcats = c("arm", "center", "tip"), annotdt = NULL,
                       catrowcol.name = catrowcol, catcol.name = catcol){
  # Originally From chrlocenrichment_asederpim.R
  # Gets categories x chromosome domain numbers, proportions. In long format oriented way - one row per
  # In: datin, data where rows are observations. Must have columns named with values of catrowcol, catcol
  #     catrowcol, name of column in datin that contains categories to use as rows of contingency table, output
  #     rowcats, values of catrowcol to tabulate
  #     catcol, name of column in datin that contains categories to use as columns of contingency table, output
  #     colcats, values of catcol to tabulate
  #     annotdt, optional one-row data.table with further information to include in output
  #     catrowcol.name, categ1name for output (might be prettier than column name) - corresponds to catrowcol
  #     catcol.name, categ2name for output (might be prettier than column name) - corresponds to catcol
  # Out: data.table with one row per combination of rowcat and column cat. (e.g., TRUE ASE and center domain)
  # ***all Ns are sum NON NA!
  #        columns (last few not in order for ease of documentation)::
  #       <any in annotdt>
  #       categ1name, catrowcol.name input - what is name of the category 1 column
  #       category1, Actual category of categ1 that this row describes (e.g. TRUE for ASE)
  #       categ2name, catcol.name input - what is name of the category 2 column
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
      data.table(categ1name = catrowcol.name,
                 category1 = y,
                 categ2name = catcol.name,
                 category2 = x,
                 total.n = datin[, sum(!is.na(get(catrowcol)) & !is.na(get(catcol)))], # total NON-NA for both now!!
                 total.thiscateg1 = datin[, sum(get(catrowcol)==y &!is.na(get(catcol)), na.rm = T)], # needs to be non-NA in BOTH to be counted
                 total.thiscateg2 = datin[, sum(get(catcol)==x & !is.na(get(catrowcol)), na.rm = T)],
                 n.thiscombo = datin[, sum(get(catcol)==x & get(catrowcol)==y, na.rm = T)])
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

modtestbinom<-function(datin, catrowcol = "combinedrp", rowcats.ref = c('conserved', 'cis', 'trans'),
                       catcol = "testDat", annotdt = NULL){
  # Does a binomial generalized linear model modeling T/F data (catcol) on another category column. Multiple times: for multiple reference level of factor
  #     so can get many comparisons
  # In: datin, data.table with all data of interest, specifically catrowcol and catcol
  #     catrowcol, character name of column that has 'explanatory' category variable
  #     rowcats.ref, values of catrowcol **to use as reference level** in tests: model will be run one time for each (data re-leveled each time)
  #     catcol, character name of column in datin that has T/F data to model (outcome variable)
  #     annotdt, optional one-row data.table with further information to include in outputs
  # Out: list of data.tables:
  #         modinfo, one-row data table describing general model fit. Columns:
  #             <any in annotDT>
  #             NullDeviance, null deviance for binomial model
  #             NullDF, DF on null deviance
  #             ResidualDeviance, residual deviance for overall binomial model
  #             ResidualDF, DF on resiudal deviance
  #             AIC, AIC for overall binomial model
  #         allmodres, data.table with coefficients/betas & stats on them for all comparisons of catrowcol vs reference values in rowcats.ref. Only non-overlapping comparisons output
  #             (i.e., results where reference category and tested category are swaps of each other not included).
  #             Columns:
  #             <any in annotDT>
  #             reference.category, Category (of catrowcol) that is set as reference for this row
  #             tested.category, Category (of catrowcol) tested against reference category here - beta is for this vs. reference category
  #             beta.estimate, Estimate (beta) from binomial glm
  #             beta.stderr, Std. Error on beta from binomial glm
  #             z.value, z value test statistic for beta from binomial glm
  #             p.value, p value (pr >|z|) from binomail glm
  #             p.bonf, Bonferroni-corrected (multi-test adjusted) p value: p.value divided by the number of tests here (each comparison included once/inverses not included)

  # --- Set up data
  datmod<-copy(datin[, c(catcol, catrowcol), with = F])
  setnames(datmod, c("y", "x"))

  # --- Run model with EACH category as 'reference' level (so can see difference in all pairs!)
  mods<-lapply(rowcats.ref[rowcats.ref%in%datmod$x], function(reflev){
    # Format data
    allvals<-as.character(datmod[, unique(x)])
    valsord<-c(reflev, allvals[-which(allvals==reflev)])
    datmod[, x:=factor(x, levels = valsord)]
    # Run model
    mod<-glm(y~x, family = "binomial", data = datmod)
    # Return
    return(mod)
  })
  names(mods)<-rowcats.ref[rowcats.ref%in%datmod$x]

  # --- Extract overall model info
  modinfo<-data.table(NullDeviance = mods[[1]]$deviance,
                      NullDF = mods[[1]]$df.null,
                      ResidualDeviance = mods[[1]]$null.deviance,
                      ResidualDF =  mods[[1]]$df.residual,
                      AIC = mods[[1]]$aic)

  # --- Extract effects for each category vs. each specified reference level
  allcoefs<-rbindlist(lapply(names(mods), function(reflev){
    out<-data.table(coef(summary(mods[[reflev]])), keep.rownames = T)
    # Remove Intercept row
    out<-out[rn!="(Intercept)", ]
    # Rename all columns as needed
    setnames(out, c("rn", "beta.estimate", "beta.stderr", "z.value", "p.value"))
    # Rename tested categories nicely
    out[, tested.category:=out[, tstrsplit(rn, "x")]$V2]
    out[, rn:=NULL]
    # Add final column, order correctly
    out[, reference.category:=reflev]
    setcolorder(out, c("reference.category", "tested.category"))

    # Return
    return(out)
  }))
  # *** narrow to UNIQUE tests - don't keep ones that are just directional inverses of each other (e.g., cis vs conserved and conserved vs cis)
  torm<-c()
  for(i in 1:(nrow(allcoefs) - 1)){
    for(j in (i+1):nrow(allcoefs)){
      if(allcoefs[j, paste(.(tested.category, reference.category), collapse = "-")]==allcoefs[i, paste(.(reference.category, tested.category), collapse = "-")]){
        torm<-c(torm, j)
      } # end if categories are reverses of each other
    }
  }
  if(length(torm)>0){
    allcoefs<-allcoefs[-torm, ]
  }

  # Add Bonferroni-adjusted p-value: for the number of among-level comparisons done here. (correct because same but inverse-direction tests removed)
  allcoefs[, p.bonf:=p.value*nrow(allcoefs)]
  allcoefs[p.bonf>1, p.bonf:=1] # don't retain p-values above 1, lol

  # --- Return
  return(list(modinfo = data.table(annotdt, modinfo),
              allmodres = data.table(annotdt, allcoefs)))
}

# --- Plotting
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

stackedbar<-function(pinhclass, mycolors, xcol = "strain", stackcol = "propName", stacknumcol = "prop", xorder = NA,
                     legendlabel = "", myxlab = "Strain", myylab = "Proportion genes", mytitle = "", mysubt = ""){
  # Modified from chrlocenrichment_asederpim.R
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

divplot<-function(onedat, xcol = "nvars.thou", ycol = "prop", xlabel = "SNVs/INDELs vs. N2 (thousands)",
                  ylabel = "", mytitle = "", mysubt = "", mysize = 1, ymincol = "lowci", ymaxcol = "highci"){
  # # Modified from chrlocenrichment_asederpim.R; which was originally One off just to make this plot. Copied from ase_de_cistransclassification.R
  plt<-ggplot(onedat, aes(eval(as.name(xcol)), eval(as.name(ycol)))) +
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
p<-arg_parser("Multipurpose script to take any per-gene 'phenotype' (categorical or quantitative) and plot and test it against ASE-derived categories (ASE/not, DE/not, regulatory pattern, inheritance mode)",
              name = "aseetc_vs_general.R", hide.opts = TRUE)

# ASE-related data inputs
p<-add_argument(p, "--resmatinput",
                help = "[[ASE-related data input]]
                Example filepath to ASE results input - output of ase_de_annotategenes_deseq2_fromemaseout.R (see that script's documentation for format details).
                **Where strain is identified in filepath, write STRAIN e.g. STRAIN_annotatedASEDEresults.txt.gz",
                type = "character")
p<-add_argument(p, "--regpats",
                help = "[[ASE-related data input]]
                Regulatory pattern classifications for all genes, strains for which classification could be made - columns strain, gene_id, regclass.
                (*_regpattern_per1plusunqalnggenestrain_<conf int>confdeaseoverlap.txt.gz output of ase_de_cistransclassifications.R)",
                type = "character")
p<-add_argument(p, "--inhmode",
                help = "[[ASE-related data input]]
                Inheritance mode classifications for all genes and strains. Columns strain, gene_id, inhmode
                (*_inheritancemode_pergenestrain.txt.gz output from f1_parental_comparisons_withinstrain.R)",
                type = "character")
p<-add_argument(p, "--aseinfo",
                help = "[[ASE-related data input]]
                File containing information on ASE strains to use to generate ASE-relevant outputs. One row per WILD strain crossed with reference.
                Columns: isotype, strain [wild strain - as in resmatinput naming], refparent [ref strain], nvars [number of variants differentiating wild strain and reference strain].
                Wild strains ordered how you want to plot them!!",
                type = "character")
p<-add_argument(p, "--informthresh",
                help = "[[ASE-related data input]]
                Gene must have this or more unique alignments in each sample to be considered informative for ASE/cis-trans analyses",
                default = 5)

# 'phenotype'-related data inputs
p<-add_argument(p, "--genelist",
                help = "[[gene-related data input]]
                OPTIONAL Filepath to list of genes to narrow to if desired [if omitted, all genes in --resmatinput will be used].
                Can be no-header list of gene_ids (same genes retained for all strains)
                OR two-column file with columns strain, gene_id (genes will be retained in strain-specific manner)",
                default = "",
                type = "character")
p<-add_argument(p, "--genesdescrip",
                help = "[[gene-related data input]]
                OPTIONAL description of the genes in genelist to pass through to outputs, plots (descriptive). Default is 'all genes' if no gene list provided",
                default = "all genes")
p<-add_argument(p, "--totest",
                help = "[[gene-related data input]]
                Path to file containing the data to test against ASE etc! Must have columns gene_id, --testcolumn value;
                if it has column strain, data will be matched in strain-specific manner; if not, data will be assigned multiply to each strain",
                type = "character")
p<-add_argument(p, "--testcolumn",
                help = "[[gene-related data input]]
                Column name in --totest file that contains the data of interest to test ASE etc against. E.g., Pi if you're looking at nucleotide diversity via a column tilted Pi in --totest file",
                type = "character")
p<-add_argument(p, "--testdatatype",
                help = "[[gene-related data input]]
                Type of data in totest file, --testcolumn value: quantitative or categorical (binary T/F!). Acceptable values: [quantitative, categorical]",
                default = "quantitative")
p<-add_argument(p, "--testdatalabel",
                help = "[[gene-related data input]]
                Description of test data to be used for plot labels and the like (presumably you should also have some sort of nod to this in output names/directories)",
                type = "character")

# Output-related arguments
p<-add_argument(p, "--outdir",
                help = "[[output-related]] Output directory. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character",
                default = "out")
p<-add_argument(p, "--baseoutname",
                help = "[[output-related]] Base name for all output files",
                type = "character",
                default = "out")

# Misc/plotting-related arguments
p<-add_argument(p, "--genomebp",
                help = "[[misc/plotting related]]
                Length of reference genome in bp. Used to get variants vs. reference per kb.
                Default is genome length from NCBI https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6/ -
                for ws235 but unlikely to have changed much (and not enough to be meaningful for our purposes)",
                default = 100286401)

# Parse
p<-parse_args(p)

#### Get in all data ####
cat("...Reading in and formatting data....\n")
# --- output directory
# Output directory
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}

# --- Set up strain info, genes to retain info
# Strain info
if(file.exists(p$aseinfo)){
  straininfo<-fread(p$aseinfo)
  if(!all(c("isotype", "strain", "refparent") %in% names(straininfo))){
    stop(paste("--aseinfo file", p$aseinfo, "does not have the appropriate column names; fix and re-try!"))
  }
}else{
  stop(paste("--aseinfo file", p$aseinfo, "does not exist; fix this filepath and re-try!"))
}

# Genes
gskeep<-genesfilter(p$genelist)
if(gskeep$gflag=="all"){
  msg<-"Retaining all genes in input ASE-related data for analysis."
}else if(gskeep$gflag=="keepsame"){
  msg<-paste("Retaining", length(gskeep$glist), "genes for analysis from file", p$genelist)
}else if(gskeep$gflag=="bystrain"){
  msg<-paste("Retaining strain-specific genes for analysis from file", p$genelist)
}

cat(paste("....", msg, "....\n"))

# --- ASE etc data
# Data CHECKS: stop before get too far [if not a check that would otherwise be implemented soon]
if(!p$testdatatype%in%c("quantitative", "categorical")){
  stop(paste("--testdatatype provided is", p$testdatatype, "- NOT ACCEPTED; type needs to be specified as 'quantitative' or 'categorical'"))
}

# Read in
dat<-freadcombdat(strains = straininfo$strain, exampaseinput = p$resmatinput, inhmode = p$inhmode, regpats = p$regpats,
                  gflag = gskeep$gflag, glist = gskeep$glist, totestf = p$totest, testcolumn = p$testcolumn, informthresh = p$informthresh)

# add 'all' strain: REPEAT each gene ***that is present in all strains ***; inform only if ALL internal strains are inform & inform = F only if NO internal strains are inform; hypdiv/lowcov/etc if that in ANY strain
alldat<-copy(dat[gene_id%in% dat[,.N, by = gene_id][N==nrow(straininfo), gene_id] ,]) # Keep only if have as many observations as strains!
setkey(alldat, gene_id)
alldat[, strain:="all"] # still have one observation per strain, but going to pool them in this fake meta-strain
## informative-ness: needs a couple steps to make sure done correctly
forinform<-alldat[, .(allstrainsinform = sum(inform)==nrow(straininfo), nostrainsinform = sum(inform==F)==nrow(straininfo)), by = gene_id]
alldat[forinform[allstrainsinform==T, gene_id], informNEW:=T]
alldat[forinform[nostrainsinform==T, gene_id], informNEW:=F]
alldat[, inform:=informNEW]
alldat[, informNEW:=NULL]
## hypdiv
alldat[, `:=`(hypdiv = ifelse(sum(hypdiv) > 0, T, F),
              lowDNACov = ifelse(sum(lowDNACov) > 0, T, F),
              highDNACov = ifelse(sum(highDNACov) > 0, T, F)),
       by = gene_id]
dat<-rbind(dat, alldat)
rm(alldat)

strainswall<-c("all", straininfo$strain)

#### Plotting and categories etc set up: come back here to update this! ####
dat[, strain:=factor(strain, levels = strainswall)]

# Regulatory patterns, inheritance mode **in order of interest** ** change if change something upstream **
regpats<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "compensating", "compensatory", "overcompensating")
inhmodes<-c("no_change", "ambiguous", "additive", paste0(straininfo[1, refparent], "_dominant"), "alt_dominant", "overdominant", "underdominant")
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
inhcols[c("additive", paste0(straininfo[1, refparent], "_dominant"), "alt_dominant", "overdominant", "underdominant")]<-rev(brewer.pal(5, "Paired"))

# Set up other color stuff
infcols<-c("gray", "blue")
names(infcols)<-c(FALSE, TRUE)
asecols<-c("gray", "red") # ASE yes/no
names(asecols)<-c(FALSE, TRUE)

# COMBINED regulatory patterns: combining all compensated together
## **add to data
cistransop<-c("compensating", "compensatory", "overcompensating")
dat[, combinedrp:=regclass]
dat[combinedrp%in%cistransop, combinedrp:="cis-trans opposing"]
### including T or F for gene is cis vs not: NA if not cis/cis enhancing or cis-trans opposing
dat[, ciscompensated:=NA]
dat[combinedrp=="cis-trans opposing", ciscompensated:=T]
dat[combinedrp%in%c("cis", "enhancing"), ciscompensated:=F]
setcolorder(dat, c("strain", "gene_id", "inform", "regclass", "combinedrp", "ciscompensated"))
## colors
combinedrps<-c("conserved", "ambiguous", "cis", "enhancing", "trans", "cis-trans opposing")
combinedrpcols<-c(na.omit(regcols[combinedrps]), "cis-trans opposing" = "#4292C6")

# --- Comparisons to do as data.table(s) so can loop through them (basically, x axis of comparison: y axis will always be test data provided in inputs)
# Two way comparison - ASE vs. not - style comparisons
twowaytests<-data.table(myname = c("ASE informative vs uninformative", "ASE vs not", "DE vs not", "DE vs not  (ASE informative genes)",
                                   "Cis-regulatory change: compensated in trans or not"),
                        shortname = c("inform", "ASE", "DE", "DEinform", "ciscompensated"), # for file etc out
                        myxname = c("ASE informative", "ASE called", "DE called", "DE called (of ASE informative genes)", "Cis change is compensated"), # for plotting
                        columnname = c("inform", "signifAtThresholds.ASE", paste0("signifAtThresholds.ParentVsParent", straininfo$refparent[1]), paste0("signifAtThresholds.ParentVsParent", straininfo$refparent[1]),
                                     "ciscompensated"), # data test: column with T/Fs for this test
                        mynarrow = c("inform==T | inform==F", "inform==T", "inform==T|inform==F", "inform==T",
                                     "inform==T & !is.na(ciscompensated)"),  # narrow data further than it already is for this test?
                        colorvec = c("asecols", "asecols", "asecols", "asecols", "asecols") # what colors to use: so far always red/black for T/F
                        )

# Multi-way comparison - regulatory pattern-style comparisons
multiwaytests<-data.table(myname = c("Regulatory pattern (all)", "Regulatory pattern (simplified)", "Inheritance mode", "Inheritance mode (ASE informative genes)"), # long/pretty description
                          shortname = c("regpatall", "regpatsimple", "inhmode", "inhmodeinformonly"), # short name for file paths etc
                          myxname = c("", "", "", ""), # filler - so matches twowaytests (sometimes that's helpful)
                          columnname = c("regclass", "combinedrp", "inhmode", "inhmode"), # column name in dat
                          mynarrow = c("inform==T", "inform==T & !is.na(combinedrp)", "inform==T | inform==F", "inform==T"), # narrow data further than it already is for this test?
                          colorvec = c("regcols", "combinedrpcols", "inhcols", "inhcols") # colors to use
                          )

# DIRECTIONALITY important comparisons ?! **** figure out how to handle this ***

#### Analyses if data is CONTINUOUS/QUANTITATIVE ####
if(p$testdatatype=="quantitative"){
  cat("....Data is quantitative (per user), doing quantitative data vs. ASE etc analyses....\n")

  # --- Mann-Whitney type data, tests
  cat("--....Doing two-way Mann-Whitney type analyses of quantitative test data vs. various gene categorizations....\n")
  # Perform tests etc
  mwout<-lapply(1:nrow(twowaytests), function(ind){
    # Test
    mymw<-rbindlist(lapply(strainswall, function(strn){
      thisdat<-dat[strain==strn & eval(parse(text = twowaytests[ind, mynarrow])), ]
      res<-wilcox.test(thisdat[get(twowaytests[ind, columnname])==T, testDat],
                       thisdat[get(twowaytests[ind, columnname])==F, testDat])
      out<-data.table(strain = strn,
                      test = twowaytests[ind, myname],
                      genes = p$genesdescrip,
                      med.T = thisdat[get(twowaytests[ind, columnname])==T, median(testDat, na.rm = T)],
                      med.F = thisdat[get(twowaytests[ind, columnname])==F, median(testDat, na.rm = T)],
                      T.greater = thisdat[get(twowaytests[ind, columnname])==T, median(testDat, na.rm = T)] >
                        thisdat[get(twowaytests[ind, columnname])==F, median(testDat, na.rm = T)],
                      ngenes.T = thisdat[, sum(get(twowaytests[ind, columnname])==T, na.rm = T)],
                      ngenes.F = thisdat[, sum(get(twowaytests[ind, columnname])==F, na.rm = T)],
                      W = res$statistic,
                      wilcox.p.value = res$p.value)

      return(out)
    }))

    # Plot
    ## Format labels
    mymw[, label:= paste(ifelse(T.greater, paste(twowaytests[ind, shortname], "> non; p ="),
                                paste(twowaytests[ind, shortname], "< non; p =")),
                         ifelse(wilcox.p.value < 0.001, paste(formattable::scientific(wilcox.p.value, digits = 1)),
                                round(wilcox.p.value, digits = 3)))]
    ## Make plot
    plt<-sinawmed(datin = dat[eval(parse(text = twowaytests[ind, mynarrow])),],
                  xcol = twowaytests[ind, columnname], ycol = "testDat", facrow = "strain", faccol = "strain",
                  colorcol = twowaytests[ind, columnname], colorvec = eval(as.name(twowaytests[ind, colorvec])),
                  boxcol = rgb(0, 0, 1, 0.6), faclabels = mymw, myxlab= twowaytests[ind, myxname],
                  myylab = p$testdatalabel, mytitle = twowaytests[ind, myname],
                  mysubt = p$genesdescrip, myscales = "fixed", facvec = strainswall)

    # Return
    return(list(mw = mymw[, .SD, .SDcols = -c(ncol(mymw))], plt = plt))
  })

  # Save test results
  write.table(rbindlist(lapply(mwout, function(x) x$mw)),
             file.path(p$outdir, paste0(p$baseoutname, "_twowaytests_mannwhitneyresults.txt")),
             sep = "\t", quote = F, row.names = F)

  # Save plots (in one PDF)
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_twowaytests_sinaplots.pdf")),
      max(8, 1.75*length(strainswall)), max(4, 1.6*length(strainswall)))
  invisible(lapply(mwout, function(x) print(x$plt)))
  invisible(dev.off())

  # --- ANOVA type tests (multi-way/multi-category)
  cat("--....Doing multiway ANOVA-type analyses of quantitative test data vs. various gene categorizations....\n")
  # Perform tests & make plots
  multwayout<-lapply(1:nrow(multiwaytests), function(ind){
    myaov<-lapply(strainswall, function(strn){
     allouts<-anovatuk(dat = dat, mystrain = strn, testinforow = multiwaytests[ind, ],
                       testagainstcol = "testDat", colorvec = eval(as.name(multiwaytests[ind, colorvec])))

      # Format for return
      anout<-data.table(strain = strn, category = multiwaytests[ind, myname],
                        genes = p$genesdescrip,
                        allouts$anout)
      tukout<-data.table(strain = strn,  category = multiwaytests[ind, myname],
                         genes = p$genesdescrip,
                         allouts$tukout)
      tuklabs<-data.table(strain = strn,  category = multiwaytests[ind, myname],
                          genes = p$genesdescrip,
                          allouts$tuklabs)
      ns<-data.table(strain = strn,  category = multiwaytests[ind, myname],
                     genes = p$genesdescrip,
                     allouts$ns)

      return(list(antab = anout, tukres = tukout, tuklabs = tuklabs, ns = ns))
    }) # end lapply over strains

    # Format test results together
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
    anplt<-sinawmed(datin = dat[eval(parse(text = multiwaytests[ind, mynarrow])) , ],
                    xcol = multiwaytests[ind, columnname], ycol = "testDat", facrow = "strain",
                    faccol = "strain", colorcol = multiwaytests[ind, columnname],
                    colorvec = eval(as.name(multiwaytests[ind, colorvec])),
                    boxcol = rgb(0, 0, 0, 0.6), faclabels = antab, myxlab = multiwaytests[ind, myname],
                    myylab = p$testdatalabel, mytitle =multiwaytests[ind, myname],
                    mysubt = p$genesdescrip, myscales = "fixed", facvec = strainswall) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ### log10 y axes
    anplt.log10<-sinawmed(datin = dat[eval(parse(text = multiwaytests[ind, mynarrow])) , ],
                    xcol = multiwaytests[ind, columnname], ycol = "testDat", facrow = "strain",
                    faccol = "strain", colorcol = multiwaytests[ind, columnname],
                    colorvec = eval(as.name(multiwaytests[ind, colorvec])),
                    boxcol = rgb(0, 0, 0, 0.6), faclabels = antab, myxlab = multiwaytests[ind, myname],
                    myylab = paste(p$testdatalabel, "(log10 scale)"), mytitle =multiwaytests[ind, myname],
                    mysubt = p$genesdescrip, myscales = "fixed", facvec = strainswall) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      scale_y_continuous(trans = "log10")

    ## Tukey significant labels
    ### Absolute y axes
    tukplt<-sinawmed(datin = dat[eval(parse(text = multiwaytests[ind, mynarrow])) , ],
                     xcol = multiwaytests[ind, columnname], ycol = "testDat", facrow = "strain",
                     faccol = "strain", colorcol = multiwaytests[ind, columnname],
                     colorvec = eval(as.name(multiwaytests[ind, colorvec])),
                     boxcol = rgb(0, 0, 0, 0.6), faclabels = F, myxlab = multiwaytests[ind, myname],
                     myylab = p$testdatalabel, mytitle =multiwaytests[ind, myname],
                     mysubt = p$genesdescrip, myscales = "fixed") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      geom_text(data = tuklabs, aes(x = eval(as.name(multiwaytests[ind, columnname])),
                                    y = Inf, label = label),
                vjust = 1.5, size = 5, fontface = "bold")
    ### log10 y axes
    tukplt.log10<-sinawmed(datin = dat[eval(parse(text = multiwaytests[ind, mynarrow])) , ],
                     xcol = multiwaytests[ind, columnname], ycol = "testDat", facrow = "strain",
                     faccol = "strain", colorcol = multiwaytests[ind, columnname],
                     colorvec = eval(as.name(multiwaytests[ind, colorvec])),
                     boxcol = rgb(0, 0, 0, 0.6), faclabels = F, myxlab = multiwaytests[ind, myname],
                     myylab = paste(p$testdatalabel, "(log10 scale)"), mytitle =multiwaytests[ind, myname],
                     mysubt = p$genesdescrip, myscales = "fixed") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      geom_text(data = tuklabs, aes(x = eval(as.name(multiwaytests[ind, columnname])),
                                    y = Inf, label = label),
                vjust = 1.5, size = 5, fontface = "bold") +
      scale_y_continuous(trans = "log10")

    # Save plots
    pdf(file.path(p$outdir, paste0(p$baseoutname, "_sinaplots_", multiwaytests[ind, shortname], ".pdf")),
        width = max(5, 1.6*length(strainswall)), height = max(5, 1.6*length(strainswall)))
    print(anplt)
    print(anplt.log10)
    print(tukplt)
    print(tukplt.log10)
    invisible(dev.off())

    # --- Return
    antab[, label:=NULL] # get rid of label column
    return(list(antab = antab, tukres = tukres, ns = ns))
  }) # end lapply over rows of multiwaytests

  # Save test results
  names(multwayout)<-multiwaytests$shortname
  ## Combine
  ansave<-rbindlist(lapply(multwayout, function(x) x$antab))
  tuksave<-rbindlist(lapply(multwayout, function(x) x$tukres))
  nsave<-rbindlist(lapply(multwayout, function(x) x$ns))
  ## Write out
  write.table(ansave, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_ANOVAresults.txt")),
              sep = "\t", quote = F, row.names = F)
  write.table(tuksave, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_TukeyHSDresults.txt")),
              sep = "\t", quote = F, row.names = F)
  write.table(nsave, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_nspercategory.txt")),
              sep = "\t", quote = F, row.names = F)

  # Summarize Tukey outs and save summaries (and/or just the significant comparisons?)
  ## Just Tukey results where ANOVA is (nominally) significant, Tukey padj is significant
  ansig<-ansave[pvalue < 0.05, ]
  write.table(ansig, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_ANOVAresults_significantonly.txt")),
              sep = "\t", quote = F, row.names = F)
  tuksig<-tuksave[paste(strain, category)%in%ansig[, paste(strain, category)] & tukey.padj < 0.05, ]
  write.table(tuksig, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_TukeyHSDresults_significantonly.txt")),
              sep = "\t", quote = F, row.names = F)
  ## Deeper level summary
  setkey(tuksave, strain, category, genes, cat1, cat2)
  tuksumm<-tuksave[,.(whichgreateravg = ifelse(mean(diff) < 0, cat2, cat1), nStrainsSigTuk = sum(tukey.padj<0.05),
                              whichstrains = paste(strain[tukey.padj<0.05], collapse = ",")),
                           by = .(category, genes, cat1, cat2)] # within site class, how many strains sig contrast for each pair
  write.table(tuksumm, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_TukeyHSDresults_countsigsummary.txt")),
              sep = "\t", quote = F, row.names = F)

  cat("....quantitative data vs. ASE etc analyses complete....\n")
} # end wrapping else if datatype is quantitative

# deal with any directionality stuff....?

#### Analyses if data is CATEGORICAL ####
if(p$testdatatype=="categorical"){
  cat("....Data is categorical (per user), doing categorical data vs. ASE etc analyses....\n")

  alltests<-rbind(twowaytests, multiwaytests)

  # --- Get n, proportion data for all
  # Get
  setkey(dat, strain)
  allprops<-rbindlist(lapply(1:nrow(alltests), function(ind){
    dat[eval(parse(text = alltests[ind, mynarrow])),
        longformnsps(datin = .SD, catrowcol = alltests[ind, columnname], rowcats = names(eval(as.name(alltests[ind, colorvec]))),
                     catcol = "testDat", colcats = c(na.omit(dat[,unique(testDat)])),
                     annotdt = data.table(testedCategory = alltests[ind, myname],
                                          genes = p$genesdescrip),
                     catrowcol.name = alltests[ind, shortname], catcol.name = p$testdatalabel),
        by = strain]
  })
  )
  # Save
  write.table(allprops, file.path(p$outdir, paste0(p$baseoutname, "_alltests_longformproportions.txt")),
              sep = "\t", quote = F, row.names = F)

  # --- Get chi-sq test results for ALL
  # Get
  chisqres<-rbindlist(lapply(1:nrow(alltests), function(ind){
    dat[eval(parse(text =  alltests[ind, mynarrow])),
        chisqfull(datin = .SD, catrowcol = alltests[ind, columnname], rowcats = names(eval(as.name(alltests[ind, colorvec]))),
                  catcol = "testDat", colcats =  c(na.omit(dat[,unique(testDat)])),
                  annotdt = data.table(testedCategory = alltests[ind, myname],
                                       genes = p$genesdescrip),
                  colcats.names = paste(p$testcolumn, c(na.omit(dat[,unique(testDat)])), sep=".")),
        by = strain]
  }))
  # Save
  write.table(chisqres, file.path(p$outdir, paste0(p$baseoutname, "_alltests_chisqtestresults_percategory.txt")),
              sep = "\t", quote = F, row.names = F)

  # Summarize ??
  chisq.summ<-unique(chisqres[, .(strain, testedCategory, genes, ChiSq, Df, pvalue)])
  write.table(chisq.summ, file.path(p$outdir, paste0(p$baseoutname, "_alltests_chisqtestresults.txt")),
              sep = "\t", quote = F, row.names = F)

  # --- Get binomial model test results for MULTI-WAY -> allows to pull out effects on individual category values [in theory]
  # IF categorical data is T/F
  if(is.logical(dat$testDat)){
    cat("---....test data is logical/binary; doing binomial glms....\n")
    # Get
    ## Add information about which reference levels to use
    rp.reflevels<-c("conserved", "cis", "trans")
    rps.reflevels<-c("conserved", "cis", "trans") # want cis-trans opposing, too, but only vs cis, trans, conserved - already covered
    im.reflevels<-c("no_change", "additive", "alt_dominant", "N2_dominant")
    multiwaytests[, reflevelvec:=c("rp.reflevels", "rps.reflevels", "im.reflevels", "im.reflevels")]
    ## Run
    allmodres<-lapply(1:nrow(multiwaytests), function(ind){
      bystrn<-lapply(strainswall, function(strn){
        modtestbinom(datin = dat[eval(parse(text =  alltests[ind, mynarrow])) & strain==strn,],
                     catrowcol = multiwaytests[ind, columnname], rowcats.ref = eval(as.name(multiwaytests[ind, reflevelvec])),
                     catcol = "testDat", annotdt = data.table(strain = strn,
                                                              testedCategory = multiwaytests[ind, myname],
                                                              genes = p$genesdescrip))
      })
      names(bystrn)<-strainswall
      return(bystrn)
    })
    ## combine
    allmodinfo<-rbindlist(lapply(allmodres, function(x) rbindlist(lapply(x, function(strn) strn$modinfo))))
    allmodcoefs<-rbindlist(lapply(allmodres, function(x) rbindlist(lapply(x, function(strn) strn$allmodres))))

    # Save
    write.table(allmodinfo, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_binomialglms_info.txt")),
                sep = "\t", quote = F, row.names = F)
    write.table(allmodcoefs, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_binomialglms_allcoefs.txt")),
                sep = "\t", quote = F, row.names = F)

    # Summarize
    ## Save all significant betas where overall chi-sq result is significant
    chisig<-chisq.summ[pvalue<0.05, ]
    coefsig<-allmodcoefs[paste(strain, testedCategory) %in% chisig[, paste(strain, testedCategory)]
                         & p.bonf < 0.05,]
    write.table(coefsig, file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_binomialglms_coefs_significantonly.txt")),
                sep = "\t", quote = F, row.names = F)
    ## Counts summaries
    setkey(allmodcoefs, strain, testedCategory, genes, reference.category, tested.category)
    coefsumm<-allmodcoefs[, .(whichgreateravg = ifelse(mean(beta.estimate) < 0, reference.category, tested.category), nStrainsSigCoef = sum(p.bonf<0.05),
                              whichstrains = paste(strain[p.bonf<0.05], collapse = ",")),
                          by = .(testedCategory, genes, reference.category, tested.category)] # within category, how many strains sig contrast for eac
    write.table(coefsumm[order(nStrainsSigCoef, decreasing = T)],
                file.path(p$outdir, paste0(p$baseoutname, "_multiwaytests_binomialglms_coefs_countsigsummary.txt")),
                sep = "\t", quote = F, row.names = F)

  }else{
    cat("---....test data IS NOT logical/binary; NOT doing binomial glms....\n")
  } # end if is.logical(testdat)

  # --- Barplots: LOTS.
  cat("--....Plotting barplots of genes in test data categories that are various other categories....\n")
  # Set up directory
  barpltdir<-file.path(p$outdir, "barplots")
  if(!dir.exists(barpltdir)){dir.create(barpltdir)}
  # Update colors to look OK with bars
  asecols<-c("gray", "black") # ASE yes/no
  names(asecols)<-c(FALSE, TRUE)
  # Test data colors
  testdatcols<-colorRampPalette(c("gray", "black"))(length(allprops[, na.omit(unique(category2))]))
  names(testdatcols)<-rev(allprops[, na.omit(unique(category2))])

  # label dt (chi-sq)
  chisq.summ[, label:=paste("ChiSq p =", ifelse(pvalue < 0.001, paste(formattable::scientific(pvalue, digits = 1)),
                                                round(pvalue, digits = 3)))]

  # Make all plots
  invisible(lapply(1:nrow(alltests), function(ind){
    thisdat<-allprops[testedCategory==alltests[ind, myname], ]
    thisdat[,strain:=factor(strain, levels = strainswall)]
    thislab<-chisq.summ[testedCategory==alltests[ind, myname], ]

    # Bar plots with testDat on x axis, segmented into other category
    pdf(file.path(barpltdir, paste0(p$baseoutname, "_barplots_", alltests[ind, shortname], "_testdatax.pdf")),
        max(8, 1.75*length(strainswall)), max(4, 1.6*length(strainswall)))
    ## Proportion within testDat category: this one has chi-sq label
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = eval(as.name(alltests[ind, colorvec])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "prop.thiscateg2",
                 xorder = thisdat[,na.omit(unique(category2))], legendlabel = alltests[ind, shortname],
                 myxlab = p$testdatalabel, myylab = "Proportion of genes in x axis category",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        geom_text(size    = 4, data    = thislab, mapping = aes(x = Inf, y = Inf, label = label), hjust   = 1.05, vjust   = 1.5)
      )
    ## Global proportions
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = eval(as.name(alltests[ind, colorvec])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "prop.total",
                 xorder = thisdat[,na.omit(unique(category2))], legendlabel = alltests[ind, shortname],
                 myxlab = p$testdatalabel, myylab = "Proportion of all genes (analyzed in this set)",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    )
    ## Numbers
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = eval(as.name(alltests[ind, colorvec])),
                 xcol = "category2", stackcol = "category1", stacknumcol = "n.thiscombo",
                 xorder = thisdat[,na.omit(unique(category2))], legendlabel = alltests[ind, shortname],
                 myxlab = p$testdatalabel, myylab = "Number genes",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    )

    invisible(dev.off())

    # Bar plots with other category on x axis, segmented into testDat
    pdf(file.path(barpltdir, paste0(p$baseoutname, "_barplots_", alltests[ind, shortname], "_", alltests[ind, shortname], "x.pdf")),
        max(8, 1.75*length(strainswall)), max(4, 1.6*length(strainswall)))
    ## Proportion within testDat category: this one has chi-sq label
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = testdatcols,
                 xcol = "category1", stackcol = "category2", stacknumcol = "prop.thiscateg1",
                 xorder = names(eval(as.name(alltests[ind, colorvec]))), legendlabel = "test data",
                 myxlab = alltests[ind, myname], myylab = "Proportion of genes in x axis category",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        geom_text(size    = 4, data    = thislab, mapping = aes(x = Inf, y = Inf, label = label), hjust   = 1.05, vjust   = 1.5)
    )
    ## Global proportions
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = testdatcols,
                 xcol = "category1", stackcol = "category2", stacknumcol = "prop.total",
                 xorder = names(eval(as.name(alltests[ind, colorvec]))), legendlabel = "test data",
                 myxlab = alltests[ind, myname], myylab = "Proportion of all genes (analyzed in this set)",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    )
    ## Numbers
    print(
      stackedbar(pinhclass = thisdat,
                 mycolors = testdatcols,
                 xcol = "category1", stackcol = "category2", stacknumcol = "n.thiscombo",
                 xorder = names(eval(as.name(alltests[ind, colorvec]))), legendlabel = "test data",
                 myxlab = alltests[ind, myname], myylab = "Number genes",
                 mytitle = alltests[ind, myname], mysubt = p$genesdescrip) +
        facet_wrap(~strain) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    )
    invisible(dev.off())

  }))

  # --- Plot proportion vs. divergence
  cat("--....Plotting proportion of genes in test data categories that are various other categories vs. divergence from N2....\n")
  # Set up directory
  divdir<-file.path(p$outdir, "divergence")
  if(!dir.exists(divdir)){dir.create(divdir)}
  # Set up # variants
  ndiv<-straininfo[, .(strain, nvars)]
  ndiv[,`:=`(nvars.thou = nvars/1e03, nvars.perkb = nvars/(p$genomebp/1e03))]
  setkey(ndiv, strain)

  invisible(lapply(1:nrow(alltests), function(ind){
    # Data set up. *** exclude 'all' metastrain!!
    pdat<-allprops[testedCategory==alltests[ind, myname] & strain!="all",
                   .(strain, testedCategory, categ1name, category1, categ2name,
                     category2, n.thiscombo, prop.thiscateg1, prop.thiscateg2,
                     low95ci.thiscateg2, high95ci.thiscateg2, low95ci.thiscateg1, high95ci.thiscateg1)]
    setkey(pdat, strain)
    pdat<-pdat[ndiv]
    pdat[,strain:=factor(strain, levels = ndiv[order(nvars.thou), strain])]
    pdat[, category1:=factor(category1, levels=names(eval(as.name(alltests[ind, colorvec]))))]

    # Label set up - include short title for non-testdata
    mylabs<-paste(alltests[ind, shortname], pdat[, unique(category1)], sep = ":\n")
    names(mylabs)<- pdat[, unique(category1)]
    testlabs<-pdat[,unique(category2)]
    names(testlabs)<-pdat[,unique(category2)]

    pdf(file.path(divdir, paste0(p$baseoutname, "_propintestvsdivergencefromN2_", alltests[ind, shortname], ".pdf")),
        7.5, max(6, length(eval(as.name(alltests[ind, colorvec]))) * 1.5))
    # MAKE PLOT: Proportion of testDat categories that are each other category
    print(
      divplot(onedat = pdat, xcol = "nvars.perkb", ycol = "prop.thiscateg2",
              xlabel = "SNVs/INDELs vs. N2, per kb",
              ylabel = "Proportion of genes in given test data/column facet\nthat are the row facet's category",
              mytitle = paste(alltests[ind, myname], "- proportion of genes testDat category that are each this category"),
              mysubt = p$genesdescrip, ymincol = "low95ci.thiscateg2", ymaxcol = "high95ci.thiscateg2") +
        facet_grid(category1~category2, scales = "free_y", labeller = labeller(category1 = as_labeller(mylabs),
                                                                               category2 = as_labeller(testlabs)))
    )

    # MAKE PLOT: proportion of ASE-relevant categories that are each testDat category
    print(
      divplot(onedat = pdat, xcol = "nvars.perkb", ycol = "prop.thiscateg1",
              xlabel = "SNVs/INDELs vs. N2, per kb",
              ylabel = "Proportion of genes in given row facett\nthat are the test data (column) facet's category",
              mytitle = paste(alltests[ind, myname], "- proportion of genes testDat category that are each this category"),
              mysubt = p$genesdescrip, ymincol = "low95ci.thiscateg1", ymaxcol = "high95ci.thiscateg1") +
        facet_grid(category1~category2, scales = "free_y", labeller = labeller(category1 = as_labeller(mylabs),
                                                                               category2 = as_labeller(testlabs)))
    )

    invisible(dev.off())
  }))

  cat("....categorical data vs. ASE etc analyses complete....\n")
} # end wrapping else if datatype is categorical

#### go back and add anything w/r/t DIRECTIONALITY? (tests of wild and N2 biased stuff differently or...) ####
# ??????

#### Script completion message & session information ####
cat(".....aseetc_vs_general.R processing complete! Session information:....\n")
sessionInfo()
