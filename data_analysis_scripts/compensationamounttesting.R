#! /usr/bin/env/ Rscript
# See if can figure out if amount of compensation is surprising from proportions of trans effects, etc: NOT polarized to cis being first
# by Avery Davis Bell, begun 2025.04.03

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
require(eulerr, lib.loc = mylibloc)
require(cowplot, lib.loc = mylibloc)

# plotting theme
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

#### Functions ####
# --- data functions
dtpci<-function(x, n, suff = "", pval = F){
  # gets prop and 95% CIs as data.table
  # If pval = T, adds binomial p-value column
  if(n<=0){
    out<-data.table(p = NaN, low95ci = NaN, high95ci = NaN)
    if(pval==T){
      out[, pvalue:=NaN]
    }
  }else{
    res<-binom.test(x, n)
    out<-data.table(p = res$estimate, low95ci = res$conf.int[1], high95ci = res$conf.int[2])
    if(pval==T){
      out[, pvalue:=res$p.value]
    }
  }
  setnames(out, paste0(names(out), suff))
  return(out)
}

comptestsone<-function(nregone, annotdt = NULL){
  # Performs tests of surprisingness of compensation: cis and trans together surprisingness (including and excluding enhancing); enhancing vs cis-trans opposing
  # In:  nregone, one-row regulatory pattern counts. Columns used: conserved, trans, cis, cis-trans opposing, enhancing (maaaybe compensatory)
  #     annotdt, optional one-row data.table to add to beginning of output
  # Out: list of data.tables:
  #     $cistransov, all results from analysis of 2x2 tables comparing cis and trans and their overlap amounts. Columns:
  #       <any in annotdt>
  #       test, description of what test done here (enhancing included vs not) - same for multiple rows
  #        propshown, in proportion columns, which proportions shown (long format!)
  #        n.trans.and.cis, Number as described - same multiple rows [enhancing included/excluded as specified]
  #        n.trans.no.cis, ""
  #        n.cis.no.trans, ""
  #        n.no.cis.trans, ""
  #        or.fet, FET odds ratio for this test - same for multiple rows
  #        pvalue.fet,  FET odds ratio for this test - same for multiple rows
  #        p,  proportion as described in 'propshown'
  #        low95ci, lower 95% binomial CI bound on proportion estimate (as described in 'prop'shown)
  #        high95ci, as above but upper bound
  #   $envsopp, one-row data.table result of binomial test cis-trans opposing proportion out of all with cis and trans effects (opp + enhancing). Columns:
  #       <any in annotdt>
  #     test, "Cis-trans opposing of all cis and trans (this + enhancing)"
  #     p, proportion
  #     low95ci, binomial 95% bounds on this
  #     high95ci, 
  #     pvalue, p-value from binomial test of this vs. 50%
  
  # --- 2x2: cis-trans overlap including enhancing
  # Ns
  ovwen<-data.table(n.trans.and.cis = nregone[, enhancing + `cis-trans opposing`],
                    n.trans.no.cis = nregone[, trans],
                    n.cis.no.trans = nregone[, cis],
                    n.no.cis.trans = nregone[,  conserved])
  # Proportions...including EXPECTED?!
  p.both.oftrans<-ovwen[, dtpci(n.trans.and.cis, n.trans.and.cis + n.trans.no.cis)]
  p.both.ofcis<-ovwen[, dtpci(n.trans.and.cis, n.trans.and.cis + n.cis.no.trans)]
  p.both.ofall<-ovwen[, dtpci(n.trans.and.cis, sum(ovwen))]
  p.both.expected<-ovwen[, dtpci(round(((n.cis.no.trans + n.trans.and.cis)/sum(ovwen) * (n.trans.no.cis + n.trans.and.cis)/sum(ovwen)) * sum(ovwen)),
                                sum(ovwen))] # compute rounded expected value as prop cis * prop trans
  # FET
  ovwen.res<-fisher.test(matrix(unlist(ovwen), nrow = 2))
  # Combine
  ovwen<-data.table(test = c("Cis-trans overlap, enhancing included"),
                    propshown = c("Overlapping of trans", "Overlapping of cis", "Overlapping of all in analysis", "expected overlap (null)"),
                    rbind(ovwen, ovwen, ovwen, ovwen), 
                    or.fet = rep(ovwen.res$estimate, 4),
                    pvalue.fet = rep(ovwen.res$p.value, 4),
                    rbind(p.both.oftrans, p.both.ofcis, p.both.ofall, p.both.expected))
  
  # --- 2x2: cis-trans overlap excluding enhancing
  ovnoen<-data.table(n.trans.and.cis = nregone[, `cis-trans opposing`],
                    n.trans.no.cis = nregone[, trans],
                    n.cis.no.trans = nregone[, cis],
                    n.no.cis.trans = nregone[,  conserved])
  # Proportions...including EXPECTED?!
  p.both.oftrans<-ovnoen[, dtpci(n.trans.and.cis, n.trans.and.cis + n.trans.no.cis)]
  p.both.ofcis<-ovnoen[, dtpci(n.trans.and.cis, n.trans.and.cis + n.cis.no.trans)]
  p.both.ofall<-ovnoen[, dtpci(n.trans.and.cis, sum(ovnoen))]
  p.both.expected<-ovnoen[, dtpci(round(((n.cis.no.trans + n.trans.and.cis)/sum(ovnoen) * (n.trans.no.cis + n.trans.and.cis)/sum(ovnoen)) * sum(ovnoen)),
                                 sum(ovnoen))] # compute rounded expected value as prop cis * prop trans
  # FET
  ovnoen.res<-fisher.test(matrix(unlist(ovnoen), nrow = 2))
  # Combine
  ovnoen<-data.table(test = c("Cis-trans overlap, enhancing excluded"),
                    propshown = c("Overlapping of trans", "Overlapping of cis", "Overlapping of all in analysis", "expected overlap (null)"),
                    rbind(ovnoen, ovnoen, ovnoen, ovnoen), 
                    or.fet = rep(ovnoen.res$estimate, 4),
                    pvalue.fet = rep(ovnoen.res$p.value, 4),
                    rbind(p.both.oftrans, p.both.ofcis, p.both.ofall, p.both.expected))
  
  # --- binomial test: enhancing vs cis-trans opposing
  envsopp<-data.table(test = "Cis-trans opposing of all cis and trans (this + enhancing)",
                      nregone[, dtpci(compensatory, compensatory + enhancing, pval = T)])
  
  # --- combine into nice format & return
  return(list(cistransov = data.table(annotdt, rbind(ovwen, ovnoen)),
              envsopp = data.table(annotdt, envsopp)))
}

movegsonce<-function(nregone, percent, annotdt = NULL){
  # Moves genes in ns around categories and re-do FET to simulate missed DE/spurious ASE. Also gets prop of ASE/DE genes this would be
  #     #     Take gene away From compensatory genes, add that gene to cis only to simulate DE was missed at that gene
  #                                           , add to conserved to simulate that ASE was spuriously called at that gene 
  #  EXCLUDES enhancing from test, uses it for calculating % of ASE/DE genes wrong here
  # In: nregone, nregone, one-row regulatory pattern counts. Columns used: conserved, trans, cis, cis-trans opposing, enhancing 
  #     percent, percent of COMPENSATORY genes to move to cis-only/conserved categories (0-100)
  #     annotdt, optional one-row data.table to add to beginning of output
  # Out: Two-row data.table, each for the downsample % specified, one for simulating DE missed and one for simulating ASE spurious. Columns:
  #       <any in annotdt>
  #       test, whether numbers were moved to simuate DE missed or ASE spurious
  #       percentCompMoved, what % compensatory genes were moved in this simulation
  #       nCompMoved, what # compensatory genes were moved in this simulation
  #       propWrongDenom, For the propWrong number, what's denominator - are we looking at % DE that were missed or what
  #       propWrong, number genes moved here divided by total number of genes of category of interest (see propWrong)
  #       or.fet, FET odds ratio for this simulation
  #       pvalue.fet, FET p-value for this simulation
  
  nrm<-round(nregone[, `cis-trans opposing`*(percent/100)])
  
  # Take gene away From compensatory genes, add that gene to cis only to simulate DE was missed at that gene
  ns<-data.table(n.both = nregone[, `cis-trans opposing`] - nrm,
                 n.trans.no.cis = nregone[, trans],
                 n.cis.no.trans = nregone[, cis] + nrm,
                 n.neither = nregone[, conserved])
  res<-fisher.test(matrix(unlist(ns), nrow = 2))
  out.demiss<-data.table(test = "DE missed: take away from compensatory, add to cis only",
                         percentCompMoved = percent,
                         nCompMoved = nrm,
                         propWrongDenom = "All with DE called: trans + cis + enhancing",
                         propWrong = nregone[, nrm/(cis + trans + enhancing)], # missed of
                         or.fet = res$estimate,
                         pvalue.fet = res$p.value)
  
  # Take gene away From compensatory genes, add to conserved to simulate that ASE was spuriously called at that gene
  ns<-data.table(n.both = nregone[, `cis-trans opposing`] - nrm,
                 n.trans.no.cis = nregone[, trans],
                 n.cis.no.trans = nregone[, cis],
                 n.neither = nregone[, conserved] + nrm)
  res<-fisher.test(matrix(unlist(ns), nrow = 2))
  out.asespur<-data.table(test = "ASE spurious: take away from compensatory, add to conserved",
                         percentCompMoved = percent,
                         nCompMoved = nrm,
                         propWrongDenom = "All with ASE called: cis + enhancing + cis-trans opposing",
                         propWrong = nregone[, nrm/(cis + enhancing + `cis-trans opposing`)], # missed of
                         or.fet = res$estimate,
                         pvalue.fet = res$p.value)
  
  # Return
  return(data.table(annotdt,
                    rbind(out.demiss, out.asespur)))
}

#### Arguments & inputs ####
# --- Command line arguments
p<-arg_parser("See if can figure out if amount of compensation is surprising from proportions of trans effects, etc", 
              name = "compensationamounttesting.R ", hide.opts = TRUE)

# Organizational & data input arguments
p<-add_argument(p, "--regpatns",
                help = "Counts of regulatory patterns within strain and geneset - *numbers_regpatterns.txt output of ase_de_cistransclassifications.R.
                **MUST be for only one 'CI' column value**",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files; should likely include strain, for example",
                type = "character",
                default = "out")
p<-add_argument(p, "--varsvsref",
                help = "Path to file containing columns strain (one row for each strain in main input), nvars (number variants vs. reference genome) [other columns optional];
                for ordering strains etc",
                type = "character")
p<-add_argument(p, "--outdir",
                help = "[[output-related]] Output directory. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character",
                default = "out")


# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# Output directory
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# Read in data
straininfo<-fread(p$varsvsref, header = T)[order(nvars), ]
regns<-fread(p$regpatns, header = T)

#### Analysis plan ####
# Is the amount of overlap between cis and trans surprising? (yes) 
#     (2x2 table you propose, including enhancing, FET)
# Is this still surprising if we exclude cases where cis and trans are in same direction? (yes) 
#     (2x2 table, excludes enhancing, FET)
# When cis and trans effects both exist at a gene, does compensatory outweigh enhancing? (obv)
#      (proportion test: compensatory/cis and trans (compensatory + enhancing) vs. enhancing/cis and trans (compensatory + enhancing)
#      BETTER: just compare compensatory/cis and trans to 0.5 with simple binomial test
# Then I’ll likely re-compute how many of which effects would have to be wrong to make results not significant 
#     (by ‘moving’ genes from compensatory category to cis-only category [as if DE were missed] 
#     and from compensatory category to conserved category [as if ASE were spurious])

#### Initial tests: 2x 2 tables and prop test for ones with cis and trans ####
# --- compute proportions, results (5plusUnqAlns & 5plusUnqAlns_exclhypdivbadcov)
props.l<-lapply(c("5plusUnqAlns", "5plusUnqAlns_exclhypdivbadcov"), function(g){
  lapply(straininfo$strain, function(strn){
    nregone <- regns[strain==strn & geneset==g, ]
    out<-comptestsone(nregone,  annotdt = data.table(strain = strn, geneset = g))
    return(out)
  })
})

# Combine
cistransov<-rbindlist(lapply(props.l, function(x) rbindlist(lapply(x, function(y) y$cistransov))))
envsopp<-rbindlist(lapply(props.l, function(x) rbindlist(lapply(x, function(y) y$envsopp))))

# Pull out FET stuff to have non-repetatively
fetclean<-unique(cistransov[, .SD,.SDcols = c(1, 2, 3, 5:10)])

# Save out
write.table(props.l, file.path(p$outdir, paste0(p$baseoutname, "_overlappingcistransprops_longformat.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(fetclean, file.path(p$outdir, paste0(p$baseoutname, "_overlappingcistrans_statres.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(envsopp, file.path(p$outdir, paste0(p$baseoutname, "_enhancingVsopposing_longformat.txt")),
            sep = "\t", quote = F, row.names = F)

#### Move genes around to see how many genes I'd have to drop to be non sig ####
#     Take gene away From compensatory genes, add that gene to cis only to simulate DE was missed at that gene
#                                           , add to conserved to simulate that ASE was spuriously called at that gene         
sims<-rbindlist(lapply(seq(1, 100), function(p){
  rbindlist(lapply(c("5plusUnqAlns", "5plusUnqAlns_exclhypdivbadcov"), function(g){
    rbindlist(lapply(straininfo$strain, function(strn){
      movegsonce(nregone =  regns[strain==strn & geneset==g, ], percent = p, annotdt = data.table(strain = strn, geneset = g))
    }))
  }))
}))

# Save
write.table(sims, gzip(file.path(p$outdir, paste0(p$baseoutname, "_simulatermcompgenes_recomputecistransoverlap.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# --- Quick plot this
sims[, strain:=factor(strain, levels = straininfo[order(nvars, strain), strain])]
plt.pcomp<-ggplot(sims, aes(percentCompMoved, -1*log10(pvalue.fet))) + 
  geom_line(aes(color = strain)) + 
  xlab("Percent compensated genes moved") +
  facet_grid(geneset~test) + myggtheme +
  geom_hline(yintercept = 0.05/nrow(straininfo), lty = "dashed") + annotate("text", x = 25, y = 10, label = "p = 0.05 after Bonferroni") +
  geom_vline(xintercept = 90, lty = "dashed") + annotate("text", x = 90, y = 250, label = "90%")

plt.demiss<-ggplot(sims[test=="DE missed: take away from compensatory, add to cis only"], aes(propWrong*100, -1*log10(pvalue.fet))) + 
  geom_line(aes(color = strain)) + 
  xlab("Percent DE calls 'missed'") +
  facet_grid(geneset~test) + myggtheme +
  geom_hline(yintercept = 0.05/nrow(straininfo), lty = "dashed") + annotate("text", x = 10, y = 10, label = "p = 0.05 after Bonferroni") +
  geom_vline(xintercept = c(30, 40), lty = "dashed") + annotate("text", x = 30, y = 250, label = "30%") + annotate("text", x = 40, y = 250, label = "40%") 

plt.asespur<-ggplot(sims[test=="ASE spurious: take away from compensatory, add to conserved"], aes(propWrong*100, -1*log10(pvalue.fet))) + 
  geom_line(aes(color = strain)) + 
  xlab("Percent ASE calls spurious") +
  facet_grid(geneset~test) + myggtheme +
  geom_hline(yintercept = 0.05/nrow(straininfo), lty = "dashed") + annotate("text", x = 10, y = 10, label = "p = 0.05 after Bonferroni") +
  geom_vline(xintercept = c(40, 50), lty = "dashed") + annotate("text", x = 50, y = 250, label = "50%") + annotate("text", x = 40, y = 250, label = "40%") 

pdf(file.path(p$outdir, paste0(p$baseoutname, "_simulatermcompgenes_plots.pdf")), 12 , 8)
print(plt.pcomp)
print(plt.demiss)
print(plt.asespur)
invisible(dev.off())

#### Plot data ####
# --- colors etc
twocols<-c("gray", "black") # ASE yes/no
names(twocols)<-c(FALSE, TRUE)

cistransov[, strain:=factor(strain, levels = straininfo[order(nvars, strain), strain])]
envsopp[, strain:=factor(strain, levels = straininfo[order(nvars, strain), strain])]

# --- Plot prop of trans and cis that are opposed: doesn't show how surprising this is but does display data. NB the 'compensated' genes are the same for both columns
#     With cis-trans as facets, strains on x axis
plt.pct.split<-ggplot(data = cistransov[propshown%in%c("Overlapping of trans", "Overlapping of cis") & geneset=="5plusUnqAlns"]) +
  geom_pointrange(aes(x = strain, y = p, ymin = low95ci, ymax = high95ci, color = strain)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +ylab("Proportion of genes with this regulation\nthat have opposing effect from other regulation") +
  xlab("") +
  facet_grid(test ~ propshown) +
  myggtheme 
pdf(file.path(p$outdir, paste0(p$baseoutname, "_propOppOfCisTrans_facets.pdf")), 10 , 8)
print(plt.pct.split)
invisible(dev.off())

# --- Plot prop of trans and cis that are opposed: doesn't show how surprising this is but does display data. NB the 'compensated' genes are the same for both columns
plt.cistrans<-ggplot(data = cistransov[propshown%in%c("Overlapping of trans", "Overlapping of cis")]) +
  geom_pointrange(aes(x = as.numeric(as.factor(propshown))*10 + as.numeric(strain), y = p, ymin = low95ci, ymax = high95ci, color = strain)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(10 + nrow(straininfo)/2, 20 + nrow(straininfo)/2), labels = c("Genes with cis effects", "Genes with trans effects")) +
  xlab("") + ylab("Proportion with other effect too") + 
  facet_grid(geneset~test) +
  myggtheme 
pdf(file.path(p$outdir, paste0(p$baseoutname, "_propOppOfCisTrans_together.pdf")), 10 , 8)
print(plt.cistrans)
invisible(dev.off())

# --- Plot proportion of all genes we'd expect to have cis & trans effect vs. proportion that have opposing cis & trans effects [and any double cis & trans effects?]
plt.expobs<-ggplot(data = cistransov[propshown%in%c("Overlapping of all in analysis", "expected overlap (null)")]) +
  geom_pointrange(aes(x = as.numeric(as.factor(propshown))*10 + as.numeric(strain), y = p, ymin = low95ci, ymax = high95ci, color = strain)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(10 + nrow(straininfo)/2, 20 + nrow(straininfo)/2), labels = c("Expected cis-trans overlap", "Observed cis-trans overlap")) +
  xlab("") + ylab("Proportion genes with cis and trans effects") + 
  facet_grid(geneset~test) +
  myggtheme 
pdf(file.path(p$outdir, paste0(p$baseoutname, "_propExpObsOverlap.pdf")), 10 , 8)
print(plt.expobs)
invisible(dev.off())

# --- Plot proportion of genes with cis and trans effects that have opposing effects (probably best in combination with another type?)
envsopp.p<-copy(envsopp)
envsopp.p[, `:=`(p=1-p)]
envsopp.p<-rbind(envsopp, envsopp.p)
envsopp.p[, opposing:=rep(c(T, F), each = nrow(envsopp))]
envsopp.p[, strain:=factor(strain, levels = straininfo[order(nvars, strain), strain])]
plt.envopp<-stackedbar(pinhclass = envsopp.p, mycolors = twocols, xcol = "strain", stackcol = "opposing", stacknumcol = "p",
                       myylab = "Proportion genes (all genes have cis and trans influences)", myleglabels = c("Enhancing", "Opposing")) +
  facet_grid(geneset~test)
pdf(file.path(p$outdir, paste0(p$baseoutname, "_propOppOfCisandTrans.pdf")), 8 , 8)
print(plt.envopp)
invisible(dev.off())

# --- Try exp vs observed where exp looks more like test statistic - sub-part of a bar?
# Proportions
plt.eobar<-ggplot(cistransov[propshown=="Overlapping of all in analysis"])

# NUMBERS!!

# --- Compare prop of genes with cis effects OVERALL vs. at genes with trans influences [opposing & all]; 
# OK new idea that may let us make the old plot but with the new data, ish. Each number below is a bar on the plot
# Section/facet 1: cis genes as compensated
# Of total genes, proportion with/without trans effects (or of genes WITHOUT cis effects?)
# Of genes WITH cis effects, proportion with/without trans effects (opposing or at all and then split out to show how its all opposing)
# Section/facet 2: trans genes as compensated
# Of total genes, proportion  with/without cis effects (or of genes WITHOUT trans effects?)
# Of genes WITH trans effects, proportion with/without cis effects  (opposing or at all and then split out to show how its all opposing)

# Format data: need p with each effect of genes without the other effect
ctclean<-unique(cistransov[propshown%in%c("Overlapping of trans", "Overlapping of cis")])
    # Add: cis of NOT trans (n.cis.no.trans/(n.cis.no.trans + n.no.cis.trans))
    #       trans of NOT cis (n.trans.no.cis/(n.trans.no.cis + n.no.cis.trans))
newct<-rbindlist(lapply(1:nrow(ctclean), function(i){
  pcnott<-dtpci(x = ctclean[i, n.cis.no.trans], n = ctclean[i, n.cis.no.trans + n.no.cis.trans])
  ptnotc<-dtpci(x = ctclean[i, n.trans.no.cis], n = ctclean[i, n.trans.no.cis + n.no.cis.trans])
  return(
  data.table(ctclean[i, .(strain, geneset, test)],
             propshown = c("cis of not trans", "trans of not cis"),
             ctclean[i, .(n.trans.and.cis, n.trans.no.cis, n.cis.no.trans, n.no.cis.trans, or.fet, pvalue.fet)],
             rbind(pcnott, ptnotc))
  )
}))
ctclean<-rbind(ctclean, newct) # 3x as long as ctclean orig

newct.p<-copy(ctclean)
newct.p[, p:=1-p]
newct.p<-rbind(copy(ctclean), newct.p)
newct.p[, othereffect:=rep(c(T, F), each = nrow(ctclean))]
newct.p[, propshown:=factor(propshown, levels = c("trans of not cis", "Overlapping of cis", "cis of not trans", "Overlapping of trans"))]

## ** I AM GETTING DUPLICATES IN HERE SOMEHOW, FIX THAT ***

# Make the plot
plt.ctvnot<-stackedbar(pinhclass = unique(newct.p[geneset=="5plusUnqAlns"]), mycolors = twocols, xcol = "propshown", stackcol = "othereffect",
                       stacknumcol = "p", myylab = "Proportion genes", myleglabels = c("No other effect", "Yes other effect"),
                       myxlab = "") +
  facet_grid(strain~test) 
pdf(file.path(p$outdir, paste0(p$baseoutname, "_propCisTransAtNotCisTrans.pdf")), 12 , 12)
print(plt.ctvnot)
invisible(dev.off())

# 
# # --- Try: eulerr showing all; cis/trans overlap; enh vs opposing breakdown
# # Practice: getting cis & trans overlap to be further characterized is the tricky bit
# test3<-list(all=c("a", "b", "c", "d", "e", "f", "g"), cis=c("a", "b", "c"), trans=c("b", "c", "d"), opp=c("c"), enh = c("b"))
# ftest3<-euler(test3)
# plot(ftest3)
# 
# # Practice: numbers
# tdat<-regns[strain=="ECA722" & geneset=="5plusUnqAns"]
# testn<-c("all" = tdat[, conserved], # n that are all only
#          "all&cis" = tdat[, cis], # n that are cis only
#          "all&trans" = tdat[, trans], # n that are trans only
#          "all&trans&cis&opposing" = tdat[, `cis-trans opposing`],
#          "all&trans&cis&enhancing" = tdat[, enhancing])
# set.seed(22)
# ftestn<-euler(testn, shape = "ellipse")
# plot(ftestn) # this is OKish; the overlap part still necessarily weird
# # what if we just take out enhancing...
# set.seed(22)
# testn2<-c("all" = tdat[, conserved], # n that are all only
#           "all&cis" = tdat[, cis], # n that are cis only
#           "all&trans" = tdat[, trans], # n that are trans only
#           "all&trans&cis&opposing" = tdat[, `cis-trans opposing`])
# ftestn2<-euler(testn2, shape = "ellipse")
# plot(ftestn2) # Better, lines up pretty nicely but not perfectly, enhancing is the 'leftover'. Possibly could give it its own color and legend somewhere?
# ftestn2.c<-euler(testn2)
# plot(ftestn2.c) # circle is WORSE: comes out over the sides
# 
# plot(ftestn2, fills = list(fill = c()),
#      legend = list(side = "right"))
# 
# plot(ftestn2, quantities = list(type = "counts") ) # label counts in segment
# 
# testn2.labs<-c(paste("all\n", tdat[, conserved + cis + trans + `cis-trans opposing` + enhancing]),
#                paste("cis\n", tdat[, cis +`cis-trans opposing` + enhancing ]),
#                paste("trans\n", tdat[, trans +`cis-trans opposing` + enhancing ]),
#                paste("compensatory\n", tdat[, `cis-trans opposing`])
#   
# ) # numbers of TOTALS
# plot(ftestn2, labels = testn2.labs) # label total n genes counts
# 
# ## Get all real numbers
# eulerns<-lapply(straininfo$strain, function(strn){
#   tdat<-regns[strain==strn & geneset=="5plusUnqAlns"]
#   testn<-c("all" = tdat[, conserved], # n that are all only
#            "all&cis" = tdat[, cis], # n that are cis only
#            "all&trans" = tdat[, trans], # n that are trans only
#            "all&trans&cis&compensatory" = tdat[, `cis-trans opposing`],
#            "all&trans&cis&enhancing" = tdat[, enhancing])
#   return(testn)
# })
# names(eulerns)<-straininfo$strain
# 
# ## Counts for possible labeling
# labcts<-lapply(straininfo$strain, function(strn){
#   tdat<-regns[strain==strn & geneset=="5plusUnqAlns"]
#   testn2.labs<-c(paste("all\n", tdat[, conserved + cis + trans + `cis-trans opposing` + enhancing]),
#                  paste("cis\n", tdat[, cis +`cis-trans opposing` + enhancing ]),
#                  paste("trans\n", tdat[, trans +`cis-trans opposing` + enhancing ]),
#                  paste("compensatory\n", tdat[, `cis-trans opposing`]),
#                  paste("enhancing\n", tdat[, enhancing])
#                  
#   )
#   return(testn2.labs)
# })
# names(labcts)<-straininfo$strain
# 
# # Seed engineering for best looking plots
# myseeds<-data.table(strain = straininfo$strain,
#                     myseed = c(22, # good
#                                4, # 4
#                                11, # good
#                                22, # good
#                                29, # good
#                                29, # good
#                                15 # good
#                                ))
# 
# ## Make & layout all these plots
# eulerfits<-lapply(names(eulerns), function(x){
#   set.seed(myseeds[strain==x, myseed])
#   return(euler(eulerns[[x]], shape = "ellipse"))
# })
# names(eulerfits)<-straininfo$strain
# 
# #     # No numbers
# eulerplts<-lapply(eulerfits, function(x){
#   plot(x, labels = list(font = c(1, 3, 3, 1, 1), fontsize = 8),
#        edges = list(col = c("black", rep("gray", 4))), # lty = c("solid", rep("dotted", 4))
#        fills = list(fill = c("white",
#                                              "lightgray",
#                                              "lightblue2",
#                                              "red4",
#                                              "orange2")))
# })
# plot_grid(plotlist = eulerplts, nrow = 3, ncol = 3, labels = names(eulerplts))
# 
# #     # With number each overlap labeled: I don't love this
# eulerplts.nlab<-lapply(eulerfits, function(x){
#   plot(x, labels = list(font = c(1, 3, 3, 1, 1), fontsize = 8),
#        quantities = list(type = "counts", 
#                          fontsize = 6), 
#        edges = list(col = c("black", rep("gray", 4))), # lty = c("solid", rep("dotted", 4))
#        fills = list(fill = c("white",
#                              "lightgray",
#                              "lightblue2",
#                              "red4",
#                              "orange2")))
# })
# plot_grid(plotlist = eulerplts.nlab, nrow = 3, ncol = 3, labels = names(eulerplts))
# 
# #     With total number in that category, regardless of overlap, labeled
# eulerplts.totlab<-lapply(names(labcts), function(x){
#   plot(eulerfits[[x]], labels =list(labels = labcts[[x]], font = c(1, 3, 3, 1, 1), fontsize = 8), 
#        edges = list(col = c("black", rep("gray", 4))), # lty = c("solid", rep("dotted", 4))
#        fills = list(fill = c("white",
#                              "lightgray",
#                              "lightblue2",
#                              "red4",
#                              "orange2")))
# })
# plot_grid(plotlist = eulerplts.totlab, nrow = 3, ncol = 3, labels = names(eulerfits)) # ** THIS IS BEST YET **
# 

#### Session info ####
cat("....compensationamuonttesting.R processing complete! Session information:....\n")
sessionInfo()