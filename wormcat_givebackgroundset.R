# (try to) Run WormCat using bespoke background sets (web and base tools don't do)
# by Avery Davis Bell, begun 2023.03.08

## How I initially installed wormcat R package:
# library("devtools")
# install_github("trinker/plotflow")
# install_github("dphiggs01/wormcat")
require(data.table)
require(argparser)
require(wormcat)

#### WormCat functions (sourced or copied and modified) ####
worm_cat_fun_giveannot<-function(file_to_process, title = "rgs", output_dir = NULL,
                                 rm_dir = FALSE, annotation_file = "whole_genome_v2_nov-11-2021.csv",
                                 input_type = "Wormbase.ID"){
  # exactly as worm_cat_fun except:
  #  takes direct path to annotation file (not relative to package etc) (won't work without full filepath)
  #   output_dir is absolute, not inside working directory
  mainDir <- getwd()
  if (is.null(output_dir)) {
    output_dir <- paste("worm-cat_", format(Sys.time(), 
                                            "%b-%d-%Y-%H:%M:%S"), sep = "")
  }
  output_dirPath <- file.path(output_dir)
  dir.create(output_dirPath)
  worm_cat_annotations <- file.path(annotation_file) 
  .worm_cat_add_categories(file_to_process, output_dirPath, 
                           worm_cat_annotations, input_type)
  .worm_cat_fisher_test(output_dirPath, worm_cat_annotations)
  
  for (i in 1:3) {
    cat_file_to_process <- sprintf("%s/rgs_fisher_cat%d.csv", 
                                   output_dir, i)
    print(sprintf("Processed %s", cat_file_to_process))
    .worm_cat_acceptable_pvalues(cat_file_to_process)
    cat_pvalue_file_to_process <- sprintf("%s/rgs_fisher_cat%d_apv.csv", 
                                          output_dir, i)
    plot_titles <- c(paste(title, "category1", sep = ":"), 
                     paste(title, "category2", sep = ":"), paste(title, 
                                                                 "category3", sep = ":"))
    .worm_cat_bubble_plot(cat_pvalue_file_to_process, plot_titles[i])
  }
  run_data <- sprintf("%s/run_data.txt", output_dir, i)
  runtime_l <- paste("runtime", format(Sys.time(), "%b-%d-%Y-%H:%M:%S"), 
                     sep = ":")
  annotation_file_l <- paste("annotation_version", annotation_file, 
                             sep = ":")
  input_type_l <- paste("input_type", input_type, sep = ":")
  cat(runtime_l, annotation_file_l, input_type_l, file = run_data, 
      sep = "\n", append = TRUE)
  files2zip <- dir(output_dirPath, full.names = TRUE)
  zip(zipfile = output_dir, files = files2zip)
  if (rm_dir == TRUE) {
    print("cleaning up")
    unlink(output_dir, TRUE)
  }
}


#### Other functions ####
makebackground<-function(backgroundlistf, gwannotf, outfile){
  # Narrows WormCat's genome-wide annotations to only include genes in backgroundlistf
  # In: backgroundlistf, path to Background gene list against which to test enrichment (e.g. expressed genes) (wormbase IDs). Column title should be 'Wormbase.ID'
  #     gwannotf, path to Genome-wide annotations from Wormcat, usually titled whole_genome_v2_nov-11-2021.csv
  #     outfile, where to write narrowed csv
  
  gwannot<-data.table(read.csv(gwannotf, header = T))
  bgs<-fread(backgroundlistf, header = T)
  
  setkey(gwannot, Wormbase.ID)
  setkey(bgs, Wormbase.ID)
  
  bgannot<-gwannot[bgs]
  write.csv(bgannot, outfile, row.names = F)
}

#### Arguments ####
p<-arg_parser("Run WormCat (http://wormcat.com/; Holdorf et al 2020 https://doi.org/10.1534/genetics.119.302919)
              using custom background geneset", 
              name = "wormcat_givebackgroundset.R", hide.opts = TRUE)

p<-add_argument(p, "--targetlist",
                help = "Gene list to test for enrichment (e.g. your upregulated genes) (wormbase IDs). Column title should be 'Wormbase.ID'",
                type = "character")
p<-add_argument(p, "--backgroundlist",
                help = "Background gene list against which to test enrichment (e.g. expressed genes) (wormbase IDs). Column title should be 'Wormbase.ID'",
                type = "character")
p<-add_argument(p, "--annots",
                help = "Genome-wide annotations from Wormcat, usually titled whole_genome_v2_nov-11-2021.csv",
                type = "character")
p<-add_argument(p, "--outdir",
                help = "Output directory path (relative or absolute); should describe this specific comparison if running multiple.",
                default = "wormcatout")
p<-add_argument(p, "--titleplot",
                help = "Title for Wormcat-generated plots",
                default = "")
p<-add_argument(p, "--wormcatcodedir",
                help = "Path to directory containing WormCat five helper function R scripts (from https://github.com/dphiggs01/Wormcat/tree/master/R)",
                default = "Wormcat/R")

p<-parse_args(p)

#### Run Wormcat ####
# Source wormcat functions directly - needed because I'm messing with worm_cat_fun; do after args parsed
source(file.path(p$wormcatcodedir, "worm_cat_acceptable_pvalues.R"))
source(file.path(p$wormcatcodedir, "worm_cat_add_categories.R"))
source(file.path(p$wormcatcodedir, "worm_cat_bubble_plot.R"))
source(file.path(p$wormcatcodedir, "worm_cat_controller.R"))
source(file.path(p$wormcatcodedir, "worm_cat_fisher_test.R"))

if(!dir.exists(p$outdir)){
  dir.create(p$outdir, recursive = T)
}

# Generate background geneset annotations
bgannotf<-file.path(p$outdir, paste0("backgroundannots_", gsub("\\s+", "", p$titleplot, fixed = F), 
                                     "_", sample.int(10e03, 1), ".csv")) # don't want spaces in file name - possible there could be some in plot title
makebackground(backgroundlistf = p$backgroundlist, gwannotf = p$annots,
               outfile = bgannotf)

# Run Wormcat
worm_cat_fun_giveannot(file_to_process = p$targetlist, title = p$titleplot, 
                       output_dir = p$outdir, annotation_file = bgannotf, 
                       input_type = "Wormbase.ID")


# Clean up: remove background geneset annotations (it's just a narrowing of the global file)
file.remove(bgannotf)