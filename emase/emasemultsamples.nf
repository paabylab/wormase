#!/usr/bin/env nextflow
/*
By Avery Davis Bell
Preps for EMASE-zero & runs EMASE-zero for multiple samples.
*/

/*
#### Set up input parameters & defaults ####
*/
// Inputs: input/output related & general
params.sampleinfo = "" // Path to tab-delimited file containing sample information. Column names (descriptions): SampleID (sample ID as in input filenames, to be used in output filenames); RefDescrip (description of reference genome(s) to use in output file names); Bwt2BasePath (base path - before .#.bt2 - for bowtie2 index files - can just put xxx if not running bowtie2); SalmonIndexDir (path to directory generated by salmon index - can just put xxx if not running salmon); PooledTranscriptLengths (path to file containing emase.pooled.transcripts.info-style transcript length information, e.g. <ref>_<strain>__transcriptlengths.txt for the strain combinations in this sample); fldMean (mean library fragment length per sample, passed to --fldMean in salmon quant), fldSD (standard deviation library fragment length per sample, passed to --fldSD in salmon quant)
params.outputdir = "" // Parent output directory. Will be created if doesn't exist.
params.refgenemappings = "" // Path to prepare-emase reference output file transcript to gene mapping (emase.gene2transcripts.tsv). Tab delimited file where the first field is the gene ID, all other fields are the transcript IDs that belong to that gene
params.fastqdir = "" // Directory containing all fastq.gz files to process. One or more per sample. **for paired end, this assumes fastqs are internal to sample-named directories here.
params.rnamap = 'salmon' // How to quantify RNA-seq for EMASE, bowtie2 mapping or salmon equivalence classes? Possible values: 'salmon' (CHOICE DEPRECATED)
params.libtype = 'SE' // SE or PE - is libary single end or paired end? So processing can be done appropriately.
params.plotinfo = "" // Path to tab-delimited file containing plotting-related sample information, passed to eqclassalnmtsummary_multsamples.R. Must include columns SampleID (sample ID as everywhere else), those passed to --groupby, --splitby, --facetby (or subsets of those used)

// Inputs: Tunable parameters for trimming
params.trimmodir = "/storage/coda1/p-apaaby3/0/shared/software/trimmomatic-0.39" // Path to trimmomatic v0.39 directory containing jar file and adapters directory (which itself contains the fasta files provided in next option, e.g. TruSeq3-SE.fa).
params.trimmofa = "TruSeq3-SE.fa" // name of fasta file within trimmodir adapters/ directory matching those used in library preparation
params.trimmoseedmism = 1 // Input to trimmomatic ILLUMINACLIP. How many of 16 bp can mismatch and still be counted as match.
params.trimmoadapclipthresh = 12 // Input to trimmomatic ILLUMINACLIP. How accurate match between adapter sequence and read must be. Each correct base adds 0.6. They recommend 7-15 (12 bases needed for 7, 25 for 15).

// Inputs: parameters for salmon
params.slibtype = "SR" // salmon --libtype option matching the library being aligned here

// Inputs: tunable parameters for emase-zero
params.emasemodel = "" // Normalization model for emase-zero (1-4)
params.emaset = 0.0001 // -t tolerance parameter for emase-zero
params.emasei = 999 // -i max # iterations parameter for emase-zero

// Inputs: parameters for eqclassalnmtsummary_multsamples.R
params.groupby = "Strain" // --groupby input of eqclassalnmtsummary_multsamples.R. Column name (in --plotinfo file) to group samples by in output plots - samples will be split by this on the same plot
params.splitby = "" // --splitby input of eqclassalnmtsummary_multsamples.R. Column name (in --plotinfo file) to split samples by, keeping all plots separate between these groups.
params.facetby = "" // --facetby input of eqclassalnmtsummary_multsamples.R. Column name (in --plotinfo file) to facet samples by, keeping them on separate plots but in the same PDF in these groups.
params.genelists = "" // --genelists input of eqclassalnmtsummary_multsamples.R. (# samples with unique alignments : # unique alignments per sample; can comma-separate to get more than 1 gene list per sample group. Sample grouping done based on by groupby, splitby, facetby inputs to this script)

// Inputs: organizational
params.alntoolsenv = '/storage/home/hcoda1/2/abell65/.conda/envs/alntools' // Path to alntools conda environment is installed in Anaconda3/2020.02 module
params.salmonenv = '/storage/home/hcoda1/2/abell65/.conda/envs/salmon' // path to conda environment where salmon is installed
params.asescriptsdir = '/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/rnaseq/ase' // path to directory containing salmonalleleeqclasses.py and eqclassalnmtsummary_multsamples.R scripts

// Housekeeping:  create output directories
// Alignment & related info
outdir = file(params.outputdir)
outdir.mkdirs()
outtrimdir = file(params.outputdir + "/triminfo")
outtrimdir.mkdirs()
outsalmdir = file(params.outputdir + "/salmonout")
outsalmdir.mkdirs()
// Quantification/EMASE. Will go in *sample-specific* directories inside here
outemasedir = file(params.outputdir + "/emase")
outemasedir.mkdirs()
// Characterization of haplotype-specific alignments/eq classes from salmon
outsalmchardir = file(params.outputdir + "/chareqclasses/persample")
outsalmchardir.mkdirs()
outsalmcharsummdir = file(params.outputdir + "/chareqclasses") // across samples directory; should have been created above

/*
#### Channel management ####
*/
//  Channel with all needed sample info as tuple.
formatChannel = Channel.fromPath(params.sampleinfo)
  .splitCsv(header: true, sep: '\t', strip: true)
  .map{row->
    return [row.SampleID, row.RefDescrip, row.Bwt2BasePath, row.SalmonIndexDir, row.PooledTranscriptLengths, row.fldMean, row.fldSD]
  }
  .set{sampleInfo}

/*
#### Processes ####
*/
if(params.libtype == 'SE'){
// do SE-specific processes
process mergeLaneFastqs{
  // Merge files across lanes so that there's one fastq per sample FOR SE DATA

  input:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(salmidx), val(trscrlns), val(meanfrag), val(sdfrag) from sampleInfo

  output:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(salmidx), val(trscrlns), val(meanfrag), val(sdfrag), path("mergedlanes.fastq") into mergedfastqs

  //when:
  //params.libtype == 'SE'

  """
  zcat ${params.fastqdir}/*${mysamp}*.fastq.gz > mergedlanes.fastq
  """
}

process trimmoIlluminaAdapters{
  // Use trimmomatic to trim Illumina adapters from merged fastqs **when they are SE**

  // save trimmomatic output summaries
  publishDir outtrimdir, mode: 'copy', pattern: "*_trimmomatic.out"

  input:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(salmidx), val(trscrlns), val(meanfrag), val(sdfrag), path(mergedfastq) from mergedfastqs

  output:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(trscrlns), path("trimmed.fastq") into trimmedfastqs // for bowtie2
  tuple val(mysamp), val(refdescrip), val(trscrlns), val(salmidx), val(meanfrag), val(sdfrag), path("trimmed.fastq") into trimmedfastqs2 // for salmon
  path("${mysamp}_trimmomatic.out") into trimlogs

  //when:
  //params.libtype == 'SE'

  """
  java -Xmx4g -jar ${params.trimmodir}/trimmomatic-0.39.jar \
  SE -threads 4 \
  ${mergedfastq} trimmed.fastq \
  ILLUMINACLIP:${params.trimmodir}/adapters/${params.trimmofa}:${params.trimmoseedmism}:30:${params.trimmoadapclipthresh} \
  2>${mysamp}_trimmomatic.out
  """
}

process salmonquantSE{
  // quantify RNA-seq data with salmon
  // Only runs when --rnamap is salmon or all, libtype is SE

  conda params.salmonenv

  // Save quantification outputs
  publishDir outsalmdir, mode: 'copy', overwrite: true

  input:
  tuple val(mysamp), val(refdescrip), val(trscrlns), val(salmidx), val(meanfrag), val(sdfrag), path(fastq) from trimmedfastqs2

  output:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path("${mysamp}_salmon_${refdescrip}") into salmouts
  tuple val(mysamp), val(refdescrip), val(trscrlns), path("${mysamp}_salmon_${refdescrip}") into salmouts2

  when:
  params.libtype == 'SE' && (params.rnamap == 'salmon' || params.rnamap == 'all')

  """
  salmon quant \
  -i ${salmidx} \
  -l ${params.slibtype} \
  -r ${fastq}
  -o ${mysamp}_salmon_${refdescrip} \
  --dumpEq \
  --fldMean ${meanfrag} \
  --fldSD ${sdfrag} \
  --rangeFactorizationBins 4 \
  --seqBias \
  --gcBias \
  -p 4
  """
}
}

if(params.libtype == 'PE'){
  // Do PE-specific processes

process trimmoIlluminaAdaptersPE{
  // Use trimmomatic to trim Illumina adapters from merged fastqs **when they are PE**
  // **currently PE assumed to be named _1 and _2.fastq.gz

  // save trimmomatic output summaries
  publishDir outtrimdir, mode: 'copy', pattern: "*_trimmomatic.out", overwrite: true
  // publishDir "${outtrimdir}/${mysamp}", mode: 'copy', pattern: "*unpaired.fastq.gz"

  input:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(salmidx), val(trscrlns), val(meanfrag), val(sdfrag) from sampleInfo

  output:
  tuple val(mysamp), val(refdescrip), val(bwtbase), val(trscrlns), path("trimmed_forward_paired.fastq"), path("trimmed_rev_paired.fastq") into trimmedfastqsPE into trimmedfastqsPE // for bowtie2
  tuple val(mysamp), val(refdescrip), val(trscrlns), val(salmidx), val(meanfrag), val(sdfrag), path("trimmed_forward_paired.fastq"), path("trimmed_rev_paired.fastq") into trimmedfastqs2PE // for salmon
  path("${mysamp}_trimmomatic.out") into trimlogsPE
  // tuple path("${mysamp}_trimmomatic.out"), path ("*unpaired.fastq.gz") into trimlogsPE

  //when:
  //params.libtype == 'PE'

  """
  java -Xmx4g -jar ${params.trimmodir}/trimmomatic-0.39.jar \
  PE -threads 4 \
  ${params.fastqdir}/${mysamp}/${mysamp}_1.fastq.gz ${params.fastqdir}/${mysamp}/${mysamp}_2.fastq.gz \
  trimmed_forward_paired.fastq trimmed_forward_unpaired.fastq \
  trimmed_rev_paired.fastq trimmed_rev_unpaired.fastq \
  ILLUMINACLIP:${params.trimmodir}/adapters/${params.trimmofa}:${params.trimmoseedmism}:30:${params.trimmoadapclipthresh}:2:True \
  2>${mysamp}_trimmomatic.out

  # Combine paired reads with unpaired for salmon. Trying it!
  # cat trimmed_forward_paired.fastq trimmed_forward_unpaired.fastq > trimmed_forward.fastq
  # cat trimmed_rev_paired.fastq trimmed_rev_unpaired.fastq > trimmed_forward.fastq

  # gzip unpaired read files for publishing - DEPRECATED
  # gzip *_unpaired.fastq
  """
}

process salmonquantPE{
  // quantify RNA-seq data with salmon
  // Only runs when --rnamap is salmon or all, libtype is PE

  conda params.salmonenv

  // Save quantification outputs
  publishDir outsalmdir, mode: 'copy', overwrite: true

  input:
  tuple val(mysamp), val(refdescrip), val(trscrlns), val(salmidx), val(meanfrag), val(sdfrag), path(fastq1), path(fastq2) from trimmedfastqs2PE

  output:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path("${mysamp}_salmon_${refdescrip}") into salmouts
  tuple val(mysamp), val(refdescrip), val(trscrlns), path("${mysamp}_salmon_${refdescrip}/aux_info/eq_classes.txt") into salmouts2

  //when:
  //params.libtype == 'PE' && (params.rnamap == 'salmon' || params.rnamap == 'all')

  """
  salmon quant \
  -i ${salmidx} \
  -l ${params.slibtype} \
  -1 ${fastq1} -2 ${fastq2} \
  -o ${mysamp}_salmon_${refdescrip} \
  --dumpEq \
  --fldMean ${meanfrag} \
  --fldSD ${sdfrag} \
  --rangeFactorizationBins 4 \
  --seqBias \
  --gcBias \
  -p 4

  # It seems that sometimes eq classes are gzipped, sometimes not? going to unzip here if so
  if [ -f ${mysamp}_salmon_${refdescrip}/aux_info/eq_classes.txt.gz ]; then
    echo "unzipping eq classes file"
    gunzip ${mysamp}_salmon_${refdescrip}/aux_info/eq_classes.txt.gz
  fi

  """
}
}

process alntoolssalm2ec{
  // Run alntools salmon2ec to get input for emase-zero
  // Only runs when --rnamap is salmon or all

  conda params.alntoolsenv

  input:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path(salmdir) from salmouts

  output:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path("*.bin") into ecsalms

  when:
  params.rnamap == 'salmon' || params.rnamap == 'all'

  """
  # doing unzip in prior script because I use this file a couple of times!
  # gunzip ${salmdir}/aux_info/eq_classes.txt.gz # alntools doesn't recognize this as zipped

  alntools salmon2ec \
  --verbose \
  ${salmdir} ${mysamp}_salm.bin
  """
}

process emasezerosalm{
  // Runs emase-zero on salmon output. emase-zero must be on path
  // only runs if --rnamap is salmon or all

  publishDir outemasedir, mode: 'copy', pattern: "*"

  input:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path(salmbin) from ecsalms

  output:
  path("*") into emasesalmout

  when:
  params.rnamap == 'salmon' || params.rnamap == 'all'

  """
  emase-zero \
  --model ${params.emasemodel} \
  -o ${mysamp}_aligned${refdescrip}_salmon_emasezero \
  -g ${params.refgenemappings} \
  -l ${trscrlns} \
  -t ${params.emaset} \
  -i ${params.emasei} \
  --verbose \
  ${salmbin}

  gzip ${mysamp}_aligned${refdescrip}_salmon_emasezero*
  """
}

process salmonalleleeqclasses{
  // Runs salmonalleleeqclasses.py to get eq, allele-specific information per transcript and gene
  // only runs if --rnamap is salmon or all

  publishDir outsalmchardir, mode: 'copy', pattern: "*"

  input:
  tuple val(mysamp), val(refdescrip), val(trscrlns), path(eqclasses) from salmouts2

  output:
  path("*") into eqinfo

  when:
  params.rnamap == 'salmon' || params.rnamap == 'all'

  """
  python ${params.asescriptsdir}/salmonalleleeqclasses.py \
  -g2t ${params.refgenemappings} \
  -eqs ${eqclasses} \
  -out ${mysamp}_aligned${refdescrip}_salmon
  """
}

process eqclassalnmtsummary{
  // Runs eqclassalnmtsummary_multsamples.R. Once all samples are through salmonalleleeqclasses
  // only runs if --rnamap is salmon or all

  publishDir outsalmcharsummdir, mode: 'copy', pattern: "*", overwrite: true

  input:
  path(alleqinfo) from eqinfo.collect()

  output:
  path("*") into eqinfomult

  when:
  params.rnamap == 'salmon' || params.rnamap == 'all'

  """
  module load r/4.2.1-tidy # not ideal to have in process but it isn't working through nf config

  Rscript ${params.asescriptsdir}/eqclassalnmtsummary_multsamples.R \
  --sampinfo ${params.plotinfo} \
  --exampinput _samp_*_genes.txt.gz \
  --outstem samplestogether \
  --groupby ${params.groupby} \
  --splitby ${params.splitby} \
  --facetby ${params.facetby} \
  --genelists ${params.genelists}
  """
}