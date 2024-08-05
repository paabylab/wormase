#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
By Avery Davis Bell
Runs Ornaments quant & optionally kallisto quant for C. elegans (or other crosses of homozygous parents where you have parental variant info)
Currently set up for paired end RNA-seq data

** first try using dsl 2/ workflow language...

*/

/*
#### Set up input parameters & defaults ####
*/

// Inputs: input/output related
params.strains = "" // Path to no-header list of strains to process - should be all strains any sample comes from. **isotypes/as in VCF
params.sampleinfo = "" // Path to tab-delimited file containing sample information. Column names (descriptions): SampleID (sample ID as in input filenames, to be used in output filenames); Strain (strain whose strain-specific ornaments/kallisto index should be used); RefDescrip (description of reference genome(s) to use in output file names)
params.out= "" // Parent output directory. Will be created if doesn't exist.
params.fastqdir = "" // Directory containing all fastq.gz files to process. One or more per sample. **for paired end, this assumes fastqs are internal to sample-named directories here.
params.vcf = "WI.20220216.hard-filter.isotype.vcf.gz" // path to VCF containing all samples of interest.
params.gtf = "c_elegans.PRJNA13758.WS283.canonical_geneset.gtf" // Path to GTF annotation file **same build/version/format/etc as used in VCF creation and transcritome FASTA
params.trfa = "c_elegans.PRJNA13758.WS283.mRNA_transcripts.fa" // Path to transcriptome fasta **same build/version/format/etc as used in VCF creation and GTF

// Inputs: run kallisto?
params.kallisto =  true // if true, kallisto run on strain-specific indexes (after creating them), if not, it isn't

// Inputs: do initial R script post ornaments
params.ornamentsrscript = true // if true, runs ornaments_initplots_data2gene.R. Only need to fill in the next options if true.
params.rscriptdir = "" // Directory containing ornaments_initplots_data2gene.R script
params.plotsampinfo = "" // path to plot script formatted sample info:  Path to sample information file containing
//			                      information on parental and F1 samples. Columns
//			                      should include: SampleID, Generation (Parental or
//			                      F1), Allele1 (NA for parental, reference strain
//			                      for F1), Allele2 (NA for parental, non-reference
//			                      strain for F1 - used as Strain for F1 in
//			                      generation/strain model), Strain (NA for F1,
//			                      strain for parental. All strains/alleles should
//                      be the way you want them to show up in outputs)
params.plotbaseoutname = "out" // Base name for all plot output files
params.plotstrains = "" // Path to strains ordered as desired in plot. **Strains, not isotypes, probably
params.tx2genef = "" // Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.
params.genegff = "" // Path to *genes only* gff3 file containing info on all gene_ids present in input counts file; includes gene location, name, biotype and other information.
params.plotexclchrs = "mtDNA" // Optional comma-separated, no space list of
//			                      chromosomes to exclude entirely from plots/downstream analyses (named as in
//			                      genegff). E.g. MtDNA cat be smart to exclude for
//			                      ASE analysis.
params.plotinclbiotype = "protein_coding" // Comma-separated list of biotypes (as in genegff) to include in plot processing

// Inputs: Tunable parameters for trimming
params.trimmodir = "" // Path to trimmomatic v0.39 directory containing jar file and adapters directory (which itself contains the fasta files provided in next option, e.g. TruSeq3-SE.fa).
params.trimmofa = "TruSeq3-PE-2.fa" // name of fasta file within trimmodir adapters/ directory matching those used in library preparation
params.trimmoseedmism = 1 // Input to trimmomatic ILLUMINACLIP. How many of 16 bp can mismatch and still be counted as match.
params.trimmoadapclipthresh = 12 // Input to trimmomatic ILLUMINACLIP. How accurate match between adapter sequence and read must be. Each correct base adds 0.6. They recommend 7-15 (12 bases needed for 7, 25 for 15).

// Inputs: organizational
params.orndir = "" // Path to ornaments directory where scripts convert_genomic_vcf_to_transcriptomic_vcf.py and reate_personalized_transcriptome.py are
params.ornconda = '' // path to conda environment set up to run ornaments & associated python scripts (as described in ornaments doc)

// Housekeeping:  create output directories
// General & trim
outdir = file(params.out)
outdir.mkdirs()
outtrimdir = file(params.out + "/triminfo")
outtrimdir.mkdirs()
// quantification
outorndir = file(params.out + "/ornaments")
outorndir.mkdirs()

if( params.kallisto ){
  outkalldir = file(params.out + "/kallisto")
  outkalldir.mkdirs()
  outkaldipdir = file(params.out + "/kallisto/diploid")
  outkaldipdir.mkdirs()
}
if( params.ornamentsrscript){
  outplotdir = file(params.out + "/ornamentsinitanalysis")
}

/*
#### Processes: per strain ####
*/

process makevcf{
  // Makes one-strain VCF where genotype is 0|1 for all variants with non-reference calls

  input:
  val(mystrain)

  output:
  tuple val(mystrain), path("$mystrain*.vcf.gz")

  """
  # Whole VCF
  bcftools view -s ${mystrain} ${params.vcf} | \
  bcftools view -i'GT="alt"' -Ov -o tmp.vcf -

  # Replace genotypes
  ## Get header for output
  bcftools view -h -s ${mystrain} ${params.vcf} | \
  tail -1 | \
  cut -f1-10 > header.tmp
  ## Create VCF
  bcftools view -H tmp.vcf | \
  cut -f1-9 | \
  awk -F : '{print \$1}' | \
  awk -F "\t" '{print \$0\"\\t\"\"0|1\"}' | \
  cat header.tmp - | gzip > ${mystrain}_01gts.vcf.gz # this NEEDS to be gzipped for ornaments!!

  # Remove unnecessary intermediate for space
  rm tmp.vcf
  """
}

process transcriptomicvcf{
  // Runs ornaments' convert_genomic_vcf_to_transcriptomic_vcf.py, which "extracts transcriptome variants located in exonic regions and transforms the genomic coordinates of variants to transcriptomic coordinates."

  conda params.ornconda

  input:
  tuple val(mystrain), path(strainvcf)

  output:
  tuple val(mystrain), path("*.vcf"), path("*.coords")

  """
  python ${params.orndir}/convert_genomic_vcf_to_transcriptomic_vcf.py \
  -v ${strainvcf} \
  -g ${params.gtf} \
  -t ${params.trfa} \
  -o ${mystrain}.transcriptome.vcf \
  -c ${mystrain}.transcriptome.coords
  """
}


process ornamentpersonalized{
  // Runs ornaments' create_personalized_transcriptome.py, which "modifies the reference transcriptome with the variant information that was prepared above, and produces an ornament personalized transcriptome."

  conda params.ornconda

  input:
  tuple val(mystrain), path(straintrvcf), path(straincoords)

  output:
  tuple val(mystrain), path(straintrvcf), path("*.fa")

  """
  python ${params.orndir}/create_personalized_transcriptome.py \
  -f ${params.trfa} \
  -v ${straintrvcf} \
  -t ornament \
  -o ${mystrain}.ornament.fa \
  -s ${mystrain}
  """
  // DON'T FORGET - POSSIBLY TEST HOW THIS WORKS IF TRY THE PHASED WAY!!
}

process ornamentsindex{
  // Builds ornaments index

  conda params.ornconda

  input:
  tuple val(mystrain), path(straintrvcf), path(ornfa)

  output:
  tuple val(mystrain), path(straintrvcf), path("*.index")

  """
  ornaments index -i ${mystrain}.ornament.index ${ornfa}
  """
}

process diptranscriptome{
  // Runs ornaments' create_personalized_transcriptome.py to get diploid transcriptome
  conda params.ornconda

  input:
  tuple val(mystrain), path(straintrvcf), path(straincoords)

  output:
  tuple val(mystrain), path("*.fa")

  """
  python ${params.orndir}/create_personalized_transcriptome.py \
  -f ${params.trfa} \
  -v ${straintrvcf} \
  -t diploid \
  -c ${straincoords} \
  -o ${mystrain}.diploid.fa \
  -s ${mystrain}
  """
}

process straintranscriptome{
 // Extracts the _R transcripts from diploid transcriptome to get strain-specific transcriptome

 input:
 tuple val(mystrain), path(dipfa)

 output:
 tuple val(mystrain), path("*_strainspec.fa")

 """
 seqkit grep -r -p "_R" ${dipfa} | \
 sed 's/_R//g' - > ${mystrain}_strainspec.fa
 """
}

process kallistoindex{
  // Runs kallisto index on personalized transcriptome from diploid (not ornaments one!)

  input:
  tuple val(mystrain), path(strainfa)

  output:
  tuple val(mystrain), path("*kallisto.index")

  """
  kallisto index \
  -i ${mystrain}.kallisto.index \
  -t 4 \
  ${strainfa}
  """
}

process kallistodipindex{
  // Runs kallisto index on DIPLOID transcriptome (to align F1s to)

  input: tuple val(mystrain), path(dipfa)

  output:
  tuple val(mystrain), path("*kallisto.index")

  """
  kallisto index \
  -i ${mystrain}.diploid.kallisto.index \
  -t 4 \
  ${dipfa}
  """
}

/*
#### Processes: per sample ####
*/

process trimmoIlluminaAdaptersPE{
  // Use trimmomatic to trim Illumina adapters from merged fastqs **when they are PE**
    // **currently PE assumed to be named _1 and _2.fastq.gz

  // save trimmomatic output summaries
  publishDir outtrimdir, mode: 'copy', pattern: "*_trimmomatic.out"
  // publishDir outtrimdir, mode: 'copy', pattern: "*unpaired.fastq.gz"

  input:
  tuple val(mysamp), val(sampstrain), val(refdescrip)

  output: // ***sampstrain needs to be FIRST to match with other channel
  tuple val(sampstrain), val(mysamp), val(refdescrip), path("trimmed_forward_paired.fastq"), path("trimmed_rev_paired.fastq"), emit: trimmedfastqs
  path("${mysamp}_trimmomatic.out"), emit: trimlogsPE
  // tuple path("${mysamp}_trimmomatic.out"), path("_unpaired.fastq.gz")

  """
  java -Xmx4g -jar ${params.trimmodir}/trimmomatic-0.39.jar \
  PE -threads 4 \
  ${params.fastqdir}/${mysamp}/${mysamp}_1.fastq.gz ${params.fastqdir}/${mysamp}/${mysamp}_2.fastq.gz \
  trimmed_forward_paired.fastq trimmed_forward_unpaired.fastq \
  trimmed_rev_paired.fastq trimmed_rev_unpaired.fastq \
  ILLUMINACLIP:${params.trimmodir}/adapters/${params.trimmofa}:${params.trimmoseedmism}:30:${params.trimmoadapclipthresh}:2:True \
  2>${mysamp}_trimmomatic.out

  # gzip unpaired read files for publishing
  # gzip *_unpaired.fastq
  """
}

process ornamentsquant{
  // Run ornaments quant on individual samples

  conda params.ornconda
  publishDir outorndir, mode: 'copy', overwrite: true

  input:
  tuple val(mystrain), path(straintrvcf), path(strainornind), val(mysamp), val(refdescrip), path(fastq1), path(fastq2)

  output:
  path("${mysamp}_ornaments_${refdescrip}")

  """
  ornaments quant \
  -i ${strainornind} \
  -o ${mysamp}_ornaments_${refdescrip} \
  --vcf ${straintrvcf} \
  --sample ${mystrain} \
  ${fastq1} ${fastq2}

  # new: gzip here!
  gzip ${mysamp}_ornaments_${refdescrip}/*.txt
  """
}

process kallistoquant{
  // Runs kallisto quant. **currently sets library as --rf-stranded, need to CHECK THIS!

  publishDir outkalldir, mode: 'copy', overwrite: true

  input:
  tuple val(mystrain), path(kalindex), val(mysamp), val(refdescrip), path(fastq1), path(fastq2)

  output:
  tuple val(mystrain), val(mysamp), path("*")
  """
  kallisto quant \
  -i ${kalindex} \
  -o ${mysamp}_kallistostrainspec_${mystrain} \
  -t 4 \
  --rf-stranded \
  ${fastq1} ${fastq2}

  # new: gzip here!
  gzip ${mysamp}_kallistostrainspec_${mystrain}/abundance.tsv
  """
}

process kallistodipquant{
  // Runs kallisto quant to diploid index. ** really only makes sense for F1 samples, but currently does for all for ease **

  publishDir outkaldipdir, mode: 'copy', overwrite: true

  input:
  tuple val(mystrain), path(kaldipindex), val(mysamp), val(refdescrip), path(fastq1), path(fastq2)

  output:
  tuple val(mystrain), val(mysamp), path("*")

  """
  kallisto quant \
  -i ${kaldipindex} \
  -o ${mysamp}_kallistodip_${refdescrip} \
  -t 4 \
  --rf-stranded \
  ${fastq1} ${fastq2}

  # new: gzip here!
  gzip ${mysamp}_kallistodip_${refdescrip}/abundance.tsv
  """
}

/*
#### Processes: combined samples ####
*/

process ornamentsinitplots{
  // Runs ornaments_initplots_data2gene.R

  publishDir outplotdir, mode: 'copy', overwrite: true

  input:
  path(ornoutdirs)

  output:
  path("*")

  """
  module load r/4.2.1-tidy # not ideal to have in process but it isn't working through nf config

  Rscript ${params.rscriptdir}/ornaments_initplots_data2gene.R \
  --sampinfo ${params.plotsampinfo} \
  --baseoutname ${params.plotbaseoutname} \
  --outdir "getwd()" \
  --allelecounts "*_SAMPID__ornaments_*/allele_counts.txt*" \
  --strains ${params.plotstrains} \
  --tx2genef ${params.tx2genef} \
  --genegff ${params.genegff} \
  --exclchrs ${params.plotexclchrs} \
  --inclbiotype ${params.plotinclbiotype}
  """
}

/*
#### Workflow ####
*/
workflow{
  // Set up channels
  // STRAIN INFO
  formatChannel = Channel.fromPath(params.strains)
    .splitCsv(header: false, sep: '\t', strip: true)
    .map{row->
      return row[0]
    }
    .set{strainInfo}
  // SAMPLE INFO
  formatChannel = Channel.fromPath(params.sampleinfo)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row->
      return [row.SampleID, row.Strain, row.RefDescrip]
    }
    .set{sampleInfo}

  // Do strain-level processes: ornaments.
  makevcf(strainInfo) | transcriptomicvcf | ornamentpersonalized | ornamentsindex

  // Do sample-level processes: ornaments
  trimmoIlluminaAdaptersPE(sampleInfo)
  // ** pre-ornaments, have to combine channels. NB collecting these before combining did NOT work
  ornamentsindex.out | combine(trimmoIlluminaAdaptersPE.out.trimmedfastqs, by: 0) | set{ forquant }
  ornamentsquant(forquant)

  // do kallisto processes if desired
  if ( params.kallisto ) {
    // strain-level processes
    diptranscriptome(transcriptomicvcf.out) | straintranscriptome | kallistoindex
    kallistodipindex(diptranscriptome.out)

    // Sample process: kallisto quant
    kallistoindex.out | combine(trimmoIlluminaAdaptersPE.out.trimmedfastqs, by: 0) | set{ forkalquant }
    kallistoquant(forkalquant)
    kallistodipindex.out | combine(trimmoIlluminaAdaptersPE.out.trimmedfastqs, by: 0) | set{ forkaldipquant }
    kallistodipquant(forkaldipquant)
  }

  // Pull ornaments outputs together, get per gene & plot (if desired)
  if ( params.ornamentsrscript ){
    ornamentsquant.out | collect | set{ allornfiles }
    ornamentsinitplots(allornfiles)
  }
}
