#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 13:41:20 2021
@author: Avery Davis Bell

Process salmon eq_classes.txt file from samon quant to a diploid (or pseudo-diploid) transcriptome
to summarize eq classes per transcript and per gene in an allele-specific manner.
Useful for identifying transcripts with allele-specific eq classes.
"""
import argparse
import gzip

#%%
def readgenes2transcriptstsv(strGns2Trs):
    """
    Processes an emase-format genes2transcripts file into dictionaries keyed by gene and transcript ID
    Input file format: each line has one gene; gene ID then all transcripts from this gene; tab-delimited

    Parameters
    ----------
    strGns2Trs : string
        Path to emase-format genes2transcripts file
        
    Returns
    -------
    hTrGeneMap : dict
        Keys are transcripts, values are gene this transcript comes from
    hGeneTrMap : dict
        Keys are gene IDs, values are gene lists of all transcripts from this gene

    """
    hTrGeneMap = {}
    hGeneTrMap = {}
    
    # Open input filestream
    istm = open(strGns2Trs, "r")
    
    # Add each transcript, gene to reciprocal dictionaries
    for strLine in istm:
        astrLine = strLine.strip().split("\t")
        hGeneTrMap[astrLine[0]] = astrLine[1:]
        for strTr in astrLine[1:]:
            hTrGeneMap[strTr] = astrLine[0]
    
    # Close input filestream
    istm.close()
    
    return([hTrGeneMap, hGeneTrMap])
#%%
def proctrline(iID, strTr, hNums, hTrs):
    """
    

    Parameters
    ----------
    iID : int
        Number of transcript being processed, to associate with transcript ID in output dicts
    strTr : string
        Transcript ID from eq class with _ separated haplotype ID, e.g. <transcript ID>_<haplotype>
    hNums : dict
        Keys are integer numbers of transcripts, as specified in eq classes. 
        Values are list: transcript this refers to, haplotype of transcript this refers to
        New iID key will be added
    hTrs : dict
        Keys are transcript IDs (no haplotype information)
        Values are list: integer number of transcript in first haplotype, integer number of transcript in second haplotype (first & second in order they were in input)
        Either new key and first haplotype number will be added, or second haplotype number will be added if transcript is already in dict

    Returns
    -------
    None - modifies in place - see inputs

    """
    
    astrTrHap = strTr.split("_")
    hNums[iID] = astrTrHap
    if astrTrHap[0] in hTrs.keys():
        hTrs[astrTrHap[0]].append(iID)
    else:
        hTrs[astrTrHap[0]] = [iID]


#%%
def readeqclasses(strEqs):
    """
    Reads in whole eq class file - transcript information as well as eq class information

    Parameters
    ----------
    strEqs : string
        Path to salmon eq classes file to read. If ends in .gz, assumes gziped and will try to use gzip to open

    Returns
    -------
    hNums : dict
        Keys are integer numbers of transcripts, as specified in eq classes. 
        Values are list: transcript this refers to, haplotype of transcript this refers to
        In order, this will match information from ambig_info.tsv!!
    hTrs : dict
        Keys are transcript IDs (no haplotype information)
        Values are list: integer number of transcript in first haplotype, integer number of transcript in second haplotype (first & second in order they were in input)
        Values are ['NA', 'NA'] in cases where only one haplotype's transcript version was present, if any
    aiHapOne: list
        List of all the transcript integer numbers from haplotype 1. NAs included in cases where only one haplotype's transcript version was present, if any
    aiHapTwo: list
        List of all the transcript integer numbers from haplotype 2. NAs included in cases where only one haplotype's transcript version was present, if any
    aaEqs : list (of lists)
        List of eq classes. Each is a list of all the integers found in that eq class.
        (first = number of transcripts in class; last = number of fragments aligning to this class; interim = transcript number IDs)
    """
    
    # Initialize
    hNums = {}
    hTrs = {}
    aaEqs = []
    
    # Open input filestream
    if strEqs.endswith(".gz"):
        istm = gzip.open(strEqs, "rb")
    else:
        istm = open(strEqs, "r")
    
    # Process file
    i = 0
    for strLine in istm:
        i += 1
        if i==1: # first line contains # transcripts
            iNTrs = int(strLine.strip())
        elif i > 2:
            if i <= (2 + iNTrs): # Process transcripts
                iID = (i - 3) # transcripts are ZERO INDEXED in eq classes file (despite no documentation to this effect)
                strTr = strLine.strip()
                proctrline(iID, strTr, hNums, hTrs) 
            else: # Process eq classes!
                astrLine = strLine.strip().split()
                aiLine = []
                for strNum in astrLine:
                    aiLine.append(int(strNum))
                aaEqs.append(aiLine)
    
    # Close input filestream
    istm.close()
    
    # Get lists with all the hap1, all the hap2 transcripts
    aiHapOne, aiHapTwo = [], []
    for strTr in hTrs.keys():
        if len(hTrs[strTr]) == 2:
            aiHapOne.append(hTrs[strTr][0])
            aiHapTwo.append(hTrs[strTr][1])
        else: # Transcript only present for one haplotype in Salmon - print record of these; propagate NAs whenever see these
            print("Transcript only present for one haplotype: " + strTr + " - will be NAed in transcript-level outputs.")
            # Update ALL records of this transcript to have NAs (keep IDs in hNums, which is also used to match with other data)
            aiHapOne.append("NA")
            aiHapTwo.append("NA")
            hTrs[strTr] = ["NA", "NA"]

    # Return
    return([hNums, hTrs, aiHapOne, aiHapTwo, aaEqs])

#%%
def proceqclass(aiEq, hNums, hTrs, aiHapOne, aiHapTwo, hTrGeneMap, hEqsPerTr):
    """
    Parses one eqclass and assigns it to each transcript in the eq class (only 1x if both haplotypes' transcript versions are included).
    Keeps track of whether both or only one haplotype represented by class; whether all transcripts in class come from the same gene

    Parameters
    ----------
    aiEq : list of ints 
        Data for one eq class. 
        First element = number of transcripts in class; last element = number of fragments aligning to this class; interim = transcript number IDs
    hNums : dict
        Output from readeqclasses()    
        Keys are integer numbers of transcripts, as specified in eq classes. 
        Values are list: transcript this refers to, haplotype of transcript this refers to
        In order, this will match information from ambig_info.tsv!!
    hTrs : dict
        Output from readeqclasses()  
        Keys are transcript IDs (no haplotype information)
        Values are list: integer number of transcript in first haplotype, integer number of transcript in second haplotype (first & second in order they were in input)
        Values are ['NA', 'NA'] in cases where only one haplotype's transcript version was present, if any
    aiHapOne : list
        Output from readeqclasses()  
        List of all the transcript integer numbers from haplotype 1. NAs included in cases where only one haplotype's transcript version was present, if any
    aiHapTwo : list
        Output from readeqclasses()  
        List of all the transcript integer numbers from haplotype 2. NAs included in cases where only one haplotype's transcript version was present, if any
    hTrGeneMap : dict
        Output from readgenes2transcriptstsv()
        Keys are transcripts, values are gene this transcript comes from
    hEqsPerTr : dict
        Dictionary with keys being transcript IDs (name, no haplotype)
        Values are lists of lists; one list per eq class
        One such list added for this eq class (see Returns) for each transcript

    Returns
    -------
    None. Modifies hEqsPerTr in place: for each transcript in eq class, appends a list describing eq class:
        [# alignments in class
        total # transcripts in class, 
        # haplotype 1 transcripts in class,
        # haplotype 2 transcripts in class,
        T/F has only one haplotype's transcripts in class,
        T/F for THIS transcript (could be different for each one in eq class), are both haplotypes for same transcript in the eq class, 
        T/F are all transcripts in this eq class from the same GENE (regardless of haplotypes)]
        str, transcript numbers in this eq class concatenated together with _ : unique haplotype class identifier (for later processing across gene)
        *** Only appends 1x even if both haplotypes' transcripts are in eq class (counting per transcript, not per haplotype)
    """
    iNts = aiEq[0] # total # transcripts in class
    iNAls = aiEq[-1] # # alignments in class]
    aiTrs = aiEq[1:len(aiEq)-1] # transcript IDs in class
    astrTrs = [str(i) for i in aiTrs] # string of transcript IDs - for output identifier
    
    # Get descriptors of this class
    iN1 = 0 # # haplotype 1 transcripts in class
    iN2 = 0 # # haplotype 2 transcripts in class
    astrGns = [] # Genes transcripts in class are from
    for iTr in aiTrs:
        if iTr in aiHapOne:
            iN1 += 1
        elif iTr in aiHapTwo:
            iN2 += 1
        # Keep track of genes transcripts are from
        strTr = hNums[iTr][0] # transcript name
        strGn = hTrGeneMap[strTr]
        if strGn not in astrGns:
            astrGns.append(strGn)
    
    bOneHap = True if iN1 == 0 or iN2 == 0 else False # T/F has only one haplotype's transcripts in class,
 
    bGene = True if len(astrGns) == 1 else False #  T/F are all transcripts in this eq class from the same GENE (regardless of haplotypes)
    
    # Add list describing class to each transcript (UNLESS it doesn't exist for both haplotypes)
    # ONLY add 1x, even if transcripts are in it from both haplotypes
    astrAddedTo = [] # keep track of which transcripts this is added to so don't add twice (i.e. for both haplotypes)
    for iTr in aiTrs:
        strTr = hNums[iTr][0] # transcript name
        if hTrs[strTr][0] != 'NA': # Check it's in both haplotypes (no NA); proceed if so
            # Figure out T/F for THIS transcript, are both haplotypes for same transcript in the eq class
            iThis = 0
            for iThisTr in hTrs[strTr]:
                if iThisTr in aiTrs:
                    iThis += 1
            bSame = True if iThis == 2 else False
            
            # Add this eq class to transcript's info IF not already added from other haplotype
            if strTr not in astrAddedTo:
                aEqInfo = [iNAls, iNts, iN1, iN2, bOneHap, bSame, bGene, "_".join(astrTrs)]
                if strTr in hEqsPerTr.keys():
                    hEqsPerTr[strTr].append(aEqInfo)
                else:
                    hEqsPerTr[strTr] = [aEqInfo]
                if bSame: # update to say this tr's been included
                    astrAddedTo.append(strTr)

#%%
def proceqclasses(hNums, hTrs, aiHapOne, aiHapTwo, aaEqs, hTrGeneMap, hGeneTrMap):
    """
    Summarizes eq classes per-transcript and per-gene, keeping track of critical information of haplotype specificity & uniqueness

    Parameters
    ----------
    hNums : dict
        Output from readeqclasses()    
        Keys are integer numbers of transcripts, as specified in eq classes. 
        Values are list: transcript this refers to, haplotype of transcript this refers to
        In order, this will match information from ambig_info.tsv!!
     hTrs : dict
        Output from readeqclasses()  
        Keys are transcript IDs (no haplotype information)
        Values are list: integer number of transcript in first haplotype, integer number of transcript in second haplotype (first & second in order they were in input)
        Values are ['NA', 'NA'] in cases where only one haplotype's transcript version was present, if any
    aiHapOne : list
        Output from readeqclasses()  
        List of all the transcript integer numbers from haplotype 1. NAs included in cases where only one haplotype's transcript version was present, if any
    aiHapTwo : list
        Output from readeqclasses()  
        List of all the transcript integer numbers from haplotype 2. NAs included in cases where only one haplotype's transcript version was present, if any
    aaEqs : list (of lists)
         Output from readeqclasses()  
        List of eq classes. Each is a list of all the integers found in that eq class.
        (first = number of transcripts in class; last = number of fragments aligning to this class; interim = transcript number IDs)
    hTrGeneMap : dict
        Output from readgenes2transcriptstsv()
        Keys are transcripts, values are gene this transcript comes from
    hGeneTrMap : dict
        Output from readgenes2transcriptstsv()
        Keys are gene IDs, values are gene lists of all transcripts from this gene

    Returns
    -------
    hTrEqSumm : dict 
        Keys are transcript IDs
        Values are lists of 12 - the following integers:
            Number of eq classes this transcript is seen in (total)
            Number of haplotype-specific eq classes (the only members are from the same haplotype)
            Number of gene-specific eq classes (the only members are from the same gene transcript is assigned to)
            Number of eq classes containing both haplotypes of this transcript
            Number of eq classes specific to haplotype and gene 
            Number of eq classes specific to haplotype, gene, and transcript: as unique as it gets (these counted as unique in the ambig.tsv salmon output file)
            Number of alignments to eq classes containing this transcript
            Number of alignments to haplotype-specific eq classes containing this transcript
            Number of alignments to gene-specific eq classes containing this transcript
            Number of alignments to eq classes containing both haplotypes of this transcript
            Number of alignments to eq classes specific to haplotype and gene containing this transcript (key for gene-level ASE!)
            Number of alignments to eq classes specific to haplotype, gene, and transcript (key for transcript-level ASE!)
            
    hGeneEqSumm : dict
        Keys are gene IDs
        Values are lists of 8 - the following integers:
            Number of unique eq classes (transcripts of) this gene present in
            Number of unique haplotype-specific eq classes (transcripts of) this gene present in
            Number of unique gene-specific eq classes (transcripts of) this gene present in
            Number of unique haplotype and gene-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique eq classes (transcripts of) this gene present in [total alignments]
            Number of alignments to unique haplotype-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique gene-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique haplotype and gene-specific eq classes (transcripts of) this gene present in (key for gene-level ASE!)

    """
    # Process all eq classes, assigning to transcript. Takes a bit of time!
    hEqsPerTr = {}
    for aiEq in aaEqs:
        proceqclass(aiEq, hNums, hTrs, aiHapOne, aiHapTwo, hTrGeneMap, hEqsPerTr)
    
    # Summarize eq classes per transcript
    hTrEqSumm = {}
    for strTr in hTrs.keys():
        # Populate output with NAs if only one haplotype was seen. 
        if hTrs[strTr][0]=='NA':  # Populate output with NAs if only one haplotype was seen. 
            hTrEqSumm[strTr] = ['NA'] * 12
        elif strTr not in hEqsPerTr.keys():  # If no eq classes, no reads - populate with 0s
            hTrEqSumm[strTr] = [0] * 12
        else:
            # Initialize counters
            iN = len(hEqsPerTr[strTr]) # Number of eq classes this transcript is seen in
            iHapSpec = 0 # Number of haplotype-specific eq classes X
            iGeneSpec = 0 # Number of gene-specific eq classes X
            iBothHapsTr = 0 # Number of eq classes with both haplotypes of this transcript X
            iHapGeneSpec = 0 # Number of eq classes specific to haplotype and gene X
            iHapGeneTrSpec = 0 # Number of eq classes specific to haplotype, gene, and transcript: as unique as it gets (these counted as unique in the ambig.tsv salmon output file) X
        
            ## Alignments in categories as above
            iAN, iAHapSpec, iAGeneSpec, iABothHapsTr, iAHapGeneSpec, iAHapGeneTrSpec = [0] * 6
        
            # Process each transcript, update counters
            for aEqInfo in hEqsPerTr[strTr]:
                iAN += aEqInfo[0]
                if aEqInfo[4]: # T/F has only one haplotype's transcripts in class
                    iHapSpec += 1
                    iAHapSpec += aEqInfo[0]
                    if aEqInfo[6]: # T/F are all transcripts in this eq class from the same GENE 
                        iHapGeneSpec += 1
                        iAHapGeneSpec += aEqInfo[0]
                        if aEqInfo[1] == 1: # only 1 transcript_haplotype in class: must be fully specific
                            iHapGeneTrSpec += 1
                            iAHapGeneTrSpec += aEqInfo[0]
                if aEqInfo[6]: # T/F are all transcripts in this eq class from the same GENE (regardless of haplotypes)
                    iGeneSpec += 1
                    iAGeneSpec += aEqInfo[0]
                if aEqInfo[5]: # T/F for THIS transcript, are both haplotypes for same transcript in the eq class
                    iBothHapsTr += 1
                    iABothHapsTr += aEqInfo[0]
            
            # Add to output
            hTrEqSumm[strTr] = [iN, iHapSpec, iGeneSpec, iBothHapsTr, iHapGeneSpec, iHapGeneTrSpec, iAN, iAHapSpec, iAGeneSpec, iABothHapsTr, iAHapGeneSpec, iAHapGeneTrSpec]
          
    
    # Summarize eq classes per gene
    hGeneEqSumm = {}
    for strGene in hGeneTrMap.keys():
        # Initialize eq class counters. Named as per transcript (keeping only relevant ones), but here identical eq classes will be collapsed across transcripts
        iN, iHapSpec, iGeneSpec, iHapGeneSpec = [0] * 4
        
        # Initialize alignment counters. Named as per transcript (keeping only relevant ones), but here identical eq classes will be collapsed across transcripts
        iAN, iAHapSpec, iAGeneSpec, iAHapGeneSpec = [0] * 4
        
        # Get all unique eq classes that map to this gene (don't double-count across transcripts)
        aaGnEqs = []
        for strTr in hGeneTrMap[strGene]:
            if strTr not in hEqsPerTr.keys(): # Happens either if NA'ed out OR if no reads. Not keeping track of which at gene level.
                continue
            for aEq in hEqsPerTr[strTr]:
                if aEq not in aaGnEqs:
                    aaGnEqs.append(aEq)
        
        # Process eq classes, updating counters
        for aEqInfo in aaGnEqs:
            iN +=1
            iAN += aEqInfo[0]
            if aEqInfo[4]: # T/F has only one haplotype's transcripts in class
                iHapSpec += 1
                iAHapSpec += aEqInfo[0]
                if aEqInfo[6]: # T/F are all transcripts in this eq class from the same GENE 
                    iHapGeneSpec += 1
                    iAHapGeneSpec += aEqInfo[0]
            if aEqInfo[6]: # T/F are all transcripts in this eq class from the same GENE (regardless of haplotypes)
                iGeneSpec += 1
                iAGeneSpec += aEqInfo[0]
        
        # Add to output
        hGeneEqSumm[strGene] = [iN, iHapSpec, iGeneSpec, iHapGeneSpec, iAN, iAHapSpec, iAGeneSpec, iAHapGeneSpec] 
        
    # Return
    return([hTrEqSumm, hGeneEqSumm])

#%%
def writeout(hTrEqSumm, hGeneEqSumm, strOutBase):
    """
    Writes transcripts & genes output files. Columns transcript_id or gene_id, then the numbers from the dicts' value lists

    Parameters
    ----------
    hTrEqSumm : dict 
        Ouptut of proceqclasses()
        Keys are transcript IDs
        Values are lists of 12 - the following integers:
            Number of eq classes this transcript is seen in (total)
            Number of haplotype-specific eq classes (the only members are from the same haplotype)
            Number of gene-specific eq classes (the only members are from the same gene transcript is assigned to)
            Number of eq classes containing both haplotypes of this transcript
            Number of eq classes specific to haplotype and gene 
            Number of eq classes specific to haplotype, gene, and transcript: as unique as it gets (these counted as unique in the ambig.tsv salmon output file)
            Number of alignments to eq classes containing this transcript
            Number of alignments to haplotype-specific eq classes containing this transcript
            Number of alignments to gene-specific eq classes containing this transcript
            Number of alignments to eq classes containing both haplotypes of this transcript
            Number of alignments to eq classes specific to haplotype and gene containing this transcript (key for gene-level ASE!)
            Number of alignments to eq classes specific to haplotype, gene, and transcript (key for transcript-level ASE!)  
    hGeneEqSumm : dict
        Output of proceqclasses()
        Keys are gene IDs
        Values are lists of 8 - the following integers:
            Number of unique eq classes (transcripts of) this gene present in
            Number of unique haplotype-specific eq classes (transcripts of) this gene present in
            Number of unique gene-specific eq classes (transcripts of) this gene present in
            Number of unique haplotype and gene-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique eq classes (transcripts of) this gene present in [total alignments]
            Number of alignments to unique haplotype-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique gene-specific eq classes (transcripts of) this gene present in
            Number of alignments to unique haplotype and gene-specific eq classes (transcripts of) this gene present in (key for gene-level ASE!)
    strOutBase : str
        Output file base. Two files will be written: <this>_eqclasses_transcripts.txt.gz
        and <this>_eqclasses_genes.txt.gz
        
    Returns
    -------
    None.

    """
    # Write out transcripts
    ostmT = gzip.open(strOutBase + "_eqclasses_transcripts.txt.gz", "wb")
    # Header
    strTHead = "\t".join(["transcript_id", "neq_all", "neq_hapspec", "neq_genespec",
                          "neq_twohaps", "neq_hapgenespec", "neq_haptranscriptspec",
                          "naln_all", "naln_hapspec", "naln_genespec",
                          "naln_twohaps", "naln_hapgenespec", "naln_haptranscriptspec"]) + "\n"
    ostmT.write(strTHead.encode("utf-8"))
    # Each one
    for strTr in sorted(hTrEqSumm.keys()):
        strWrite = strTr + "\t" + "\t".join([str(i) for i in hTrEqSumm[strTr]]) + "\n"
        ostmT.write(strWrite.encode("utf-8"))
    ostmT.close()
    
    # Write out genes
    ostmG = gzip.open(strOutBase + "_eqclasses_genes.txt.gz", "wb")
    # Header
    strGHead = "\t".join(["gene_id",
                          "neq_all", "neq_hapspec", "neq_genespec", "neq_hapgenespec",
                          "naln_all", "naln_hapspec", "naln_genespec", "naln_hapgenespec",]) + "\n"
    ostmG.write(strGHead.encode("utf-8"))
    # Each one
    for strG in sorted(hGeneEqSumm.keys()):
        strWrite = strG + "\t" + "\t".join([str(i) for i in hGeneEqSumm[strG]]) + "\n"
        ostmG.write(strWrite.encode("utf-8"))
    ostmG.close()  

#%%
def main():
    # Parse arguments
    argp = argparse.ArgumentParser(prog = "salmonalleleeqclasses.py",
                                   description = "Processes salmon eq_classes.txt file from salmon quant to a diploid (or pseudo-diploid) transcriptome \
                                       to summarize eq classes per transcript and per gene in an allele-specific manner \
                                    Useful for identifing transcripts with allele-specific eq classes.")
    argp.add_argument("-g2t", dest = "strGns2Trs", required = True, type = str,
                      metavar = "emase-format genes2transcripts mapping file",
                      help = "Path to emase-format genes2transcripts mapping file for reference genome. \
                          Each line has one gene; gene ID then all transcripts from this gene; tab-delimited")
    argp.add_argument("-eqs", dest = "strEqs", required = True, type = str, 
                      metavar = "salmon eq_classes.txt file",
                      help = "Path to salmon eq classes file to read. If ends in .gz, assumes gziped and will try to use gzip to open\
                          careful, gzip seems to be behaving weirdly with these specific files.")
    argp.add_argument("-out", dest = "strOutBase",
                      metavar = "output filestem/base", 
                      help = "Output file base. Two files will be written: <this>_eqclasses_transcripts.txt.gz \
                          and <this>_eqclasses_genes.txt.gz")
    args = argp.parse_args()
    
    # Run
    print("...Reading in gene-transcript mappings...")
    hTrGeneMap, hGeneTrMap = readgenes2transcriptstsv(args.strGns2Trs)
    print("...Gene-transcript mappings read in...")
    print("...Reading in eq classes...")
    hNums, hTrs, aiHapOne, aiHapTwo, aaEqs = readeqclasses(args.strEqs)
    print("...Eq classes read in...")
    print("...Processing eq classes...")
    hTrEqSumm, hGeneEqSumm = proceqclasses(hNums, hTrs, aiHapOne, aiHapTwo, aaEqs, hTrGeneMap, hGeneTrMap)
    print("...eq classes processed...")
    print("...writing output...")
    writeout(hTrEqSumm, hGeneEqSumm, args.strOutBase)
    print("...output written...")
    
    ### *******try encoding thing in reverse for gzip open before giving up on that***********

if __name__=="__main__":
	main()

# FIX GZIP ISSUE.....doesn't work for input; some weird encoding






