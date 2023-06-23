"""
This module contains the functions to rename and retrieve sequences from the T. congo fasta files generated from the TransDecoder output
"""

#rename fasta file with protein coding sequences inside
def rename_seq(filename):
    
    from Bio import SeqIO
    import re

    """
    This function generates a new name for sequences in the fasta file passed.

    Args:
        filename (fasta file): fasta file containing sequences to rename

    Returns:
        fasta file: filename with new sequence names

    Examples:
        > rename_seq('Tcongo_epimastigote.adj.fas.transdecoder.pep')

    will return value example
    """
    
    #creating a new file to input the data into
    newfile = open(f"named_{filename}", "a")
    
    #setting seqID count to 0
    n = 0
    
    #compiling regular expression rules to identify life cycle stage
    pattern = re.compile(r'(?<![\w])TCep(?![a-z])'
                         r'|(?<![\w])TCbln(?![a-z])'
                         r'|(?<![\w])Tcbl(?![a-z])'
                         r'|(?<![\w])TCmcn(?![a-z])'
                         r'|(?<![\w])TCmc(?![a-z])'
                         r'|(?<![\w])TCpc(?![a-z])', flags=re.IGNORECASE)
    
    for record in SeqIO.parse(filename, "fasta"):
        # this regular expression searches for and extracts the T congo growth stage
        ##it keeps a count of the seqID for each sequence for the lifecyle stage
        real_name = pattern.search(record.description)
        if real_name:
            if 'TCep' in str(real_name) or 'Tcep' in str(real_name):
                n = n + 1
                newfile.write(f">T.congo_epimastigote {n}")
            elif 'TCbln' in str(real_name):
                n += 1
                newfile.write(f">T.congo_nor_bloodstream {n}")
            elif 'TCmcn' in str(real_name):
                n = n + 1
                newfile.write(f">T.congo_nor_metacyclic {n}")
            elif 'TCpc' in str(real_name):
                n = n + 1
                newfile.write(f">T.congo_procyclic {n}")
            elif 'TCbl' in str(real_name):
                n = n + 1
                newfile.write(f">T.congo_reg_bloodstream {n}")
            elif 'TCmc' in str(real_name):
                n = n + 1
                newfile.write(f">T.congo_reg_metacyclic {n}")

        # regular expression that searches for and extracts orf type of the sequence
        orf = re.search("ORF type:.*? ", record.description)
        if orf:
            newfile.write(f" {orf.group().replace(':', '-').replace('ORF type', 'ORF_type')}")

        # regular expression that searches for and extracts the length of the sequence
        le = re.search("len:.*? ", record.description)
        if le:
            newfile.write(f" {le.group().replace(':', '-')}")

        # regular expression that searches for and extracts the score of the CDS
        score = re.search("score=.*? ", record.description)
        if score:
            newfile.write(f" {score.group().replace('=', '-')}")

        # regular expression that searches for and extracts the bases of the CDS
        bases = re.search("[0-9]+-[0-9]+[(](.)[)]", record.description)
        if bases:
            newfile.write(f" bases;{bases.group().replace('(',' ').rstrip (')')}")

        file = newfile.write(f"\n{record.seq}\n")
        
    return file


#retrieve sequence by seqID from rename_seq() output file
def retrieve_seq(named_fasta_file, seqID):
    """
    This function searches for and outputs to a new file the sequence assigned to the seqID.

    Args:
        named_fasta_file (fasta file): fasta file containing named sequences with seqID number

    Returns:
        fasta file: filename with the sequence and header corresponding to the seqID number

    Examples:
        > retrieve_seq('named_Tcongo_nor_bloodstream.adj.fas.transdecoder.pep', 201)

    will return value example
    """
    
    from Bio import SeqIO
    import re
    import pytest
    
    #compiling regular expression rules to identify life cycle stage and seqID
    pattern = re.compile(r'T.congo_epimastigote\s[0-9]+\s'
                         r'|T.congo_nor_bloodstream\s[0-9]+\s'
                         r'|T.congo_reg_bloodstream\s[0-9]+\s'
                         r'|T.congo_nor_metacyclic\s[0-9]+\s'
                         r'|T.congo_reg_metacyclic\s[0-9]+\s'
                         r'|T.congo_procyclic\s[0-9]+\s', 
                        flags=re.IGNORECASE)
    
    #warning flag set up
    error_flag = False
    
    #searches through each fasta sequence header 
    for seq in SeqIO.parse(named_fasta_file, 'fasta'):
        
        #searches for the seqID passed by the user
        number = re.search(f" {seqID} ", seq.description)
        if number:
            #if the number is found in the fasta header it proceeds to the next loop which searches for the life cycle stage
            
            #number found so no error to report
            error_flag = True
            stage = pattern.search(seq.description)
            if stage:
                #if one of the life cycle stages below is found in the sequence header, then open new file with life cycle name and seqID number in the name
                if 'T.congo_epimastigote' in str(stage):
                    newfile = open(f"Tcongo_epimastigote_{seqID}.fas", "a")
                
                elif 'T.congo_nor_bloodstream' in str(stage):
                    newfile = open(f"Tcongo_nor_bloodstream_{seqID}.fas", "a")
                
                elif 'T.congo_nor_metacyclic' in str(stage):
                    newfile = open(f"Tcongo_nor_metacyclic_{seqID}.fas", "a")
                
                elif 'T.congo_procyclic' in str(stage):
                    newfile = open(f"Tcongo_procyclic_{seqID}.fas", "a")
                
                elif 'T.congo_reg_bloodstream' in str(stage):
                    newfile = open(f"Tcongo_reg_bloodstream_{seqID}.fas", "a")
                
                elif 'T.congo_reg_metacyclic' in str(stage):
                    newfile = open(f"Tcongo_reg_metacyclic_{seqID}.fas", "a")
        
            #append the sequence header and the coding sequence to the new file
            newfile = newfile.write(f"{seq.description}\n{seq.seq}")
    
    if error_flag == True:
        return newfile
    elif error_flag == False:
        raise ValueError(f"Sequence ID {seqID} not found in {named_fasta_file}")
        return