import pandas as pd #need to compare to flumut
import math #need math.ceil for aa position

import pandas as pd #need to compare to flumut
import math #need math.ceil for aa position

CODON = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',

    # All codons containing N -> ambiguous
    'NAA':'ambiguous','NAC':'ambiguous','NAG':'ambiguous','NAT':'ambiguous',
    'NCA':'ambiguous','NCC':'ambiguous','NCG':'ambiguous','NCT':'ambiguous',
    'NGA':'ambiguous','NGC':'ambiguous','NGG':'ambiguous','NGT':'ambiguous',
    'NTA':'ambiguous','NTC':'ambiguous','NTG':'ambiguous','NTT':'ambiguous',

    'ANA':'ambiguous','ANC':'ambiguous','ANG':'ambiguous','ANT':'ambiguous',
    'CNA':'ambiguous','CNC':'ambiguous','CNG':'ambiguous','CNT':'ambiguous',
    'GNA':'ambiguous','GNC':'ambiguous','GNG':'ambiguous','GNT':'ambiguous',
    'TNA':'ambiguous','TNC':'ambiguous','TNG':'ambiguous','TNT':'ambiguous',

    'AAN':'ambiguous','ACN':'ambiguous','AGN':'ambiguous','ATN':'ambiguous',
    'CAN':'ambiguous','CCN':'ambiguous','CGN':'ambiguous','CTN':'ambiguous',
    'GAN':'ambiguous','GCN':'ambiguous','GGN':'ambiguous','GTN':'ambiguous',
    'TAN':'ambiguous','TCN':'ambiguous','TGN':'ambiguous','TTN':'ambiguous',

    'NNN':'ambiguous'
}

def nuc2aa(mut1,gappy_ref,n_before_orf,n_after_orf,known_prior=True): 

    '''
    mut1: string in format such as 505_G_A (Position_BeforeNucleotide_AfterNucleotide)
        position should be ALIGNMENT/COLUMN position, NOT nucleotide position (see note 3). It includes any gaps before
        BeforeNucleotide is whatever the nucleotide is before mutation, AferNucleotide is after mutation. If Before_Nucleotide
        is unknown, see variable known_prior
    gappy_ref: string, should be the reference sequence aligned to the rest of the sequences while maintaining length of alignment, e.g.
        mafft --add refseq.fasta --keeplength alignment.fasta > ref_and_alignment.fasta
    n_before_orf: interger, how many nucleotides before the open reading frame. ==21 nucleotides for HA
    n_after_orf: interger, last nucleotide position of the open reading frame. == n_before_orf + len(orf) == 1728 for HA 
    known_prior: set to False if the prior nucleotide is unknown, and set mut1 in format Position_AfterNucleotide

    note1: if nucleotide is ambiguous, it will return as ambiguous
    note2: any U should be converted to T
    note3: this function also works with the "true" nucleotide position rather than column position i.e. in mut1 variable
    However, the input reference sequence must not contain gaps in this case (raw sequence, not aligned)
    
    2026/10/19
    Deanna Tait
    '''
    
    #just in case mutation and/or sequence aren't capitalized
    mut1 = mut1.upper()
    gappy_ref = gappy_ref.upper()
    
    #Get nucleotide position
    column_pos = int(mut1.split(sep = '_')[0]) #column position, aligned to all amino acids
    cut_ref = gappy_ref[0:column_pos].replace('-','') #remove all trailing and gaps. THIS LEAVES THE TRUE NUCLEOTIDE POSITION IN GSGD, MARRIES COLUMN WITH GSGD
    nuc_pos = len(cut_ref) #final nucleotide position

    #Get amino acid position
    orf = gappy_ref.replace('-','')[n_before_orf:n_after_orf] #make open reading frame 
    len_to_mut = cut_ref[n_before_orf:]
    aa_pos = math.ceil(len(len_to_mut)/3)
    
    #Check type
    if nuc_pos <= (n_before_orf): #outside of open reading frame
        return("non-coding before")
    elif nuc_pos > n_after_orf: #outside of open reading frame
        return("non-coding after")

    if known_prior==True:
        if mut1.split(sep = '_')[2] not in ["A","C","G","T","*"] or mut1.split(sep = '_')[1] not in ["A","C","G","T","*"]:
            return ("invalid")
        elif mut1.split(sep = '_')[1] == "*":
            return("insertion")
        elif mut1.split(sep = '_')[2] == "*":
            return("deletion")
    if known_prior==False:
        if mut1.split(sep = '_')[1] not in ["A","C","G","T","*"]:
            return ("invalid")
        elif mut1.split(sep = '_')[1] == "*":
            return("deletion")        

    
    #SNP
    orf_pos = nuc_pos-n_before_orf #position in open reading frame

    if known_prior==True: #mut1 provided in the format Position_BeforeNucleotide_AfterNucleotide
        if orf_pos % 3 == 1: #first position
            #print full aa mutation
            return(CODON[mut1.split(sep = '_')[1]+str(orf[orf_pos])+str(orf[orf_pos+1])]+ 
                  str(aa_pos)+#AA position
                  CODON[mut1.split(sep = '_')[2]+str(orf[orf_pos])+str(orf[orf_pos+1])])
        
        elif orf_pos % 3 == 2: #second position
            #print full aa mutation
            return(CODON[orf[orf_pos-2]+mut1.split(sep = '_')[1]+orf[orf_pos]]+
                  str(aa_pos)+#AA position
                  CODON[orf[orf_pos-2]+mut1.split(sep = '_')[2]+orf[orf_pos]])
            
        elif orf_pos % 3 == 0: #third position
            #print full aa mutation
            return(CODON[orf[orf_pos-3]+orf[orf_pos-2]+mut1.split(sep = '_')[1]]+
                  str(aa_pos)+#AA position
                  CODON[orf[orf_pos-3]+orf[orf_pos-2]+mut1.split(sep = '_')[2]])
            
    if known_prior==False: #mut1 provided in the format Position_AfterNucleotide
        if orf_pos % 3 == 1: #first position
            return(CODON[str(orf[orf_pos-1])+str(orf[orf_pos])+str(orf[orf_pos+1])]+ 
                  str(aa_pos)+#AA position
                  CODON[mut1.split(sep="_")[1].upper()+str(orf[orf_pos])+str(orf[orf_pos+1])])
        elif orf_pos % 3 == 2: #second position
            return(CODON[orf[orf_pos-2]+orf[orf_pos-1]+orf[orf_pos]]+
                  str(aa_pos)+#AA position
                  CODON[orf[orf_pos-2]+mut1.split(sep="_")[1].upper()+orf[orf_pos]])
            
        elif orf_pos % 3 == 0: #third position
            return(CODON[orf[orf_pos-3]+orf[orf_pos-2]+orf[orf_pos-1]]+
                  str(aa_pos)+#AA position
                  CODON[orf[orf_pos-3]+orf[orf_pos-2]+mut1.split(sep="_")[1].upper()])
    
def get_type(df):

    '''
    Takes mutation in format "AminoacidPositionAminoacid" e.g. I67G. Must be single-letter format. 
    df: pandas dataframe with at least columns "AA mutation" and "type" and at least one row
    It may also be specifically one of the strings: [insertion, non-coding before, unknkown after, invalid, deletion]
    If it does not follow any of those formats it will be misclassified as a non-synonymous mutation
    
    2026/08/19
    Author: Deanna Tait
    '''
    
    for i in range(len(df)):
    	if df.loc[i,"AA mutation"] == "insertion": 
    		df.loc[i,"type"] = "insertion"
    	elif df.loc[i,"AA mutation"] == "non-coding before":
    		df.loc[i,"type"] = "non-coding before"
    	elif df.loc[i,"AA mutation"] == "unknown after":
    		df.loc[i,"type"] = "unknown after"
    	elif df.loc[i,"AA mutation"] == "invalid":
    		df.loc[i,"type"] = "invalid"
    	elif df.loc[i,"AA mutation"][0] == "deletion":
    		df.loc[i,"type"] = "deletion"
    	elif df.loc[i,"AA mutation"][0] == df.loc[i,"AA mutation"][-1]:
    		df.loc[i,"type"] = "synonymous"
    	else:
    		df.loc[i,"type"] = "non-synonymous"
    return(df)

def check_ifin_flumut(df,flumut_set):

    '''
    Check if a mutation is in FluMut tool database (https://izsvenezie-virology.github.io/FluMut/docs/markers)
    df: pandas dataframe with at least columns "AA mutation" "In FluMut" and "In FluMut EX" and at least one row
        AA mutation should be a string, other two columns set as False
    flumut_set: set of mutations (string)
    format of both df["AA mutation"] and flumut_set should be AminoacidPositionAminoacid e.g. I67K

    Note: Position of HA is different, -16 because signal peptide removed. Explained here:
    https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/numbering.html
    They subtract 16 and use an American Wigeon from 2021, and after checking that the 
    sequence of my own american wigeon from 2021 matched the first 17 amino acids of GSGD,
    Subtract 16 as well. Any mutations in the signal peptide region cannot be detected. 
    This does not change the position as listed by the AA mutation column

    2026/08/19
    Deanna Tait
    '''
    

    flumut_ex = list(flumut_set)
    for i in range(len(flumut_ex)):
        flumut_ex[i] = flumut_ex[i][1:]  
    #Check if mutation is in either flumut or the flumut after list
    for i in range(len(df)):
        if df.loc[i,"AA mutation"] in flumut_set:
            df.loc[i,"In FluMut"] = True
            df.loc[i,"In FluMut EX"] = True
        elif df.loc[i,"AA mutation"][1:] in flumut_ex:
            df.loc[i,"In FluMut EX"] = True
    if is_ha==True:
        for i in range(len(df)):
            if df.loc[i,"AA mutation"][1].isnumeric()==True:
                if df.loc[i,"AA mutation"][0]+str(int(df.loc[i,"AA mutation"][1:-1])-16)+df.loc[i,"AA mutation"][-1] in flumut_set:
                    df.loc[i,"In FluMut"] = True
                    df.loc[i,"In FluMut EX"] = True
                elif str(int(df.loc[i,"AA mutation"][1:-1])-16)+df.loc[i,"AA mutation"][-1] in flumut_ex:
                    df.loc[i,"In FluMut EX"] = True
    return(df)

def main():
    pass
if __name__ == "__main__":
    main()