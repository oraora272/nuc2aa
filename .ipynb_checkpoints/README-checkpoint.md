# nuc2aa
Bioinformatics tool to convert nucleotide mutations to amino acid mutations.

Example: Enter the nucleotide mutation 679_C_A and a reference sequence; the function will return the resulting amino acid mutation L219I.

## Usage:

mut="3_a_t"
ref="aaatgtttttt"
n_bef=2
n_aft=6

nuc2aa(mut,ref,n_bef,n_aft,known_prior=True)

Will return:

'M1L'

## Anticipated updates:

- Allow ambiguous nucleotides

- Return position in case of insertion/deletion

- Automatic detection of start and stop codon, with nucleotides before/after