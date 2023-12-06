# GAMMA_DB_Maker
GAMMA expects a gene database with all of the genes in the forward reading frame with stop codons. This script will take in a nucleotide multifasta of genes and put them into the correct reading frame for translating the gene with GAMMA. It will remove any nonstandard (i.e., NOT "A", "C", "T", or "G" bases), consider all 6 potential reading frames (3 forward and 3 reverse) to find a start codon and no nonsense mutations in the gene, and will remove any genes with the same sequence.

***Note: The script expects your DB name to have a ".fasta" extension, so you may encounter errors if it is ".fna" or something other than ".fasta"***

Usage:

> python3 GAMMA_DB_Maker.py My_Genes_DB.fasta

It will generate the following output fasta files:
  - My_Genes_DB_rb.fasta: Fasta file only with standard bases (any genes with nonstandard bases are removed).
  - My_Genes_DB_nsb.fasta: Fasta file with genes that have nonstandard bases, these are excluded from the final output.
  - My_Genes_DB_nonstop.fasta: Fasta file with genes that do not have a stop codon at the end of the gene.
  - My_Genes_DB_nonstandard: Fasta file with genes that the script could not find one or more of a standard start or stop codon, and/or had nonsense mutations from any of the six reading frames. This is usually an indication of an incomplete/noncoding gene sequence.
  - My_Genes_DB_correct.fasta: Fasta file with genes in the correct reading frame with start and stop codons and no nonsense mutations.
  - My_Genes_DB_correct.fasta: Fasta file with genes in the correct reading frame with start and stop codons and no nonsense mutations, and any duplicates removed.

You may want to still use/edit/add genes from any of these outputs back into your final DB, but GAMMA will point out these changes (i.e., if you have a partial gene without a stop codon at the end it will be labeled as a "nonstop mutation" even if there is an exact match in your assembly to the gene, or a truncation will be identified in the gene if it's an exact match if there's a truncation in the gene in your database). Of course, if you're only interested in seeing the nucleotide overlap (without translating your gene, you can use GAMMA-S, which is just for the sequence.
