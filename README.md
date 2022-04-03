# bioresear
A bioinformatics workflow manager that allows to standardise genome data formats and perform protein sequencing

To start using the functionalities start by using file.id()

--> Main Functionalities <--
.id()
-> extracts the informations about the sequence imported from the file

.dna()
-> extracts only the dna sequence imported from the file

.compl()
-> prints the dna sequence imported and compute the complementary dna sequence (base change)

.freq()
-> computes information about the frequency of the dna nucleotides (A, G, C, T)

.rna()
-> prints the complementary sequence and computes the rna sequence from the complementary dna sequence

.gc()
-> computes the GC content of the dna sequence (in %)

.base()
-> computes and prints the rna in 3 base pairs

.open_r()
-> computes the open reading frame of the dna and complementary dna sequences

.prot_seq()
-> computes the protein bases extracted from the open reading frame dna sequences

.protein()
-> computes the proteins sequences in the dna sequence imported for blast evaluation while saving the file as .txt to the directory set as "direct_protein"


--> Additional Functionalities <--
1) Possibility to export protein outputs.
Nedded to create a folder called "Exported Proteins" when set up the whole program
