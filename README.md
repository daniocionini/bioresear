# bioresear
A bioinformatics workflow manager that allows to standardise genome data formats and perform protein sequencing

--> Main Functionalities <--
.id()
-> extracts informations about the sequence imported

.dna()
-> extracts the dna sequence imported

.compl()
-> extracts the dna sequence imported and compute the complementary dna sequence (base change)

.freq()
-> computes information about the frequency of the dna nucleotides

.rna()
-> computes the rna sequence from the complementary dna sequence

.gc()
-> computes the GC content of the dna sequence

.base()
-> computes the rna for 3 base pairs

.open_r()
-> computes the open reading frame of the dna and complementary dna sequences

.prot_seq()
-> computes the protein bases extracted from the open reading frame dna sequences

.protein()
-> computes the proteins present in the dna sequence imported for blast evaluation

.blast()
-> 



--> Additional Functionalities <--
1) Possibility to export protein outputs.
Nedded to create a folder called "Exported Proteins" when set up the whole program

2) In order to use .blast() you need to save the proteins_file.txt into directory as stated into the documentation

