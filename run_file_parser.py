from structures import *
import sys
import io
import os
from io import StringIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pypdb
import untangle
import requests

class bioresear:

    def __init__(self, input, file_type):
        self.input = input
        self.type = file_type
        self.doc = (open("D:\GitHub\easy_seq_file_parser\documentation.txt", "r")).read()


    def id(self):
        if self.type == 'fasta':
            id = txt_file[:txt_file.index("\n")]
            print(id)
        else:
            id = txt_file[:txt_file.index("other;")]
            print(id)


    def dna(self):
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            print(dseq)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            print(dseq)



    def compl(self):
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            print(dseq)
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            print(cdna)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            print(dseq)
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            print(cdna)



    def freq(self):
        all_freq = {}
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for i in dseq:
                if i in all_freq:
                    all_freq[i] += 1
                else:
                    all_freq[i] = 1
            print(all_freq)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for i in dseq:
                if i in all_freq:
                    all_freq[i] += 1
                else:
                    all_freq[i] = 1
            print(all_freq)


    def rna(self):
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            rseq = cdna.replace('T', 'U')
            print(cdna)
            print(rseq)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            rseq = cdna.replace('T', 'U')
            print(cdna)
            print(rseq)


    def gc(self):
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            c_cnt = dseq.count("C")
            g_cnt = dseq.count("G")
            tot = dseq.count("")
            unk_dna = c_cnt + g_cnt
            gc_count = (unk_dna * 100) / tot
            gc_count = round(gc_count, 2)
            print(gc_count, '%')
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            c_cnt = dseq.count("C")
            g_cnt = dseq.count("G")
            tot = dseq.count("")
            unk_dna = c_cnt + g_cnt
            gc_count = (unk_dna * 100) / tot
            gc_count = round(gc_count, 2)
            print(gc_count, '%')


    def base(self):
        n = 3
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            rseq = cdna.replace('T', 'U')
            base = [rseq[i:i+n] for i in range(0, len(rseq), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            print(base)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            rseq = cdna.replace('T', 'U')
            base = [rseq[i:i+n] for i in range(0, len(rseq), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            print(base)


    def open_r(self):
        n = 3
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            print(uno_frame)
            print(due_frame)
            print(tre_frame,'\n')
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            print(n_uno_frame)
            print(n_due_frame)
            print(n_tre_frame)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            print(uno_frame)
            print(due_frame)
            print(tre_frame, '\n')
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            print(n_uno_frame)
            print(n_due_frame)
            print(n_tre_frame)

    def prot_seq(self):
        n = 3
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            for pair, prot in dna_codons.items():
                uno_frame = uno_frame.replace(pair.upper(), prot.lower())
            transl_1 = uno_frame.upper()
            print(transl_1)
            for pair, prot in dna_codons.items():
                due_frame = due_frame.replace(pair.upper(), prot.lower())
            transl_2 = due_frame.upper()
            print(transl_2)
            for pair, prot in dna_codons.items():
                tre_frame = tre_frame.replace(pair.upper(), prot.lower())
            transl_3 = tre_frame.upper()
            print(transl_3)
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            for pair, prot in dna_codons.items():
                n_uno_frame = n_uno_frame.replace(pair.upper(), prot.lower())
            n_transl_1 = n_uno_frame.upper()
            print(n_transl_1)
            for pair, prot in dna_codons.items():
                n_due_frame = n_due_frame.replace(pair.upper(), prot.lower())
            n_transl_2 = n_due_frame.upper()
            print(n_transl_2)
            for pair, prot in dna_codons.items():
                n_tre_frame = n_tre_frame.replace(pair.upper(), prot.lower())
            n_transl_3 = n_tre_frame.upper()
            print(n_transl_3)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            for pair, prot in dna_codons.items():
                uno_frame = uno_frame.replace(pair.upper(), prot.lower())
            transl_1 = uno_frame.upper()
            print(transl_1)
            for pair, prot in dna_codons.items():
                due_frame = due_frame.replace(pair.upper(), prot.lower())
            transl_2 = due_frame.upper()
            print(transl_2)
            for pair, prot in dna_codons.items():
                tre_frame = tre_frame.replace(pair.upper(), prot.lower())
            transl_3 = tre_frame.upper()
            print(transl_3)
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            for pair, prot in dna_codons.items():
                n_uno_frame = n_uno_frame.replace(pair.upper(), prot.lower())
            n_transl_1 = n_uno_frame.upper()
            print(n_transl_1)
            for pair, prot in dna_codons.items():
                n_due_frame = n_due_frame.replace(pair.upper(), prot.lower())
            n_transl_2 = n_due_frame.upper()
            print(n_transl_2)
            for pair, prot in dna_codons.items():
                n_tre_frame = n_tre_frame.replace(pair.upper(), prot.lower())
            n_transl_3 = n_tre_frame.upper()
            print(n_transl_3)


    def protein(self):
        old_stdout = sys.stdout
        new_stdout = io.StringIO()
        sys.stdout = new_stdout
        n = 3
        if self.type == 'fasta':
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            for pair, prot in dna_codons.items():
                uno_frame = uno_frame.replace(pair.upper(), prot.lower())
            transl_1 = uno_frame.upper()
            transl_1 = transl_1.replace(" ", "")
            for pair, prot in dna_codons.items():
                due_frame = due_frame.replace(pair.upper(), prot.lower())
            transl_2 = due_frame.upper()
            transl_2 = transl_2.replace(" ", "")
            for pair, prot in dna_codons.items():
                tre_frame = tre_frame.replace(pair.upper(), prot.lower())
            transl_3 = tre_frame.upper()
            transl_3 = transl_3.replace(" ", "")
            current_prot = []
            proteins = []
            for aa in transl_1:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            current_prot = []
            proteins = []
            for aa in transl_2:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            for aa in transl_3:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            dseq = txt_file[txt_file.index("\n"):]
            dseq = dseq.replace('\n','')
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            for pair, prot in dna_codons.items():
                n_uno_frame = n_uno_frame.replace(pair.upper(), prot.lower())
            n_transl_1 = n_uno_frame.upper()
            n_transl_1 = n_transl_1.replace(" ", "")
            for pair, prot in dna_codons.items():
                n_due_frame = n_due_frame.replace(pair.upper(), prot.lower())
            n_transl_2 = n_due_frame.upper()
            n_transl_2 = n_transl_2.replace(" ", "")
            for pair, prot in dna_codons.items():
                n_tre_frame = n_tre_frame.replace(pair.upper(), prot.lower())
            n_transl_3 = n_tre_frame.upper()
            n_transl_3 = n_transl_3.replace(" ", "")
            for aa in n_transl_1:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            current_prot = []
            proteins = []
            for aa in n_transl_2:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            for aa in n_transl_3:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            output = new_stdout.getvalue()
            sys.stdout = old_stdout
            output = output.replace("[", "")
            output = output.replace("]", "")
            output = output.replace("' ", "")
            output = output.replace("'", "")
            output = output.replace(", ", "\n")
            print(output)
            with open(direct_protein, "w") as f:
                for line in output:
                    f.write(line)
        else:
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            tre_frame = base3
            for pair, prot in dna_codons.items():
                uno_frame = uno_frame.replace(pair.upper(), prot.lower())
            transl_1 = uno_frame.upper()
            transl_1 = transl_1.replace(" ", "")
            for pair, prot in dna_codons.items():
                due_frame = due_frame.replace(pair.upper(), prot.lower())
            transl_2 = due_frame.upper()
            transl_2 = transl_2.replace(" ", "")
            for pair, prot in dna_codons.items():
                tre_frame = tre_frame.replace(pair.upper(), prot.lower())
            transl_3 = tre_frame.upper()
            transl_3 = transl_3.replace(" ", "")
            current_prot = []
            proteins = []
            for aa in transl_1:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            current_prot = []
            proteins = []
            for aa in transl_2:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            for aa in transl_3:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            dseq = txt_file[txt_file.index("other;"):]
            dseq = dseq.replace('other;','' )
            dseq = ''.join([c for c in dseq if c not in list_of_numbers])
            dseq = dseq.replace('/', '')
            dseq = dseq.replace(' ', '')
            dseq = dseq.replace('\n', '')
            dseq = dseq.upper()
            for nuc, comp in {"T":"A", "C":"G", "G":"C", "A":"T"}.items():
                dseq = dseq.replace(nuc.upper(),comp.lower())
            cdna = dseq.upper()
            cdna =''.join(reversed(cdna))
            base = [cdna[i:i+n] for i in range(0, len(cdna), n)]
            base = str(base)
            base = base.replace("[", "")
            base = base.replace("]", "")
            base = base.replace("'", "")
            base = base.replace(",", "")
            n_uno_frame = base
            cdna2 = cdna[1:]
            base2 = [cdna2[i:i+n] for i in range(0, len(cdna2), n)]
            base2 = str(base2)
            base2 = base2.replace("[", "")
            base2 = base2.replace("]", "")
            base2 = base2.replace("'", "")
            base2 = base2.replace(",", "")
            n_due_frame = base2
            cdna3 = cdna[2:]
            base3 = [cdna3[i:i+n] for i in range(0, len(cdna3), n)]
            base3 = str(base3)
            base3 = base3.replace("[", "")
            base3 = base3.replace("]", "")
            base3 = base3.replace("'", "")
            base3 = base3.replace(",", "")
            n_tre_frame = base3
            for pair, prot in dna_codons.items():
                n_uno_frame = n_uno_frame.replace(pair.upper(), prot.lower())
            n_transl_1 = n_uno_frame.upper()
            n_transl_1 = n_transl_1.replace(" ", "")
            for pair, prot in dna_codons.items():
                n_due_frame = n_due_frame.replace(pair.upper(), prot.lower())
            n_transl_2 = n_due_frame.upper()
            n_transl_2 = n_transl_2.replace(" ", "")
            for pair, prot in dna_codons.items():
                n_tre_frame = n_tre_frame.replace(pair.upper(), prot.lower())
            n_transl_3 = n_tre_frame.upper()
            n_transl_3 = n_transl_3.replace(" ", "")
            for aa in n_transl_1:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            current_prot = []
            proteins = []
            for aa in n_transl_2:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            for aa in n_transl_3:
                if aa == '_':
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append(" ")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            print(proteins)
            output = new_stdout.getvalue()
            sys.stdout = old_stdout
            output = output.replace("[", "")
            output = output.replace("]", "")
            output = output.replace("' ", "")
            output = output.replace("'", "")
            output = output.replace(", ", "\n")
            print(output)
            with open(direct_protein, "w") as f:
                for line in output:
                    f.write(line)

                    
#-------------------------------------------- PRINT MAIN
opener = open("YOUR_PATH", "r")
direct_protein = "YOUR_PATH\protein_list.txt"
txt_file = opener.read()
file = bioresear(txt_file, 'fasta') #file type
#--------------------------------------------
