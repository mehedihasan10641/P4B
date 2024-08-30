# This is a sample Python script.

# Calculating AT content

# dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
#
# x = dna.count("A")
# y = dna.count("T")
# len_dna = len(dna)
# at_content = (x+y)/len_dna
# print(at_content)

# Complementing DNA
# dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
# w = dna.replace("A","t")
# x = w.replace("T","a")
# y = x.replace("G","c")
# z = y.replace("C","g")
# print(z.upper())

# Restriction fragment lengths
# dna = "ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT"
# cut_pos = dna.find("AATTC")
# frag1 = dna[:cut_pos]
# frag2 = dna[cut_pos:]
# print("1st fragment" , frag1 ,"\nLength" , len(frag1) , "\n2nd fragment" , frag2, "\nLength" , len(frag2))

# Splicing out introns, part one and part two
# dna = "ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"
# exon1 = dna[:63]
# exon2 = dna[90:]
# exon = exon1 + exon2
# print("Coding region is", exon)
# intron = dna[63:90]
# print("Percentage of coding seq", (len(exon)/len(dna))*100)
# print(exon1.upper()+intron.lower()+exon2.upper())


# Splitting genomic DNA
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\reading_files\exercises\genomic_dna.txt")
# dna = file.read().rstrip("\n")
# exon1 = dna[:63]
# exon2 = dna[90:]
# exon = exon1 + exon2
# intron = dna[63:90]
# with open("coding.txt","w") as file1, open("non_coding.txt","w") as file2:
#     file1.write(exon)
#     file2.write(intron)

# Writing a fasta file
# file = open("fasta_file.txt","a")
# heading1 = ">ABC123"
# heading2 = ">DEF456"
# heading3 = ">HIJ789"
# seq1 = "ATCGTACGATCGATCGATCGCTAGACGTATCG\n"
# seq2 = "actgatcgacgatcgatcgatcacgact\n"
# seq3 = "ACTGAC-ACTGT--ACTGTA----CATGTG\n"
# file.write(heading1 + "\n" + seq1)
# file.write(heading2 + "\n" + seq2.upper())
# file.write(heading3 + "\n" + seq3.replace("-",""))
# file.close()

# Writing multiple FASTA files

# heading1 = ">ABC123"
# heading2 = ">DEF456"
# heading3 = ">HIJ789"
# seq1 = "ATCGTACGATCGATCGATCGCTAGACGTATCG\n"
# seq2 = "actgatcgacgatcgatcgatcacgact\n"
# seq3 = "ACTGAC-ACTGT--ACTGTA----CATGTG\n"
# with open("ABC123.fasta","w") as file1, open("DEF456.fasta","w") as file2, open("HIJ789.fasta","w") as file3:
#     file1.write(heading1 + "\n" + seq1)
#     file2.write(heading2 + "\n" + seq2.upper())
#     file3.write(heading3 + "\n" + seq3.replace("-", ""))


# Processing DNA in a file
# seq = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\lists_and_loops\exercises\input.txt")
# for line in seq:
#     trim = open("trimmed.txt","a")
#     clean_seq = line[14:]
#     trim.write(clean_seq)
#     print(line + "Sequence Length is",len(line))

# Multiple exons from genomic DNA
# file1 = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\lists_and_loops\exercises\genomic_dna.txt")
# file2 =  open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\lists_and_loops\exercises\exons.txt")
# dna = file1.read()
# for line in file2:
#     position = line.rstrip().split(",")
#     start_pos = int(position[0])
#     stop_pos = int(position[1])
#     exon_seq = dna[start_pos:stop_pos]
#     print(exon_seq)
#     file = open("exon.txt","a")
#     file.write(exon_seq)


# Percentage of amino acid residues, part one
# def my_function(protein,amino):
#     protein_len = len(protein)
#     amino_count = protein.count(amino.upper())
#     protein_per = (amino_count/protein_len)*100
#     return round(protein_per)
#
# assert my_function("MSRSLLLRFLLFLLLLPPLP", "M") == 5
# assert my_function("MSRSLLLRFLLFLLLLPPLP", "r") == 10
# assert my_function("MSRSLLLRFLLFLLLLPPLP", "L") == 50
# assert my_function("MSRSLLLRFLLFLLLLPPLP", "Y") == 0


# Percentage of amino acid residues, part two
# def my_function(protein,amino = ["A","I","L","M","F","W","Y","V"]):
#     protein_len = len(protein)
#     amino_count = 0
#     for acid in amino:
#         amino_num = protein.count(acid)
#         amino_count += amino_num
#     protein_per = (amino_count/protein_len)*100
#     return round(protein_per)
#
# assert my_function("MSRSLLLRFLLFLLLLPPLP", ["M"]) == 5
# assert my_function("MSRSLLLRFLLFLLLLPPLP", ['M', 'L']) == 55
# assert my_function("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L']) == 70
# assert my_function("MSRSLLLRFLLFLLLLPPLP") == 65


# Several species
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\conditional_tests\exercises\data.csv")
# for line in file:
#     data = line.rstrip("\n").split(",")
#     species_name = data[0]
#     seq = data[1]
#     gene = data[2]
#     express = data[3]
#     if species_name == "Drosophila melanogaster" or species_name == "Drosophila simulans":
#         print(gene)

# Length range
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\conditional_tests\exercises\data.csv")
# for line in file:
#     data = line.rstrip("\n").split(",")
#     species_name = data[0]
#     seq = data[1]
#     gene = data[2]
#     express = data[3]
#     if len(seq) >= 90 and len(seq) <= 110:
#         print(gene)

# AT content

# def my_function(dna, sig_figs = 2):
#     length = len(dna)
#     dna = dna.upper()
#     a_count = dna.count("A")
#     t_count = dna.count("T")
#     at_content = (a_count + t_count)/length
#     return round(at_content,sig_figs)
#
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\conditional_tests\exercises\data.csv")
# for line in file:
#     data = line.rstrip("\n").split(",")
#     seq = data[1]
#     gene = data[2]
#     express = int(data[3])
#     if my_function(seq) < 0.5 and express > 200:
#         print(gene)

# Complex condition
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\conditional_tests\exercises\data.csv")
# for line in file:
#     data = line.rstrip("\n").split(",")
#     species_name = data[0]
#     seq = data[1]
#     gene = data[2]
#     express = int(data[3])
#     if (gene.startswith("k") or gene.startswith("h")) and not species_name == "Drosophila melanogaster":
#         print(gene)

# High Low Medium
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\conditional_tests\exercises\data.csv")
# for line in file:
#     data = line.rstrip("\n").split(",")
#     species_name = data[0]
#     seq = data[1]
#     gene = data[2]
#     express = int(data[3])
#     if my_function(seq) > 0.65:
#         print("AT content is high.")
#     elif my_function(seq) <0.65 and my_function(seq) >0.45:
#         print("AT content is medium")
#     else:
#         print("AT content is low")

# Accession names
#access = ["xkn59438", "yhdck2", "eihd39d9", "chdsye847", "hedle3455", "xjhd53e", "45da", "de37dp"]
#import re
# for x in access:
#     if re.search(r"5",x):
#         print(x)

# for x in access:
#     if re.search(r"d|e",x):
#         print(x)

# for x in access:
#     if re.search(r"de",x):
#         print(x)

# for x in access:
#     if re.search(r"d.e",x):
#         print(x)

# for x in access:
#     if re.search(r"(d.*e|e.*d)",x):
#         print(x)

# for x in access:
#     if x.startswith("x") or x.startswith("y"):
#         print(x)

# for x in access:
#     if (x.startswith("x") or x.startswith("y")) and x.endswith("e"):
#         print(x)

# for x in access:
#     if re.search(r"[12345789]{3,100}",x):
#         print(x)

# for x in access:
#     if re.search(r"d[arp]$",x):
#         print(x)


# Double digest
# file = open(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\regular_expressions\exercises\dna.txt")
# seq = file.read().rstrip("\n")
#
# import re
# x = re.finditer(r"A[ATCG]TAAT",seq)
# cut_pos = [0]
# for match in x:
#     pos = match.start() + 3
#     cut_pos.append(pos)
# y = re.finditer(r"GC[AG][AT]TG",seq)
# for match in y:
#     pos = match.start() + 4
#     cut_pos.append(pos)
# cut_pos.append(len(seq))
# cut_pos.sort()
# for i in range(1,len(cut_pos)):
#     frag = seq[cut_pos[i-1]:cut_pos[i]]
#     print("Fragment length is",len(frag))

# DNA Transation
# gencode = {
# 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
# 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
# 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
# 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
# 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
# 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
# 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
# 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
# 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
# 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
# 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
# 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
# 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
# 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
# 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
# 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
#
# def translate_dna(dna):
#     protein = ""
#     for x in range(0,len(dna)-2,3):
#         start_pos = x
#         stop_pos = x + 3
#         trineucleotide = dna[start_pos:stop_pos].upper()
#         decode = gencode.get(trineucleotide,"X")
#         protein += decode
#     return protein.replace("_","")
#
# assert(translate_dna("ATGTTCGGT")) == "MFG"
# assert(translate_dna("ATCGATCGATCGTTGCTTATCGATCAG")) == "IDRSLLIDQ"
# assert(translate_dna("actgatcgtagcttgcttacgtatcgtat")) == "TDRSLLTYR"
# assert(translate_dna("ACGATCGATCGTNACGTACGATCGTACTCG")) == "TIDRXVRSYS"

# import os
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\100_199")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\200_299")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\300_399")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\400_499")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\500_599")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\600_699")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\700_799")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\800_899")
# os.mkdir(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\Above_900")
#
# for x in os.listdir(r"C:\Users\KHAN GADGET\Desktop\Onedrive\p4b_exercises\exercises and examples\working_with_the_filesystem\exercises\Dna"):
#     dna_file = open(r"C:/Users/KHAN GADGET/Desktop/Onedrive/p4b_exercises/exercises and examples/working_with_the_filesystem/exercises/Dna/" + x)
#     for seq in dna_file:
#         if len(seq) >= 100 and len(seq) <= 199:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\100_199\100_199.txt","a")
#             file.write(seq)
#         elif len(seq) >= 200 and len(seq) <= 299:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\200_299\200_299.txt","a")
#             file.write(seq)
#         elif len(seq) >= 300 and len(seq) <= 399:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\300_399\300_399.txt", "a")
#             file.write(seq)
#         elif len(seq) >= 400 and len(seq) <= 499:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\400_499\400_499.txt","a")
#             file.write(seq)
#         elif len(seq) >= 500 and len(seq) <= 599:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\500_599\500_599.txt","a")
#             file.write(seq)
#         elif len(seq) >= 600 and len(seq) <= 699:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\600_699\600_699.txt","a")
#             file.write(seq)
#         elif len(seq) >= 700 and len(seq) <= 799:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\700_799\700_799.txt","a")
#             file.write(seq)
#         elif len(seq) >= 800 and len(seq) <= 899:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\800_899\800_899.txt","a")
#             file.write(seq)
#         else:
#             file = open(r"C:\Users\KHAN GADGET\PycharmProjects\P4B\Above_900\Above_900.txt","a")
#             file.write(seq)


# dna = "ACTGTA"
kmer_length = 5
cutoff_number = 5
dict = {}
seqs = ""
import os
for x in os.listdir(r"E:\Documents\Python\p4b_exercises\exercises and examples\working_with_the_filesystem\exercises\Dna"):
    dna_file = open(r"E:/Documents/Python/p4b_exercises/exercises and examples/working_with_the_filesystem/exercises/Dna/" + x)
    for dna in dna_file:
        seq = dna.rstrip("\n")
        seqs = seqs + seq
print(seqs)

for start in range(0,len(seqs)-(kmer_length-1)):
    neucleotide = seqs[start:start+ kmer_length]
    count = seqs.count(neucleotide)
    dict[neucleotide] = count
print(dict)



