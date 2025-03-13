#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 13:36:00 2025

@author: expeditoolimi
"""
list_1 = [1,2,3,4,5] # a list of integers
print('List of Integers:', list_1)

list_2 = ['A', 'T', 'G', 'C'] # a list of strings
print('List of Strings:', list_2)

list_3 = [[1,2,3], ['A', 'T', 'G', 'C']] # a list of lists
print('List of Lists:', list_3)

new_list=list_1 + list_2
print('new_list:')

list = [1,2,3]
print('Here is a list:', list)
list.append(4)
print('Here is a list with a new item:', list)
list.remove(2)
print('Here is a list with the second item removed:', list)

list = list_2
print('Here is a list:', list)
list.append("B")
print('Here is a list with a new item:', list)
list.remove("C")
print('Here is a list with the second item removed:', list)


lst = [10, 20, 30, 40, 50] 
print(lst[:]) 
print(lst[:3])

list = ['a', 'b', 'c', 'd', 'e']
print(list[3])# retrieve the item in '3' position
print(list[1:4]) # retrieve items in positions 1-3
print(list[0:4:2]) # retrieve items in positions 0-3, skipping every other item
print(list[:4]) # retreive items 0-3
print(list[4])# retrieve the item in '3' position

dna_sequences = ['ACTGGT', 'TGGCAG', 'CGTTAA']

for sequence in dna_sequences:
  print(sequence)


dna_sequences = ['ACTGGT', 'TGGCAG', 'CGTTAA']
gc_content = [(seq.count('G') + seq.count('C')) / len(seq) for seq in dna_sequences]
print('The GC content of each sequence is:', gc_content)

#more complex example
seque=['GTGTCAGCAGCCGCGGTAAGACGGGGGGGGCAAGTGTTCTTCGGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTGAAAGCCGCCAAAAACTGGCGGAATGCTCTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGTAGGGGTGAAATCCGGAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGAGGTGCGAAAGCATGGGGAGCGAACAGGATTAGAAACCCGAGTAGTCC','GTGCCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTGGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGAAACCCCAGTAGTCC']

for item in seque:
  print(seque)

gc_content = [(seq.count('G') + seq.count('C')) / len(seq) for seq in seque]
print('The GC content of each sequence is:', gc_content)

#working with doctionaries 
patient_info = {'Patient 1':{'Height':170, 'Weight':60}, 'Patient 2':{'Height':160, 'Weight':50}}
print(patient_info) # prints contents of dictionary

dna_to_rna = {'A':'A', 'T':'U', 'C':'C', 'G':'G'}
print(dna_to_rna['T']) # retrieve value associated with the key 'T'


dna_to_rna = {'A':'A', 'T':'U', 'C':'C', 'G':'G'}
print(dna_to_rna)

dna_to_rna.pop('T') # remove value by key
print(dna_to_rna)


dna_to_rna['T'] = 'U' # add key-value pair to dictionary
print(dna_to_rna)


gene_annotations = [
    ('GeneA',1000,2000,'protein_coding'),
    ('GeneB',3000,4000,'non_coding'),
    ('GeneC',5000,6000,'protein_coding')]

#gene_dictionary = [{'gene_name':name, 'start_position':start, 'end_position':end, 'gene_type',gene_type} for name, start, end, gene_type in gene_annotations]



#writing conditionals

x=15
if x>10: 
  print('x is greater than 10')


x=5
if x>10: 
  print('x is greater than 10')
else:
  print ('x is less than or equal to 10')


#Nested conditonal
x=5
y=6
if x==y:
  print('x and y are equal')
else:
  if x<y:
    print('x is less than y') # in this case, this second print statement is executed.
  else:
    print('x is greater than y')



#Control Flow With Loops
#In Python, loops enable you to repeat a block of code multiple times, making it easier to process data or perform repetitive tasks. There are two types of loops you should be familiar with: for loops and while loops. A for loop iterates over a sequence (such as a list, tuple, string, or range) and executes a code block for each element in that sequence. The syntax of a for loop is as follows:
# List of nucleotide sequences
dna_sequences = ['ACTGGT', 'TGGCAG', 'CGTTAA']

# Loop through each sequence in the list
for sequence in dna_sequences:
    length=len(sequence)
    # Print the length of the current sequence
    print(f"The length of the sequence '{sequence}' is {length} bases")



#the word sequence above is nothing speacial
# List of nucleotide sequences
dna_sequences = ['ACTGGT', 'TGGCAG', 'CGTTAA']

# Loop through each sequence in the list
for file in dna_sequences:
    length=len(file)
    # Print the length of the current sequence
    print(f"The length of the sequence '{file}' is {length} bases")




def print_sequences(sequences):
  """
  Prints the length of each nucleotide sequence in a list of sequences.

  Args:
    sequences: A list of strings, where each string represents a nucleotide sequence.
  """
  for sequence in sequences:
    length = len(sequence)
    print(f"Length of sequence: {length}")

# Example usage:
dna_sequences = ['ACTGGT', 'TGGCAG', 'CGTTAA']
print_sequences(dna_sequences)




#examople from LLM
def print_sequence_lengths(sequences):
  """
  Prints the length of each nucleotide sequence in a list of sequences.

  Args:
    sequences: A list of strings, where each string represents a nucleotide sequence.
  """
  for file in sequences:
    length = len(file)
    print(f"Length of sequence: {length}")

# Example usage:
nucleotide_sequences = ["ATCG", "GGCC", "AATT", "CGTACGT"]
print_sequence_lengths(nucleotide_sequences)



#code from deepseek
# List of nucleotide sequences
sequences = ["ATCG", "GATTACA", "TTAGGG", "ACGTACGT", "CGTA"]

# Loop through each sequence in the list
for file in sequences:
    # Print the length of the current sequence
    print(f"The length of the sequence '{sequence}' is {len(sequence)} nucleotides")




# Example usage:


def print_sequence_lengths(sequences):
  """
  Prints the length of each nucleotide sequence in a list of sequences.

  Args:
    sequences: A list of strings, where each string represents a nucleotide sequence.
  """
  for sequence in sequences:
    length = len(sequence)
    print(f"Length of sequence: {length}")

sequences = ["ATCG", "GATTACA", "TTAGGG", "ACGTACGT", "CGTA"]
print_sequence_lengths(sequence)





#WHILE LOOP command
#Whereas for loops iterates over a sequence and executes a code block for each element in that sequence, while loops repeatedly execute a block of code as long as a specified condition is true. The syntax of a while loop is as follows:
while condition:
# code to execute as long as the condition is True\
#As previously mentioned, for loops are commonly used in bioinformatics for iterating over sequences or datasets. While loops, on the other hand, are useful in scenarios where the number of iterations is not known beforehand or when iterating until a specific condition is met.
#In bioinformatics, a while loop might be used to simulate a scenario where a biological process continues until a certain condition is met. Let's consider a simplified example where we simulate a mutation in a DNA sequence until a specific mutation pattern is achieved:

import random
dna_seq = 'TGGATCCATGCA'
target_mutation = 'TGGATCCATGCT'
mutations = 0


while dna_seq != target_mutation:
  position = random.randint(0, len(dna_seq)-1)
  mutated_base = random.choice('ATCG'.replace(dna_seq[position],''))
  dna_seq = dna_seq[:position]+mutated_base+dna_seq[position+1:]
  mutations +=1

print(f'Target mutation achieved after {mutations} mutations.')
print('Final DNA sequence:', dna_seq)


import random 
dna_seq="GTGTCAGCAGCCGCGGTAAGACGGG"
target_mutation="GTGTCAGCAGCCGCGGTAAGTCGGG"
mutations = 0

while dna_seq != target_mutation:
  position = random.randint(0, len(dna_seq)-1)
  mutated_base = random.choice('ATCG'.replace(dna_seq[position],''))
  dna_seq = dna_seq[:position]+mutated_base+dna_seq[position+1:]
  mutations +=1

print(f'Target mutation achieved after {mutations} mutations.')
print('Final DNA sequence:', dna_seq)
#


dna_seq = 'ACGTGTCAGTGGGAC

def dna_to_rna(any_dna_sequence):
  rna_sequence = dna_sequence.repace('T':'U')
  return rna_sequence

rna_seq = dna_to_rna(dna_seq)
print(rna_seq) # the code will return the rna sequence ACGUGUCAGUGGGAC