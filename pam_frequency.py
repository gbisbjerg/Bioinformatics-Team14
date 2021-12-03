import pandas as pd
from Trie import *


translation = {
    "A" : ["A"], 
    "C" : ["C"],
    "G" : ["G"],
    "T" : ["T"],
    "U" : ["U"],
    "R" : ["G", "A"],
    "Y" : ["C" , "T"],
    "K" : ["G", "T"], 
    "M" : ["A", "C"], 
    "S" : ["G", "C"], 
    "W" : ["A", "T"], 
    "B" : ["G", "T", "C"], 
    "D" : ["G", "A", "T"], 
    "H" : ["A", "C", "T"],  
    "V" : ["G", "C", "A"], 
    "N" : ["A", "G", "C", "T"]
}


def process_pam_csv(filename): 
    df = pd.read_csv (filename)
    return df 


def create_possibilities(filename): 
    df = process_pam_csv(filename)
    list_of_lists =  []
    all_lists = []
    for ind in df.index:
        string_ = df['sequence'][ind]
        list_ = [char for char in string_]

        for c in list_: 
            options = translation[c]
            if len(list_of_lists) == 0: 
                for o in options:
                    list_of_lists.append(o)
            else: 
                list_of_lists2 = []
                j = 0 
                for j in range(0, len(list_of_lists)): 
                    for o in options:
                        new_item = list_of_lists[j] + o
                        list_of_lists2.append(new_item)
                    j+=1 
                list_of_lists = list_of_lists2.copy()

        all_lists.append(list_of_lists)
        list_of_lists = []

    df['possibilities'] = all_lists
    return df 


def build_pam_trie(df):
   trie = Trie()
   for i in df.index:
      pam = df['sequence'][i]
      sequences = df['possibilities'][i]
      for sequence in sequences:
         trie.insert(pam, sequence)
   return trie


def process_genome_txt(file_name):
   genome = ""
   f = open(file_name)
   lines = f.readlines()
   for line in lines:
      genome += line[ : -1]
      f.close()
   return genome


def generate_pam_positions(trie):
   pam_positions = {}
   genome = process_genome_txt("sample1_mutation.txt")
   for i in range(len(genome)):
      trie.search(genome, i, pam_positions)
   return pam_positions


def build_pam_map(df):
   pam_map = {}
   max_len = 0
   for i in df.index:
      pam = df['sequence'][i]
      sequences = df['possibilities'][i]
      for sequence in sequences:
         max_len = max(max_len, len(sequence))
         if sequence in pam_map:
            pam_map[sequence].append(pam)
         else:
            pam_map[sequence] = [pam]
   return pam_map, max_len


def generate_pam_positions_naive(pam_map, max_len):
   pam_positions = {}
   genome = process_genome_txt("sample1_mutation.txt")
   for start in range(len(genome)):
      for offset in range(1, min(max_len, len(genome) - start)):
         sequence = genome[start : start + offset]
         if sequence in pam_map:
            pams = pam_map[sequence]
            for pam in pams:
               if pam in pam_positions:
                  pam_positions[pam].append(start)
               else:
                  pam_positions[pam] = [start]
   return pam_positions


def main():
   df = create_possibilities("pam_raw.csv")
   pam_trie = build_pam_trie(df)
   pam_positions = generate_pam_positions(pam_trie)
   print(pam_positions)
   pam_map, max_len = build_pam_map(df)
   pam_positions = generate_pam_positions_naive(pam_map, max_len)
   print(pam_positions)


if __name__ == '__main__': 
    main()
