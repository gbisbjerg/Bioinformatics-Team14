import pandas as pd
from PamTrie import *
import time

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
    print(df)
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


def generate_pam_positions_trie(trie, genome):
   pam_positions = {}
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


def generate_pam_positions_hashmap(pam_map, max_len, genome):
   pam_positions = {}
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


def pick_best_pam(genome, pam_positions, index_of_mutation):

   smallest_num_pams = len(genome)
   min_ = len(genome)
   closest_pam = ""

   for pam in pam_positions.keys(): 
      if len(pam_positions[pam]) < smallest_num_pams: 
         smallest_num_pams =  len(pam_positions[pam])
         positions_list = pam_positions[pam]

         for p in positions_list:
            distance = abs(index_of_mutation[0] - p)
            if distance < min_: 
               min_ =  distance
               closest_pam = pam 

   if pam == "": 
      return pam_positions.keys[0]

   print("\nIdeal PAM= {}".format(closest_pam))
   print("\tOccurs {} times in the sequence".format(smallest_num_pams))
   print("\tIs {} neucleotide(s) away from the mutation".format(min_))

   return closest_pam


def MutationLocation(source_file, reference_file):
    source_file = open(source_file, "r")
    reference_file = open(reference_file, "r")
    source = source_file.read()
    reference = reference_file.read()

    mutation_index = []
    for i, nucleotide in enumerate(source):
        if (reference[i] != nucleotide):
            mutation_index.append(i)

    return mutation_index

def main():
   df = create_possibilities("pam_raw.csv")
   mutated_file = "sample_2k_mutation.txt"
   regular_file = "sample_2k.txt"
   genome = process_genome_txt(mutated_file)
   mutated_genome = process_genome_txt(regular_file)
   index_of_mutation = MutationLocation(mutated_file, regular_file)
   start1 = time.time()
   pam_trie = build_pam_trie(df)
   start2 = time.time()
   pam_positions = generate_pam_positions_trie(pam_trie, genome)
   end1 = time.time()
   start3 = time.time()
   pam_map, max_len = build_pam_map(df)
   start4 = time.time()
   pam_positions = generate_pam_positions_hashmap(pam_map, max_len, genome)
   end2 = time.time()
   print("Hashmap (no construction) " + str(round(end2 - start4, 9)))
   print("Trie (no construction) " + str(round(end1 - start2, 9)))
   print("Hashmap (with construction) " + str(round(end2 - start3, 9)))
   print("Trie (with construction) " + str(round(end1 - start1, 9)))

   best_pam = pick_best_pam(genome, pam_positions, index_of_mutation)
   pam_info = df.loc[df['sequence'] == best_pam]
   print("Additional PAM Info:\n{}".format(pam_info))



if __name__ == '__main__': 
    main()
