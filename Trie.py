class TrieNode:

   def __init__(self):
      self.keys = set()
      self.children = {"A": None, "G": None, "C": None, "T": None}


class Trie:

   def __init__(self):
      self.root = TrieNode()

   def insert(self, key, value):
      curr = self.root
      for i in range(len(value)):
         char = value[i]
         if curr.children[char] is None:
            curr.children[char] = TrieNode()
         next = curr.children[char]
         if i == len(value) - 1:
            next.keys.add(key)
         curr = next

   def search(self, genome, index, positions):
      curr = self.root
      i = index
      while curr is not None and i < len(genome):
         char = genome[i]
         for key in curr.keys:
            if key in positions:
               positions[key].append(index)
            else:
               positions[key] = [index]
         curr = curr.children[char]
         i += 1
