import sys

'''
The function takes in two files and returns an array of all of the point mutations
No allignment occures so files are expected to be of the same length and starting point
'''
def MutationLocation(source_file, reference_file):
    source = source_file.read()
    reference = reference_file.read()

    mutation_index = []
    for i, nucleotide in enumerate(source):
        if (reference[i] != nucleotide):
            mutation_index.append(i)

    return mutation_index

def main():
    source_file = open(sys.argv[1], "r")
    reference_file = open(sys.argv[2], "r")

    mutation_index = MutationLocation(source_file, reference_file)

    source_file.close()
    reference_file.close()

if __name__ == "__main__":
    main()
