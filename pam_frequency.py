
import pandas as pd 
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
        string_ = df['squence'][ind]
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

def main():
    df = create_possibilities("pam_raw.csv")
    print("df output", df)
            

if __name__ == '__main__': 
    main()













