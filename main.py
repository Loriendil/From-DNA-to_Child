from MolecularBiology import *

def main():
    single_molecule_dna = DNA_Molecule()
    single_molecule_rna = RNA_Molecule(single_molecule_dna)
    print("RNA:")
    print(single_molecule_rna.strand)

if __name__ == "__main__":
    main()