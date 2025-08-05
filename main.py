from MolecularBiology import *

def main():
    single_molecule_dna = DnaMolecule()
    single_molecule_rna = RnaMolecule(single_molecule_dna)
    print(single_molecule_rna.strand.name)
    print(single_molecule_rna.strand.content)

if __name__ == "__main__":
    main()