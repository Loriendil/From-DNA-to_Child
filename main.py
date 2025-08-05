from MolecularBiology import *

def main():
    single_molecule_dna = DnaMolecule()
    print(single_molecule_dna.dna_3_to_5_strand.name)
    print(single_molecule_dna.dna_3_to_5_strand.content)
    print(single_molecule_dna.dna_5_to_3_strand.name)
    print(single_molecule_dna.dna_5_to_3_strand.content)
    single_molecule_rna = RnaMolecule(single_molecule_dna)
    print(single_molecule_rna.strand.name)
    print(single_molecule_rna.strand.content)

if __name__ == "__main__":
    main()