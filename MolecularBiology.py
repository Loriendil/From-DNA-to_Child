import random

# Генетический код для трансляции
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

class DNA_Molecule:
    """
    Класс молекулы ДНК
    Входные:
    Параметр dna_length: (целочисленное) длина ДНК состоит из
    длины кодона в нуклеотидах (для человека 3)
    умноженной на число 61 кодонов (для человека 61),
    которыми кодируются все аминокислоты (для человека 20).
    Приватные:
    Параметр dna_matrix: (строковое) матричная цепь ДНК (направление 5'-3').
    Параметр dna_complement: (строковое) комплементарная цепь ДНК (направление 3'-5').
    """
    def __init__(self, epigenetic_marks=None):
        self.codon_length = 3
        self.amount_of_codon = 61
        self.size = 1
        self.dna_length = self.codon_length * self.amount_of_codon * self.size
        self.epigenetic_marks = epigenetic_marks or {}
        # генерирую матричную цепь ДНК
        self.dna_matrix_strand = self.__set_dna_matrix_strand(self.dna_length)
        # генерирую комплементарную цепь ДНК
        self.dna_complement_strand = self.__set_dna_complement_strand(self.dna_matrix_strand)

    @staticmethod
    def __set_dna_matrix_strand(dna_length:int):
        dna_nucleotides = ['A', 'T', 'G', 'C']
        return ''.join(random.choices(dna_nucleotides, k=dna_length))

    @staticmethod
    def __set_dna_complement_strand(dna_matrix:str):
        # Словари для быстрых преобразований
        dna_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(dna_complement[nt] for nt in dna_matrix)

    def __set_genetic_mutations(self):
        pass

class RNA_Molecule:
    """
    Класс молекулы РНК
    Входной:
    Параметр complimentary_strand комплементарная цепь ДНК
    """
    def __init__(self, mol_dna:DNA_Molecule):
        """ Транскрибирует ДНК (5'->3') в РНК (5'->3').
        Из сообщения нейросети:
            Читаем ДНК в обратном порядке (от 3' к 5').
            Одновременно заменяем нуклеотиды.
            Объединяем в строку РНК 5'->3'.
        """
        # Читаем ДНК с конца (3'->5') и заменяем нуклеотиды.
        self.mol_dna = mol_dna
        dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        self.strand = ''.join(dna_to_rna[nt] for nt in reversed(mol_dna.dna_complement_strand))

class Protein:
    def __init__(self, mol_rna:RNA_Molecule):
        pass