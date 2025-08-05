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

class Strand:
    def __init__(self, name:str, content:str):
        self.name:str = name
        self.content:str = content

class DnaMolecule:
    """
    Класс молекулы ДНК
    Входные:
    Параметр dna_length: (целочисленное) длина ДНК состоит из
    длины кодона в нуклеотидах (для человека 3)
    умноженной на число 61 кодонов (для человека 61),
    которыми кодируются все аминокислоты (для человека 20).
    Приватные:
    Параметр dna_5_to_3_strand: (строковое) цепь ДНК (направление 5'-3').
    Параметр dna_3_to_5_strand: (строковое) цепь ДНК (направление 3'-5').
    """
    def __init__(self, epigenetic_marks=None):
        self.codon_length = 3
        self.amount_of_codon = 61
        self.size = 1
        self.dna_length = self.codon_length * self.amount_of_codon * self.size
        self.epigenetic_marks = epigenetic_marks or {}

        # генерирую цепь ДНК с направлением 5'-3'
        self.dna_5_to_3_strand = Strand("5'-3' direction of DNA",self.__set_dna_5_to_3_strand(self.dna_length))
        # генерирую цепь ДНК с направлением 3'-5'
        self.dna_3_to_5_strand = Strand("3'-5' direction of DNA", self.__set_dna_3_to_5_strand(self.dna_5_to_3_strand.content))

    @staticmethod
    def __set_dna_5_to_3_strand(dna_length:int):
        dna_nucleotides = ['A', 'T', 'G', 'C']
        return ''.join(random.choices(dna_nucleotides, k=dna_length))

    @staticmethod
    def __set_dna_3_to_5_strand(dna_sequence:str):
        # Словари для быстрых преобразований
        dna_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(dna_complement[nt] for nt in dna_sequence)

    def __set_genetic_mutations(self):
        pass

class RnaMolecule:
    """
    Класс молекулы РНК
    Входной:
    Параметр экземпляр класса DnaMolecule, который является программной репрезентацией одной молекулы ДНК
    """
    def __init__(self, mol_dna:DnaMolecule):
        """ Транскрипция ДНК (5'->3') в РНК (5'->3').
        Тут реализуемся очень упрощённая работа РНК-полимеразы II, благодаря которой получается пре-мРНК:
            1. Читаем ДНК в обратном порядке (от 3' к 5').
            2. Одновременно заменяем нуклеотиды.
            3. Объединяем в строку РНК 5'->3'.
        За кадром работы конструктора остаётся моделирование сборки ферментного и белкового комплекса на регулярных
        участках ДНК и начала процесса транскрипции (инициация), удлинение полинуклеотидной цепи РНК (элонгация),
        окончание транскрипции (терминация).
        """
        # Читаем ДНК с конца (3'->5') и заменяем нуклеотиды.
        self.mol_dna = mol_dna
        dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        self.strand = Strand("", ''.join(dna_to_rna[nt] for nt in reversed(self.mol_dna.dna_3_to_5_strand.content)))

        # Результатом работы конструктора до вызова метода splicing() пре-мРНК.

    def splicing(self, sequence:Strand)->Strand:
        pass

class Protein:
    def __init__(self, mol_rna:RnaMolecule):
        pass