import hashlib, time, uuid, os
from collections import defaultdict
from typing import Tuple, Optional, Iterable, List
import re, random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Strand:
    def __init__(self, name: str, content: str):
        self.name: str = name
        self.content: str = content


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

    def __init__(self, epigenetic_marks=None,
                 length: int | None = None,
                 run_id: str | None = None,
                 gc: float = 0.41
                 ):
        self.codon_length = 3
        self.amount_of_codon = 61
        self.size = 1

        self.dna_length = length or (self.codon_length * self.amount_of_codon * self.size)
        self.epigenetic_marks = epigenetic_marks or {}

        # 1) Получаем семя. Если хотите воспроизводимость — передайте run_id (строка).
        self.run_id = run_id or str(uuid.uuid4())
        # seed = self.derive_seed(deterministic_key=self.run_id)
        #
        # # 2) Генерируем 5'→3' с заданным GC (чуть ближе к «человеческому» профилю).
        # stream = self.byte_stream(seed)
        # seq5 = self.bytes_to_dna(stream, self.dna_length, self.weights_for_gc(gc))
        #
        # # 3) Вторая нить — обратный комплемент (а не просто комплемент без разворота).
        # self.dna_5_to_3_strand = Strand("5'-3' direction of DNA", seq5)
        # self.dna_3_to_5_strand = Strand("3'-5' direction of DNA", self.reverse_complement(seq5))

        # 1. Загружаем FASTA-файл с реальной последовательностью нуклеотидов в память
        in_fa = r'E:\local_repository\Python-repository\From-DNA-to_Child\split\chr1.fa'
        real_dna: str = ''
        with open(in_fa) as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                seq = Seq(str(record.seq).upper())
                real_dna = str(seq)

        # 2. Получаем семя. Если хотите воспроизводимость — передайте run_id (строка).
        self.run_id = run_id or str(uuid.uuid4())

        # 3. Строим модель
        order = 2  # Используем цепь 2-го порядка (зависимость от 2 предыдущих нуклеотидов)
        model = self.build_markov_chain(real_dna, order)

        # 4. Выбираем начальный контекст (например, случайный из модели)
        seed_context = random.choice(list(model.keys()))

        # 5. Генерируем правдоподобную последовательность
        generator = self.markov_dna_generator(seed_context, model, length=1000)
        plausible_dna = ''.join(generator)

        # 6. Вторая нить — обратный комплемент (а не просто комплемент без разворота).
        self.dna_5_to_3_strand = Strand("5'-3' direction of DNA", plausible_dna)
        self.dna_3_to_5_strand = Strand("3'-5' direction of DNA", self.reverse_complement(plausible_dna))

    @staticmethod
    def build_markov_chain(real_dna_sequence: str, order: int = 1) -> dict:
        """
        Строит модель Марковской цепи на основе реальной последовательности ДНК.
        Возвращает словарь, где ключ - это кортеж из n предыдущих символов (контекст),
        а значение - это список нуклеотидов, которые следуют за этим контекстом.
        """
        model = defaultdict(list)
        for i in range(order, len(real_dna_sequence)):
            context = tuple(real_dna_sequence[i - order:i])
            next_char = real_dna_sequence[i]
            model[context].append(next_char)
        return model

    @staticmethod
    def markov_dna_generator(seed_context: tuple, markov_model: dict, length: int) -> Iterable[str]:
        """
        Генератор последовательности ДНК на основе модели Марковской цепи.
        """
        context = seed_context
        for _ in range(length):
            # Получаем возможные следующие символы для текущего контекста
            next_options = markov_model.get(context, [])
            if not next_options:
                # Если контекст неизвестен (не был в обучающей выборке), выбираем случайно
                next_char = random.choice('ATGC')
            else:
                # Выбираем случайный следующий символ из возможных
                next_char = random.choice(next_options)
            yield next_char
            # Обновляем контекст: убираем первый символ и добавляем новый
            context = context[1:] + (next_char,)

    @staticmethod
    def derive_seed(seed: bytes | None = None, *, deterministic_key: str | None = None) -> bytes:
        if seed is not None:
            return hashlib.sha256(seed).digest()
        if deterministic_key is not None:
            return hashlib.sha256(deterministic_key.encode()).digest()
        payload = f"{time.time_ns()}|{os.getpid()}|{uuid.getnode()}|{uuid.uuid4()}".encode()
        return hashlib.sha256(payload).digest()

    @staticmethod
    def byte_stream(seed: bytes) -> Iterable[int]:
        """
        SHA-256 в режиме счётчика: на каждом шаге хешируем seed||counter и отдаём байты.
        """
        counter = 0
        while True:
            block = hashlib.sha256(seed + counter.to_bytes(8, 'big')).digest()
            counter += 1
            for b in block:
                yield b

    @staticmethod
    def weights_for_gc(gc: float = 0.41) -> tuple[float, float, float, float]:
        """
        Доли A,T,G,C при заданном GC-содержании (для человека ~0.41).
        """
        gc = min(max(gc, 0.0), 1.0)
        at = 1.0 - gc
        return (at / 2, at / 2, gc / 2, gc / 2)  # A,T,G,C

    @staticmethod
    def bytes_to_dna(biter: Iterable[int], length: int,
                     weights: tuple[float, float, float, float]) -> str:
        """
        Преобразует поток байт в ДНК с заданными частотами A/T/G/C.
        Без 'random': просто пороговое сравнение в [0..255].
        """
        a, t, g, c = weights
        tA = int(round(a * 256))
        tT = tA + int(round(t * 256))
        tG = tT + int(round(g * 256))
        out = []
        it = iter(biter)
        for _ in range(length):
            x = next(it)
            if x < tA:
                out.append('A')
            elif x < tT:
                out.append('T')
            elif x < tG:
                out.append('G')
            else:
                out.append('C')
        return ''.join(out)

    @staticmethod
    def reverse_complement(dna: str) -> str:
        return dna.translate(str.maketrans('ATGC', 'TACG'))[::-1]

    @staticmethod
    def __set_genetic_mutations(self):
        pass


class RnaMolecule:
    """
    Класс молекулы РНК
    Входной:
    Параметр экземпляр класса DnaMolecule, который является программной репрезентацией одной молекулы ДНК
    """

    def __init__(self, mol_dna: DnaMolecule, gene_source: str = "3'-5' direction of DNA"):
        pass
        # """ Транскрипция ДНК (5'->3') в РНК (5'->3').
        # Тут реализуемся очень упрощённая работа РНК-полимеразы II, благодаря которой получается пре-мРНК:
        #     1. Читаем ДНК в обратном порядке (от 3' к 5').
        #     2. Одновременно заменяем нуклеотиды.
        #     3. Объединяем в строку РНК 5'->3'.
        # За кадром работы конструктора остаётся моделирование сборки ферментного и белкового комплекса на регулярных
        # участках ДНК и начала процесса транскрипции (инициация), удлинение полинуклеотидной цепи РНК (элонгация),
        # окончание транскрипции (терминация).
        # """
        # # Читаем ДНК с конца (3'->5') и заменяем нуклеотиды.
        # self.name: str = 'mature mRNA'
        # self.mol_dna = mol_dna
        # dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        #
        # if gene_source == "3'-5' direction of DNA":
        #     # шаблонная (antisense) нить -> комплементарная пре-мРНК
        #     first_rna_transcript = Strand(
        #         "from template (3'-5' DNA)",
        #         ''.join(dna_to_rna[nt] for nt in self.mol_dna.dna_3_to_5_strand.content)
        #     )
        # else:
        #     # Для кодирующей цепи: реверс + комплементарная замена
        #     reversed_strand = self.mol_dna.dna_5_to_3_strand.content[::-1]
        #     first_rna_transcript = Strand("from coding (5'-3' DNA)",
        #                                   ''.join(dna_to_rna[nt] for nt in reversed_strand))
        #
        # # Кэпирование, полиаденилирование пропускаем, так как химические процессы в чистом виде нет цели симулировать.
        # # Результатом работы конструктора до вызова метода splicing() пре-мРНК.
        # pre_mrna: Strand = self.validate_rna_transcript(first_rna_transcript)
        # self.strand, self.log = self.splicing(pre_mrna)

    @staticmethod
    def validate_rna_transcript(raw: Strand) -> Strand:
        content = raw.content.upper()
        # Удаляем недопустимые символы
        clean_content = re.sub(r"[^AUGC]", "", content)
        return Strand(raw.name, clean_content)

    @staticmethod
    def is_py_tract(seq: str, min_len: int = 12, frac: float = 0.8) -> bool:
        if len(seq) < min_len:
            return False
        cu = sum(1 for c in seq if c in "CU")
        return (cu / len(seq)) >= frac

    @staticmethod
    def find_branch_a(intron_seq: str, min_off: int = 18, max_off: int = 40) -> Optional[int]:
        """
        Ищем A в мотиве ~YNYURAY (7-мер) в окне [min_off..max_off] от 3'-конца интрона.
        Возвращаем индекс 'A' внутри интрона или None.
        """
        L = len(intron_seq)
        # регион для сканирования (относительно 3'-конца)
        start = max(0, L - max_off - 7)  # 7-мер
        end = max(0, L - min_off)  # end-exclusive

        def is_Y(b):
            return b in "CU"

        def is_R(b):
            return b in "AG"

        best = None
        for i in range(start, end - 6):
            w = intron_seq[i:i + 7]
            if (is_Y(w[0]) and w[1] in "ACGU" and is_Y(w[2]) and
                    w[3] == "U" and is_R(w[4]) and w[5] == "A" and is_Y(w[6])):
                best = i + 5  # позиция A
                break
        return best

    @staticmethod
    def splicing(sequence: Strand,
                 min_intron_len: int = 50,  # чуть поднять минимум
                 max_intron_len: int = 500_000,
                 branch_min_offset: int = 18,
                 branch_max_offset: int = 40,
                 py_tract_min_len: int = 12,
                 py_tract_threshold: float = 0.8,
                 allow_noncanonical: bool = False) -> Tuple[str, list]:
        """
        # --- Параметры биологической эвристики (можно настроить) ---
        min_intron_len = 30         # минимальная длина интрона (в нуклеотидах) — меньше обычно редкость
        max_intron_len = 500000     # максимум (защита от бесконечных поисков)
        branch_min_offset = 18      # минимальное расстояние от branch-A до 3'AG
        branch_max_offset = 40      # максимальное расстояние от branch-A до 3'AG
        py_tract_min_len = 5        # минимальная длина рассматриваемого Py-tract окна
        py_tract_threshold = 0.6    # доля C/U в окне для пометки как Py-tract
        allow_noncanonical = True   # разрешить проверку GC-AG и AT-AC (U12) как опцию
        """
        """
            Возвращает (mature_seq, log_list).
            log_list: список кортежей (intron_start, intron_end, donor, acceptor, score, comment)
              где intron_start,intron_end — координаты в исходной последовательности (0-based, end exclusive),
              donor — 5' два нуклеотида интрона, acceptor — 3' два нуклеотида интрона,
              score — оценка доверия 0..1, comment — текстовое примечание.
            """
        # Канонические пары для сайтов сплайсинга
        canon_pairs = [("GU", "AG")]
        if allow_noncanonical:
            canon_pairs += [("GC", "AG"), ("AU", "AC")]

        seq = sequence.content
        mature_parts = []  # Собранные экзоны
        log = []  # Лог интронов
        pos = 0  # Текущая позиция в последовательности

        while pos < len(seq):
            found_intron = False
            best_candidate = None

            # Поиск донорного сайта
            for donor, acceptor in canon_pairs:
                # Проверяем только если достаточно места для донора
                if pos + 2 > len(seq):
                    continue

                # Проверяем совпадение с донорным сайтом
                if seq[pos:pos + 2] == donor:
                    # Корректный расчет границ поиска акцептора
                    acceptor_start = pos + min_intron_len
                    acceptor_end = min(pos + max_intron_len + 1, len(seq))

                    # Поиск акцептора в допустимых пределах
                    j = acceptor_start
                    while j < acceptor_end - 1:  # -1 чтобы осталось место для двух нуклеотидов
                        if seq[j:j + 2] == acceptor:
                            intron_start = pos
                            intron_end = j + 2
                            intron_len = intron_end - intron_start

                            # Проверка минимальной длины интрона
                            if intron_len < min_intron_len:
                                j += 1
                                continue

                            intron_seq = seq[intron_start:intron_end]

                            # Проверка биологических критериев
                            branch_idx = RnaMolecule.find_branch_a(
                                intron_seq, branch_min_offset, branch_max_offset
                            )
                            py_region = seq[max(j - 20, intron_start):j]
                            py_ok = RnaMolecule.is_py_tract(
                                py_region, py_tract_min_len, py_tract_threshold
                            )

                            # Расчет confidence score
                            score = 0.6 if (donor, acceptor) == ("GU", "AG") else 0.3
                            if branch_idx is not None:
                                score += 0.25
                            if py_ok:
                                score += 0.2
                            score = min(score, 1.0)

                            # Выбор лучшего кандидата
                            if score >= 0.8 and (best_candidate is None or score > best_candidate[0]):
                                best_candidate = (
                                    score, intron_start, intron_end,
                                    donor, acceptor,
                                    f"branch={'ok' if branch_idx is not None else 'no'};py_ok={py_ok}"
                                )

                        j += 1

            # Обработка найденного интрона
            if best_candidate:
                score, start, end, d, a, comment = best_candidate

                # Сохраняем экзон перед интроном
                if start > pos:
                    mature_parts.append(seq[pos:start])

                # Логируем интрон
                log.append((start, end, d, a, round(score, 3), comment))

                # Перемещаем позицию за интрон
                pos = end
                found_intron = True
            else:
                # Не найден интрон - двигаемся дальше
                pos += 1

        # Добавляем последний экзон
        if mature_parts:
            # Добавляем остаток после последнего интрона
            if pos < len(seq):
                mature_parts.append(seq[pos:])
        else:
            # Если не было ни одного интрона - вся последовательность является экзоном
            mature_parts.append(seq)

        return ''.join(mature_parts), log


class Protein:
    def __init__(self, mol_rna: RnaMolecule):
        pass
        # self.mol_rna = mol_rna
        # # Генетический код для трансляции
        # # https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
        # genetic_code = {
        #     'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        #     'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        #     'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        #     'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        #     'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        #     'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        #     'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        #     'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        #     'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        #     'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        #     'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        #     'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        #     'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        #     'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        #     'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        #     'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        # }
        # start_codon = "AUG"
        # start_pos = self.mol_rna.strand.find(start_codon)
        # if start_pos == -1:
        #     self.sequence = ""
        #     return
        #
        # protein_seq = []
        # for i in range(start_pos, len(self.mol_rna.strand), 3):
        #     codon = self.mol_rna.strand[i:i + 3]
        #     if codon in genetic_code:
        #         aa = genetic_code[codon]
        #         if aa == '*':  # Стоп-кодон
        #             break
        #         protein_seq.append(aa)
        #
        # self.sequence = ''.join(protein_seq)