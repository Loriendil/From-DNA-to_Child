from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import gzip
import hashlib
import os, re
import random

@dataclass
class Strand:
        name: str
        content: str

class Utility:
    ALPHABET = "ACGT"
    C2I = {c: i for i, c in enumerate(ALPHABET)}  # 'A'->0, 'C'->1, 'G'->2, 'T'->3
    I2C = {i: c for c, i in C2I.items()}

    # IUPAC-коды и комплементарность (в верхнем регистре)
    _IUPAC = "ACGTRYSWKMBDHVN"
    _RC = "TGCAYRWSMKVHDBN"  # A<->T, C<->G; R<->Y; S<->S; W<->W; K<->M; B<->V; D<->H; N<->N
    RC_TABLE = str.maketrans(_IUPAC, _RC)

    @staticmethod
    def open_text_auto(path: str):
        """Открываем как текст; поддерживаем .gz автоматически."""
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "rt")

    @staticmethod
    def sha256_fasta_stream(path: str, *, include_headers: bool = False) -> bytes:
        """
        Стриминговый SHA-256 по FASTA: по умолчанию хэш только по нуклеотидам,
        заголовки строк '>' пропускаем; пробелы/переводы строк игнорируем.
        Возвращает байтовый digest (32 байта).
        """
        h = hashlib.sha256()
        with Utility.open_text_auto(path) as f:
            for line in f:
                if not include_headers and line.startswith(">"):
                    continue
                s = line.strip().upper()
                if s:
                    h.update(s.encode("ascii"))
        return h.digest()

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Обратный комплемент с поддержкой IUPAC-кодов."""
        return seq.upper().translate(Utility.RC_TABLE)[::-1]

    @staticmethod
    def weights_for_gc(gc: float = 0.41) -> Tuple[float, float, float, float]:
        """
        Целевые доли A,T,G,C при заданном GC (для человека ~0.41).
        Возвращаем кортеж в порядке A,C,G,T? -> У нас порядок A,C,Г,T (см. C2I).
        ВНИМАНИЕ: C2I = {A:0, C:1, G:2, T:3}
        """
        gc = min(max(gc, 0.0), 1.0)
        at = 1.0 - gc
        # вернём в порядке A, C, G, T
        return (at / 2, gc / 2, gc / 2, at / 2)

    @staticmethod
    def kmer_to_index(kmer: str) -> int:
        """Кодируем k-мер в целочисленный индекс по основанию 4."""
        x = 0
        for ch in kmer:
            x = x * 4 + Utility.C2I[ch]
        return x

    @staticmethod
    def index_to_kmer(idx: int, k: int) -> str:
        """Обратное преобразование индекса -> k-мер."""
        out = []
        for _ in range(k):
            out.append(Utility.I2C[idx % 4])
            idx //= 4
        return "".join(reversed(out))

class BuildingMarkovModel:
    ALPHABET = "ACGT"
    C2I = {c: i for i, c in enumerate(ALPHABET)}  # 'A'->0, 'C'->1, 'G'->2, 'T'->3
    I2C = {i: c for c, i in C2I.items()}

    # IUPAC-коды и комплементарность (в верхнем регистре)
    _IUPAC = "ACGTRYSWKMBDHVN"
    _RC = "TGCAYRWSMKVHDBN"  # A<->T, C<->G; R<->Y; S<->S; W<->W; K<->M; B<->V; D<->H; N<->N
    RC_TABLE = str.maketrans(_IUPAC, _RC)

    @staticmethod
    def markov_counts_from_fasta(path: str, k: int = 2) -> Dict[int, List[int]]:
        """
        Строит словарь counts: k-мер (индекс) -> [cA, cC, cG, cT].
        Стримингово читает FASTA; окна, пересекающие не-ACGT символы (напр. 'N'), обнуляются.
        """
        if k <= 0:
            raise ValueError("k должен быть >= 1")
        counts: Dict[int, List[int]] = {}
        # для сдвига окна: base = 4^(k-1)
        base = 4 ** (k - 1)
        prefix_idx = 0  # индекс k-мера для первых k букв
        have = 0  # сколько валидных (ACGT) букв накопили подряд

        with Utility.open_text_auto(path) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                for ch in line.strip().upper():
                    val = Utility.C2I.get(ch)
                    if val is None:
                        # встретили N или иной символ — разрываем окно
                        have = 0
                        prefix_idx = 0
                        continue
                    if have < k:
                        # накапливаем первые k букв
                        prefix_idx = prefix_idx * 4 + val
                        have += 1
                    else:
                        # есть валидный префикс длины k, обновляем переход на val
                        row = counts.get(prefix_idx)
                        if row is None:
                            row = [0, 0, 0, 0]
                            counts[prefix_idx] = row
                        row[val] += 1
                        # сдвигаем окно: выкидываем старший символ, добавляем val
                        prefix_idx = ((prefix_idx % base) * 4) + val
        return counts

    @staticmethod
    def row_probs(row: List[int],
                   alpha: float,
                   gc_weights: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
        """
        Нормируем строку счётчиков в вероятности с псевдосчётчиками:
        p ~ count + alpha * gc_weight.
        Это аккуратно переносит глобальный GC-профиль на «редкие» контексты.
        """
        a, c, g, t = row
        gwA, gwC, gwG, gwT = gc_weights
        rA = a + alpha * gwA
        rC = c + alpha * gwC
        rG = g + alpha * gwG
        rT = t + alpha * gwT
        s = rA + rC + rG + rT
        if s == 0:  # полностью пустая строка (не встретилась ни разу)
            return (gwA, gwC, gwG, gwT)
        return (rA / s, rC / s, rG / s, rT / s)

    @staticmethod
    def sample_idx_by_probs(probs: Tuple[float, float, float, float], rng: random.Random) -> int:
        """Сэмплируем индекс 0..3 по распределению probs (кумулятивно)."""
        r = rng.random()
        cum = 0.0
        for i, p in enumerate(probs):
            cum += p
            if r <= cum:
                return i
        return 3  # защита от потерь точности float

class DnaMolecule:
    """
        Синтетическая молекула ДНК, сгенерированная на основе:
          • референтной хромосомы (FASTA),
          • Марковской модели порядка k,
          • детерминированного сида (SHA-256 от FASTA + run_id),
          • и глобального GC-содержания для мягкого сглаживания.
        """

    def __init__(
            self,
            fasta_path: str,
            length: int,
            *,
            k: int = 2,
            gc: float = 0.41,
            run_id: Optional[str] = None,
            smoothing_alpha: float = 0.5,
            seed_include_headers: bool = False,
    ) -> None:
        if length < max(16, k + 1):
            raise ValueError("length слишком мал для генерации последовательности.")
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA не найден: {fasta_path}")

        self.k = k
        self.length = length
        self.gc = gc
        self.run_id = run_id

        # 1) Стриминговый SHA-256 от FASTA
        digest = Utility.sha256_fasta_stream(fasta_path, include_headers=seed_include_headers)

        # 2) Сид: смесь (FASTA_digest + run_id?) -> 64-битное целое
        material = digest + (run_id.encode("utf-8") if run_id else b"")
        seed64 = int.from_bytes(hashlib.sha256(material).digest()[:8], "big")
        rng = random.Random(seed64)

        # 3) Модель Маркова: k-мер -> счётчики A/C/G/T
        counts = BuildingMarkovModel.markov_counts_from_fasta(fasta_path, k=k)

        # 4) Выбор стартового k-мера:
        #    берём самый «массовый» (макс. сумма переходов),
        #    чтобы старт не попадал в «мёртвую» зону.
        if counts:
            start_idx = max(counts.keys(), key=lambda i: sum(counts[i]))
            seed_kmer = Utility.index_to_kmer(start_idx, k)
        else:
            # Если хромосома слишком короткая/пустая: старт — по GC весам
            gw = Utility.weights_for_gc(gc)
            seed_kmer = "".join(
                Utility.I2C[BuildingMarkovModel.sample_idx_by_probs(gw, rng)] for _ in range(k)
            )

        # 5) Генерация длиной `length`, включая seed_kmer в начало
        gw = Utility.weights_for_gc(gc)
        alpha = smoothing_alpha
        base = 4 ** (k - 1)

        seq_chars: List[str] = list(seed_kmer)
        prefix_idx = Utility.kmer_to_index(seed_kmer)

        while len(seq_chars) < length:
            row = counts.get(prefix_idx, [0, 0, 0, 0])
            probs = BuildingMarkovModel.row_probs(row, alpha=alpha, gc_weights=gw)
            nxt = BuildingMarkovModel.sample_idx_by_probs(probs, rng)
            seq_chars.append(Utility.I2C[nxt])
            prefix_idx = ((prefix_idx % base) * 4) + nxt

        seq5to3 = "".join(seq_chars)

        # 6) Вторая нить — обратный комплемент (биологически корректно)
        self.dna_5_to_3_strand = Strand("5'-3' direction of DNA", seq5to3)
        self.dna_3_to_5_strand = Strand("3'-5' direction of DNA", Utility.reverse_complement(seq5to3))

    # Оставляю публичные методы, если они вам нужны:
    @staticmethod
    def derive_seed(seed: Optional[bytes] = None, *, deterministic_key: Optional[str] = None) -> bytes:
        """Совместимость с вашим API — может пригодиться в других сценариях."""
        if seed is not None:
            return hashlib.sha256(seed).digest()
        if deterministic_key is not None:
            return hashlib.sha256(deterministic_key.encode()).digest()
        # Неденерминированный fallback (обычно лучше избегать)
        payload = f"{os.getpid()}|{os.urandom(16).hex()}".encode()
        return hashlib.sha256(payload).digest()

    @staticmethod
    def __set_genetic_mutations(self):
        pass


class RnaMolecule:
    """
    Класс молекулы РНК
    Входной:
    Параметр экземпляр класса DnaMolecule, который является программной репрезентацией одной молекулы ДНК
    """

    def __init__(self, mol_dna: DnaMolecule, gene_source: str = "3'-5' direction of DNA")->None:
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
        self.name: str = 'mature mRNA'
        self.mol_dna = mol_dna
        dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}

        if gene_source == "3'-5' direction of DNA":
            # шаблонная (antisense) нить -> комплементарная пре-мРНК
            first_rna_transcript = Strand(
                "from template (3'-5' DNA)",
                ''.join(dna_to_rna[nt] for nt in self.mol_dna.dna_3_to_5_strand.content)
            )
        else:
            # Для кодирующей цепи: реверс + комплементарная замена
            reversed_strand = self.mol_dna.dna_5_to_3_strand.content[::-1]
            first_rna_transcript = Strand("from coding (5'-3' DNA)",
                                          ''.join(dna_to_rna[nt] for nt in reversed_strand))

        # Кэпирование, полиаденилирование пропускаем, так как химические процессы в чистом виде нет цели симулировать.
        # Результатом работы конструктора до вызова метода splicing() пре-мРНК.
        pre_mrna: Strand = self.validate_rna_transcript(first_rna_transcript)
        self.strand, self.log = self.splicing(pre_mrna)

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