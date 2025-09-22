from MolecularBiology import *
from pathlib import Path

def main():
    fasta = r"E:\local_repository\Python-repository\From-DNA-to_Child\split\chr1.fa"  # <-- ваш путь (можно .gz)
    L = 100_000
    K = 2

    # 1) Просто уникальная «человеко‑подобная» ДНК (GC ~ 41 %)
    dna = DnaMolecule(fasta_path=fasta, length=L, k=K, gc=0.41, run_id= 'char_06CQ3B77EPH3BFKNNAK5KCDFY0-MY')
    print(f"ID: {dna.run_id}")  # может быть None, если не передавали
    print(f"5'->3' length: {len(dna.dna_5_to_3_strand.content)}")
    print(f"3'->5' length: {len(dna.dna_3_to_5_strand.content)}")


    # 2) Повторить тот же запуск:
    same_dna = DnaMolecule(fasta_path=fasta, length=L, k=K, gc=0.41, run_id=dna.run_id)
    print(f'ID: {same_dna.run_id}')
    assert dna.dna_5_to_3_strand.content == same_dna.dna_5_to_3_strand.content
    print("Повтор воспроизведён: последовательности идентичны.")

    # 3) пишем результат в файл
    rna = RnaMolecule(dna)
    text = ("ID: " + dna.run_id + "\n" + "3'-5' direction of DNA\n" +
            dna.dna_5_to_3_strand.content + "\nRNA strand\n" + rna.strand)
    path = Path("out.txt")
    path.write_text(text, encoding="utf-8")  # создаст/перезапишет файл

if __name__ == "__main__":
    main()