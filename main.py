from MolecularBiology import *

def main():
    # 1) Просто уникальная «человеко‑подобная» ДНК (GC ~ 41 %)
    dna = DnaMolecule(length=100_000, gc=0.41)  # run_id создастся сам
    print(f'ID: {dna.run_id}')  # сохраните, чтобы воспроизвести
    print(f"Name DNA: {dna.dna_5_to_3_strand.name}")
    print(f"Name DNA: {dna.dna_3_to_5_strand.name}")
    #print(f"5'-3' strand DNA: {dna.dna_5_to_3_strand.content}")
    #print(f"3'-5' strand DNA: {dna.dna_3_to_5_strand.content}")
    # 2) Повторить тот же запуск:
    same_dna = DnaMolecule(length=100_000, gc=0.41, run_id=dna.run_id)
    print(f'ID: {same_dna.run_id}')  # сохраните, чтобы воспроизвести
    print(f"Name DNA: {same_dna.dna_5_to_3_strand.name}")
    print(f"Name DNA: {same_dna.dna_3_to_5_strand.name}")
    #print(f"5'-3' strand DNA: {same_dna.dna_5_to_3_strand.content}")
    #print(f"3'-5' strand DNA: {same_dna.dna_3_to_5_strand.content}")

if __name__ == "__main__":
    main()