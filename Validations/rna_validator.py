#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple human mRNA sanity validator.

Checks:
- first AUG and Kozak context
- longest ORF length
- trailing poly(A) tail length and internal poly(A) islands
- PAS (AAUAAA variants) 10–30 nt upstream of tail
- stop codon densities per frame

Usage: python rna_validator.py /path/to/file.txt
The script will try to extract the first RNA-like block (ACGU-only).
"""

import sys, re
from collections import defaultdict

stop_codons = {"UAA","UAG","UGA"}

def read_text(p):
    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        return f.read()

def extract_rna_blocks(text):
    return [(m.start(), m.end(), m.group(0)) for m in re.finditer(r'[ACGU]{30,}', text)]

def first_aug(s): 
    return s.find("AUG")

def orf_from(s, start_idx):
    if start_idx < 0: return None
    i = start_idx
    while i+3 <= len(s):
        cod = s[i:i+3]
        if cod in stop_codons:
            return {"orf_start": start_idx, "orf_end": i+3, "stop": cod, "length_nt": i+3-start_idx}
        i += 3
    return {"orf_start": start_idx, "orf_end": len(s), "stop": None, "length_nt": len(s)-start_idx}

def kozak(s, start_idx):
    minus3 = s[start_idx-3] if start_idx-3 >= 0 else None
    plus4 = s[start_idx+4] if start_idx+4 < len(s) else None
    return {"minus3": minus3, "plus4": plus4, "ok": (minus3 in ("A","G")) and (plus4 == "G")}

def trailing_polyA(s):
    m = re.search(r'A+$', s)
    return (len(m.group(0)), m.start()) if m else (0, None)

def internal_polyA_runs(s, min_len=8, tail_start=None):
    runs = []
    for m in re.finditer(rf'A{{{min_len},}}', s):
        if tail_start is not None and m.start() >= tail_start:
            continue
        runs.append((m.start(), m.end()-m.start()))
    return runs

def pas_hits(s, cleavage_idx, window=(10,30)):
    motifs = {"AAUAAA","AUUAAA","AGUAAA","AAUAUA","AAUACA","AAUAAG"}
    start = max(0, cleavage_idx - window[1] - 6)
    end = max(0, cleavage_idx - window[0])
    hits = []
    region = s[start:end]
    for i in range(len(region)-5):
        kmer = region[i:i+6]
        if kmer in motifs:
            hits.append((start+i, kmer, cleavage_idx - (start+i)))
    return hits

def stops_in_frame(s, frame):
    c=0
    for i in range(frame, len(s)-2, 3):
        if s[i:i+3] in stop_codons: c+=1
    return c

def window(s, center, flank=30):
    a = max(0, center-flank)
    b = min(len(s), center+flank)
    return s[a:b]

def main(p):
    text = read_text(p)
    blocks = extract_rna_blocks(text)
    if not blocks:
        print("No RNA-like blocks (ACGU) found.")
        return 1
    start, end, seq = blocks[0]
    print(f"RNA block at bytes {start}-{end}, length {len(seq)} nt")
    fa = first_aug(seq)
    print(f"First AUG at: {fa}")
    if fa != -1:
        kz = kozak(seq, fa)
        print(f"Kozak: -3={kz['minus3']}, +4={kz['plus4']}, match={kz['ok']}")
        orf = orf_from(seq, fa)
        print(f"First ORF: length={orf['length_nt']} nt, stop={orf['stop']}, stop_pos={orf['orf_end']-3 if orf['stop'] else None}")
        print("First AUG context (±30):", window(seq, fa, 30))
        if orf["stop"]:
            print("First ORF stop context (±30):", window(seq, orf["orf_end"]-3, 30))
    # Longest ORF
    best=None
    for i in range(len(seq)-2):
        if seq[i:i+3]=="AUG":
            o=orf_from(seq,i)
            if o and o["stop"] and (not best or o["length_nt"]>best["length_nt"]):
                best=o
    if best:
        print(f"Longest ORF: {best['length_nt']} nt ({best['length_nt']//3} aa), stop={best['stop']}, start={best['orf_start']}, end={best['orf_end']}")
        print("Longest ORF start snippet:", seq[best['orf_start']:best['orf_start']+60])
        print("Longest ORF end snippet:", seq[best['orf_end']-60:best['orf_end']])
    # PolyA
    tail_len, tail_start = trailing_polyA(seq)
    print(f"Trailing poly(A) length: {tail_len} nt, tail_start={tail_start}")
    runs = internal_polyA_runs(seq, min_len=8, tail_start=tail_start)
    print(f"Internal poly(A) runs (len≥8): {len(runs)}. Examples (pos,len) first 10:", runs[:10])
    # PAS
    if tail_start is not None:
        pas = pas_hits(seq, tail_start, window=(10,30))
        print("PAS hits 10–30 nt upstream of tail (pos,motif,distance):", pas[:10])
    # Stops density
    for f in range(3):
        c = stops_in_frame(seq,f)
        dens = round(c/(len(seq)/1000.0),2)
        print(f"Stops frame {f}: {c} ({dens} per kb)")
    # Tail context
    print("Tail last 120 nt:", seq[-120:])
    return 0

if __name__=="__main__":
    if len(sys.argv)<2:
        print("Usage: python rna_validator.py /path/to/file.txt")
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
