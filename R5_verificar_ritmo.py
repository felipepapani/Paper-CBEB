# -*- coding: utf-8 -*-
"""
R-5 (prioridade nº 1 do parecer CBEB) — Verificação de ritmo em [-30 s, 0)
===========================================================================
Pergunta que decide o conteúdo do artigo (R3/F3.1):
   Nas janelas dos últimos 30 s antes do vfon, quantos dos 20 registros
   estão em RITMO SINUSAL, quantos em TAQUICARDIA VENTRICULAR (VT/VFL) e
   quantos em outro ritmo?

Se a maioria estiver em VT, o "estado pré-fibrilatório" descrito no artigo
é, na verdade, taquicardia ventricular — e o texto precisa ser reformulado
como "caracterização topológica da transição VT->FV" (o que continua
publicável e é honesto).

COMO A SDDB ANOTA RITMO
-----------------------
As anotações de ritmo ficam em ann.aux_note, como strings iniciadas por '('
colocadas NO INÍCIO do episódio; o ritmo vale até a próxima anotação de ritmo.
Rótulos (WFDB / MIT):
   (N ou (NSR = ritmo sinusal normal
   (VT        = taquicardia ventricular
   (VFL       = flutter ventricular
   (VF/(VFIB  = fibrilação ventricular
   (AFIB      = fibrilação atrial   (PM = marcapasso, etc.)

Este script NÃO reprocessa sinal — só lê os .atr e os headers que você já usa.

REQUISITOS: wfdb (o mesmo do seu pipeline). Ajuste apenas DATA_DIR.
"""

import os
import glob
import numpy as np
import pandas as pd
import wfdb

# ---------------------------------------------------------------------------
# AJUSTE AQUI: pasta onde estão os .hea/.dat/.atr da SDDB (a que seu código já usa)
# ---------------------------------------------------------------------------
DATA_DIR = "/content/drive/MyDrive/TCC/matrizes"   # <-- troque se necessário
FS = 250                                            # Hz (SDDB)
JANELA_S = 30                                       # últimos 30 s antes do vfon

# Rótulos de ritmo (sem o parêntese inicial)
SINUSAL = {"N", "NSR"}
TAQUI_V = {"VT", "VFL"}          # taquicardia/flutter ventricular = o que R3 suspeita
FIB_V   = {"VF", "VFIB"}         # fibrilação já instalada
OUTROS_CONHECIDOS = {"AFIB", "PM", "B", "BII", "SBR", "SVTA", "NOD",
                     "VER", "HGEA", "ASYS", "NOISE", "BI"}


def limpar_rotulo(aux):
    """Extrai o rótulo de ritmo de uma entrada aux_note ('(VT\\x00' -> 'VT')."""
    if aux is None:
        return None
    s = str(aux).replace("\x00", "").strip()
    if not s:
        return None
    if s.startswith("("):
        s = s[1:]
    return s.strip() or None


def vfon_do_header(rec_base):
    """Lê o vfon (em amostras) do comentário do header. Retorna None se ausente."""
    hea = rec_base + ".hea"
    if not os.path.exists(hea):
        return None
    with open(hea, errors="ignore") as f:
        for linha in f:
            if "vfon" in linha.lower():
                # formato típico: "#vfon: HH:MM:SS"  ou "# VF onset: HH:MM:SS"
                try:
                    hh, mm, ss = map(int, linha.split(":", 1)[1].strip().split(":"))
                    return (hh * 3600 + mm * 60 + ss) * FS
                except Exception:
                    pass
    return None


def ritmo_no_intervalo(rec_base, ini_amostra, fim_amostra):
    """
    Devolve o conjunto de ritmos VIGENTES no intervalo [ini, fim] em amostras.
    Como o ritmo vale até a próxima anotação, considera-se:
      - toda anotação de ritmo cujo instante cai dentro do intervalo, e
      - a última anotação de ritmo ANTES do início (ritmo herdado).
    """
    ann = wfdb.rdann(rec_base, "atr")
    # pares (amostra, rotulo) só das entradas que têm rótulo de ritmo
    ritmos = []
    for amostra, aux in zip(ann.sample, ann.aux_note):
        rot = limpar_rotulo(aux)
        if rot:
            ritmos.append((int(amostra), rot))
    if not ritmos:
        return set(), "sem anotacao de ritmo"

    vigentes = []
    # ritmo herdado: último rótulo antes de ini
    herdado = None
    for amostra, rot in ritmos:
        if amostra <= ini_amostra:
            herdado = rot
        elif ini_amostra < amostra <= fim_amostra:
            vigentes.append(rot)
    if herdado is not None:
        vigentes.insert(0, herdado)

    return set(vigentes), ";".join(dict.fromkeys(vigentes)) or "(nenhum no intervalo)"


def classificar(conj):
    """Classifica o conjunto de ritmos do intervalo numa categoria única."""
    if conj & TAQUI_V:
        return "TV"           # taquicardia/flutter ventricular presente
    if conj & FIB_V:
        return "FV_ja"        # já fibrilando dentro dos 30 s
    if conj & SINUSAL and not (conj & (TAQUI_V | FIB_V)):
        return "sinusal"
    if conj & {"AFIB"}:
        return "FA"
    if conj & {"PM"}:
        return "marcapasso"
    if not conj:
        return "indefinido"
    return "outro"


def descobrir_registros(pasta):
    """Lista os record-base (sem extensão) que têm .hea nesta pasta ou subpastas."""
    bases = set()
    for hea in glob.glob(os.path.join(pasta, "**", "*.hea"), recursive=True):
        bases.add(os.path.splitext(hea)[0])
    return sorted(bases, key=lambda p: os.path.basename(p))


# ---------------------------------------------------------------------------
# EXECUÇÃO
# ---------------------------------------------------------------------------
def main():
    registros = descobrir_registros(DATA_DIR)
    if not registros:
        print(f"Nenhum .hea encontrado em {DATA_DIR}. Ajuste DATA_DIR.")
        return

    linhas = []
    for base in registros:
        rec = os.path.basename(base)
        vfon = vfon_do_header(base)
        if vfon is None:
            linhas.append({"record": rec, "vfon_s": None, "categoria": "sem_vfon",
                           "ritmos_no_intervalo": "-"})
            continue

        ini = max(0, vfon - JANELA_S * FS)   # -30 s
        fim = vfon                            # 0 (instante do vfon)
        try:
            conj, detalhe = ritmo_no_intervalo(base, ini, fim)
        except FileNotFoundError:
            linhas.append({"record": rec, "vfon_s": vfon // FS, "categoria": "sem_.atr",
                           "ritmos_no_intervalo": "-"})
            continue

        linhas.append({
            "record": rec,
            "vfon_s": vfon // FS,
            "categoria": classificar(conj),
            "ritmos_no_intervalo": detalhe,
        })

    df = pd.DataFrame(linhas).sort_values("record").reset_index(drop=True)

    print("=" * 74)
    print("  R-5 — RITMO EM [-30 s, 0) POR REGISTRO".center(74))
    print("=" * 74)
    print(df.to_string(index=False))

    print("\n" + "=" * 74)
    print("  CONTAGEM POR CATEGORIA".center(74))
    print("=" * 74)
    cont = df["categoria"].value_counts()
    print(cont.to_string())

    n_sinusal = int((df.categoria == "sinusal").sum())
    n_tv = int((df.categoria == "TV").sum())
    n_fv = int((df.categoria == "FV_ja").sum())
    n_val = int(df.categoria.isin(["sinusal", "TV", "FV_ja", "FA",
                                   "marcapasso", "outro"]).sum())

    print("\n" + "=" * 74)
    print("  VEREDITO".center(74))
    print("=" * 74)
    print(f"  Registros avaliados        : {n_val}")
    print(f"  Sinusal em [-30 s, 0)      : {n_sinusal}")
    print(f"  Taquicardia ventricular    : {n_tv}")
    print(f"  Já em FV nesse intervalo   : {n_fv}")
    if n_tv + n_fv > n_sinusal:
        print("\n  >> A MAIORIA NÃO está em ritmo sinusal nos últimos 30 s.")
        print("     O objeto descrito como 'estado pré-fibrilatório' é, em grande")
        print("     parte, TAQUICARDIA VENTRICULAR. O texto do artigo precisa ser")
        print("     reformulado como caracterização da transição VT->FV (honesto e")
        print("     ainda publicável).")
    else:
        print("\n  >> A maioria ESTÁ em ritmo sinusal até perto do evento.")
        print("     O enquadramento 'pré-fibrilatório' se sustenta; reporte a")
        print("     contagem no texto como evidência contra a hipótese de VT (R3).")
    print("=" * 74)

    # salva CSV para você anexar como suplementar / citar no texto
    saida = os.path.join(os.path.dirname(DATA_DIR) or ".", "R5_ritmo_por_registro.csv")
    try:
        df.to_csv(saida, index=False)
        print(f"\nSalvo: {saida}")
    except Exception as e:
        print(f"\n(Não consegui salvar CSV: {e})")


if __name__ == "__main__":
    main()
