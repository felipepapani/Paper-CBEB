# -*- coding: utf-8 -*-
"""
R-4 (bloqueante do parecer CBEB) — Controle do contraste de domínio
====================================================================
Objeção de R2/F2.1: o contraste amplitude vs. VFC compara 20 pacientes
(amplitude) contra 7 (VFC). A diferença de Cohen's d pode vir da amostra,
não do domínio. Este script ELIMINA esse confundimento recomputando o d
far/near do braço de AMPLITUDE restrito EXATAMENTE aos 7 pacientes do
braço de VFC.

Se o d de amplitude nesses 7 permanecer ≈ 0,19 (como nos 20), o contraste
ganha uma perna real: mesma amostra, mesma pergunta, domínios diferentes.
Se subir muito, você precisa saber disso ANTES do revisor.

NÃO reprocessa sinal — só filtra os manifests CSV que você já tem.
Reusa a MESMA função d_per_patient da célula 111 (faixa longe/perto, por record).
"""

import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Reutiliza csv_path e FEATURES se já existirem no notebook; senão, define.
# ---------------------------------------------------------------------------
try:
    csv_path
except NameError:
    DRIVE_TCC = globals().get("DRIVE_TCC", "/content/drive/MyDrive/TCC")
    CSV_DIR = globals().get("CSV_DIR", os.path.join(DRIVE_TCC, "csvs"))
    def csv_path(nome):
        if not nome.endswith(".csv"):
            nome += ".csv"
        for base in (CSV_DIR, DRIVE_TCC):
            p = os.path.join(base, nome)
            if os.path.exists(p):
                return p
        return os.path.join(CSV_DIR, nome)

FEATURES = globals().get("FEATURES", ['RR','DET','L','L_max','DIV','L_entr','LAM',
        'TT','V_max','V_entr','W','W_max','W_entr','DET_RR','LAM_DET'])

P_AMP = 10   # mesmo manifesto de amplitude usado na Fig. 4 (troque p/ 40 se for o caso)


def cohen_d_pareado(far, near):
    """
    Cohen's d PAREADO (d_z): far e near são medidos no MESMO paciente.
    d_z = média(diferença) / desvio-padrão(diferença).
    É o estimador correto para o desenho far/near (R2/F2.2 pede declarar isto).
    """
    dif = np.asarray(far) - np.asarray(near)
    if len(dif) < 2 or dif.std(ddof=1) == 0:
        return np.nan
    return dif.mean() / dif.std(ddof=1)


def d_por_paciente(df_in, features, col_band="faixa", registros=None, pareado=True):
    """
    Cohen's d (far vs near), UM valor por paciente por faixa.
    - registros: se fornecido, restringe a esse subconjunto (o filtro do R-4).
    - pareado=True  -> d_z (recomendado; mesmo paciente nas duas faixas)
      pareado=False -> d_s (entre grupos), como estava na célula 111.
    """
    rows = []
    for feat in features:
        if feat not in df_in.columns:
            continue
        sub = df_in[df_in[col_band].isin(["longe", "perto"])].copy()
        if registros is not None:
            sub = sub[sub["record"].astype(str).isin([str(r) for r in registros])]
        agg = sub.groupby(["record", col_band])[feat].mean().reset_index()
        counts = agg.groupby("record")[col_band].nunique()
        completos = counts[counts == 2].index         # tem far E near
        if len(completos) < 3:
            continue
        ok = agg[agg["record"].isin(completos)].sort_values("record")
        far = ok[ok[col_band] == "longe"][feat].values
        near = ok[ok[col_band] == "perto"][feat].values
        if pareado:
            d = cohen_d_pareado(far, near)
        else:
            sp = np.sqrt(((len(far)-1)*far.std(ddof=1)**2 + (len(near)-1)*near.std(ddof=1)**2)
                         / (len(far)+len(near)-2))
            d = (far.mean()-near.mean())/sp if sp > 0 else np.nan
        rows.append({"feature": feat, "n_patients": len(completos),
                     "cohen_d": round(abs(d), 3) if np.isfinite(d) else np.nan})
    return pd.DataFrame(rows).sort_values("cohen_d", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# EXECUÇÃO
# ---------------------------------------------------------------------------
def main():
    # 1. carrega os dois manifests
    df_amp = pd.read_csv(csv_path(f"manifesto_p{P_AMP:02d}.csv"))
    df_hrv = pd.read_csv(csv_path("manifesto_hrv_rqa.csv"))
    if "aprovada" in df_hrv.columns:
        df_hrv = df_hrv[df_hrv.aprovada]

    # 2. descobre os 7 registros do braço de VFC (os que aparecem no manifesto HRV)
    sete = sorted(df_hrv["record"].astype(str).unique(), key=lambda x: int(x))
    print("=" * 72)
    print(f"  Registros do braço de VFC (n={len(sete)}): {sete}")
    print("=" * 72)

    # 3. os três d far/near:
    #    (a) amplitude nos 20 (como no paper)   (b) amplitude nos MESMOS 7   (c) VFC nos 7
    d_amp20_z = d_por_paciente(df_amp, FEATURES, pareado=True)
    d_amp7_z  = d_por_paciente(df_amp, FEATURES, registros=sete, pareado=True)
    d_hrv7_z  = d_por_paciente(df_hrv, FEATURES, pareado=True)
    # versão d_s (entre grupos) só para amplitude-7, para comparar com o d do paper
    d_amp7_s  = d_por_paciente(df_amp, FEATURES, registros=sete, pareado=False)

    comp = (d_amp20_z[["feature", "cohen_d"]].rename(columns={"cohen_d": "AMP_20_dz"})
            .merge(d_amp7_z[["feature", "cohen_d"]].rename(columns={"cohen_d": "AMP_7_dz"}), on="feature", how="outer")
            .merge(d_amp7_s[["feature", "cohen_d"]].rename(columns={"cohen_d": "AMP_7_ds"}), on="feature", how="outer")
            .merge(d_hrv7_z[["feature", "cohen_d"]].rename(columns={"cohen_d": "HRV_7_dz"}), on="feature", how="outer"))
    comp = comp.sort_values("HRV_7_dz", ascending=False).reset_index(drop=True)

    print("\n  Cohen's d por paciente — far vs near")
    print("  AMP_20 = amplitude, 20 pac | AMP_7 = amplitude nos 7 do braço VFC | HRV_7 = VFC")
    print("  dz = pareado (recomendado) · ds = entre grupos (como na Fig. 4 atual)\n")
    print(comp.to_string(index=False))

    amp20 = d_amp20_z.cohen_d.max()
    amp7  = d_amp7_z.cohen_d.max()
    hrv7  = d_hrv7_z.cohen_d.max()
    print("\n" + "=" * 72)
    print("  VEREDITO (maior d de cada braço, pareado d_z)".center(72))
    print("=" * 72)
    print(f"  Amplitude, 20 pacientes : {amp20:.3f}")
    print(f"  Amplitude, MESMOS 7     : {amp7:.3f}   <-- o número que o R-4 pede")
    print(f"  VFC, 7 pacientes        : {hrv7:.3f}")
    if amp7 < 0.2 <= hrv7:
        print("\n  >> CONTRASTE SOBREVIVE ao controle de amostra: mesmos 7 pacientes,")
        print("     amplitude continua sem separar (d<0,2) e VFC separa. Esta é a")
        print("     frase forte para responder R2/F2.1 no texto.")
    elif amp7 >= 0.2:
        print("\n  >> ATENÇÃO: nos mesmos 7, a amplitude também separa (d>=0,2). Parte do")
        print("     contraste vinha da amostra, não do domínio. O texto precisa dizer isto")
        print("     honestamente — e o enquadramento do título precisa ser suavizado.")
    print("=" * 72)

    saida = os.path.join(os.path.dirname(csv_path("x.csv")), "R4_contraste_mesmos7.csv")
    try:
        comp.to_csv(saida, index=False)
        print(f"\nSalvo: {saida}")
    except Exception as e:
        print(f"\n(Não consegui salvar CSV: {e})")

    return comp


if __name__ == "__main__":
    main()
