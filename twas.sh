#!/usr/bin/env bash
set -euo pipefail

DISC="/Users/an/Desktop/Discovery data"
LDREF="/Users/an/Desktop/LDREF/1000G.EUR."
FUSION="/Users/an/Desktop/fusion_twas-master/FUSION.assoc_test.R"
WEIGHTS_DIR="/Users/an/Desktop/twas/WEIGHTS/GTExv8.ALL.Brain_Nucleus_accumbens_basal_ganglia"

OUTSUM="/Users/an/Desktop/FUSION_SUMSTATS_BY_PHENO"
OUTTWAS="/Users/an/Desktop/TWAS_RESULTS_BY_PHENO"

mkdir -p "$OUTSUM" "$OUTTWAS"

phenoname() {
  case "$1" in
    1) echo "Nicotine_Composite_Score" ;;
    2) echo "Alcohol_Consumption_Composite_Score" ;;
    3) echo "Alcohol_Dependence_Composite_Score" ;;
    4) echo "Illicit_Drug_Composite_Score" ;;
    5) echo "Behavioral_Disinhibition_Composite_Score" ;;
    *) echo "UNKNOWN_PHENO_$1" ;;
  esac
}

convert_assoc_to_fusion() {
  local infile="$1"
  local outfile="$2"

  python3 - "$infile" "$outfile" <<'PY'
import sys, math
import pandas as pd

infile, outfile = sys.argv[1], sys.argv[2]
df = pd.read_csv(infile, sep="\t", dtype=str)

# drop repeated header rows if present
df = df[~((df.iloc[:,0].str.lower()=="chr") & (df.iloc[:,1].str.lower()=="rs"))].copy()
df.columns = [c.strip() for c in df.columns]

need = ["chr","rs","ps","allele1","allele0"]
miss = [c for c in need if c not in df.columns]
if miss:
  raise SystemExit(f"{infile}: missing {miss}. Have: {list(df.columns)}")

# choose p-value column
pcol = None
for cand in ["p_wald","p_lrt","p_score","p"]:
  if cand in df.columns:
    pcol = cand
    break
if pcol is None:
  raise SystemExit(f"{infile}: no p-value column. Have: {list(df.columns)}")

CHR = pd.to_numeric(df["chr"], errors="coerce")
BP  = pd.to_numeric(df["ps"], errors="coerce")
SNP = df["rs"].astype(str)
A1  = df["allele1"].astype(str)
A2  = df["allele0"].astype(str)
P   = pd.to_numeric(df[pcol], errors="coerce")

beta = pd.to_numeric(df["beta"], errors="coerce") if "beta" in df.columns else None
se   = pd.to_numeric(df["se"],   errors="coerce") if "se"   in df.columns else None

if beta is not None and se is not None:
  Z = beta / se
else:
  # fallback: unsigned Z from p (direction unknown)
  p = P.copy()
  p[p<=0] = 1e-300
  p[p>1]  = 1.0
  try:
    import mpmath as mp
    def z_from_p(pv):
      # two-sided p -> |z|
      # invnorm(1-p/2)
      return float(mp.sqrt(2) * mp.erfinv(1 - pv))
    Z = p.apply(z_from_p).abs()
  except Exception:
    # crude fallback
    Z = p.apply(lambda pv: math.sqrt(max(0.0, -2.0*math.log(max(pv,1e-300)))))

  Z = pd.to_numeric(Z, errors="coerce")

out = pd.DataFrame({"CHR":CHR,"SNP":SNP,"BP":BP,"A1":A1,"A2":A2,"Z":Z,"P":P})
out = out.dropna(subset=["CHR","BP","Z","P"])
out = out[(out["CHR"]>=1) & (out["CHR"]<=22)]
out = out[(out["P"]>0) & (out["P"]<=1)]
out["CHR"] = out["CHR"].astype(int)
out["BP"]  = out["BP"].astype(int)

out.to_csv(outfile, sep="\t", index=False)
PY
}

for p in 1 2 3 4 5; do
  PNAME="$(phenoname "$p")"
  echo "=============================="
  echo "PHENO $p => $PNAME"
  echo "=============================="

  SUMDIR="$OUTSUM/$PNAME"
  TWASDIR="$OUTTWAS/$PNAME"
  mkdir -p "$SUMDIR" "$TWASDIR"

  # convert chr1-22 -> FUSION.sumstats
  for c in $(seq 1 22); do
    INFILE="$DISC/chr${c}_pheno${p}_lmm.assoc.txt"
    OUTFILE="$SUMDIR/${PNAME}.chr${c}.FUSION.sumstats"
    if [[ ! -f "$INFILE" ]]; then
      echo "WARN missing: $INFILE"
      continue
    fi
    echo "  -> converting chr$c"
    convert_assoc_to_fusion "$INFILE" "$OUTFILE"
  done

  # run TWAS chr1-22
  for c in $(seq 1 22); do
    GW="$SUMDIR/${PNAME}.chr${c}.FUSION.sumstats"
    if [[ ! -f "$GW" ]]; then
      echo "WARN missing GW: $GW"
      continue
    fi
    echo "  -> TWAS chr$c"
    Rscript "$FUSION" \
      --sumstats "$GW" \
      --weights "$WEIGHTS_DIR" \
      --weights_dir "$WEIGHTS_DIR" \
      --ref_ld_chr "$LDREF" \
      --chr "$c" \
      --out "$TWASDIR/${PNAME}_NAcc_chr${c}" \
      > "$TWASDIR/chr${c}.log" 2>&1
  done

  # merge dats
  echo -e "ID\tCHR\tP0\tP1\tTWAS.Z\tTWAS.P" > "$TWASDIR/${PNAME}_NAcc_ALLCHR.tsv"
  for c in $(seq 1 22); do
    DAT="$TWASDIR/${PNAME}_NAcc_chr${c}.dat"
    [[ -f "$DAT" ]] && tail -n +2 "$DAT" >> "$TWASDIR/${PNAME}_NAcc_ALLCHR.tsv" || true
  done

  echo "DONE $PNAME"
done

echo "=============================="
echo "ALL DONE"
echo "Sumstats root: $OUTSUM"
echo "TWAS root:     $OUTTWAS"
