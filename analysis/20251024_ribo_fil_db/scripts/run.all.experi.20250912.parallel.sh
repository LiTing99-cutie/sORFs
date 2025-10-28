#!/usr/bin/env bash
set -euo pipefail
meta_file=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
J=${J:-4}

gen_list() { cut -d',' -f1 "$meta_file" | sed 's/-/_/g'; }

run_one() {
  sample="$1"
  bash -c '
    set -euo pipefail
    mkdir -p ../processed/db_search_merge logs
    log="logs/'"$sample"'.log"
    {
      echo "***Processing '"$sample"' at $(date +%F" "%T)"
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate base
      unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
      eval "$(perl -Mlocal::lib=--deactivate-all)" 2>/dev/null || true

      msf_closed=$(find ../processed/db_search_closed -type f -name psm.tsv -path "*'"$sample"'*" | head -n1 || true)
      msf_open=$(  find ../processed/db_search_open   -type f -name psm.tsv -path "*'"$sample"'*" | head -n1 || true)
      pfind_closed=$(ls ../processed/pFind_res_20251026/closed/*"'"$sample"'"*.spectra 2>/dev/null | head -n1 || true)
      pfind_open=$(  ls ../processed/pFind_res_20251026/open/*"'"$sample"'"*.spectra   2>/dev/null | head -n1 || true)
      if [[ -z "$msf_closed" || -z "$msf_open" || -z "$pfind_closed" || -z "$pfind_open" ]]; then
        echo "SKIP '"$sample"': missing inputs"; exit 0; fi

      out_merge="../processed/db_search_merge/'"$sample"'.tsv"
      python merge.db.search.res.20250909.py --msf_closed "$msf_closed" --msf_open "$msf_open" \
        --pfind_closed "$pfind_closed" --pfind_open "$pfind_open" --out "$out_merge"
      [[ -s "$out_merge" ]] || { echo "ERROR: merge output empty"; exit 1; }

      python fold_by_peptide_il.20250909.py -i "$out_merge" -o "../processed/db_search_merge/'"$sample"'.peptide.tsv"
      bash protein.map.20250909.uni.sh "'"$sample"'"
      echo "***Done '"$sample"' at $(date +%F" "%T)"
    } >"$log" 2>&1
  '
}

export -f run_one
gen_list | xargs -n1 -P "$J" -I{} bash -c 'run_one "$@"' _ {}
