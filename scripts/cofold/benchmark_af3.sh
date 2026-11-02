#!/usr/bin/env bash
# AF3 benchmark: positive controls, decoys, yeast case, and length scaling.
# One queue per GPU, run in parallel. Inputs are regenerated into /tmp/pe_big/cofold/inputs.
# Usage: benchmark_af3.sh [queue ...]   (queues: q0 q1 q2 q3; default all)
set -uo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"  # repo root
PY=${AF3_ENV:-/tmp/af3_env}/bin/python
MK="$PY $P/scripts/cofold/make_af3_input.py"
RUN=$P/scripts/cofold/run_af3.sh
IN=/tmp/pe_big/cofold/inputs; mkdir -p $IN $P/logs/cofold
CFF='SMILES:Cn1cnc2c1c(=O)n(C)c(=O)n2C'        # caffeine (decoy)
AZM_SMI='SMILES:CC(=O)Nc1nnc(s1)S(N)(=O)=O'      # acetazolamide via SMILES

job() {  # job <gpu> <name> <make args...>
  local g=$1 n=$2; shift 2
  $MK --name "$n" --out $IN/$n.json "$@" || { echo "make failed: $n" >&2; return; }
  $RUN -i $IN/$n.json -j bench/$n -g $g -m grow > $P/logs/cofold/bench_$n.out 2>&1
  echo "$(date -Is) done $n rc=$? $(tail -1 $P/results/cofold/bench/$n/run_af3.log)"
}
q0() { job 0 cdk2_atp_mg       --acc P24941 --ligand CCD:ATP --ligand CCD:MG
       job 0 cdk2_atp_mg_notmpl --acc P24941 --ligand CCD:ATP --ligand CCD:MG --no-templates
       job 0 cdk2_decoy_azm    --acc P24941 --ligand CCD:AZM
       job 0 cdk2_decoy_cff    --acc P24941 --ligand "$CFF"; }
q1() { job 1 ca2_zn_azm        --acc P00918 --ligand CCD:ZN --ligand CCD:AZM
       job 1 ca2_zn_azmsmiles  --acc P00918 --ligand CCD:ZN --ligand "$AZM_SMI"
       job 1 ca2_zn_decoy_atp  --acc P00918 --ligand CCD:ZN --ligand CCD:ATP
       job 1 ca2_zn_decoy_cff  --acc P00918 --ligand CCD:ZN --ligand "$CFF"; }
q2() { job 2 yADH1_nad_zn      --acc P00330 --ligand CCD:NAD --ligand CCD:ZN --ligand CCD:ZN
       job 2 yADH1_decoy_cff   --acc P00330 --ligand CCD:ZN --ligand CCD:ZN --ligand "$CFF"
       job 2 scale_cdk2x4_atp  --acc P24941 --copies 4 --ligand CCD:ATP --seeds 1,2
       job 2 scale_cdk2x5_atp  --acc P24941 --copies 5 --ligand CCD:ATP --seeds 1,2; }
q3() { job 3 pygl_plp          --acc P06737 --ligand CCD:PLP
       job 3 scale_pyglx2_plp  --acc P06737 --copies 2 --ligand CCD:PLP --ligand CCD:PLP --seeds 1,2
       job 3 scale_pyglx3_plp  --acc P06737 --copies 3 --ligand CCD:PLP --ligand CCD:PLP --ligand CCD:PLP --seeds 1,2; }
QS=${*:-q0 q1 q2 q3}
for q in $QS; do $q > $P/logs/cofold/bench_$q.log 2>&1 & done
wait
echo all done
