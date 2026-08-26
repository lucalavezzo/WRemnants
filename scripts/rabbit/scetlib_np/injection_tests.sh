#!/bin/bash
# SCETlib NP ParamModel injection / closure test suite.
#
# Four blocks:
#   (1) 10 Asimov stat toys, with priors.
#   (2) 10 Asimov stat toys, no priors.
#   (3) Pseudo-data closure with 10 shifted POI start points, no priors.
#   (4) Pseudo-data closure with 10 shifted POI start points, with priors.
#
# (1)/(2) test pull distributions under data noise.
# (3) tests pure optimizer convergence from a non-trivial start.
# (4) tests how well the data can override a mis-centered prior
#     (xparamdefault shifts BOTH start and prior mean together).
#
# Run inside the singularity with setup.sh sourced. Each fit dumps its
# log + fitresults.hdf5 into a subdir of $OUT_BASE.

set -e
set -u
set -o pipefail

WREM_BASE=${WREM_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}

# ---- inputs (override via env vars) ------------------------------------
FIT_HDF5=${FIT_HDF5:-/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260528_debug_SCETlibPOIModel/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_excludeSCETlibNP/ZMassDilepton.hdf5}
BTGRID_DIR=${BTGRID_DIR:-/ceph/submit/data/user/l/lavezzo/zstuff/Z_COM13_CT18Z_N3p0LL_btgrid_fineall/}
OUT_BASE=${OUT_BASE:-/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/scetlib_np_injection_$(date +%Y%m%d_%H%M%S)}
N_TOYS=${N_TOYS:-10}
MAXITER=${MAXITER:-50}

mkdir -p "$OUT_BASE"
echo "Output base: $OUT_BASE"
echo "FIT_HDF5:    $FIT_HDF5"
echo "BTGRID_DIR:  $BTGRID_DIR"
echo "N_TOYS:      $N_TOYS"
echo "MAXITER:     $MAXITER"

# ---- common args -------------------------------------------------------
# Freeze the 5 ParamModel POUs without theorist priors so only
# lambda2, lambda4, lambda2_nu actually float. The '^scetlibNP.*'
# regex also freezes the old discrete-template NP systs so they
# don't double-count with our continuous λ POUs.
COMMON_FREEZE=(
    --freezeParameters
        lambda_inf lambda6 lambda_inf_nu lambda4_nu
        '^scetlibNP.*'
)

# ALL model knobs are --paramModel spec tokens (the model does not read
# SCETLIB_NP_* env vars anymore): inputs as key=value, Gaussian priors via
# priors=1 (rabbit#133 dropped the --paramModelPriors CLA), start shifts via
# xparam_default=name=val,... . Tokens must live INSIDE the --paramModel
# spec: --paramModel is action=append, so repeating the option would
# composite a second model instead of overriding.
PM_ARGS=(
    --paramModel wremnants.postprocessing.scetlib_np.SCETlibNPParamModel
        "btgrid_dir=$BTGRID_DIR"
)
PM_ARGS_PRIORS=("${PM_ARGS[@]}" priors=1)

COMMON_ARGS=(
    "$FIT_HDF5"
    -v 4
    "${COMMON_FREEZE[@]}"
    --noHessian --noEDM
    --minimizerMaxiter "$MAXITER"
)

# ---- runner ------------------------------------------------------------
run_fit () {
    local name="$1"; shift
    local outdir="$OUT_BASE/$name"
    mkdir -p "$outdir"
    echo
    echo "================================================================"
    echo "  $name"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  out: $outdir"
    echo "  extra: $*"
    echo "================================================================"
    local t0=$SECONDS
    rabbit_fit.py "${COMMON_ARGS[@]}" "$@" -o "$outdir/" 2>&1 | tee "$outdir/log.txt"
    echo "  wall: $((SECONDS - t0))s" | tee -a "$outdir/log.txt"
}

# ---- (1) and (2): Asimov toys ------------------------------------------
echo
echo "###### Block (1): $N_TOYS Asimov toys WITH priors ######"
run_fit "01_asimov_priors" -t "$N_TOYS" "${PM_ARGS_PRIORS[@]}"

echo
echo "###### Block (2): $N_TOYS Asimov toys WITHOUT priors ######"
run_fit "02_asimov_nopriors" -t "$N_TOYS" "${PM_ARGS[@]}"

# ---- (3) and (4): Shift scan with pseudoData ---------------------------
# Each entry: "<env_var_value> <slug>".
# Shifts are kept modest (within a few prior σ) so the integral stays
# in its physical region. Extend the list as you want broader coverage.
SHIFTS=(
    "lambda2=0.5             shift_lambda2_p010"
    "lambda2=0.3             shift_lambda2_m010"
    "lambda2=0.6             shift_lambda2_p020"
    "lambda4=0.55            shift_lambda4_p015"
    "lambda4=0.25            shift_lambda4_m015"
    "lambda2_nu=0.20         shift_lambda2nu_p005"
    "lambda2_nu=0.10         shift_lambda2nu_m005"
    "lambda2=0.55,lambda4=0.55    shift_lam2lam4_p"
    "lambda2=0.25,lambda4=0.25    shift_lam2lam4_m"
    "lambda2=0.6,lambda4=0.55,lambda2_nu=0.18   shift_multi"
)

echo
echo "###### Block (3): pseudoData closure, ${#SHIFTS[@]} shifts, NO priors ######"
for entry in "${SHIFTS[@]}"; do
    SHIFT=$(echo "$entry" | awk '{print $1}')
    SLUG=$(echo "$entry"  | awk '{print $2}')
    run_fit "03_inject_nopriors__${SLUG}" \
        -t 0 --pseudoData nominal "${PM_ARGS[@]}" "xparam_default=$SHIFT"
done

echo
echo "###### Block (4): pseudoData closure, ${#SHIFTS[@]} shifts, WITH priors ######"
for entry in "${SHIFTS[@]}"; do
    SHIFT=$(echo "$entry" | awk '{print $1}')
    SLUG=$(echo "$entry"  | awk '{print $2}')
    run_fit "04_inject_priors__${SLUG}" \
        -t 0 --pseudoData nominal "${PM_ARGS_PRIORS[@]}" "xparam_default=$SHIFT"
done

# ---- (5): multi-param random + corner-case shifts, all 8 POUs floating ----
# This block tests the optimizer over a broad swath of the 8-D ParamModel
# parameter space, with deliberate coverage of the (negative lambda4,
# small lambda_inf) corner where the tanh-based factorisation is most
# sensitive. The generator script enforces:
#   - hard physical bounds on positive denominators
#   - tanh_2 stability: lambda4 + (1/3) * lambda2^3 / lambda_inf^2 > 0
#
# All 8 ParamModel POUs are left floating. The COMMON_FREEZE entries
# above freeze 4 POUs + discrete-NP systs; we OVERRIDE that here by
# passing a fresh `--freezeParameters` AFTER COMMON_ARGS. rabbit's
# argparse uses nargs="+", so the last `--freezeParameters` wins and
# only the '^scetlibNP.*' regex (discrete NP systs) survives.

TEST5_SHIFTS_FILE="$OUT_BASE/test5_shifts.txt"

echo
echo "###### Block (5): pseudoData multi-shift, all 8 POUs floating ######"
echo "Generating test 5 shift list -> $TEST5_SHIFTS_FILE"
python3 "$WREM_BASE/scripts/rabbit/scetlib_np/gen_test5_shifts.py" \
    --output "$TEST5_SHIFTS_FILE"

N_T5=$(wc -l < "$TEST5_SHIFTS_FILE")
echo "Block (5): $N_T5 entries from $TEST5_SHIFTS_FILE"

while IFS= read -r entry; do
    [ -z "$entry" ] && continue
    SHIFT=$(echo "$entry" | awk '{print $1}')
    SLUG=$(echo "$entry"  | awk '{print $2}')
    run_fit "05_test5__${SLUG}" \
        -t 1 --pseudoData nominal "${PM_ARGS[@]}" "xparam_default=$SHIFT" \
        --freezeParameters '^scetlibNP.*'
done < "$TEST5_SHIFTS_FILE"

echo
echo "================================================================"
echo "All done. Results under $OUT_BASE"
echo "Per-test logs in <subdir>/log.txt; postfits in <subdir>/fitresults*.hdf5"
echo "================================================================"
