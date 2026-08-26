#!/usr/bin/env bash
# =============================================================================
# SCETlib NP comparison suite: NEW param model vs OLD NP nuisance variations.
#
# Six tests (run as stages; see STAGES below):
#   1. asimov_new   : Asimov fit, new SCETlibNPParamModel, + Hessian/Impacts
#                     (two passes: noHessian fit -> straight-through GN cov pass)
#   2. asimov_old   : Asimov fit, OLD scetlibNP nuisance variations, no param
#                     model (datacard remade without --excludeNuisances), 1-pass
#                     Hessian + Impacts
#   3. plot_impacts : overlay impacts on pdfAlphaS, new (#1) vs old (#2)
#   4. plot_cov     : corr/cov matrices for #1 and #2
#   5. toys         : 5 toys, t5 start, --unblind -> table (lambdas + pdfAlphaS)
#   6. toys_priors  : as #5 + Gaussian priors (priors=1 spec token) centred on t5
#
# What they answer: how the uncertainties change with the new model (Asimov,
# #1 vs #2); how the lambdas + pdfAlphaS sit in toy postfits (#5); whether
# alphaS is pulled spuriously by NP/noise (#5); and how Gaussian constraints
# shift those conclusions (#6).
#
# RUN INSIDE the wmassdev container with setup.sh sourced. Detached launch:
#   IMG=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
#   OUT=/ceph/.../alphaS/scetlib_np_suite_$(date +%Y%m%d_%H%M%S)
#   APPTAINER_BIND="/tmp,/run,/home,/work,/ceph" setsid nohup singularity exec --cleanenv "$IMG" \
#     bash -c "source <WREM>/setup.sh && OUT=$OUT bash <WREM>/scripts/rabbit/scetlib_np/compare_suite.sh all" \
#     </dev/null > "$OUT/suite.log" 2>&1 &
#
# Stage subset:  ... scetlib_np_compare_suite.sh datacard asimov_old   (etc.)
# =============================================================================
set -uo pipefail

# ---------- fixed inputs -----------------------------------------------------
WREM=/home/submit/lavezzo/alphaS/main/WRemnants
# new-model datacard (scetlibNP EXCLUDED; the param model replaces them):
FIT_NEW=/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260528_debug_SCETlibPOIModel/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_excludeSCETlibNP/ZMassDilepton.hdf5
# histmaker output, to re-make the OLD-variations datacard (scetlibNP kept):
HISTMAKER=/scratch/submit/cms/areimers/alphas/histmaker/AlphaS/Theorymodels/mz_dilepton_scetlib_dyturbo_LatticeNPLambda4Bugfix_FranksVals_CT18Z_N3p0LL_N2LO_Corr_maxFiles_m1.hdf5
BTGRID=/ceph/submit/data/user/l/lavezzo/zstuff/Z_COM13_CT18Z_N3p0LL_btgrid_fineall/
MODEL=wremnants.postprocessing.scetlib_np.SCETlibNPParamModel
# model inputs as key=value spec tokens; ALL model knobs are spec tokens
# (the model does not read SCETLIB_NP_* env vars anymore). Truth start =
# no xparam_default token.
PM_INPUTS=("btgrid_dir=$BTGRID")
FREEZE=(lambda_inf lambda6 lambda_inf_nu lambda4_nu '^scetlibNP.*')
FITVAR="ptll-yll-cosThetaStarll_quantile-phiStarll_quantile"
T5="lambda2=0.6,lambda4=0.55,lambda2_nu=0.18,delta_lambda2=0.1"
MAXITER=50

# ---------- knobs (env-overridable) ------------------------------------------
OUT=${OUT:?set OUT to the suite parent dir}
THREADS=${THREADS:-64}
NTOYS=${NTOYS:-5}
TOY_CONCURRENT=${TOY_CONCURRENT:-2}

mkdir -p "$OUT"
export OMP_NUM_THREADS="$THREADS"
export TF_NUM_INTRAOP_THREADS="$THREADS"
export TF_NUM_INTEROP_THREADS=2

log() { echo "[suite $(date +%H:%M:%S)] $*"; }

DC_NEW_DIR="$OUT/asimov_newmodel"          # #1 fit + cov/
DC_OLD_PARENT="$OUT/datacard_oldvar"       # setupRabbit output parent
ASIMOV_OLD_DIR="$OUT/asimov_oldvar"        # #2 fit
DC_NEW_PRIORS_DIR="$OUT/asimov_newmodel_priors"   # #1b: new model + Gaussian priors
PLOTS="${PLOTS:-$OUT/plots}"          # override to send plots elsewhere (e.g. public_html)
TOYS_DIR="$OUT/toys_t5_unblind"
TOYS_PRIORS_DIR="$OUT/toys_t5_unblind_priors"

# resolve the remade old-var datacard (setupRabbit names the subdir from fitvar+postfix)
oldvar_datacard() { ls "$DC_OLD_PARENT"/ZMassDilepton_*withSCETlibNP/ZMassDilepton.hdf5 2>/dev/null | head -1; }

# =============================================================================
# Stage 0: remake the OLD-variations datacard (drop --excludeNuisances)
# =============================================================================
stage_datacard() {
  log "STAGE datacard: setupRabbit (keeping scetlibNP nuisances)"
  mkdir -p "$DC_OLD_PARENT"
  python3 "$WREM/scripts/rabbit/setupRabbit.py" \
    -i "$HISTMAKER" \
    --fitvar "$FITVAR" \
    -o "$DC_OLD_PARENT/" \
    --noi alphaS --pdfUncFromCorr --npUnc LatticeNoConstraintsFranks \
    --axlim ptll 0j 44j -p 'LatticeNPLambda4Bugfix_FranksVals' --pseudoData nominal \
    --postfix withSCETlibNP
  local dc; dc=$(oldvar_datacard)
  log "datacard -> ${dc:-<NOT FOUND>}"
  [[ -n "$dc" ]]
}

# =============================================================================
# Stage 1: Asimov, NEW param model (truth start) + Hessian/Impacts (two-pass)
# =============================================================================
stage_asimov_new() {
  log "STAGE asimov_new (pass A: fit, --noHessian, truth start)"
  mkdir -p "$DC_NEW_DIR"
  rabbit_fit.py "$FIT_NEW" -v 4 \
    --paramModel "$MODEL" "${PM_INPUTS[@]}" \
    --freezeParameters "${FREEZE[@]}" \
    --minimizerMaxiter "$MAXITER" \
    -t 0 --pseudoData nominal \
    --noHessian --noEDM \
    -o "$DC_NEW_DIR/"

  log "STAGE asimov_new (pass B: straight-through GN cov + impacts)"
  cd "$WREM"
  rabbit_fit.py "$FIT_NEW" -v 4 \
      --paramModel "$MODEL" "${PM_INPUTS[@]}" hessian_straightthrough=1 hessian_gn=1 \
      --freezeParameters "${FREEZE[@]}" \
      --externalPostfit "$DC_NEW_DIR/fitresults.hdf5" \
      --externalPostfitResult nominal --noFit \
      -t 0 --pseudoData nominal \
      --doImpacts \
      -o "$DC_NEW_DIR/cov/"
}

# =============================================================================
# Stage 1b: Asimov, NEW param model + Gaussian priors on lambda (truth start)
#           + Hessian/Impacts (two-pass). Priors centred at truth (xparamdefault
#           = truth since no xparam_default shift). The priors=1 spec
#           token (rabbit#133 dropped the --paramModelPriors CLA) is passed in
#           BOTH passes so the prior penalty enters the Hessian.
# =============================================================================
stage_asimov_new_priors() {
  log "STAGE asimov_new_priors (pass A: fit, --noHessian, priors=1, truth start)"
  mkdir -p "$DC_NEW_PRIORS_DIR"
  rabbit_fit.py "$FIT_NEW" -v 4 \
    --paramModel "$MODEL" "${PM_INPUTS[@]}" priors=1 \
    --freezeParameters "${FREEZE[@]}" \
    --minimizerMaxiter "$MAXITER" \
    -t 0 --pseudoData nominal \
    --noHessian --noEDM \
    -o "$DC_NEW_PRIORS_DIR/"

  log "STAGE asimov_new_priors (pass B: straight-through GN cov + impacts, priors=1)"
  cd "$WREM"
  rabbit_fit.py "$FIT_NEW" -v 4 \
      --paramModel "$MODEL" "${PM_INPUTS[@]}" priors=1 hessian_straightthrough=1 hessian_gn=1 \
      --freezeParameters "${FREEZE[@]}" \
      --externalPostfit "$DC_NEW_PRIORS_DIR/fitresults.hdf5" \
      --externalPostfitResult nominal --noFit \
      -t 0 --pseudoData nominal \
      --doImpacts \
      -o "$DC_NEW_PRIORS_DIR/cov/"
}

# =============================================================================
# Stage 2: Asimov, OLD variations, no param model, 1-pass Hessian + Impacts
# =============================================================================
stage_asimov_old() {
  local dc; dc=$(oldvar_datacard)
  [[ -n "$dc" ]] || { log "asimov_old: old-var datacard missing; run 'datacard' first"; return 1; }
  log "STAGE asimov_old: 1-pass Asimov fit + impacts on $dc"
  mkdir -p "$ASIMOV_OLD_DIR"
  rabbit_fit.py "$dc" -v 4 \
    --minimizerMaxiter "$MAXITER" \
    -t 0 --pseudoData nominal \
    --doImpacts \
    -o "$ASIMOV_OLD_DIR/"
}

# =============================================================================
# Stage 3: overlay impacts on pdfAlphaS, new vs old   (VALIDATE CLI before run)
# =============================================================================
stage_plot_impacts() {
  # Canonical alphaS impacts options (mirrors WRemnantsHelpers/workflows/
  # pullsAndImpacts.sh) + the new-vs-old overlay via --referenceFile.
  # NB: the canonical script also runs a 2nd pass with --impactType global,
  # which needs --globalImpacts at fit time; we only have traditional
  # (--doImpacts), and --globalImpacts OOMs for the param model, so skip it.
  log "STAGE plot_impacts (canonical alphaS impacts, new vs old)"
  mkdir -p "$PLOTS"
  local new="$DC_NEW_DIR/cov/fitresults.hdf5"
  local old="$ASIMOV_OLD_DIR/fitresults.hdf5"
  local styles="$WREM/wremnants/utilities/styles/styles.py"
  rabbit_plot_pulls_and_impacts.py "$new" --result nominal \
    --referenceFile "$old" --refResult nominal \
    --name "new param model" --refName "old NP variations" \
    --poi pdfAlphaS --config "$styles" --scaleImpacts 2.0 --showNumbers \
    --oneSidedImpacts --grouping alphaS -n 50 \
    --impactTitle '<i>α</i><sub>S</sub> in 10<sup>-3</sup>' \
    --title CMS --subtitle Preliminary --otherExtensions pdf png \
    -o "$PLOTS/impacts_new_vs_old/"
}

# =============================================================================
# Stage 3b: impacts with Gaussian priors on the NEW model. Two overlays:
#   (a) new+priors  vs old variations  (both constrained -> apples-to-apples)
#   (b) new+priors  vs new unconstrained (isolates the prior's effect)
# =============================================================================
stage_plot_impacts_priors() {
  log "STAGE plot_impacts_priors (new+priors vs old, and vs new-unconstrained)"
  mkdir -p "$PLOTS"
  local newp="$DC_NEW_PRIORS_DIR/cov/fitresults.hdf5"
  local new="$DC_NEW_DIR/cov/fitresults.hdf5"
  local old="$ASIMOV_OLD_DIR/fitresults.hdf5"
  local styles="$WREM/wremnants/utilities/styles/styles.py"
  local common=(--poi pdfAlphaS --config "$styles" --scaleImpacts 2.0 --showNumbers
                --oneSidedImpacts --grouping alphaS -n 50
                --impactTitle '<i>α</i><sub>S</sub> in 10<sup>-3</sup>'
                --title CMS --subtitle Preliminary --otherExtensions pdf png)
  rabbit_plot_pulls_and_impacts.py "$newp" --result nominal \
    --referenceFile "$old" --refResult nominal \
    --name "new model + priors" --refName "old NP variations" \
    "${common[@]}" -o "$PLOTS/impacts_newpriors_vs_old/"
  rabbit_plot_pulls_and_impacts.py "$newp" --result nominal \
    --referenceFile "$new" --refResult nominal \
    --name "new model + priors" --refName "new model (free)" \
    "${common[@]}" -o "$PLOTS/impacts_newpriors_vs_newfree/"
}

# =============================================================================
# Stage 4: corr/cov matrices for each fit              (VALIDATE CLI before run)
# =============================================================================
stage_plot_cov() {
  log "STAGE plot_cov (corr + cov, each fit; param subset)"
  mkdir -p "$PLOTS/cov_newmodel" "$PLOTS/cov_oldvar"
  # --params takes EXACT param names (not regexes); an unmatched name crashes
  # rabbit_plot_cov (list[boolmask] bug). Old-var scetlibNP names discovered
  # from the fit (4 of them for this datacard).
  local new_params=(pdfAlphaS lambda2 lambda4 lambda2_nu delta_lambda2)
  local old_params=(pdfAlphaS chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 scetlibNPgammaLambda2 chargeVgenNP0scetlibNPZdelta_lambda2)
  local annot=(--showNumbers --scaleTextSize 0.5)
  for corr in "" "--correlation"; do
    rabbit_plot_cov.py "$DC_NEW_DIR/cov/fitresults.hdf5" --result nominal \
      --params "${new_params[@]}" "${annot[@]}" $corr -o "$PLOTS/cov_newmodel/"
    rabbit_plot_cov.py "$ASIMOV_OLD_DIR/fitresults.hdf5" --result nominal \
      --params "${old_params[@]}" "${annot[@]}" $corr -o "$PLOTS/cov_oldvar/"
  done
}

# =============================================================================
# Stage 5/6: toy ensembles (t5 start, unblind) -> tables
# =============================================================================
stage_toys() {
  log "STAGE toys (t5, --unblind, $NTOYS toys)"
  python3 "$WREM/scripts/rabbit/scetlib_np/toys_suite.py" \
    --outbase "$TOYS_DIR" --ntoys "$NTOYS" --max-concurrent "$TOY_CONCURRENT" \
    --threads "$THREADS" --unblind
  python3 "$WREM/scripts/rabbit/scetlib_np/toys_make_table.py" \
    --indir "$TOYS_DIR" --outdir "$TOYS_DIR" --name toys_postfit_table
}

stage_toys_priors() {
  log "STAGE toys_priors (t5, --unblind, priors=1, $NTOYS toys)"
  python3 "$WREM/scripts/rabbit/scetlib_np/toys_suite.py" \
    --outbase "$TOYS_PRIORS_DIR" --ntoys "$NTOYS" --max-concurrent "$TOY_CONCURRENT" \
    --threads "$THREADS" --unblind --priors
  python3 "$WREM/scripts/rabbit/scetlib_np/toys_make_table.py" \
    --indir "$TOYS_PRIORS_DIR" --outdir "$TOYS_PRIORS_DIR" --name toys_postfit_table
}

# ---------- dispatcher -------------------------------------------------------
ALL=(datacard asimov_new asimov_new_priors asimov_old plot_impacts plot_impacts_priors plot_cov toys toys_priors)
run_stage() {
  case "$1" in
    datacard)            stage_datacard ;;
    asimov_new)          stage_asimov_new ;;
    asimov_new_priors)   stage_asimov_new_priors ;;
    asimov_old)          stage_asimov_old ;;
    plot_impacts)        stage_plot_impacts ;;
    plot_impacts_priors) stage_plot_impacts_priors ;;
    plot_cov)            stage_plot_cov ;;
    toys)                stage_toys ;;
    toys_priors)         stage_toys_priors ;;
    *) log "unknown stage '$1'"; return 2 ;;
  esac
}

STAGES=("$@")
[[ ${#STAGES[@]} -eq 0 || "${STAGES[0]}" == "all" ]] && STAGES=("${ALL[@]}")
log "OUT=$OUT  THREADS=$THREADS  stages: ${STAGES[*]}"
for s in "${STAGES[@]}"; do
  log ">>> begin $s"
  if run_stage "$s"; then log "<<< done $s"; else log "!!! FAILED $s (continuing)"; fi
done
log "suite finished."
