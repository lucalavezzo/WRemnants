# `scetlib_np` — SCETlib non-perturbative param model

The fit-time NP model (`SCETlibNPParamModel`) that reconstructs σ_gen from a
cached SCETlib bt-grid, folds it through the response matrix R, and applies
`rnorm = σ_reco(λ)/σ_reco(λ_c)` in the rabbit fit — plus the tools to validate it.

Run everything inside the wmass singularity with the venv + `setup.sh` sourced.

## User entry points

### Validation / agreement — the two you'll normally use
| command | what it checks |
|---|---|
| `python -m …scetlib_np.validate_agreement --reference card --datacard <hdf5> [--outdir <dir>]` | model σ_reco **and** σ_gen at λ_central vs the **datacard itself** (`norm[signal]` + `N_gen`). No external inputs. |
| `python -m …scetlib_np.validate_agreement --reference histmaker --datacard <hdf5> --histmaker <hdf5> [--plot-out <path>] [--gen-histmaker <hdf5>] [--variation lambda21.0 …]` | the same vs an external **histmaker** `nominal` / gen MC, plus the reco λ-**variations** (`Corr[var]/Corr[pdf0]`). |
| `python -m …scetlib_np.sigma_gen_at_lambda --theory-corr <pkl.lz4> [--datacard <hdf5>] [--lambdas lambda2=0.5 …] [--fitresult <hdf5>] [--plot <path>]` | σ_gen at an **arbitrary λ tune** (cardless, gen-only) vs the official **TheoryCorrection**. Distinct: no datacard/response needed, any λ. |

### Analysis / inspection
| command | what it does |
|---|---|
| `python -m …scetlib_np.fitresult_lambdas <fitresults.hdf5> …` | read fitted NP λ out of a rabbit fitresults (table / form-factor curves / toys). |
| `python -m …scetlib_np.np_function_plots …` | plot the NP form factors: CS γ_ν^NP(b_T) and TMD F_eff(b_T, y). |
| `python -m …scetlib_np.point_to_binned …` | convert a POINT-spectrum SCETlib pickle → a binned `{hist: hist.Hist}`. |
| `python -m …scetlib_np.lambda_central <hdf5>` | inspect the central NP (λ) tune carried by a datacard/correction. |
| `python -m …scetlib_np.response_matrix <hdf5>` | load / inspect the (reco × gen) response matrix R. |

### Developer validation, smoke & timing (`…scetlib_np.validation.<name>`)
Not everyday tools — deeper cross-checks of the bt-grid factorization:
`native_validation` (native-binning vs SCETlib spectrum-mode ref),
`resum_validation` (resummed σ_gen from two sources), `export_spectrum`
(export σ onto a SCETlib run's grid), `gen_level_smoke` (gen_level=1 fold-free),
`factorized_parity` (legacy vs factorized reconstruction), `truth_start_grid`
(random-truth POI recovery), `damping_wall_dispatch` (NPDampingWall dispatch),
`timing` (one-σ_gen cost).

## Library modules (imported, not run)
`param_model.py` (the rabbit fit model) · `sigma_gen.py` (datacard-free σ_gen core) ·
`btgrid_{cache,integrate,tf}.py` (bt-grid load + Hankel/TF integration) ·
`response_matrix.py` (R) · `params.py`, `lambda_central.py` (axes / central λ) ·
`np_damping_wall.py` (fit-time physical-NP regularizer) ·
`validation/agreement.py` (shared reference loaders/aligners for the CLIs) ·
`param_model_diagnostics.py` (**Layer 0**: the numpy in-fit `run_reco_guard`, the
postfit pathology detectors, and `run_card_diagnostics` — the impl behind
`validate_agreement --reference card`) · `validation_plots.py`, `plot_output.py`
(shared plotting/array helpers + the one save entry point).
