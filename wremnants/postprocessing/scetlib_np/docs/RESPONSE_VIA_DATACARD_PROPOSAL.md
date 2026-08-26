# Proposal: carry the SCETlib-NP response matrix R through the setupRabbit datacard

**Status:** design agreed and decisions finalized (2026-06-22); not yet implemented.
The rabbit-side worktree is set up (see "Implementation plan → A"). Pick-up doc for a
fresh session.
**Context:** WRemnants PR #701 (`scetlib-np-param-model`). `/home/submit/lavezzo/alphaS/main/WRemnants`
is checked out on that branch. This is review comment **#6**.

---

## Goal

Today the `SCETlibNPParamModel` reads the reco×gen response matrix **R** (and the
gen-total **N_gen**) from a *separate* file — historically a hand-made skim in
`wremnants-data` (`.../scetlib_np/mz_dilepton_unfolding_R_skim.hdf5`), an obvious
staleness/freshness risk.

We want R to live **inside the setupRabbit datacard** (the fit input the
ParamModel already receives as `indata`), so the fit is self-contained and R is
always consistent with the run that produced the datacard — the same philosophy
as the λ_central metadata propagation already done in this PR
(`histmaker_tools._add_scetlib_np_lambda_central` → `lambda_central`).

## Decisions (agreed with Luca)

1. **Extraction happens in setupRabbit; the ParamModel reads ONLY the datacard.**
   - setupRabbit runs the existing extraction (`response_matrix.load_R`) on its
     input histmaker output and writes the *finished* arrays into the datacard.
   - The ParamModel reads R **only** from `indata` (the datacard). **No optional
     external-file override** — this was considered and explicitly walked back,
     because supporting an external source would force the extraction logic to
     live at fit time again. One source, one path.
   - ⇒ `unfolding_hdf5_path` is **removed** from the ParamModel (see "Cleanup").
2. **Storage = a generic `auxiliary` group in the datacard**, name `"scetlib_np"`,
   holding the response arrays as named datasets (R, N_gen, axis names/edges).
   This mirrors how `external_terms` is stored. The arrays ARE "x/y/z, all those
   things" — just self-described keys rather than literal `x/y/z`.
3. **`N_gen` is REQUIRED everywhere — the σ_gen-proxy fallback is removed.**
   Previously, when the gen-total hist was absent, the ParamModel substituted
   `σ_gen(λ_central)` for N_gen (`param_model.py` ≈ 762/886). That path is
   **circular** — the central reco prediction collapses to `σ_reco(λ_c) = Σ_g R_raw`
   (the row-sum of R), so the absolute-normalization closure becomes trivially true
   and tests nothing. Nothing relies on it (no tests/callers feed a gen-total-less
   input; the reco-validation closure ran on the true N_gen, UL component). So:
   - **setupRabbit raises** at datacard build if `load_R` returns `N_gen is None`
     (do not write R without N_gen).
   - **param_model raises** if `indata.auxiliary["scetlib_np"]` has no `N_gen`
     (defensive guard against a hand-made card).
   - ⇒ The auxiliary contract is **R + N_gen + axes, all required** — no optional
     `N_gen` branch in the reader, no degraded mode anywhere.
   - `response_matrix.load_R` itself stays permissive (it may still return
     `N_gen=None`); enforcement lives in the callers. Cleaner layering than baking
     policy into the extractor. The unrelated `safe_N_gen` per-bin zero guard
     (`param_model.py` ≈ 888) stays — it handles individual empty gen bins, a
     different concern.
4. **Fail loud on the target analysis.** For the Z-dilepton setup the ParamModel
   targets, if `load_R` fails (wrong sample key, missing axes, gen-shape mismatch,
   missing N_gen), setupRabbit **raises** at datacard-build time rather than
   silently skipping — so the problem surfaces immediately, not later as a confusing
   "datacard has no scetlib_np auxiliary" fit-time error. Silent skip is reserved
   only for analyses that are *not* the target (no response hist present at all).
5. **Multi-`--inputFile` guard.** setupRabbit can run over several input files
   (e.g. the wwidth / decorMassWidth multi-file branch). Extract from whichever
   input contains the response hist; **raise if more than one does** (avoids the
   `add_auxiliary` dup-name throw). Extraction happens **once, after the per-file
   loop, before `writer.write(...)`** — not per loop iteration.

## Mechanism — mirror `external_terms` (the precedent)

`external_terms` is the existing template for "stash extra structured arrays in
the fit hdf5 and expose them on `FitInputData`". Verified present on the chosen
base (upstream `origin/main`, see A):

- **Write** (rabbit `rabbit/tensorwriter.py`):
  - `self.external_terms = []` in `__init__` (line 74).
  - `add_external_likelihood_term(...)` (line 1326) appends to `self.external_terms`.
  - At `write()` (lines 2147-2173): `f.create_group("external_terms")`, one
    subgroup per term, numeric arrays via `h5pyutils_write.writeFlatInChunks`,
    the `params` name list via a vlen-str dataset.
- **Read** (rabbit `rabbit/inputdata.py`):
  - `FitInputData.__init__` line 125: `self.metadata = pickle_load_h5py(f["meta"])`
    (this is how the ParamModel already gets λ_central), and
  - line 200: `self.external_terms = read_external_terms_from_h5(f.get("external_terms"))`
    (`read_external_terms_from_h5` lives in `rabbit/rabbit/external_likelihood.py`
    line 34).

We replicate this for `auxiliary`.

**Shape is handled for free.** `writeFlatInChunks` already stores the array's
`original_shape` as an HDF5 attr (`rabbit/h5pyutils_write.py:59`) and `maketensor`
recovers it on read (`rabbit/h5pyutils_read.py:6-27`). So a multi-dim R (and the
1-D edge arrays) round-trip losslessly **as long as the reader uses `maketensor`** —
there is no manual flatten/reshape convention to copy. (This was an earlier "risk";
it is a non-issue.) Store R/N_gen as **float64** so the ~0.14% closure is preserved
bit-for-bit.

---

## Implementation plan

### A. rabbit (submodule — separate PR + submodule bump in #701)

`rabbit` is the submodule at `./rabbit`, currently on branch `luca-dev` with an
**uncommitted, in-progress** composite-param-model change (`param_models/*`). To
avoid colliding with that effort, the auxiliary work lives in an **isolated
worktree**:

```
# already created:
/tmp/rabbit-aux   branch `auxiliary-group`   based off upstream origin/main
                  (WMass/rabbit), HEAD 2299ee9 (incl. merged PR #143)
```

Base rationale: the auxiliary feature touches only `tensorwriter.py`,
`external_likelihood.py`, `inputdata.py` — **zero file overlap** with the in-progress
`param_models/*` edits, and `external_terms` (the template) already exists on
upstream `origin/main`, so it is a clean PR base. NB the "local `main` is 238 behind"
git message is a **stale local ref**: the real distance from the #701 pin (`6bfa4ff`)
to `origin/main` is only **23 behind / 11 ahead**. The 11 fork-only commits (incl.
`b2dd17f "changes for NP scetlib model"`, `843fbbd` gaussian-constraints-on-params)
are NOT on `origin/main` and are needed by the #701 ParamModel *fit* — so:
- the rabbit-side **unit test** (round-trip below) runs standalone in `/tmp/rabbit-aux`;
- the **end-to-end fit test** (B/C) runs from the WRemnants tree against an
  *integration* rabbit (merge `auxiliary-group` into the #701 rabbit line, then bump
  the submodule pointer).

Changes (in `/tmp/rabbit-aux`):

1. **`rabbit/rabbit/tensorwriter.py`**
   - Add `self.auxiliary = []` in `__init__` (next to `self.external_terms`, line 74).
   - Add method, mirroring `add_external_likelihood_term`:
     ```python
     def add_auxiliary(self, name, datasets):
         """Store a named bundle of arbitrary arrays in the output, under the
         'auxiliary' hdf5 group. `datasets` is a dict {str: np.ndarray|str|list}.
         Not used by the fit itself; available to ParamModels via
         FitInputData.auxiliary[name]."""
         if any(a["name"] == name for a in self.auxiliary):
             raise RuntimeError(f"auxiliary '{name}' already added")
         self.auxiliary.append({"name": name, "datasets": dict(datasets)})
     ```
   - In `write()`, after the `external_terms` block (≈ line 2173), add an
     `auxiliary` group, mirroring it. Per-key dispatch:
     - float/numeric arrays (R, N_gen, edges) → `writeFlatInChunks` (shape recovered
       on read via `maketensor`'s `original_shape`);
     - string lists (axis names) → vlen-str dataset like `external_terms`'s `params`
       (lines 2151-2157).
2. **`rabbit/rabbit/external_likelihood.py`** (or a new small module): add
   `read_auxiliary_from_h5(group)` returning `{name: {key: ndarray|list[str]}}`,
   mirroring `read_external_terms_from_h5`. Decode numeric datasets with
   `maketensor` (→ correct shape automatically), vlen-str datasets to a list of str.
   Return `{}` when `group is None`.
3. **`rabbit/rabbit/inputdata.py`**: in `FitInputData.__init__`, add
   `self.auxiliary = read_auxiliary_from_h5(f.get("auxiliary"))` (next to line 200).
   Default to `{}` when the group is absent (back-compatible read).
4. **Unit test**: `add_auxiliary` → `write` → `read_auxiliary_from_h5` round-trip,
   asserting a multi-dim float64 array and a str list survive bit-for-bit
   (incl. shape).
5. rabbit linting: `isort`/`black`/`flake8`/`autoflake`/`pylint` (see rabbit
   CLAUDE.md). Run before committing.

### B. WRemnants `scripts/rabbit/setupRabbit.py`

- The datacard meta is written at **line 3589** (`meta = {...}`) and serialized at
  **line 3597** (`writer.write(...)`), **after** the per-input-file loop
  (`for i, ifile in enumerate(args.inputFile)`, line 3413). setupRabbit reads its
  input via `Datagroups(inputFile)` (line 1178).
- **Once, after the loop and before `writer.write(...)`**, identify the input that
  carries the unfolding response and stash R:
  ```python
  from wremnants.postprocessing.scetlib_np import response_matrix as fz_R

  # Find the input(s) that contain the response hist. The presence of
  # `nominal_prefsr_yieldsUnfolding` under the Z sample IS the guard signal —
  # no new CLI flag needed (mirrors the λ_central guard intent).
  resp_inputs = [f for f in args.inputFile if fz_R.has_response(f)]  # small helper, see D
  if len(resp_inputs) > 1:
      raise RuntimeError(
          f"Multiple inputs carry the SCETlib-NP response hist {resp_inputs}; "
          "expected exactly one (the Z dilepton --unfolding run).")
  if resp_inputs:
      R_info = fz_R.load_R(resp_inputs[0])      # reuse the existing extraction
      if R_info.get("N_gen") is None:           # Decision 3: N_gen is required
          raise RuntimeError(
              f"{resp_inputs[0]}: response present but gen-total (N_gen) missing; "
              "rebuild the histmaker with the prefsr xnorm hist.")
      writer.add_auxiliary("scetlib_np", {
          "R":         R_info["R"],             # ndarray (reco..., gen...), float64
          "N_gen":     R_info["N_gen"],         # ndarray (gen...), float64  (required)
          "reco_axes": [n for n, _ in R_info["reco_axes"]],   # ordered names
          "gen_axes":  [n for n, _ in R_info["gen_axes"]],
          # one edges dataset per axis (variable length):
          **{f"edges__{n}": e for n, e in R_info["reco_axes"] + R_info["gen_axes"]},
      })
  # else: not an unfolding run — leave the card without the aux group (a non-target
  # analysis). The ParamModel will raise clearly if such a card is then fed to it.
  ```
  - **Fail loud (Decision 4):** `load_R` raises `KeyError` (wrong sample key /
    missing hist / no `output`) and `ValueError` (missing axes / gen-shape mismatch).
    Do **not** wrap the `load_R`/`add_auxiliary` block in a blanket `except: pass` —
    once `has_response` says the response IS there, any extraction failure should
    propagate and abort the build. The silent path is only "no response hist at all"
    (the `if resp_inputs:` guard), which is the correct no-op for non-target analyses.
  - `load_R`'s return keys are confirmed (see D): `R, N_gen, reco_axes, gen_axes,
    reco_shape, gen_shape, source`. Only `R`, `N_gen`, and the two `(name, edges)`
    axis lists are consumed downstream; `reco_shape`/`gen_shape`/`source` need not be
    stored (the ParamModel recomputes shapes from the cropped R).

### C. WRemnants `wremnants/postprocessing/scetlib_np/param_model.py`

- **Read R from `indata.auxiliary["scetlib_np"]`** instead of a file path. Replace
  the current `fz_R.load_R(unfolding_hdf5_path)` call (line 740) with a small reader
  that reconstructs the dict shape `load_R` returns (`R`, `N_gen`, and `reco_axes` /
  `gen_axes` as ordered `[(name, edges), ...]` lists) from the auxiliary datasets.
  Reassemble the axis lists from the explicit `reco_axes`/`gen_axes` name lists +
  the `edges__<name>` datasets (do **not** rely on HDF5 key-iteration order).
- **Remove `unfolding_hdf5_path`** from the constructor signature (line 472), its
  docstring (≈ 502), and the `is None` required-check (line 581). Clear error if the
  auxiliary is missing:
  ```python
  aux = getattr(indata, "auxiliary", {}) or {}
  if "scetlib_np" not in aux:
      raise ValueError(
          "SCETlibNPParamModel: datacard has no 'scetlib_np' auxiliary (the "
          "response matrix R). Rebuild the datacard with the setupRabbit that "
          "embeds it from a mz_dilepton --unfolding input.")
  ```
- **Remove the σ_gen-proxy fallback (Decision 3):** delete the
  `else: self._N_gen_flat = None` branch (lines 761-768) and the proxy substitution
  `N_gen = self._N_gen_flat if ... else gen_flat` (line 886); require N_gen and raise
  if the auxiliary lacks it. Strip the "(If the gen-total hist is absent…)" caveat
  from the Step-4 docstring (lines 161-163) and the WARNING print.
- This **supersedes the interim edit** currently in the tree (see Cleanup).

### D. WRemnants `wremnants/postprocessing/scetlib_np/response_matrix.py`

- `load_R` stays as the single extraction implementation, now called from
  setupRabbit (B) rather than the ParamModel. No behavior change to the SUM-vs-UL
  helicity logic — it is the load-bearing subtlety. Keep `load_R` permissive about
  `N_gen=None` (callers enforce; Decision 3).
- Add a small **`has_response(path)`** helper (open the file, check the Z sample's
  `output` for `nominal_prefsr_yieldsUnfolding`) so setupRabbit can cheaply detect
  the target input without a full `load_R`. Optional small refactor so the
  ParamModel-side "rebuild dict from arrays" reader can share post-processing — not
  required.

---

## Auxiliary contents (what the ParamModel needs)

From `load_R`'s return, stored in `auxiliary["scetlib_np"]` (all required):
- **R** — reco×gen, **float64**, with the ptVGen-overflow gen bin appended.
- **N_gen** — gen-total (UL helicity), same gen binning, **float64**. Required
  (Decision 3).
- **reco_axes / gen_axes** — ordered axis-name lists.
- **edges__<name>** — one 1-D float64 edge array per reco and gen axis.

The ParamModel rebins the btgrid onto the gen edges and labels the reco output, so
it never needs the original hist. R is dense reco×gen —
`(ptll·yll·cosThetaStarll_quantile·phiStarll_quantile) × (ptVGen·absYVGen)`,
~MB scale — negligible in a ~1 GB datacard.

## Test

Histmaker output that already contains the needed hists (provided by Luca):
```
/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260411_histmaker_dilepton_unfolding/mz_dilepton.hdf5
```
Plan:
1. **rabbit unit test** (in `/tmp/rabbit-aux`, standalone): `add_auxiliary` /
   `read_auxiliary_from_h5` round-trip — multi-dim float64 array + str list survive
   bit-for-bit incl. shape.
2. Run setupRabbit on the above (or a datacard built from it) and confirm the
   `auxiliary/scetlib_np` group is written and matches `load_R` on the same file.
3. Run a ParamModel fit reading R from the datacard (against the *integration*
   rabbit, A) and check the λ_central closure is unchanged vs the old external-skim
   path — should be bit-for-bit, since `load_R` is reused and float64 round-trips
   losslessly.

## Cleanup / current working-tree state (important for the new session)

- `/home/...WRemnants` is on `scetlib-np-param-model`, tracking
  `origin/scetlib-np-param-model`. There is an **uncommitted batch** of review
  fixes #1-5 and #7 plus an **interim #6 edit** that made `unfolding_hdf5_path`
  *required* and dropped the wremnants-data skim default
  (`param_model.py` ≈ 472/502/581). **That interim edit must be reverted/replaced**
  by step C (remove the arg entirely; read from the datacard).
- `./rabbit` is on `luca-dev` with uncommitted composite-param-model edits — left
  untouched. The auxiliary work is in the `/tmp/rabbit-aux` worktree (A).
- Nothing is committed/pushed yet for the WRemnants batch — the plan was to `black` +
  commit + push once the review pass is done.
- Backups from the earlier branch consolidation: branches
  `backup/local-main-preupdate`, `backup/pr701-preupdate`, and
  `~/alphaS/HOME_WIP_BACKUP_*.tgz`.

## Resolved (was "Open questions / risks")

- **N_gen key name / presence** — confirmed: key is `N_gen`, can be `None` from
  `load_R`. Now **required** (Decision 3): both callers raise if absent; writer skips
  no key (N_gen always written for the target).
- **`writeFlatInChunks` flattening** — non-issue: `original_shape` attr +
  `maketensor` recover shape automatically. Reader just uses `maketensor`.
- **Which setupRabbit path builds the ParamModel datacard** — resolved by the guard
  signal: presence of `nominal_prefsr_yieldsUnfolding` under the Z sample (the
  `has_response` helper, D). The closure suite's `--excludeSCETlibNP` card must be
  built from a `--unfolding` run; otherwise the guard skips and the ParamModel errors
  clearly (acceptable, by design).
- **rabbit base** — upstream `origin/main` worktree (A); the "238 behind" was a stale
  local ref (real: 23 behind / 11 ahead of the pin).
- rabbit `auxiliary` is intentionally **generic** (any ParamModel can stash arrays) —
  keep it un-specialized; "scetlib_np" is just our entry name.
