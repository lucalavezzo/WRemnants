"""Record the NP anchor -- the lambda values the theory correction was made at.

The scetlib_ad ParamModel returns ``sigma_gen(p) / sigma_gen(p_anchor)``, which
multiplies the card's templates. At the fit start ``p == p_anchor``, so the
ratio is exactly 1 and the yields ARE the card's templates -- whatever NP tune
those templates were built at. A cache built at a different tune than the card
therefore looks perfect in every prefit check and is wrong only in the
derivatives, which perturb around the wrong origin. That failure is invisible
without a recorded anchor to compare against, which is why this exists.

Write side (histmaker): :func:`build_lambda_central_meta` parses the runcard out
of the correction pkl and the histmaker stores it in its output metadata; rabbit
propagates that into the datacard.
Read side (fit): ``response.np_anchor_from_meta`` and the ``check_anchor`` guard
in ``param_model``.

Ported to scetlib_ad on 2026-09-08. It previously lived in the scetlib_np
package, which is NOT tracked on this branch, so a histmaker run from here
silently recorded no anchor and the guard was inert. The port deliberately does
NOT depend on a parameter registry: it records EVERY numeric key the
Nonperturbative section carries, so a new lambda cannot go unrecorded. The
reader keeps only what it knows and the guard compares only
``LAMBDA_CENTRAL_KEYS``, so extra keys are harmless.
"""

import os
import pickle

import lz4.frame

from wremnants.postprocessing.scetlib_ad.response import NP_ANCHOR_META_KEY
from wremnants.utilities import common as wrem_common

META_KEY = NP_ANCHOR_META_KEY

EFF_MODEL_KEY = "np_model"
GNU_MODEL_KEY = "np_model_nu"


def _find_nonperturbative(corr_dict):
    """[(basename, Nonperturbative dict)] for every basename in the pkl.

    A correction pkl carries several basenames (resummed SCETlib file,
    fixed-order singular file, gen hist, ...) and the resummed and singular
    files can hold DIFFERENT runcards, so keep all of them and let
    :func:`_select_resummed` choose.
    """
    meta = corr_dict.get("file_meta_data")
    if not isinstance(meta, dict):
        raise KeyError("Correction pkl has no 'file_meta_data' entry.")
    out = []
    for basename, file_meta in meta.items():
        if not isinstance(file_meta, dict):
            continue
        cfg = file_meta.get("config")
        if not isinstance(cfg, dict):
            continue
        npert = cfg.get("Nonperturbative")
        if isinstance(npert, dict):
            out.append((basename, npert))
    return out


def _select_resummed(sections):
    """The resummed prediction's runcard, which is the central NP tune.

    A scetlib_dyturbo correction is built from a resummed SCETlib file plus a
    fixed-order *singular* file subtracted in the matching, and only the
    resummed file's runcard is central. ``make_theory_corr.py`` distinguishes
    them by the ``sing`` substring in the filename.
    """
    resummed = [item for item in sections if "sing" not in item[0]]
    if len(resummed) == 1:
        return resummed[0]
    if not resummed:
        raise KeyError(
            "No resummed (non-'sing') basename carries a Nonperturbative "
            f"section; basenames seen: {[bn for bn, _ in sections]}."
        )
    raise KeyError(
        "Multiple resummed basenames carry a Nonperturbative section "
        f"({[bn for bn, _ in resummed]}); cannot pick the central runcard."
    )


def _parse_section(npert):
    """Split one Nonperturbative dict into the eff / gnu groups.

    Every numeric entry is recorded, split on the ``_nu`` suffix: the gamma_nu
    (CS) parameters carry it and the F_eff (TMD) ones do not. That reproduces
    the split the scetlib_np writer produced key-for-key, without needing its
    parameter registry to enumerate them.
    """
    eff_params, gnu_params = {}, {}
    if EFF_MODEL_KEY in npert:
        eff_params[EFF_MODEL_KEY] = npert[EFF_MODEL_KEY]
    if GNU_MODEL_KEY in npert:
        gnu_params[GNU_MODEL_KEY] = npert[GNU_MODEL_KEY]
    for key, value in npert.items():
        if key in (EFF_MODEL_KEY, GNU_MODEL_KEY):
            continue
        try:
            val = float(value)
        except (TypeError, ValueError):
            continue  # a non-numeric setting, not a lambda
        (gnu_params if key.endswith("_nu") else eff_params)[key] = val
    return eff_params, gnu_params


def extract_lambda_central(corr_dict, tag, proc):
    """``{tag, basename, eff_params, gnu_params}`` from a loaded correction pkl."""
    sections = _find_nonperturbative(corr_dict)
    if not sections:
        raise KeyError(
            f"No Nonperturbative section in correction pkl for tag={tag!r}, "
            f"proc={proc!r}."
        )
    basename, npert = _select_resummed(sections)
    eff_params, gnu_params = _parse_section(npert)
    return dict(
        tag=tag, basename=basename, eff_params=eff_params, gnu_params=gnu_params
    )


def _correction_pkl_path(tag, proc, data_dir=None):
    data_dir = data_dir if data_dir is not None else wrem_common.data_dir
    return os.path.join(data_dir, "TheoryCorrections", f"{tag}_Corr{proc}.pkl.lz4")


def build_lambda_central_meta(theory_corr_tags, procs=("Z", "W"), data_dir=None):
    """``{proc: lambda_central}`` for the histmaker's output metadata.

    Reads the CENTRAL correction pkl (``theory_corr_tags[0]``; the rest are
    pdfvars / pdfas) for each proc. Procs whose pkl is absent or carries no
    Nonperturbative section are skipped -- most analyses have no SCETlib NP
    correction. This is the only place the upstream pkl is read for the anchor.
    """
    if not theory_corr_tags:
        return {}
    tag = theory_corr_tags[0]
    out = {}
    for proc in procs:
        path = _correction_pkl_path(tag, proc, data_dir=data_dir)
        if not os.path.exists(path):
            continue
        try:
            with lz4.frame.open(path, "rb") as fh:
                corr_dict = pickle.load(fh)
            out[proc] = extract_lambda_central(corr_dict, tag, proc)
        except KeyError:
            continue  # pkl present, not an NP correction
    return out
