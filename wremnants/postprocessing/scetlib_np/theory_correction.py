"""Resolve the param model's PDF-dependent inputs from the theory correction.

The param model reconstructs a matched σ_gen from three PDF-dependent inputs:

  * ``btgrid_dir``  — resummed perturbative bT-space integrand (NP-off), on ``/scratch``.
  * ``nnlo_sing``   — fixed-order SINGULAR (SCETlib), subtracted from
  * ``dyturbo_fo``  — fixed-order FULL, to form the nonsingular σ_ns; either a DYTurbo
    ``…scetlibmatch.txt`` (``-g scetlib_dyturbo``) or an NNLOjet stitched-export base
    name (``-g scetlib_nnlojet``) — the key keeps its historical name.

Two of these — ``nnlo_sing`` and ``dyturbo_fo`` — plus the ``pdf_set`` are read straight
from the theory-correction pkl the datacard picked (:func:`read_corr_inputs`): they are
recorded in the correction's own build metadata (``meta_data["command"]`` corrFiles and
``file_meta_data[...].config.QCD.pdf_set``) and NOWHERE else, so reading them from the
correction is the only way to get them without a bigger refactor.

KNOWN, ACCEPTED RISK (for now): this reads the CorrZ pkl at fit time. If that pkl were
rebuilt AFTER the datacard was made, the fit would read the new inputs against a nominal
template baked from the old ones. The correct fix is to freeze the actual FO/singular
pieces (or σ_ns itself) into the datacard — deferred to the broader theory-corrections
refactor. Until then this is the lightweight choice, and it is at least consistent with
whatever the correction currently records.

The ``btgrid`` is the one input NOT recorded in any correction (it is a separate
resummed-integrand production ``make_theory_corr`` never sees), so a small
``pdf_set → btgrid`` map (:data:`PDF_BTGRID`) supplies it, guarded by
:func:`check_btgrid_pdf`.
"""

import os
import pickle
import re
import shlex

import lz4.frame

from wremnants.utilities import common as wrem_common
from wremnants.utilities.data_paths import getDataPath


def _btgrid(subdir):
    """Absolute btgrid dir under the scratch scetlib_np area."""
    base = getDataPath(fallback="/scratch/submit/cms/wmass/NanoAOD")
    return os.path.join(os.path.dirname(base), "scetlib_np", subdir)


# The ONE code map: the btgrid is not recorded in any correction (separate production),
# so pin it per PDF here. base.conf pdf_set is cross-checked in check_btgrid_pdf.
PDF_BTGRID = {
    "CT18ZNNLO": _btgrid("Z_COM13_CT18Z_N3p0LL_btgrid_fineall"),
    "MSHT20nnlo_as118": _btgrid("Z_COM13_MSHT20_N3p0LL_btgrid_fineall"),
    "MSHT20an3lo_as118": _btgrid("Z_COM13_MSHT20aN3LO_N4p0LL_btgrid_fineall"),
}


def _correction_pkl_path(tag, proc, data_dir=None):
    data_dir = data_dir if data_dir is not None else wrem_common.data_dir
    return os.path.join(data_dir, "TheoryCorrections", f"{tag}_Corr{proc}.pkl.lz4")


def _pdf_set_from_file_meta(corr):
    """The single ``pdf_set`` shared by the correction's SCETlib inputs."""
    seen = set()
    for _bn, fm in (corr.get("file_meta_data") or {}).items():
        cfg = (fm or {}).get("config") if isinstance(fm, dict) else None
        pdf = (
            ((cfg or {}).get("QCD") or {}).get("pdf_set")
            if isinstance(cfg, dict)
            else None
        )
        if pdf:
            seen.add(pdf)
    if len(seen) != 1:
        raise ValueError(
            f"expected exactly one pdf_set across the correction's SCETlib inputs, "
            f"got {sorted(seen)}"
        )
    return seen.pop()


def _corrfiles_from_command(command):
    """The ``-c/--corrFiles`` argument list from a stored make_theory_corr command."""
    toks = shlex.split(command or "")
    files, i = [], 0
    while i < len(toks):
        if toks[i] in ("-c", "--corrFiles"):
            i += 1
            while i < len(toks) and not toks[i].startswith("-"):
                files.append(toks[i])
                i += 1
        else:
            i += 1
    return files


def read_corr_inputs(tag, proc="Z", data_dir=None):
    """Read ``{pdf_set, nnlo_sing, dyturbo_fo}`` from the correction pkl for ``tag``.

    Opens ``{tag}_Corr{proc}.pkl.lz4`` (in wremnants-data by default), takes ``pdf_set``
    from its ``file_meta_data`` config and the singular / DYTurbo paths from its stored
    build command. See the module docstring for the accepted fit-time-read caveat.
    """
    path = _correction_pkl_path(tag, proc, data_dir)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"theory correction pkl for tag={tag!r} not found: {path}"
        )
    with lz4.frame.open(path, "rb") as f:
        corr = pickle.load(f)
    # Local import: sigma_gen pulls in TensorFlow, which this light module must not
    # require at import time (the FO-flavour predicate is shared with the σ_ns reader).
    from wremnants.postprocessing.scetlib_np.sigma_gen import is_nnlojet_export

    pdf_set = _pdf_set_from_file_meta(corr)
    nnlo_sing = dyturbo_fo = None
    for c in _corrfiles_from_command((corr.get("meta_data") or {}).get("command", "")):
        b = os.path.basename(c)
        if c.endswith(".txt") or "{scale}" in c:
            dyturbo_fo = c  # DYTurbo FO (-g scetlib_dyturbo)
        elif b.endswith(".pkl") and "sing" in b:
            nnlo_sing = c
        elif is_nnlojet_export(c):
            dyturbo_fo = c  # NNLOjet stitched-export base name (-g scetlib_nnlojet)
    if nnlo_sing is None or dyturbo_fo is None:
        raise ValueError(
            f"could not identify nnlo_sing / FO-full in the correction's build command "
            f"for tag={tag!r} (pkl {path}); corrFiles seen do not match the expected "
            f"'…sing….pkl' + ('…{{scale}}….txt' | NNLOjet stitched-export base name)."
        )
    return dict(pdf_set=pdf_set, nnlo_sing=nnlo_sing, dyturbo_fo=dyturbo_fo)


def btgrid_for(pdf_set):
    """The btgrid dir mapped to ``pdf_set`` (raises if unmapped — add to PDF_BTGRID)."""
    if pdf_set not in PDF_BTGRID:
        raise KeyError(
            f"no btgrid mapped for pdf_set={pdf_set!r}; add it to PDF_BTGRID in "
            f"theory_correction.py. Known: {sorted(PDF_BTGRID)}."
        )
    return PDF_BTGRID[pdf_set]


def _pdf_set_from_conf(conf_path):
    """Best-effort ``[QCD] pdf_set`` from a SCETlib ``base.conf`` (regex, comment-safe)."""
    try:
        with open(conf_path) as f:
            for line in f:
                m = re.match(r"\s*pdf_set\s*=\s*([^#\s]+)", line)
                if m:
                    return m.group(1).strip()
    except OSError:
        return None
    return None


def check_btgrid_pdf(btgrid_dir, pdf_set):
    """Guard: the btgrid's ``base.conf`` pdf_set must match the correction's ``pdf_set``.

    Catches a btgrid labelled for one PDF but generated with another (the historical
    pdf_set-mislabel). Missing/unreadable conf → warning; genuine mismatch → raise.
    """
    found = _pdf_set_from_conf(os.path.join(btgrid_dir, "base.conf"))
    if found is None:
        print(
            f"[theory_correction] WARNING: no pdf_set in {btgrid_dir}/base.conf; "
            f"skipping PDF-consistency guard.",
            flush=True,
        )
        return
    if found != pdf_set:
        raise ValueError(
            f"btgrid PDF mismatch: {btgrid_dir}/base.conf has pdf_set={found!r} but the "
            f"correction expects {pdf_set!r}. Fix the btgrid or PDF_BTGRID — do NOT fit "
            f"with a mislabelled grid."
        )
