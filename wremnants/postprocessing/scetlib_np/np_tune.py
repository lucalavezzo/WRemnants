"""Resolve :class:`~params.NPTune` objects from their sources, for every CLI.

The single ingestor. :mod:`params` owns the tune TYPE (forms + λ, validated
together); this module owns where a tune comes FROM — a fitresults' postfit λ, a
card's propagated λ_central, a theory-correction runcard, explicit ``--lambdas``
— and the argparse surface that selects among them. It is kept out of
:mod:`params` because the readers pull in hdf5/pickle/rabbit, which the
lightweight tools must not pay for.

Three tunes, and the distinction matters:

  * **base** — the card's λ_central. What ``SigmaGenModel`` is CONSTRUCTED at:
    the σ_gen(λ_c) normalization reference and the positivity guard. Its forms
    are **card-locked** (:func:`resolve_base_tune` never applies ``--np-model``),
    mirroring ``param_model``, where the R denominator must stay
    histmaker-consistent. Stamping a CLI form onto the base is what used to
    produce a tune whose λ dict lacked that form's parameters.
  * **num** / **den** — the tunes actually EVALUATED. Each is the base plus its
    own λ, with its own optional form override; unset λ therefore agree between
    the two sides by construction.

Per-side argparse via :func:`add_tune_args`, which emits ``--num-*`` / ``--den-*``
and nothing else. The historical unprefixed (``--fitresult``, ``--lambdas``,
``--np-model``, ``--np-model-nu``, ``--result``, ``--label``) and ``--ratio-*``
(``--ratio-to``, ``--ratio-lambdas``, ``--ratio-result``, ``--ratio-label``)
spellings were kept as aliases for a while and were REMOVED on 2026-08-07 (Luca:
one spelling per flag). Consequence to know about: a command recorded in an older
``.log`` sidecar or study logbook in those spellings will no longer re-run —
prefix its tune flags with ``--num-`` (or ``--den-`` for the ``--ratio-*`` ones).

A ``--*-fitresult`` may name the run DIRECTORY instead of the hdf5 inside it
(:func:`resolve_fitresult_path`); every CLI normalizes its own paths through
:func:`normalize_fitresult_args` right after parsing, so nothing downstream ever
sees a directory.
"""

import os

from wremnants.postprocessing.scetlib_np.params import (
    ALL_PARAMS,
    EFF_MODEL_KEY,
    GNU_MODEL_KEY,
    NPTune,
    parse_lambda_overrides,
)

_SIDE_WORD = {"num": "numerator", "den": "denominator"}


def add_tune_args(parser, side="num", *, group=None, models=True, label=True):
    """Add one side's tune flags (``--{side}-fitresult`` / ``-lambdas`` / …).

    ``side`` is ``"num"`` or ``"den"``. One spelling per flag — the historical
    unprefixed / ``--ratio-*`` aliases are gone (see the module docstring).
    Returns the argument group. ``models=False`` omits the form overrides (for a
    CLI that must not let the form move); ``label=False`` omits the legend label.
    """
    if side not in _SIDE_WORD:
        raise ValueError(f"add_tune_args: side must be num/den, got {side!r}")
    word = _SIDE_WORD[side]
    g = group or parser.add_argument_group(f"{word} λ tune")

    def add(name, *args, **kw):
        g.add_argument(f"--{side}-{name}", *args, **kw)

    add(
        "fitresult",
        dest=f"{side}_fitresult",
        default=None,
        help=f"fitresults hdf5 whose POSTFIT λ define the {word} tune",
    )
    add(
        "result",
        dest=f"{side}_result",
        default=None,
        help=f"result group suffix for the {word} fitresults",
    )
    add(
        "lambdas",
        dest=f"{side}_lambdas",
        default=None,
        help=f"explicit λ 'name=val,...' for the {word}, applied on top of its "
        "fitresults (or of the base tune when there is none). A λ the chosen "
        "forms ignore is a hard error, not a silent no-op",
    )
    if models:
        add(
            "np-model",
            dest=f"{side}_np_model",
            default=None,
            help=f"F_eff form to EVALUATE the {word} at (default: the fit's "
            "numerator form, else the card form). Never changes the "
            "construction/denominator form. λ the new form needs are filled "
            "from the registry defaults and reported",
        )
        add(
            "np-model-nu",
            dest=f"{side}_np_model_nu",
            default=None,
            help=f"γ_ν^NP form to EVALUATE the {word} at (see --{side}-np-model)",
        )
    if label:
        add(
            "label",
            dest=f"{side}_label",
            default=None,
            help=f"legend label for the {word} curve",
        )
    return g


def _get(args, side, name, default=None):
    return getattr(args, f"{side}_{name}", default) or default


# The file a rabbit fit writes into its output directory. A run directory holds
# this plus the follow-up passes' own subdirectories (``cov/``, ``saturated/``),
# so the top level is unambiguously "the fit".
FITRESULT_NAME = "fitresults.hdf5"


def resolve_fitresult_path(path, what="--fitresult"):
    """Let a fitresults flag name the run DIRECTORY rather than the hdf5 in it.

    ``…/<run>`` -> ``…/<run>/fitresults.hdf5``, announced on stdout. Anything that
    is not a directory (a file, ``None``) is returned untouched, so this is safe to
    apply to every path unconditionally.

    Deliberately resolves to the FIT pass at the top level, never to
    ``cov/fitresults.hdf5``: the covariance pass is a separate run whose λ are only
    the fit's by construction (``--externalPostfit --noFit``), and silently reading
    a different file than the one asked for is how a stale cov pass would get
    plotted as the fit. Point the flag at ``<run>/cov`` to get that one.
    """
    if not path or not os.path.isdir(path):
        return path
    cand = os.path.join(path, FITRESULT_NAME)
    if not os.path.exists(cand):
        raise ValueError(
            f"{what}: {path} is a directory but has no {FITRESULT_NAME} in it"
        )
    print(f"[np_tune] {what}: {path} is a directory -> {cand}")
    return cand


def normalize_fitresult_args(args, sides=("num", "den"), extra=("meta_from",)):
    """Rewrite the fitresults paths on ``args`` in place, run dir -> hdf5.

    Call once right after ``parse_args``, before anything reads a path. ``extra``
    names non-side path attributes to normalize too (``--meta-from``, whose λ come
    from the same kind of file). Raises ``ValueError`` (catch it and ``p.error``) on
    a directory with no fitresults in it."""
    for side in sides:
        attr = f"{side}_fitresult"
        if getattr(args, attr, None):
            setattr(
                args,
                attr,
                resolve_fitresult_path(getattr(args, attr), f"--{side}-fitresult"),
            )
    for attr in extra:
        if getattr(args, attr, None):
            setattr(
                args,
                attr,
                resolve_fitresult_path(
                    getattr(args, attr), "--" + attr.replace("_", "-")
                ),
            )
    return args


def side_requested(args, side):
    """True if anything selects a tune for ``side`` (so a ratio was asked for)."""
    return any(
        _get(args, side, n) for n in ("fitresult", "lambdas", "np_model", "np_model_nu")
    )


def resolve_base_tune(args, *, what="base tune"):
    """The CONSTRUCTION tune: the card λ_central, with its own card-locked forms.

    Priority (unchanged from the old ``resolve_base_lambda``): ``--meta-from``
    hdf5 metadata > the ``--theory-corr`` file's Nonperturbative runcard > the
    canonical FranksVals default. Deliberately ignores every ``--*-np-model``:
    the construction point is the card's own prediction, and a form the card did
    not use has no λ here. Raises if the card itself names a form SCETlib cannot
    produce (``*_sigmoid``), which would mean a corrupt card.
    """
    from wremnants.postprocessing.scetlib_np.validation.agreement import (
        resolve_base_lambda,
    )

    lc = resolve_base_lambda(args)
    tune = NPTune.create(
        lc["eff_params"][EFF_MODEL_KEY],
        lc["gnu_params"][GNU_MODEL_KEY],
        base={
            k: v
            for k, v in {**lc["eff_params"], **lc["gnu_params"]}.items()
            if k in ALL_PARAMS
        },
        what=what,
    )
    if tune.np_model.endswith("_sigmoid"):
        raise ValueError(
            f"{what}: card form np_model={tune.np_model!r} is a turn-off variant, "
            "which SCETlib cannot produce — it is an evaluation-only form. A card "
            "claiming it is corrupt."
        )
    if tune.filled:
        print(
            f"[np_tune] {what}: λ {', '.join(tune.filled)} absent from the card, "
            "taken from the registry defaults"
        )
    return tune


def resolve_side_tune(args, side, base, *, which="postfit", inherit=None, quiet=False):
    """The EVALUATION tune for one side: base λ, then its fitresults, then its λ.

    λ precedence, lowest first: the base tune, the fitresults' ``which`` values,
    the side's ``--lambdas``. Note λ fall back to the BASE, never to the other
    side — that is what makes ``--num-lambdas lambda_inf=2 --den-lambdas
    lambda_inf=1`` a clean two-point scan rather than a delta on the numerator.

    Forms: the side's explicit ``--{side}-np-model(-nu)``, else its own
    fitresults' numerator forms (``np_model_(nu_)fit``), else ``inherit``'s forms
    (pass the numerator tune when resolving the denominator), else the base's.
    Forms inherit across sides where λ do not, so that overriding ONE form on the
    denominator leaves the other where the numerator had it. Without that, a
    ``--den-np-model`` against a numerator whose forms came from a fitresults
    would silently drop the *other* form back to the card's.

    Returns ``(tune, provenance)``; ``provenance`` is a short human string naming
    where the forms and λ came from, which callers print so a plot is never
    ambiguous about what it drew.
    """
    fr = _get(args, side, "fitresult")
    src = []
    values = {}
    eff_form = gnu_form = None

    if fr:
        from wremnants.postprocessing.scetlib_np import lambda_central as lc
        from wremnants.postprocessing.scetlib_np.fitresult_lambdas import (
            _flat_values,
        )

        values.update(_flat_values(fr, which=which, result=_get(args, side, "result")))
        eff_form, gnu_form = lc.read_np_models(fr)
        src.append(f"{which} λ from {fr}")

    lam = parse_lambda_overrides(_get(args, side, "lambdas"), what=f"--{side}-lambdas")
    if lam:
        values.update(lam)
        src.append(
            "--%s-lambdas %s" % (side, ",".join(f"{k}={v:g}" for k, v in lam.items()))
        )

    fallback = inherit or base
    np_model = _get(args, side, "np_model") or eff_form or fallback.np_model
    np_model_nu = _get(args, side, "np_model_nu") or gnu_form or fallback.np_model_nu
    if _get(args, side, "np_model") or _get(args, side, "np_model_nu"):
        src.append(f"forms overridden -> {np_model}/{np_model_nu}")
    elif not fr and inherit is not None:
        src.append(f"forms inherited -> {np_model}/{np_model_nu}")
    if not src:
        src.append("base tune only")

    tune = NPTune.create(
        np_model,
        np_model_nu,
        values,
        base=base.values,
        what=f"--{side}-lambdas",
    )
    prov = f"[{_SIDE_WORD[side]}] " + "; ".join(src)
    if tune.filled:
        prov += (
            f"; λ {', '.join(tune.filled)} not supplied by that source -> "
            "registry defaults"
        )
    if not quiet:
        print(prov)
        print(f"    -> {tune.describe()}")
    return tune, prov


def tune_diff(a, b):
    """``{name: (a, b)}`` for the λ / forms that DIFFER between two tunes.

    What identifies two curves in a legend: both sides are the same base plus
    their own overrides, so spelling out the full tune says nothing.
    """
    out = {}
    for k in (EFF_MODEL_KEY, GNU_MODEL_KEY):
        va = a.np_model if k == EFF_MODEL_KEY else a.np_model_nu
        vb = b.np_model if k == EFF_MODEL_KEY else b.np_model_nu
        if va != vb:
            out[k] = (va, vb)
    for p in ALL_PARAMS:
        va, vb = a.values.get(p), b.values.get(p)
        if va != vb:
            out[p] = (va, vb)
    return out


def tune_labels(a, b, fallback="base tune", max_show=3):
    """``(label_a, label_b)`` naming only what differs between the two tunes."""
    d = tune_diff(a, b)
    if not d:
        return fallback, fallback

    def fmt(v):
        return v if isinstance(v, str) else f"{v:g}" if v is not None else "-"

    keys = list(d)
    tail = f" (+{len(keys) - max_show} more)" if len(keys) > max_show else ""
    return (
        ", ".join(f"{k}={fmt(d[k][0])}" for k in keys[:max_show]) + tail,
        ", ".join(f"{k}={fmt(d[k][1])}" for k in keys[:max_show]) + tail,
    )


__all__ = [
    "add_tune_args",
    "normalize_fitresult_args",
    "resolve_base_tune",
    "resolve_fitresult_path",
    "resolve_side_tune",
    "side_requested",
    "tune_diff",
    "tune_labels",
]
