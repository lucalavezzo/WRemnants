"""One save entry point for every SCETlib-NP plot.

Wraps the two ``wums`` output helpers so each plot is written in BOTH formats and
with a provenance sidecar, instead of a bare ``fig.savefig`` (single format, no
log):

* :func:`wums.plot_tools.save_pdf_and_png` -> ``{basename}.pdf`` + ``{basename}.png``
* :func:`wums.output_tools.write_index_and_log` -> ``{basename}.log`` (the exact
  command line, the parsed args, a timestamp, and the git hash/diff) plus the
  ``index.php`` web-gallery template. The gallery globs ``*.png`` and links each
  plot to its same-basename ``.log``/``.pdf``, so on the webdir "click a plot ->
  read the command that made it" just works.

Import is cheap (only ``os`` at module load); ``wums`` is imported lazily inside
:func:`save_plot`, so packaged modules can import this without paying for it.
"""

import os


def save_plot(outdir, basename, fig=None, args=None, meta_info=None, dpi=None):
    """Write ``fig`` to ``{outdir}/{basename}.{pdf,png}`` + a provenance log/index.

    Parameters
    ----------
    outdir : str
        Output directory (created if missing). Falsy -> current directory.
    basename : str
        Filename stem, WITHOUT extension.
    fig : matplotlib Figure, optional
        Figure to save; if ``None`` the current pyplot figure is used.
    args : argparse.Namespace, optional
        The script's parsed args, recorded in the ``.log``. The command line is
        captured from ``sys.argv`` regardless, so ``None`` (e.g. env-var-driven
        scripts) still logs the invocation.
    meta_info : dict, optional
        Extra ``key -> value`` entries appended to the ``.log`` (e.g. config that
        does not live in ``args``, like env-var settings or fit inputs).
    dpi : int, optional
        PNG resolution. The PDF is vector, so this only affects the raster
        output. ``None`` keeps the matplotlib default.
    """
    from wums import output_tools, plot_tools

    outdir = outdir or "."
    os.makedirs(outdir, exist_ok=True)
    if dpi is not None and fig is not None:
        fig.set_dpi(dpi)
    plot_tools.save_pdf_and_png(outdir, basename, fig=fig)
    output_tools.write_index_and_log(
        outdir,
        basename,
        analysis_meta_info=meta_info or {},
        args=args,
    )


def split_outpath(out_path, default_name=None):
    """Split a user ``--out``-style path into ``(outdir, basename)``.

    Accepts a file path (``dir/name.png``), a bare name, or -- when
    ``default_name`` is given -- a directory (trailing slash or no suffix), in
    which case ``default_name`` is appended. The returned ``basename`` has no
    extension, ready for :func:`save_plot`.
    """
    if default_name is not None and (
        out_path.endswith(("/", os.sep))
        or os.path.isdir(out_path)
        or not os.path.splitext(out_path)[1]
    ):
        out_path = os.path.join(out_path, default_name)
    outdir = os.path.dirname(out_path) or "."
    basename = os.path.splitext(os.path.basename(out_path))[0]
    return outdir, basename
