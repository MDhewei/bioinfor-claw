"""
Shared publication-quality plot style for all bioinfor-claw skills.

Usage — add these two lines near the top of every plotting script:

    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '_shared'))
    from plot_style import init_style

Then, after parsing args (or at the start of your main function):

    init_style()                          # defaults: Arial, 12pt base, 300 DPI
    init_style(font_size=14)              # override base font size
    init_style(font_family='Helvetica')   # override font family

If the script exposes --font-family / --font-size as CLI flags:

    init_style(
        font_family=getattr(args, 'font_family', None),
        font_size=getattr(args, 'font_size', None),
    )

The module also provides:
  - `PALETTE` — a colorblind-safe qualitative palette
  - `save_fig()` — a wrapper around savefig with consistent defaults
  - `smart_labels()` — automatic repositioning of overlapping text labels
"""

from __future__ import annotations

import sys
import warnings
from typing import Optional, Sequence

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


# ──────────────────────────────────────────────────────────────────────────────
# Font resolution: prefer Arial → Helvetica → Liberation Sans → DejaVu Sans
# ──────────────────────────────────────────────────────────────────────────────
_FONT_PREFERENCE = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']


def _best_available_font(requested: Optional[str] = None) -> str:
    """Return the best available font family from the preference chain.
    If `requested` is given and available, use it; otherwise walk the chain."""
    available = {f.name for f in fm.fontManager.ttflist}
    if requested and requested in available:
        return requested
    if requested:
        # Warn but don't crash — fall through to preference chain
        warnings.warn(
            f"[plot_style] Requested font '{requested}' not found on this system. "
            f"Falling back through preference chain: {_FONT_PREFERENCE}",
            stacklevel=3,
        )
    for font in _FONT_PREFERENCE:
        if font in available:
            return font
    return 'sans-serif'  # matplotlib ultimate fallback


# ──────────────────────────────────────────────────────────────────────────────
# Default values — "compact, publication-ready, no further refinement needed"
# ──────────────────────────────────────────────────────────────────────────────
DEFAULT_FONT_SIZE   = 12      # base text size (pt)
DEFAULT_DPI         = 300
DEFAULT_LINE_WIDTH  = 0.8     # axis spine / tick width
DEFAULT_TICK_WIDTH  = 0.6
DEFAULT_TICK_LENGTH = 3.5

# Colorblind-safe qualitative palette (8 colors from Okabe-Ito + Paul Tol)
PALETTE = [
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermilion
    '#CC79A7',  # pink
    '#999999',  # grey
]


def init_style(
    font_family: Optional[str] = None,
    font_size:   Optional[float] = None,
    dpi:         Optional[int] = None,
) -> str:
    """Apply publication-quality matplotlib rcParams globally.

    Parameters
    ----------
    font_family : str or None
        Font family name.  None → auto-resolve via _FONT_PREFERENCE.
    font_size : float or None
        Base font size in pt.  None → DEFAULT_FONT_SIZE (12).
    dpi : int or None
        Figure DPI.  None → DEFAULT_DPI (300).

    Returns
    -------
    str  — the resolved font family name (useful for logging).
    """
    family = _best_available_font(font_family)
    size   = float(font_size) if font_size else DEFAULT_FONT_SIZE
    _dpi   = int(dpi) if dpi else DEFAULT_DPI

    # ── Font ──────────────────────────────────────────────────────────────
    plt.rcParams.update({
        'font.family':     'sans-serif',
        'font.sans-serif': [family, 'Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size':       size,
        'axes.titlesize':  size + 2,   # slightly larger titles
        'axes.labelsize':  size + 1,   # axis labels
        'xtick.labelsize': size - 1,
        'ytick.labelsize': size - 1,
        'legend.fontsize': size - 1,
        'figure.titlesize': size + 3,
    })

    # ── Layout & spines ───────────────────────────────────────────────────
    plt.rcParams.update({
        'axes.linewidth':    DEFAULT_LINE_WIDTH,
        'xtick.major.width': DEFAULT_TICK_WIDTH,
        'ytick.major.width': DEFAULT_TICK_WIDTH,
        'xtick.major.size':  DEFAULT_TICK_LENGTH,
        'ytick.major.size':  DEFAULT_TICK_LENGTH,
        'xtick.minor.size':  DEFAULT_TICK_LENGTH * 0.6,
        'ytick.minor.size':  DEFAULT_TICK_LENGTH * 0.6,
        'xtick.direction':   'out',
        'ytick.direction':   'out',
        'axes.spines.top':   False,
        'axes.spines.right': False,
    })

    # ── Figure defaults ───────────────────────────────────────────────────
    # NOTE: We deliberately do NOT enable figure.constrained_layout.use here.
    # Many existing scripts call plt.tight_layout() or fig.tight_layout(),
    # and constrained_layout + tight_layout conflict on certain matplotlib
    # versions, causing colorbar / legend positioning to break.
    # Instead we rely on savefig.bbox='tight' for a similar cropping effect.
    plt.rcParams.update({
        'figure.dpi':        _dpi,
        'savefig.dpi':       _dpi,
        'savefig.bbox':      'tight',
        'savefig.pad_inches': 0.15,
    })

    # ── Color cycle ───────────────────────────────────────────────────────
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=PALETTE)

    # ── Math text: use same sans-serif font for math mode too ─────────
    plt.rcParams.update({
        'mathtext.fontset': 'custom',
        'mathtext.rm':      family,
        'mathtext.it':      f'{family}:italic',
        'mathtext.bf':      f'{family}:bold',
    })

    return family


# ──────────────────────────────────────────────────────────────────────────────
# Text adjustment — prevent overlapping labels
# ──────────────────────────────────────────────────────────────────────────────
_HAS_ADJUST_TEXT = False
try:
    from adjustText import adjust_text as _adjust_text
    _HAS_ADJUST_TEXT = True
except ImportError:
    pass


def smart_labels(
    texts: list,
    ax=None,
    arrowprops: Optional[dict] = None,
    **kwargs,
):
    """Reposition overlapping text labels using adjustText (if available).

    Parameters
    ----------
    texts : list of matplotlib.text.Text
        The text objects to adjust (returned by ax.text() calls).
    ax : matplotlib Axes or None
        The axes to adjust on.  None → current axes.
    arrowprops : dict or None
        Arrow style connecting labels to their original positions.
        None → thin grey arrow.
    **kwargs
        Extra keyword arguments passed to adjust_text().

    Returns
    -------
    bool — True if adjustment was applied, False if adjustText not available.

    Usage
    -----
        from plot_style import init_style, smart_labels

        init_style()
        fig, ax = plt.subplots()
        texts = [ax.text(x, y, label) for x, y, label in my_data]
        smart_labels(texts, ax=ax)
    """
    if not _HAS_ADJUST_TEXT or not texts:
        return False
    if ax is None:
        ax = plt.gca()
    if arrowprops is None:
        arrowprops = dict(arrowstyle='->', color='#999999', alpha=0.4, lw=0.6)
    _adjust_text(texts, ax=ax, arrowprops=arrowprops, **kwargs)
    return True


def save_fig(
    fig,
    path: str,
    dpi: Optional[int] = None,
    transparent: bool = False,
    close: bool = True,
):
    """Save a figure with consistent defaults, then close it.

    Parameters
    ----------
    fig : matplotlib Figure
    path : str or Path
    dpi : int or None   (None → uses rcParams savefig.dpi)
    transparent : bool  (True → transparent background, for web/slides)
    close : bool        (True → plt.close(fig) after saving)
    """
    kwargs = {'bbox_inches': 'tight', 'pad_inches': 0.15}
    if dpi is not None:
        kwargs['dpi'] = dpi
    if transparent:
        kwargs['transparent'] = True
    fig.savefig(str(path), **kwargs)
    if close:
        plt.close(fig)
    print(f"Saved: {path}")
