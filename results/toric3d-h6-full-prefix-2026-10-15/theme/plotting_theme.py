"""Reusable styling helpers; all data selection and fitting stay with callers."""
from pathlib import Path

import matplotlib.pyplot as plt

STYLE_PATH = Path(__file__).with_name('styles') / 'paper-no-grid.mplstyle'
PALETTE = ('#0072B2', '#E69F00', '#009E73', '#CC79A7', '#7F7F7F', '#D55E00')


def paper_theme():
    """Use as a context manager so the theme does not alter unrelated plots."""
    return plt.style.context(str(STYLE_PATH))


def style_axes(axis, *, inset=False):
    """Also remove pre-existing grids, including minor grids and inset grids."""
    axis.grid(False, which='both', axis='both')
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color('0.25')
        spine.set_linewidth(0.8)
    if inset:
        axis.tick_params(labelsize=8, length=3, pad=2)
        axis.xaxis.label.set_fontsize(9)
        axis.yaxis.label.set_fontsize(9)
        axis.xaxis.labelpad = 1
        axis.yaxis.labelpad = 2
    return axis


def save_figure(figure, stem):
    """Write a raster preview and editable vector formats at the same dimensions."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in ('png', 'pdf', 'svg'):
        path = stem.with_suffix('.' + extension)
        figure.savefig(path)
        paths.append(path)
    return paths
