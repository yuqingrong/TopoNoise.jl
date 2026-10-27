# Paper theme without grids

`paper-no-grid.mplstyle` follows the supplied scientific-figure reference:
white background, thin boxed axes, blue/orange/green/pink series, compact
typography, circular markers and capped error bars supplied by the plotting
code. Gridlines are disabled for both main axes and insets.

For any Matplotlib script:

```python
with plt.style.context('/path/to/TopoNoise.jl/examples/styles/paper-no-grid.mplstyle'):
    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=uncertainty, marker='o')
```

Scripts in `examples` can also import `paper_theme`, `style_axes`, and
`save_figure` from `plotting_theme`. Call `style_axes(ax, inset=True)` for
compact inset labels. `save_figure(fig, stem)` exports PNG, PDF, and SVG.
Use a consistent size-to-color mapping across related figures. Color alone
does not indicate whether a fit is accepted.

The style controls appearance only. Axes ranges, uncertainty definitions,
fit-window shading, statistical labels, and candidate versus accepted
threshold status must be chosen explicitly by the calling script. A shaded
fit window must not be described as a threshold confidence interval.
