from collections import OrderedDict

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def compact_etable_plot(df, title=None, size=2.2, col_wrap=8,
                        xticks=np.arange(0, 7), yticks=np.arange(-2, 3),
                        xlim=(0.0, 6.0), ylim=(-2.5, 2.5),
                        marker='', ms=0, lw=1,
                        aspect=0.7, scores=None):
    if not scores:
        scores = OrderedDict(
            fa_atr='green',
            fa_rep='red',
            fa_sol='blue',
            fa_elec='orange',
            lk_ball_wtd='purple',
            score='black',
        )
    elif isinstance(scores,str):
        colors = 'black red green blue orange purple'.split()
        tmp = OrderedDict()
        for i, k in enumerate(scores.split()):
            tmp[k] = colors[i]
        scores = tmp
    sns.set(style="ticks")
    grid = sns.FacetGrid(df, col='t', size=size, col_wrap=col_wrap, aspect=aspect)
    grid.map(plt.axhline, y=0, ls=":", c='.5')
    for score, col in list(scores.items()):
        grid.map(plt.plot, 'dist', score, marker=marker, lw=lw, ms=ms, color=col)
    grid.set(xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim)
    grid.fig.tight_layout(h_pad=1, w_pad=0)
    if title:
        grid.fig.suptitle(title)
    return grid
