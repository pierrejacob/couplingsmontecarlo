# coding: utf8
import numpy as np

from couplings_core import maximal_coupling_two_bernoullis

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p, q = 0.3, 0.4
nsamples = int(1e4)
X, Y, gama = maximal_coupling_two_bernoullis(p, q, nsamples)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax.set_xlim((-0.5, 1.5))
ax.set_ylim((-0.5, 1.5))

ax.axvline(x=0.5, ymin=0, ymax=1, ls='--', c='k')
ax.axhline(y=0.5, xmin=0, xmax=1, ls='--', c='k')

ax.set_xticks(range(2))
ax.set_xlabel(r'$X$')
ax.set_yticks(range(2))
ax.set_ylabel(r'$Y$', rotation=0)

for x in range(2):
    for y in range(2):
        text = ax.text(x, y, np.round(gama[x, y], decimals=2),
                       ha="center", va="center",
                       color="k")

ax_barx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_barx.set_ylim((0, 1))
ax_barx.tick_params(axis="x", labelbottom=True)

p_hat = np.mean(X)
ax_barx.bar((0, 1), height=(1 - p_hat, p_hat),
            width=0.5, color='C0', alpha=0.7)
ax_barx.axhline(p, 0., 0.9, ls='--', c='k')

ax_bary = fig.add_subplot(gs[1, 1], sharey=ax)
ax_bary.set_xlim((0, 1))
ax_bary.tick_params(axis="y", labelleft=False)

q_hat = np.mean(Y)
ax_bary.barh((0, 1), height=0.5, width=(1 - q_hat, q_hat),
             color='C1', alpha=0.7)
ax_bary.axvline(q, 0., 0.9, ls='--', c='k')

plt.savefig('plots/maximal_coupling_mixture_two_bernoullis.pdf')
