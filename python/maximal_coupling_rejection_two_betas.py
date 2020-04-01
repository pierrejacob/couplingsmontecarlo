# coding: utf8
import numpy as np
from scipy.stats import beta

from couplings_core import maximal_coupling_rejection

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p, q = beta(5, 5), beta(2, 1.5)
p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

nb_points = int(1e4)
# r_states = [None] * nb_points  # random
r_states = range(nb_points)
res = np.array([maximal_coupling_rejection(p, q, random_state=r)
                for r in r_states])

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histy.tick_params(axis="y", labelleft=False)

eps = 0.05
support = [0.0 - eps, 1.0 + eps]
x = np.linspace(0.0, 1.0, int(1e2))

ax.scatter(res[:, 0], res[:, 1],
           s=10, alpha=0.3,
           c='black')
ax.set_xlim(support)
ax.set_ylim(support)

ax_histx.plot(x, p.pdf(x), c='black')
ax_histx.hist(res[:, 0], density=True,
              color=p_col)

ax_histy.plot(q.pdf(x), x, c='black')
ax_histy.hist(res[:, 1], density=True,
              color=q_col, orientation='horizontal')

plt.savefig('plots/maximal_coupling_rejection_two_betas.pdf')
