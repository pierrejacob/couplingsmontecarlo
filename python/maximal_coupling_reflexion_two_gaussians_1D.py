# coding: utf8
import numpy as np
from scipy.stats import norm

from couplings_core import maximal_coupling_reflexion_two_gaussians

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


mu1, mu2 = np.array([-1]), np.array([2])
cov = np.array([0.7**2], ndmin=2)

nsamples = int(1e4)
p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

X, Y = maximal_coupling_reflexion_two_gaussians(mu1, mu2, cov, nsamples)

min_x, max_x = np.min(X), np.max(X)
min_y, max_y = np.min(Y), np.max(Y)

x = np.linspace(min_x - 1, max_x + 1, 200)
y = np.linspace(min_y - 1, max_y + 1, 200)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

ax.scatter(X, Y, s=1, c='C0')

ax.set_xlim((min_x, max_x))
ax.set_ylim((min_y, max_y))

ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

ax_histx.tick_params(axis="x", labelbottom=False)

ax_histx.plot(x,
              norm.pdf(x, loc=mu1, scale=np.sqrt(cov[0, 0])),
              c='black')
ax_histx.hist(X, density=True, bins=50,
              color=p_col)

ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

ax_histy.tick_params(axis="y", labelleft=False)

ax_histy.plot(norm.pdf(y, loc=mu2, scale=np.sqrt(cov[0, 0])),
              y, c='black')
ax_histy.hist(Y, density=True, bins=50,
              color=q_col, orientation='horizontal')

plt.savefig('plots/maximal_coupling_reflexion_two_gaussians_1D.pdf')
