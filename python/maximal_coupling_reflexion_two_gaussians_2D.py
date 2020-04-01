# coding: utf8
import numpy as np
from scipy.stats import norm

from couplings_core import maximal_coupling_reflexion_two_gaussians

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


mu1, mu2 = np.array([-1, 0]), np.array([1, -1])
cov = np.array([[0.8, -0.3], [-0.3, 0.2]], ndmin=2)

nsamples = int(5e4)

X, Y = maximal_coupling_reflexion_two_gaussians(mu1, mu2, cov, nsamples)

min_x = min(np.min(X[:, 0]), np.min(Y[:, 0]))
max_x = max(np.max(X[:, 0]), np.max(Y[:, 0]))

min_y = min(np.min(X[:, 1]), np.min(Y[:, 1]))
max_y = max(np.max(X[:, 1]), np.max(Y[:, 1]))


x = np.linspace(min_x - 1, max_x + 1, 200)
y = np.linspace(min_y - 1, max_y + 1, 200)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

ax.scatter(X[:, 0], X[:, 1], s=1, c='C0', alpha=0.1)
ax.scatter(Y[:, 0], Y[:, 1], s=1, c='C1', alpha=0.1)

ax.set_xlim((min_x, max_x))
ax.set_ylim((min_y, max_y))


ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

ax_histx.tick_params(axis="x", labelbottom=False)

ax_histx.plot(x,
              norm.pdf(x, loc=mu1[0], scale=np.sqrt(cov[0, 0])),
              c='black')
ax_histx.hist(X[:, 0], density=True, bins=50,
              color='C0', alpha=0.7)

ax_histx.plot(x,
              norm.pdf(x, loc=mu2[0], scale=np.sqrt(cov[0, 0])),
              c='black')
ax_histx.hist(Y[:, 0], density=True, bins=50,
              color='C1', alpha=0.7)

ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

ax_histy.tick_params(axis="y", labelleft=False)

ax_histy.plot(norm.pdf(y, loc=mu1[1], scale=np.sqrt(cov[1, 1])),
              y, c='black')
ax_histy.hist(X[:, 1], density=True, bins=50,
              color='C0', alpha=0.7, orientation='horizontal')

ax_histy.plot(norm.pdf(y, loc=mu2[1], scale=np.sqrt(cov[1, 1])),
              y, c='black')
ax_histy.hist(Y[:, 1], density=True, bins=50,
              color='C1', alpha=0.7, orientation='horizontal')

plt.savefig('plots/maximal_coupling_reflexion_two_gaussians_2D.pdf')
