# coding: utf8
import numpy as np
from scipy.stats import expon, uniform

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

# inverse cdf transform
target = expon(loc=0.0, scale=1.0)

nsamples = int(1e4)
u = uniform.rvs(size=nsamples)
x = -np.log(u)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

ax.scatter(x, u,
           s=10, alpha=0.3,
           c='black')

ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histx.tick_params(axis="x", labelbottom=False)

pts = np.linspace(0.0, np.max(x) + 1.0, int(1e2))
ax_histx.plot(pts, target.pdf(pts), c='black')
ax_histx.hist(x, density=True,
              color=q_col)

ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
ax_histy.tick_params(axis="y", labelleft=False)

pts = np.linspace(0.0, 1.0, int(1e2))
ax_histy.plot(uniform.pdf(pts), pts, c='black')
ax_histy.hist(u, density=True,
              color=p_col, orientation='horizontal')

plt.savefig('plots/inverse_cdf_exponential.pdf')
