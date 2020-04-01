import numpy as np
from scipy.stats import beta
from utils import check_random_state

from config_plot import set_backend
set_backend()

import matplotlib.pyplot as plt


p, q = beta(5, 5), beta(10, 1.5)
eps = 0.05
support = (0.0 - eps, 1.0 + eps)


def T_opt(x, dist_p, dist_q):
    return dist_q.ppf(dist_p.cdf(x))


seed = 0
rng = check_random_state(seed)
nsamples = 10000

X = np.sort(p.rvs(size=nsamples, random_state=rng))
Y = T_opt(X, p, q)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(2, 7), height_ratios=(7, 2),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.1, hspace=0.1)

ax = fig.add_subplot(gs[0, 1])
ax.set_xlim(support)
ax.set_ylim(support)

pts = np.linspace(*support, 200)
p_pdf = p.pdf(pts)
q_pdf = q.pdf(pts)

ax.scatter(X, np.zeros_like(X), s=1, c='C0', alpha=0.1)
ax.scatter(np.zeros_like(Y), Y, s=1, c='C1', alpha=0.1)
ax.plot(pts, T_opt(pts, p, q),
        c='k', label=r'$F_Y^{-1}\circ F_X$')

plt.legend(loc='center right')

ax_histx = fig.add_subplot(gs[1, 1], sharex=ax)
ax_histx.tick_params(axis="x", bottom=False, labelbottom=False)
ax_histx.invert_yaxis()

ax_histx.plot(pts, p_pdf,
              c='black')
ax_histx.hist(X, density=True, bins=50,
              color='C0', alpha=0.7)

ax_histy = fig.add_subplot(gs[0, 0], sharey=ax)
ax_histy.tick_params(axis="y", left=False, labelleft=False)
ax_histy.invert_xaxis()

ax_histy.plot(q_pdf, pts,
              c='black')
ax_histy.hist(Y, density=True, bins=50,
              color='C1', alpha=0.7, orientation='horizontal')

plt.savefig('plots/common_random_number_coupling_two_betas.pdf')
