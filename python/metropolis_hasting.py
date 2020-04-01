# coding: utf8
import numpy as np
from scipy.stats import norm, multivariate_normal, uniform

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

mean = np.zeros(2)
cov = np.array([[1.0, 0.8],
                [0.8, 1.0]])
target = multivariate_normal(mean, cov)

nb_iterations = int(1e4)

chain = np.zeros((nb_iterations, 2))
proposal = np.zeros((nb_iterations, 2))

chain[0] = [5.0, -5.0]
proposal[0] = chain[0]

naccepts = 0
for it in range(1, nb_iterations):
    proposal[it] = norm.rvs(loc=chain[it - 1])  # cov = identity
    u = uniform.rvs()
    if np.log(u) < target.logpdf(proposal[it]) - target.logpdf(chain[it - 1]):
        naccepts += 1
        chain[it] = proposal[it]
    else:
        chain[it] = chain[it - 1]

print("accept rate", naccepts / nb_iterations)

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(10, 7), height_ratios=(7, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax1 = fig.add_subplot(gs[0, 0])
ax1.tick_params(axis="x", labelbottom=False)
ax1.plot(chain[:100, 0], '.-',
         c=q_col, label=r'chain$_x$')
ax1.scatter(x=range(100), y=proposal[:100, 0],
            c=[p_col], label=r'proposal$_x$')

ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
ax2.tick_params(axis="x", labelbottom=False)
ax2.tick_params(axis="y", labelleft=False)
pts = np.linspace(-4, 4, 100)
ax2.plot(norm.pdf(pts), pts, c=q_col)
ax2.hist(chain[100:, 0], density=True,
         color=q_col, orientation='horizontal')

ax3 = fig.add_subplot(gs[1, 0], sharey=ax1)
ax3.plot(chain[:100, 0], '.-',
         c=q_col, label=r'chain$_y$')
ax3.scatter(x=range(100), y=proposal[:100, 1],
            c=[p_col], label=r'proposal$_y$')

ax4 = fig.add_subplot(gs[1, 1], sharey=ax3)
ax4.tick_params(axis="y", labelleft=False)
pts = np.linspace(-4, 4, 100)
ax4.plot(norm.pdf(pts), pts, c=q_col)
ax4.hist(chain[100:, 1], density=True,
         color=q_col, orientation='horizontal')

plt.savefig('plots/metropolis_1.pdf')

fig = plt.figure(figsize=(8, 8))

gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

burn_in = 0
ax.scatter(chain[burn_in:, 0], chain[burn_in:, 1], s=1)

x = np.linspace(np.min(chain[burn_in:, 0]), np.max(chain[burn_in:, 0]), 100)
y = np.linspace(np.min(chain[burn_in:, 1]), np.max(chain[burn_in:, 1]), 100)
X, Y = np.meshgrid(x, y)
Z = target.pdf(np.array((X.ravel(), Y.ravel())).T).reshape(X.shape)

ax.contour(X, Y, Z, levels=np.max(Z) * np.array([0.05, 0.1, 0.2, 0.5, 0.9]))

ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histx.tick_params(axis="x", labelbottom=False)

ax_histx.plot(x, norm.pdf(x), c='black')
ax_histx.hist(chain[burn_in:, 0], density=True,
              color=q_col)

ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
ax_histy.tick_params(axis="y", labelleft=False)

ax_histy.plot(norm.pdf(y), y, c='black')
ax_histy.hist(chain[burn_in:, 1], density=True,
              color=p_col, orientation='horizontal')

plt.savefig('plots/metropolis_2.pdf')
