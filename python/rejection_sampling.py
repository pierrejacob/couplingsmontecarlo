# coding: utf8
import numpy as np
from scipy.stats import beta, multivariate_normal, uniform
from scipy.optimize import fminbound

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

# rejection sampling
nsamples = int(1e3)

mean = np.zeros(2)
cov = np.array([[1.0, 0.8],
                [0.8, 1.0]])
gauss = multivariate_normal(mean, cov)


def accept_x_plus_y_leq_1(rv):
    # constraint: x + y <= 1
    x = rv.rvs()
    while np.sum(x) > 1:
        x = rv.rvs()
    return x


x = np.array([accept_x_plus_y_leq_1(gauss) for _ in range(nsamples)])

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(x[:, 0], x[:, 1])
ax.plot((-4, 4), (5, -3))

eps = 5e-1
box = (-4 + eps, 4 - eps)
ax.set_xlim(box)
ax.set_ylim(box)

plt.tight_layout()
plt.savefig('plots/rejection_gaussian.pdf')

# or...
# sample from a beta(3, 2) by sampling from uniforms

target = beta(3, 2)
mode = fminbound(lambda x: -target.logpdf(x), x1=0.0, x2=1.0)

nsamples = int(1e3)
u_x = uniform.rvs(loc=0, scale=1, size=nsamples)
u_y = uniform.rvs(loc=0, scale=target.pdf(mode), size=nsamples)
accepts = np.log(u_y) < target.logpdf(u_x)

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(u_x[accepts], u_y[accepts], c=q_col)
ax.scatter(u_x[~accepts], u_y[~accepts], c=p_col)

pts = np.linspace(0.0, 1.0, int(1e2))
ax.plot(pts, target.pdf(pts), c='black')

plt.tight_layout()
plt.savefig('plots/rejection_beta.pdf')
