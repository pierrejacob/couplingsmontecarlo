# coding: utf8
import numpy as np
from scipy.stats import multivariate_normal

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt
from couplings_core import optimal_transport_two_gaussians

mu_1 = np.zeros(2)
cov_1 = np.array([[2.5, 0.8],
                  [0.8, 0.5]])
mu_2 = 3.0 * np.ones(2)
cov_2 = np.array([[2, 0.25],
                  [0.25, 0.5]])

nsamples = int(1e4)

X, Y = optimal_transport_two_gaussians(mu_1, cov_1,
                                       mu_2, cov_2,
                                       nsamples)
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(X[:, 0], X[:, 1], s=0.1)
ax.scatter(Y[:, 0], Y[:, 1], s=0.1)

min_x = min(np.min(X[:, 0]), np.min(Y[:, 0]))
max_x = max(np.max(X[:, 0]), np.max(Y[:, 0]))
x_mesh = np.linspace(min_x - 1, max_x + 1, 200)

min_y = min(np.min(X[:, 1]), np.min(Y[:, 1]))
max_y = max(np.max(X[:, 1]), np.max(Y[:, 1]))
y_mesh = np.linspace(min_y - 1, max_y + 1, 200)

X_mesh, Y_mesh = np.meshgrid(x_mesh, y_mesh)
levels = np.array([0.05, 0.1, 0.2, 0.5, 0.9])

for m, cov in zip((mu_1, mu_2), (cov_1, cov_2)):
    target = multivariate_normal(m, cov)
    Z = target.pdf(np.array((X_mesh.ravel(),
                             Y_mesh.ravel())).T).reshape(X_mesh.shape)
    ax.contour(X_mesh, Y_mesh, Z,
               levels=np.max(Z) * levels)

plt.savefig('plots/optimal_transport_two_gaussians_2D.pdf')
