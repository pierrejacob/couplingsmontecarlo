# coding: utf8
import numpy as np
from scipy.stats import gamma
from utils import check_random_state

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)
seed = 0
rng = check_random_state(seed)

p = gamma(4, 2)
q = gamma(2, 1)

nsamples = int(1e4)
x_q = q.rvs(size=nsamples, random_state=rng)
log_w = p.logpdf(x_q) - q.logpdf(x_q)

w = np.exp(log_w - np.max(log_w))
w /= np.sum(w)  # normalized weights

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 8))

ax1.set_title(r'$q$')
pts = np.linspace(0.0, np.max(x_q), 100)
ax1.plot(pts, q.pdf(pts), c=q_col)
ax1.hist(x_q, density=True, bins=20,
         color=q_col)

ax2.set_title(r'mimic $p$ samples, using reweighted $q$ samples')
pts = np.linspace(0.0, np.max(x_q), 100)
ax2.plot(pts, p.pdf(pts), c=p_col)
ax2.hist(x_q, density=True, weights=w, bins=20,
         color=p_col)

plt.tight_layout()
plt.savefig('plots/importance_sampling_two_gammas.pdf')
