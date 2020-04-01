# coding: utf8
import numpy as np
from scipy.stats import beta

from config_plot import set_backend
set_backend()
import matplotlib.pyplot as plt


p, q = beta(5, 5), beta(2, 1.5)
p_col, q_col = (1, 0.1, 0.3), (0.3, 0.1, 1)

x = np.linspace(start=0.0, stop=1.0, num=int(1e3))
p_x, q_x = p.pdf(x), q.pdf(x)

fig, (ax, ax1) = plt.subplots(2, 1,
                              figsize=(12, 8),
                              sharex=True)

# Plot P, Q, and min(P, Q)
ax.set_title(r'Densities $p - q$')
ax.plot(x, p_x, x, q_x, color='black')
ax.fill_between(x, p_x, 0,
                facecolor=p_col, interpolate=True, alpha=0.75,
                label=r'$p$')
ax.fill_between(x, q_x, 0,
                facecolor=q_col, interpolate=True, alpha=0.75,
                label=r'$q$')

ax.fill_between(x, np.minimum(p_x, q_x), 0,
                hatch='x',
                label=r'$\min(p,q)$')

ax.legend(loc='upper left')
ax.set_frame_on(False)
ax.tick_params(axis='both', labelsize=0, length=0)

# Plot P - Q
ax1.set_title(r'$p - q$')
ax1.plot(x, p_x - q_x,
         color='black',
         label=r'$p-q$')
ax1.plot(x, np.sign(p_x - q_x),
         label=r'$sign(p-q)$',
         linestyle='dashed')
ax1.fill_between(x, p_x - q_x, 0,
                 facecolor=(0.9, 0.8, 0.2), interpolate=True, alpha=1)

ax1.legend(loc='upper left')
ax1.set_frame_on(False)
ax1.tick_params(axis='both', labelsize=0, length=0)

plt.tight_layout()
plt.savefig('plots/total_variation_distance.pdf')
