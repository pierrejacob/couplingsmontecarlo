# coding: utf8
import numpy as np
from scipy.stats import uniform, beta
from scipy.integrate import quad
from scipy.optimize import brentq

from utils import check_random_state


class p_distribution(object):

    def __init__(self):
        self.supp = np.array([0.0, 1.0])
        self.mix_components = (0.25, 0.5, 0.25)
        self.distrib = (beta(40, 10), beta(20, 20), beta(10, 45))

    def pdf(self, x):
        return np.dot(self.mix_components, [d.pdf(x) for d in self.distrib])

    def cdf(self, x):
        return np.dot(self.mix_components, [d.cdf(x) for d in self.distrib])

    def ppf(self, y):

        y = np.array(y, ndmin=1)

        def inv_cdf(z):
            return brentq(lambda x: self.cdf(x) - z, *self.supp)

        return np.array([inv_cdf(z) for z in y])

    def rvs(self, size=1, random_state=None):
        rng = check_random_state(random_state)
        unif = uniform.rvs(size=size, random_state=rng)

        return self.ppf(unif)


class q_distribution(object):

    def __init__(self):
        self.supp = np.array([0.0, 1.0])
        self.norm_const = quad(self.unnormalized_pdf,
                               self.supp[0],
                               self.supp[1])[0]

    def unnormalized_pdf(self, x):
        return np.exp(beta.logpdf(x, 4, 4) + np.cos(5 * np.pi * x))

    def pdf(self, x):
        x = np.array(x, ndmin=1)
        y = np.zeros_like(x)
        indic = np.logical_and(self.supp[0] <= x, x <= self.supp[1])
        if np.any(indic):
            y[indic] = self.unnormalized_pdf(x[indic]) / self.norm_const
        return y

    def cdf(self, x):
        x = np.array(x, ndmin=1)
        return np.array([quad(self.pdf, a=0.0, b=y)[0] for y in x])

    def ppf(self, y):

        def inv_cdf(z):
            return brentq(lambda x: self.cdf(x) - z, *self.supp)

        return np.array([inv_cdf(z) for z in y])

    def rvs(self, size=1, random_state=None):
        rng = check_random_state(random_state)
        unif = uniform.rvs(size=size, random_state=rng)

        return self.ppf(unif)


if __name__ == '__main__':

    import pickle as pkl

    from config_plot import set_backend
    set_backend()
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    file_name = 'plots/common_random_number_coupling_mixture.pkl'

    try:
        with open(file_name, 'rb') as f:
            data = pkl.load(f)

            seed = data.get('seed')
            nsamples = data.get('nsamples')
            print('nsamples = {}'.format(nsamples))
            p = data.get('p_distribution').get('class')
            p_samples = data.get('p_distribution').get('samples')
            q = data.get('q_distribution').get('class')
            q_samples = data.get('q_distribution').get('samples')

    except FileNotFoundError:
        print('No previous samples were saved under {}'.format(file_name))
        print('Generating samples may take a while')

        p, q = p_distribution(), q_distribution()

        nsamples = int(1e4)
        seed = 0
        rng = check_random_state(seed)
        unif_01 = rng.rand(nsamples)

        p_samples, q_samples = p.ppf(unif_01), q.ppf(unif_01)
        data = {'seed': seed,
                'nsamples': nsamples,
                'p_distribution': {'class': p, 'samples': p_samples},
                'q_distribution': {'class': q, 'samples': q_samples}}

        with open(file_name, 'wb') as f:
            pkl.dump(data, f)

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(3, 1,
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[1, 0])
    step = 5
    p_samp = p_samples[::step]
    q_samp = q_samples[::step]
    lines = np.c_[q_samp,
                  np.zeros_like(q_samp),
                  p_samp,
                  np.ones_like(p_samp)].reshape(nsamples // step, 2, 2)
    ax.add_collection(LineCollection(lines, linewidths=0.01, colors='k'))
    ax.scatter(p_samples, np.zeros_like(p_samples), s=0.1)
    ax.scatter(q_samples, np.ones_like(q_samples), s=0.1)

    eps = 5e-2

    ax_hist_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_hist_top.tick_params(axis="x", labelbottom=False)
    x_pts = np.linspace(q.supp[0] - eps, q.supp[1] + eps, 200)
    ax_hist_top.plot(x_pts, q.pdf(x_pts), c='black')
    ax_hist_top.hist(q_samples, density=True, bins=50,
                     color='C1', alpha=0.7)

    ax_hist_bottom = fig.add_subplot(gs[2, 0], sharex=ax)
    x_pts = np.linspace(p.supp[0] - eps, p.supp[1] + eps, 200)
    ax_hist_bottom.plot(x_pts, p.pdf(x_pts), c='black')
    ax_hist_bottom.hist(p_samples, density=True, bins=50,
                        color='C0', alpha=0.7)
    ax_hist_bottom.invert_yaxis()

    plt.savefig('plots/common_random_number_coupling_mixture_arrangement.pdf')

    def T_opt(x, dist_p, dist_q):
        return dist_q.ppf(dist_p.cdf(x))

    X = p_samples
    Y = q_samples  # T_opt(X, p, q)

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=(2, 7), height_ratios=(7, 2),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.15, hspace=0.1)

    ax = fig.add_subplot(gs[0, 1])

    ax.scatter(X, np.zeros_like(X), s=1, c='C0', alpha=0.1)
    ax.scatter(np.zeros_like(Y), Y, s=1, c='C1', alpha=0.1)

    eps = 5e-2
    x_pts = np.linspace(p.supp[0] - eps, p.supp[1] + eps, 200)
    ax.plot(x_pts, T_opt(x_pts, p, q),
            c='k', label=r'$F_Y^{-1}\circ F_X$')
    plt.legend(loc='center right')

    ax_histx = fig.add_subplot(gs[1, 1], sharex=ax)
    ax_histx.tick_params(axis="x", bottom=False, labelbottom=False,
                         top=True, labeltop=False)
    x_pts = np.linspace(p.supp[0] - eps, p.supp[1] + eps, 200)
    ax_histx.plot(x_pts, p.pdf(x_pts),
                  c='black')
    ax_histx.hist(X, density=True, bins=50,
                  color='C0', alpha=0.7)
    ax_histx.invert_yaxis()

    ax_histy = fig.add_subplot(gs[0, 0], sharey=ax)
    ax_histy.tick_params(axis="y",
                         left=False, labelleft=False,
                         right=True, labelright=False)
    x_pts = np.linspace(q.supp[0] - eps, q.supp[1] + eps, 200)
    ax_histy.plot(q.pdf(x_pts), x_pts,
                  c='black')
    ax_histy.hist(Y, density=True, bins=50,
                  color='C1', alpha=0.7, orientation='horizontal')
    ax_histy.invert_xaxis()

    plt.savefig('plots/common_random_number_coupling_mixture.pdf')
