# coding: utf8
import numpy as np
import scipy.linalg as la
from scipy.stats import multivariate_normal

from utils import check_random_state


def maximal_coupling_rejection(p, q, nsamples=1, random_state=None):
    """Rejection sampler to obtain pairs (X, Y) such that :math:`X \\sim p`, :math:`Y \\sim q` and the event :math:`\\{X=Y\\}` occurs with maximal probability.

    :param p:
        distribution of X

    :param q:
        distribution of Y

    Both parameter distributions ``p, q`` must have the methods

        - ``.rvs()`` to generate samples
        - ``.logpdf()`` to evaluate the log probability density function

    .. seealso::

        `scipy.stats module <https://docs.scipy.org/doc/scipy/reference/stats.html>`_
    """

    rng = check_random_state(random_state)

    for _ in range(nsamples):
        X = p.rvs(random_state=rng)
        if np.log(rng.rand()) < (q.logpdf(X) - p.logpdf(X)):
            return X, X

        Y = q.rvs(random_state=rng)
        while np.log(rng.rand()) < (p.logpdf(Y) - q.logpdf(Y)):
            Y = q.rvs(random_state=rng)

    return X, Y


def maximal_coupling_two_bernoullis(p, q, nsamples=1, random_state=None):
    """Rejection sampler to obtain pairs (X, Y) such that :math:`X \\sim p`, :math:`Y \\sim q` and the event :math:`\\{X=Y\\}` occurs with maximal probability.

    .. seealso::

        `scipy.stats module <https://docs.scipy.org/doc/scipy/reference/stats.html>`_
    """

    rng = check_random_state(random_state)

    if p == q:
        X = rng.rand(nsamples) < p
        return X, X

    ber_p, ber_q = np.array((1 - p, p)), np.array((1 - q, q))

    min_pq = np.minimum(ber_p, ber_q)
    int_min_pq = np.sum(min_pq)

    ber_p_tld = (ber_p - min_pq) / (1.0 - int_min_pq)
    ber_q_tld = (ber_q - min_pq) / (1.0 - int_min_pq)

    gam = np.diag(min_pq) + (1.0 - int_min_pq) * np.outer(ber_p_tld, ber_q_tld)

    sample = rng.choice(4, size=nsamples, p=gam.ravel())
    X = np.array([0, 0, 1, 1])
    Y = np.array([0, 1, 0, 1])

    return X[sample], Y[sample], gam


def maximal_coupling_reflexion_two_gaussians(mu_1, mu_2, cov, nsamples=1, random_state=None):

    rng = check_random_state(random_state)

    L = la.cholesky(cov, lower=True)
    z = la.solve_triangular(L, mu_1 - mu_2,
                            lower=True, unit_diagonal=False)
    e = z / la.norm(z)

    dim = len(cov)

    X = rng.randn(nsamples, dim)
    Y = X.copy()

    gauss = multivariate_normal(mean=np.zeros(dim))  # gaussian N(0, I)
    accept = np.log(rng.rand(nsamples)) < gauss.logpdf(X + z) - gauss.logpdf(X)

    Y[accept] += z
    Y[~accept] -= np.outer(X[~accept].dot(e), 2.0 * e)  # Householder reflexion

    return X.dot(L.T) + mu_1, Y.dot(L.T) + mu_2


def optimal_transport_two_gaussians(mu_1, cov_1, mu_2, cov_2, nsamples=1, random_state=None):
    """
    .. seealso::

        Remark 2.31 of Peyré and Cuturi (2020)
        http://arxiv.org/pdf/1803.00567.pdf
    """
    rng = check_random_state(random_state)

    dim = len(cov_1)

    sqrt_cov_1 = la.sqrtm(cov_1)
    inv_sqrt_cov_1 = la.inv(sqrt_cov_1)

    scaling_matrix = la.sqrtm(sqrt_cov_1.dot(cov_2).dot(sqrt_cov_1))
    scaling_matrix = inv_sqrt_cov_1.dot(scaling_matrix).dot(inv_sqrt_cov_1)

    X = rng.randn(nsamples, dim)
    X = X.dot(sqrt_cov_1) + mu_1

    Y = mu_2 + (X - mu_1).dot(scaling_matrix)

    return X, Y


def common_random_number_coupling(p, q, nsamples=1, random_state=None):
    """
    :param p:
        distribution of X with methods

        - ``.rvs()`` to generate samples from distribution :math:`p`
        - ``.cdf()`` to evaluate the cumulative distribution function

    :param q:
        distribution of Y with method

        - ``.ppf()`` to evaluate the inverse cumulative distribution function

    .. seealso::

        `scipy.stats module <https://docs.scipy.org/doc/scipy/reference/stats.html>`_
    """

    rng = check_random_state(random_state)

    X = p.rvs(size=nsamples, random_state=rng)
    Y = q.ppf(p.cdf(X))

    return X, Y
