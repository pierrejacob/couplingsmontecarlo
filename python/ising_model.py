# coding: utf8
import numpy as np
from utils import check_random_state


# Gibbs sampling on Ising model
class IsingModel(object):
    """Naive implementation of Gibbs sampler for the 2D Ising model with joint probability proportional to

    .. math::

        \\exp( \\beta \\sum_{\\left\\{ (i, j), (k, l) \\right\\}} s_{ij} s_{kl}),

    where the sum runs over all pairs of neighoring states with periodic boundary conditions and :math:`s_{ij} \\in \\left\\{ -1, 1 \\right\\}`.

    :param beta:
        Inverse temperature parameter :math:`\\beta > 0`
    :type beta:
        float

    :param state_matrix:
        Matrix :math:`S` encoding the current states of the configuration.

        .. math::

            S_{ij} =
            \\begin{cases}
                \\text{True} & \\text{if} s_{ij} = 1\\
                \\text{False} & \\text{if} s_{ij} = -1
            \\end{cases}

    :type state_matrix:
        boolean array
    """

    def __init__(self, beta, state_matrix):
        self.beta = beta
        self.state_matrix = state_matrix
        self.shape = state_matrix.shape

        sum_ngbr_state = np.array([-4, -2, 0, 2, 4])
        self.cond_proba = 1.0 / (1.0 + np.exp(-2 * beta * sum_ngbr_state))

    def neighbors(self, i, j):
        """Indices of neighbors with periodic boundary conditions"""
        n_rows, n_cols = self.shape
        idx_row = (i, i, (i + n_rows - 1) % n_rows, (i + 1) % n_rows)
        # i, i, i - 1, (i + 1) % n_rows
        idx_col = ((j + n_cols - 1) % n_cols, (j + 1) % n_cols, j, j)
        # j - 1, (j + 1) % n_cols, j, j
        return idx_row, idx_col

    def sample_gibbs_sweep(self, nb_sweeps=1, random_sate=None):

        rng = check_random_state(random_sate)
        n_rows, n_cols = self.shape

        for _ in range(nb_sweeps):
            for i in range(n_rows):
                for j in range(n_cols):
                    sum_ngbr_01 = sum(self.state_matrix[self.neighbors(i, j)])
                    cond_proba = self.cond_proba[sum_ngbr_01]
                    self.state_matrix[i, j] = rng.rand() < cond_proba

        return self.state_matrix, rng


if __name__ == '__main__':

    from config_plot import set_backend
    set_backend()
    import matplotlib.pyplot as plt

    seed = 0
    rng = check_random_state(seed)

    beta = 0.42
    n_rows, n_cols = 30, 30
    s_mat = rng.rand(n_rows, n_cols) < 0.5

    ising = IsingModel(beta, s_mat)

    fig, axes = plt.subplots(nrows=1, ncols=3,
                             figsize=(15, 5),
                             sharex=False, sharey=False)

    nb_sweeps = (0, 2, 100)

    for swp, ax in zip(nb_sweeps, axes):
        s_mat, rng = ising.sample_gibbs_sweep(swp, random_sate=rng)
        ax.imshow(ising.state_matrix, cmap=plt.cm.gray)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    plt.tight_layout()
    plt.savefig('plots/ising_model.pdf')
