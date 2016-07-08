import pickle

import numpy as np
import scipy

__all__ = ['_normalize_eigensystem', '_solve_msm_eigensystem', 'verbosedump',
           'verboseload', 'check_random_state']


def _normalize_eigensystem(u, lv, rv):
    """Normalize the eigenvectors of a reversible Markov state model according
    to our preferred scheme.
    """
    # first normalize the stationary distribution separately
    lv[:, 0] = lv[:, 0] / np.sum(lv[:, 0])

    for i in range(1, lv.shape[1]):
        # the remaining left eigenvectors to satisfy
        # <\phi_i, \phi_i>_{\mu^{-1}} = 1
        lv[:, i] = lv[:, i] / np.sqrt(np.dot(lv[:, i], lv[:, i] / lv[:, 0]))

    for i in range(rv.shape[1]):
        # the right eigenvectors to satisfy <\phi_i, \psi_j> = \delta_{ij}
        rv[:, i] = rv[:, i] / np.dot(lv[:, i], rv[:, i])

    return u, lv, rv


def _solve_msm_eigensystem(transmat, k):
    """Find the dominant eigenpairs of an MSM transition matrix
    Parameters
    ----------
    transmat : np.ndarray, shape=(n_states, n_states)
        The transition matrix
    k : int
        The number of eigenpairs to find.
    Notes
    -----
    Normalize the left (:math:`\phi`) and right (:math:``\psi``) eigenfunctions
    according to the following criteria.
      * The first left eigenvector, \phi_1, _is_ the stationary
        distribution, and thus should be normalized to sum to 1.
      * The left-right eigenpairs should be biorthonormal:
        <\phi_i, \psi_j> = \delta_{ij}
      * The left eigenvectors should satisfy
        <\phi_i, \phi_i>_{\mu^{-1}} = 1
      * The right eigenvectors should satisfy <\psi_i, \psi_i>_{\mu} = 1
    Returns
    -------
    eigvals : np.ndarray, shape=(k,)
        The largest `k` eigenvalues
    lv : np.ndarray, shape=(n_states, k)
        The normalized left eigenvectors (:math:`\phi`) of ``transmat``
    rv :  np.ndarray, shape=(n_states, k)
        The normalized right eigenvectors (:math:`\psi`) of ``transmat``
    """
    u, lv, rv = scipy.linalg.eig(transmat, left=True, right=True)
    order = np.argsort(-np.real(u))
    u = np.real_if_close(u[order[:k]])
    lv = np.real_if_close(lv[:, order[:k]])
    rv = np.real_if_close(rv[:, order[:k]])
    return _normalize_eigensystem(u, lv, rv)


def verbosedump(value, fn):
    """Verbose wrapper around dump"""
    print('Saving "%s"... (%s)' % (fn, type(value)))

    with open(fn, 'wb') as f:
        pickle.dump(value, f)


def verboseload(fn):
    """Verbose wrapper around pickle.load
    """
    print('loading "%s"...' % fn)

    with open(fn, 'rb') as f:
        return pickle.load(f)


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance
    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (int, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
