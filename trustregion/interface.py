# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from ._trs import trsapp, trsbox, trslin

__all__ = ['solve']

ZERO_THRESH = 1e-14


def _to_array(X, lbl):
    try:
        return np.asarray_chkfinite(X)
    except ValueError:
        raise ValueError('%s contains Nan/Inf values' % lbl)


def check_dimension(x, dim, lbl, check_real=True):
    if check_real and not np.all(np.isreal(x)):
        raise ValueError('%s must be a real array' % (lbl))
    if len(x.shape) != dim:
        raise ValueError('%s must be %g-dimensional (%s.shape=%s)' % (lbl, dim, lbl, str(x.shape)))
    return


def solve(g, H, delta, sl=None, su=None, verbose_output=False):
    """
    Approximately solve the trust-region subproblem:
    min_{s}     g*s + 0.5*s*H*s
    subject to  ||s||_2 <= delta
                (sl <= s <= su)

    The matrix H must be symmetric and the bounds must satisfy sl <= 0 <= su.
    If the problem is linear, H can be None.

    Outputs are:
    - s: the (approximate) minimizer of the subproblem (numpy.ndarray)
    In addition, if verbose_output=True, outputs include:
    - gnew: the gradient of the objective at the new point, gnew = g + H*s (numpy.ndarray)
    - crvmin: the minimum curvature of H found (float), with convention:
        crvmin = 0 if ||s||_2 = delta
        crvmin = -1 if every step of the search hit the sl/su bounds
        crvmin > 0 is the smallest curvature of H found

    :param g: gradient of model (n-dimensional vector; as numpy.ndarray)
    :param delta: trust-region radius (must be >= 0; as float)
    :param H: Hessian of model (if not None, n*n symmetric real matrix; as numpy.ndarray)
    :param sl: lower bounds on step (n-dimensional vector; as numpy.ndarray)
    :param su: upper bounds on step (n-dimensional vector; as numpy.ndarray)
    :param verbose_output: whether to return full output or just solution s (bool)
    :return: s[, gnew, crvmin]
    """
    # Convert to desired types
    g = _to_array(g, 'g')
    H = None if H is None else _to_array(H, 'H')
    try:
        delta = float(delta)
    except:
        raise ValueError('delta must be a float')
    sl = None if sl is None else _to_array(sl, 'sl')
    su = None if su is None else _to_array(su, 'su')
    if verbose_output not in [True, False]:
        raise ValueError('verbose_output must be a bool')
    # Check inputs
    check_dimension(g, 1, 'g')
    n = len(g)
    if H is not None:
        check_dimension(H, 2, 'H')
        if H.shape != (n,n):
            raise ValueError("H must be square with same dimension as g (got H.shape = %s and len(g)=%g)" % (str(H.shape), n))
    if (sl is not None and su is None) or (sl is None and su is not None):
        raise ValueError("Must specify none or both of sl and su")
    if sl is not None:
        check_dimension(sl, 1, 'sl')
        if sl.shape != (n,):
            raise ValueError("sl must have same length as g (got len(sl) = %g and len(g) = %g)" % (len(sl), n))
    if su is not None:
        check_dimension(su, 1, 'su')
        if su.shape != (n,):
            raise ValueError("su must have same length as g (got len(su) = %g and len(g) = %g)" % (len(su), n))
    # Classify type of problem
    linear = (H is None) or (np.allclose(H, 0.0))
    box_constrained = (sl is not None)
    # Types checked, now validate values
    if delta < 0.0:
        raise ValueError("delta must be non-negative (got delta = %g)" % delta)
    if not linear:
        if not np.allclose(H, H.T):
            raise ValueError("H must be symmetric")
    if box_constrained:
        # assert np.all(su >= sl), "Must have sl <= su (in all components)"
        if np.any(sl > 0.0):
            raise ValueError("Must have sl <= 0 (in all components)")
        if np.any(su < 0.0):
            raise ValueError("Must have su >= 0 (in all components)")

    # Solve (using different routines, depending
    if delta <= ZERO_THRESH or (box_constrained and np.min(su-sl) <= ZERO_THRESH):
        # Trivial if delta=0 or sl=su
        s = np.zeros((n,))
        crvmin = 0.0 if delta <= ZERO_THRESH else -1.0
    elif linear and not box_constrained:
        # Minimise linear function inside a ball - analytic solution
        s = -g * (delta / np.linalg.norm(g))
        crvmin = 0.0
    elif linear:
        # Linear with box constraints - TRSLIN
        s, crvmin = trslin_interface(g, sl, su, delta)
    elif not box_constrained:
        # Quadratic without box constraints - TRSAPP
        # min_{s} g*(xopt+s) + 0.5*(xopt+s)*H*(xopt+s)
        # s.t. ||s||_2 <= delta
        s, crvmin = trsapp_interface(g, H, delta)
    else:
        # Quadratic with box constraints - TRSBOX
        trs_inside_box = np.max(sl) <= -delta and np.min(su) >= delta
        if trs_inside_box:
            # No way the box constraints can be active, so ignore
            s, crvmin = trsapp_interface(g, H, delta)
        else:
            s, _, crvmin = trsbox_interface(g, H, sl, su, delta)
    # Ensure box constraints are exactly satisfied (in case of small rounding errors)
    if box_constrained:
        s = np.minimum(np.maximum(sl, s), su)
    # Return relevant output
    if verbose_output:
        gnew = g.copy() if linear else g + H.dot(s)  # gradient of objective at solution
        return s, gnew, crvmin
    else:
        return s


def trslin_interface(g, sl, su, delta):
    # Basic interface to Fortran routine TRSLIN
    # min_{s} g*s
    # s.t.    ||s||_2 <= delta
    #         sl <= s <= su
    s, crvmin = trslin(g, sl, su, delta)
    return s, crvmin


def trsapp_interface(g, H, delta):
    # Basic interface to Fortran routine TRSAPP
    # min_{s} g*(xopt+s) + 0.5*(xopt+s)*H*(xopt+s)
    # s.t. ||s||_2 <= delta
    xopt = np.zeros((len(g),))
    s, crvmin = trsapp(xopt, g, H, delta)
    return s, crvmin


def trsbox_interface(g, H, sl, su, delta):
    # Basic interface to Fortran routine TRSBOX
    # min_{s} g*s + 0.5*s*H*s
    # s.t.    ||s||_2 <= delta
    #         sl <= xopt+s <= su
    xopt = np.zeros((len(g),))
    s, gnew, crvmin = trsbox(xopt, g, H, sl, su, delta)
    return s, gnew, crvmin
