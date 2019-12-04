"""
Use trustregion in an optimization algorithm: Newton's method,
made globally convergent using trust-regions, with box constraints.
"""
# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import trustregion


def minimize(f, gradf, hessf, x0, xl, xu, delta0=1.0, deltamin=1e-6, gtol=1e-6, maxiter=100, verbose=False):
    """
    Solve the nonconvex optimization problem
      min_{x} f(x) subject to xl <= x <= xu
    using Newton's method, made globally convergent with trustregion.

    Simple implementation of trust-region methods, based on Algorithm 4.1 from
      Nocedal & Wright, Numerical Optimization, 2nd edn (2006)

    :param f: objective function, f : np.ndarray -> float
    :param gradf: gradient of objective, gradf : np.ndarray -> np.ndarray
    :param hessf: Hessian of objective, hessf : np.ndarray -> np.ndarray
    :param x0: starting point of solver, np.ndarray
    :param xl: lower bounds on x0, np.ndarray
    :param xu: upper bounds on x0, np.ndarray
    :param delta0: initial trust-region radius, float
    :param deltamin: final trust-region radius, float
    :param gtol: terminate when ||gradf(x)|| <= gtol, float
    :param maxiter: terminate after maxiter iterations
    :param verbose: whether to print information at each iteration, bool
    :return: solution x, number of iterations
    """
    xk = np.maximum(np.minimum(xu, x0), xl)  # current iterate, project x0 into box
    deltak = delta0
    if verbose:
        print("{0:^10}{1:^10}{2:^15}{3:^15}".format("k", "f(xk)", "||gradf(xk)||", "xk"))
        np.set_printoptions(precision=4, suppress=True)
    k = -1
    while k < maxiter:
        k += 1
        # Evaluate objective
        fk = f(xk)
        gk = gradf(xk)
        Hk = hessf(xk)
        # With box constraints, we have a new criticality measure
        # See Theorem 12.1.6 of Conn, Gould & Toint, Trust-Region Methods (2000)
        # Note: if unbounded, this reduces to np.linalg.norm(gk)
        crit_measure = abs(np.dot(gk, trustregion.solve(gk, None, 1.0, sl=xl - xk, su=xu - xk)))
        if verbose:
            print("{0:^10}{1:^10.4f}{2:^15.2e}{3:^15}".format(k, fk, crit_measure, str(xk)))
        # Check termination
        if crit_measure <= gtol or deltak <= deltamin:
            break  # quit loop
        # Step calculation
        sk = trustregion.solve(gk, Hk, deltak, sl=xl-xk, su=xu-xk)
        model_value = fk + np.dot(gk, sk) + 0.5 * np.dot(sk, Hk.dot(sk))  # mk(sk)
        rhok = (fk - f(xk + sk)) / (fk - model_value)
        # Update trust-region radius
        if rhok < 0.25:
            deltak = 0.25 * deltak
        elif rhok > 0.75 and abs(np.linalg.norm(sk) - deltak) < 1e-10:
            deltak = min(2 * deltak, 1e10)
        else:
            deltak = deltak
        # Update iterate
        if rhok > 0.01:
            xk = xk + sk
        else:
            xk = xk
    return xk, k


if __name__ == '__main__':
    print("Using Newton's method to minimize the 2D Rosenbrock function (with box constraints)")
    # See https://en.wikipedia.org/wiki/Rosenbrock_function
    f = lambda x: (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2
    gradf = lambda x: np.array([400.0 * x[0] ** 3 - 400.0 * x[0] * x[1] + 2 * x[0] - 2.0, 200.0 * (x[1] - x[0] ** 2)])
    hessf = lambda x: np.array([[1200.0 * x[0] ** 2 - 400.0 * x[1] + 2.0, -400.0 * x[0]], [-400.0 * x[0], 200.0]])
    x0 = np.array([-1.2, 1.0])  # (note: infeasible x0)
    xl = np.array([0.0, 0.0])
    xu = np.array([0.8, 0.6])
    xmin, niters = minimize(f, gradf, hessf, x0, xl, xu, maxiter=200, verbose=True)
    print("Found solution xmin = %s after %g iterations" % (str(xmin), niters))
    print("Done")
