"""
My Python implementation of TRSBOX from DFO-LS
"""

# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
from math import sqrt
import numpy as np
import trustregion


__all__ = ['trsbox', 'trsbox_geometry']

ZERO_THRESH = 1e-14

def sumsq(x):
    return np.dot(x,x)

def trsbox(xopt, g, hess, sl, su, delta):
    n = xopt.size
    assert xopt.shape == (n,), "xopt has wrong shape (should be vector)"
    assert g.shape == (n,), "g and xopt have incompatible sizes"
    assert hess.shape == (n,n), "hess and xopt have incompatible sizes"
    assert sl.shape == (n,), "sl and xopt have incompatible sizes"
    assert su.shape == (n,), "su and xopt have incompatible sizes"
    assert delta > 0.0, "delta must be strictly positive"
    # Assume g and hess have full quadratic model for objective
    # i.e. skip straight to label 8 in DFBOLS version

    # The sign of G(I) gives the sign of the change to the I-th variable
    # that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
    # or not to fix the I-th variable at one of its bounds initially, with
    # NACT being set to the number of fixed variables. D and GNEW are also
    # set for the first iteration. DELSQ is the upper bound on the sum of
    # squares of the free variables. QRED is the reduction in Q so far.

    iterc = 0
    nact = 0  # number of fixed variables

    xbdi = np.zeros((n,), dtype=int)  # fix x_i at bounds? [values -1, 0, 1]
    xbdi[(xopt <= sl) & (g >= 0.0)] = -1
    xbdi[(xopt >= su) & (g <= 0.0)] = 1

    d = np.zeros((n,))
    s = np.zeros((n,))
    gnew = g.copy()
    qred = 0.0
    delsq = delta ** 2
    crvmin = -1.0
    beta = 0.0  # label 20

    need_alt_trust_step = False  # will either quit main CG loop to finish, or do alternative step
    MAX_LOOP_ITERS = 100 * n ** 2  # avoid infinite loops
    # while True:  # main CG loop [label 30]
    for ii in range(MAX_LOOP_ITERS):
        # print('Iteration %g, s =' % ii, s)
        # print('xbdi =', xbdi)
        # print('gnew =', gnew)
        # print('beta =', beta)
        s[xbdi != 0] = 0.0
        if beta == 0.0:
            s[xbdi == 0] = -gnew[xbdi == 0]
        else:
            s[xbdi == 0] = beta * s[xbdi == 0] - gnew[xbdi == 0]
        stepsq = sumsq(s)
        # print('After setting')
        # print('s =', s)
        # print('stepsq =', stepsq)

        if stepsq == 0.0:
            need_alt_trust_step = False
            break  # break and quit

        if beta == 0.0:
            gredsq = stepsq
            itermax = iterc + n - nact

        if iterc == 0:
            gredsq0 = gredsq

        # Exit conditions
        if gredsq <= min(1.0e-6 * gredsq0, 1.0e-18) or gredsq * delsq <= min(1.0e-6 * qred ** 2, 1.0e-18):  # DFBOLS
            need_alt_trust_step = False
            break  # break and quit

        # Multiply the search direction by the second derivative matrix of Q and
        # calculate some scalars for the choice of steplength. Then set BLEN to
        # the length of the the step to the trust region boundary and STPLEN to
        # the steplength, ignoring the simple bounds.

        hs = hess.dot(s)
        # print('hs =', hs)

        # label 50
        ds = np.dot(s[xbdi == 0], d[xbdi == 0])
        shs = np.dot(s[xbdi == 0], hs[xbdi == 0])
        resid = delsq - sumsq(d[xbdi == 0])
        if resid <= 0.0:
            need_alt_trust_step = True
            break  # break and calculate alt step instead

        # print('shs = %.15g' % shs)
        # print('resid = %.15g' % resid)
        # print('ds = %.15g' % ds)
        # print('stepsq = %.15g' % stepsq)
        temp = sqrt(stepsq * resid + ds ** 2)
        blen = (resid / (temp + ds) if ds >= 0.0 else (temp - ds) / stepsq)
        stplen = (blen if shs <= 0.0 else min(blen, gredsq / shs))

        # Exit condition
        if stplen <= 1.0e-30:  # DFBOLS
            need_alt_trust_step = False
            break  # break and quit

        # Reduce STPLEN if necessary in order to preserve the simple bounds,
        # letting IACT be the index of the new constrained variable.
        # print('Setting iact')
        iact = None
        for i in range(n):
            if s[i] != 0.0:
                temp = (su[i] - xopt[i] - d[i] if s[i] > 0.0 else sl[i] - xopt[i] - d[i]) / s[i]
                # print('i =%g, temp = %g, stplen = %g' % (i, temp, stplen))
                if temp < stplen:
                    stplen = temp
                    iact = i

        # Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
        sdec = 0.0
        if stplen > 0.0:
            iterc += 1
            temp = shs / stepsq
            if iact is None and temp > 0.0:
                crvmin = min(crvmin, temp) if crvmin != -1.0 else temp
            ggsav = gredsq
            gnew += stplen * hs
            d += stplen * s
            gredsq = sumsq(gnew[xbdi == 0])
            sdec = max(stplen * (ggsav - 0.5 * stplen * shs), 0.0)
            qred += sdec

        # Restart the conjugate gradient method if it has hit a new bound.
        if iact is not None:
            nact += 1
            xbdi[iact] = (1 if s[iact] >= 0.0 else -1)
            delsq = delsq - d[iact] ** 2
            if delsq <= 0.0:
                need_alt_trust_step = True
                # print('Alt step: delsq = %g' % delsq)
                break  # break and calculate alt step instead
            beta = 0.0  # label 20
            continue  # restart loop (new CG iteration)

        # If STPLEN is less than BLEN, then either apply another conjugate
        # gradient iteration or RETURN.
        if stplen >= blen:
            need_alt_trust_step = True
            # print('Alt step: stplen = %g, blen = %g' % (stplen, blen))
            # print('stplen - blen = %g' % (stplen-blen))
            break  # break and calculate alt step instead

        # Exit condition
        if iterc == itermax or sdec <= 1.0e-6 * qred:  # DFBOLS
            need_alt_trust_step = False
            break  # break and quit

        beta = gredsq / ggsav
        continue  # new CG iteration
    # end of CG loop

    # either done or need to take and alternative step
    if need_alt_trust_step:
        # print('Fi?ished main loop, alt step now')
        crvmin = 0.0
        d, gnew = alt_trust_step(n, xopt, hess, sl, su, d, xbdi, nact, gnew, qred)
        return d, gnew, crvmin
    else:
        # print('Finished main loop, no alt step')
        return d_within_bounds(d, xopt, sl, su, xbdi), gnew, crvmin


# Alternative Trust Region Step (label 100 of TRSBOX in BOBYQA, where crvmin=0)
def alt_trust_step(n, xopt, hess, sl, su, d, xbdi, nact, gnew, qred):
    MAX_LOOP_ITERS = 100 * n ** 2  # avoid infinite loops
    # while True:  # label 100 here
    for ii in range(MAX_LOOP_ITERS):
        if nact >= n - 1:
            return d_within_bounds(d, xopt, sl, su, xbdi), gnew

        # Prepare for the alternative iteration by calculating some scalars
        # and by multiplying the reduced D by the second derivative matrix of
        # Q, where S holds the reduced D in the call of GGMULT.
        s = np.zeros((n,))
        s[xbdi == 0] = d[xbdi == 0]
        dredsq = sumsq(d[xbdi == 0])
        dredg = np.dot(d[xbdi == 0], gnew[xbdi == 0])
        gredsq = sumsq(gnew[xbdi == 0])

        # Label 210 (crvmin = 0, itcsav = iterc)
        hs = hess.dot(s)

        hred = hs.copy()
        # quit 210 by goto 120

        # Let the search direction S be a linear combination of the reduced D
        # and the reduced G that is orthogonal to the reduced D.
        restart_alt_loop = False  # once the below loop finishes, quit unless need to go again
        # while True:  # label 120
        for jj in range(MAX_LOOP_ITERS):
            temp = gredsq * dredsq - dredg ** 2
            if temp <= 1.0e-4 * qred ** 2:
                restart_alt_loop = False
                break  # quit inner label 120 loop and return results
            temp = sqrt(temp)
            s = np.zeros((n,))
            s[xbdi == 0] = (dredg * d[xbdi == 0] - dredsq * gnew[xbdi == 0]) / temp
            sredg = -temp

            # By considering the simple bounds on the variables, calculate an upper
            # bound on the tangent of half the angle of the alternative iteration,
            # namely ANGBD, except that, if already a free variable has reached a
            # bound, there is a branch back to label 100 after fixing that variable.
            free_variable_reached_bound = False
            angbd = 1.0
            iact = None
            for i in range(n):
                if xbdi[i] == 0:
                    tempa = xopt[i] + d[i] - sl[i]
                    tempb = su[i] - xopt[i] - d[i]
                    if tempa <= 0.0:
                        nact += 1
                        xbdi[i] = -1
                        free_variable_reached_bound = True
                        break  # skip the rest of this for loop
                    elif tempb <= 0.0:
                        nact += 1
                        xbdi[i] = 1
                        free_variable_reached_bound = True
                        break  # skip the rest of this for loop
                    ssq = d[i] ** 2 + s[i] ** 2
                    temp = ssq - (xopt[i] - sl[i]) ** 2
                    if temp > 0.0:
                        temp = sqrt(temp) - s[i]
                        if angbd * temp > tempa:
                            angbd = tempa / temp
                            iact = i
                            xsav = -1
                    temp = ssq - (su[i] - xopt[i]) ** 2
                    if temp > 0.0:
                        temp = sqrt(temp) + s[i]
                        if angbd * temp > tempb:
                            angbd = tempb / temp
                            iact = i
                            xsav = 1
            # End for loop
            if free_variable_reached_bound:  # deal with break conditions above
                restart_alt_loop = True
                break  # quit inner label 120 loop and restart alt iteration loop (label 100)

            # Label 210 (crvmin = 0, itcsav < iterc since iterc+=1 earlier)
            hs = hess.dot(s)

            # Label 150
            # Calculate HHD and some curvatures for the alternative iteration.
            shs = np.sum(s[xbdi == 0] * hs[xbdi == 0])
            dhs = np.sum(d[xbdi == 0] * hs[xbdi == 0])
            dhd = np.sum(d[xbdi == 0] * hred[xbdi == 0])

            # Seek the greatest reduction in Q for a range of equally spaced values
            # of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
            # the alternative iteration.
            redmax = 0.0
            isav = -1
            redsav = 0.0
            temp = 0.0  # force scope outside i loop below since needed later
            iu = int(17 * angbd + 3.1)
            for i in range(iu):  # i = 0, ..., iu-1
                angt = angbd * float(i + 1) / float(iu)
                sth = 2.0 * angt / (1.0 + angt ** 2)
                temp = shs + angt * (angt * dhd - 2.0 * dhs)
                rednew = sth * (angt * dredg - sredg - 0.5 * sth * temp)
                if rednew > redmax:
                    redmax = rednew
                    isav = i
                    rdprev = redsav
                elif i == isav + 1:
                    rdnext = rednew
                redsav = rednew

            # Return if the reduction is zero. Otherwise, set the sine and cosine
            # of the angle of the alternative iteration, and calculate SDEC.
            if isav == -1:
                restart_alt_loop = False
                break  # quit inner label 120 loop and return results

            if isav < iu - 1:
                temp = (rdnext - rdprev) / (2.0 * redmax - rdprev - rdnext)
                angt = angbd * (float(isav + 1) + 0.5 * temp) / float(iu)

            cth = (1.0 - angt ** 2) / (1.0 + angt ** 2)
            sth = 2.0 * angt / (1.0 + angt ** 2)
            temp = shs + angt * (angt * dhd - 2.0 * dhs)
            sdec = sth * (angt * dredg - sredg - 0.5 * sth * temp)

            if sdec <= 0.0:
                restart_alt_loop = False
                break  # quit inner label 120 loop and return results

            # Update GNEW, D and HRED. If the angle of the alternative iteration
            # is restricted by a bound on a free variable, that variable is fixed
            # at the bound.
            gnew += (cth - 1.0) * hred + sth * hs
            d[xbdi == 0] = cth * d[xbdi == 0] + sth * s[xbdi == 0]
            dredg = np.dot(d[xbdi == 0], gnew[xbdi == 0])
            gredsq = sumsq(gnew[xbdi == 0])
            hred = cth * hred + sth * hs

            qred += sdec
            if iact is not None and isav == iu - 1:
                nact += 1
                xbdi[iact] = xsav
                restart_alt_loop = True
                break  # quit inner label 120 loop and restart alt iteration loop (label 100)

            if (sdec <= 0.01 * qred):
                restart_alt_loop = False
                break  # quit inner label 120 loop and return results
            continue  # back to inner label 120 loop

        # End inner label 120 loop

        if restart_alt_loop:
            continue
        else:
            break  # end outer loop and quit

    # End while True (label 100)
    return d_within_bounds(d, xopt, sl, su, xbdi), gnew


def d_within_bounds(d, xopt, sl, su, xbdi):
    # Used in TRSBOX, force d to be within bounds
    # In Fortran code, is at label 190
    xnew = np.maximum(np.minimum(xopt + d, su), sl)
    xnew[xbdi == -1] = sl[xbdi == -1]
    xnew[xbdi == 1] = su[xbdi == 1]
    d = xnew - xopt
    return d


if __name__ == '__main__':
    # A test problem where Python got better results than Fortran (before bugfix)
    delta = 0.48
    xopt = np.array([0.230233016244113, 0.0474202712953557])
    g = np.array([-87.7045737137454, -38.6055462975271])
    sl = np.array([-0.8, -2.7]) - xopt
    su = np.array([2.2, 0.2]) - xopt
    H = np.array([[943.577752759425, 433.953396751178], [433.953396751178, 200]])

    pred = lambda s: -np.dot(g, s) - 0.5*np.dot(s, H.dot(s))

    print('==============================')
    print('Python')
    print('==============================')
    s1, gnew1, crvmin1 = trsbox(np.zeros((2,)), g, H, sl, su, delta)
    print('')
    print('s =', s1)
    print('gnew =', gnew1)
    print('crvmin = ', crvmin1)
    print('pred =', pred(s1))
    print('xnew  =', xopt + s1)
    print('')
    print('==============================')
    print('Fortran')
    print('==============================')
    s2, gnew2, crvmin2 = trustregion.solve(g, H, delta, sl=sl, su=su, verbose_output=True)
    print('')
    print('s =', s2)
    print('gnew =', gnew2)
    print('crvmin = ', crvmin2)
    print('pred =', pred(s2))
    print('xnew  =', xopt + s2)
    print('')
    print('Done')