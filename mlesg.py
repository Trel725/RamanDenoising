import numpy as np
from scipy.signal import savgol_filter
import numba

const1 = np.array([-0.00000000000020511637148964194179,
                   0.00000000017456738606819003079,
                   -0.00000006108370350381461555,
                   0.000011326343084825261969,
                   -0.001196766893248055099,
                   0.072278932710701945807,
                   -2.3738144356931862866,
                   36.543962848162728108])

const2 = np.array([0.000000049662053296461794872,
                   -0.000045735657830153646274,
                   0.013867582814278392803,
                   -1.7581457324626106331,
                   88.680082559339525284])


def get_minmax(snr):
    if snr <= 200:
        p = const1
        minm = np.round(np.polyval(p, snr))
        p2 = const2
        maxm = np.round(np.polyval(p2, snr))
    else:
        minm = 2
        maxm = 5
    return minm, maxm


def calculate_m(y, wvnumbers, peak_locs, g_sigma,
                snr=None,
                minm=None,
                maxm=None,
                mu=None):
    g_sigma = g_sigma * (wvnumbers[1] - wvnumbers[0])

    if mu is None:
        mu = 0

    if snr is not None:
        minm, maxm = get_minmax(snr)

    else:
        filtered = savgol_filter(y - mu, 9, 3)
        noise = y - mu - filtered
        snr = (y - mu).max() / noise.std()
        minm, maxm = get_minmax(snr)

    x = y - mu
    n = len(wvnumbers)
    G = np.zeros((peak_locs.size, n))

    for i in range(peak_locs.size):
        nom = -(wvnumbers - peak_locs[i])**2
        denom = 2 * g_sigma**2
        G[i, :] = 1 - np.exp(nom / denom)

    print(G.shape)
    m = G.min(axis=0)
    return m * (maxm - minm) + minm


@numba.njit
def main_job(x, n, m, j, xe, xdash, sigma, p, lmbd):
    for i in range(n):
        if m[i] > j:
            mle_range = np.arange(xe[i] - 3 * sigma,
                                  xe[i] + 3 * sigma,
                                  6 * sigma / 100)
            le_mle_range = np.zeros(mle_range.size)
            for k in range(mle_range.size):
                tmp1 = np.abs(mle_range[k] - xdash[i])**p
                limit1 = lmbd * tmp1
                tmp2_nom = (x[i] - mle_range[k])**2
                tmp2_denom = 2 * sigma**2
                limit2 = tmp2_nom / tmp2_denom
                le_mle_range[k] = limit1 + limit2

            b = le_mle_range.argmin()
            xe[i] = mle_range[b]
    return xe


def MLESG_core(y, m, v=5, q=7, lmbd=1.8, p=0.4, mu=None):
    if mu is None:
        mu = 0

    m = np.round(m).astype(int)
    x = y - mu
    n = x.size
    temp = savgol_filter(x, 9, 3)
    noise = x - temp
    sigma = noise.std()

    for j in range(m.max()):
        if j < min(m):
            v = 5
            q = 7
        elif (j >= m.max() - round(m.max() / 5)):
            q += 4
            lmbd *= 10
        else:
            v = 3
            q = 5

        if j == 0:
            xe = x
        xdash = savgol_filter(xe, q, v)

        xe = main_job(x, n, m, j, xe, xdash, sigma, p, lmbd)

    return xe


def MLESG(y, wvnumbers, peak_locs, g_sigma=10):
    m = calculate_m(y, wvnumbers, peak_locs, g_sigma)
    xe = MLESG_core(y, m)
    return xe
