# -*- coding: utf-8 -*-
# signal.py
# author : Antoine Passemiers

import numpy as np
from scipy.signal import savgol_filter


def moving_average(signal, wsize, overlapping=True):
    """Computes moving average on a signal.

    If `overlapping` is True, then the size of the input
    signal is preserved. Otherwise, only the centers
    of non-overlapping windows are returned.

    Parameters:
        signal (:obj:`np.ndarray`): Input signal
        wsize (int): Size of the running window.

    Returns:
        :obj:`np.ndarray`: Filtered signal.
    """
    means = np.cumsum(signal, dtype=np.float)
    offset = wsize // 2
    diff = means[wsize:] - means[:-wsize]
    means[offset:offset+len(diff)] = diff[:] / wsize

    # Keep borders as such
    means[:offset] = signal[:offset]
    means[offset+len(diff):] = signal[offset+len(diff):]

    # If specified, only return centers of non-overlapping
    # windows
    if not overlapping:
        offset = (wsize // 2) if wsize < len(signal) else 0
        means = means[offset::wsize]
    return means


def savitzky_golay_filter(signal, wsize):
    """Order-2 Savitzky-Golay filter.

    Parameters:
        signal (:obj:`np.ndarray`): Input signal.
        wsize (int): Filter window size.

    Returns:
        :obj:`np.ndarray`: Filtered signal.
    """
    return savgol_filter(signal, wsize, 2, deriv=0)
