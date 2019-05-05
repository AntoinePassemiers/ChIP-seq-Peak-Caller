# -*- coding: utf-8 -*-
# signal.py
# author : Antoine Passemiers

import numpy as np


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
    means[wsize:] = means[wsize:] - means[:-wsize]
    means[wsize-1:] /= wsize
    if not overlapping:
        offset = (wsize // 2) if wsize < len(signal) else 0
        means = means[offset::wsize]
    return means
