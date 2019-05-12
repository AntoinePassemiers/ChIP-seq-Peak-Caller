# -*- coding: utf-8 -*-
# peak_calling.py
# author : Antoine Passemiers

from chipseq.chromosome import CalledPeaks
from chipseq.signal import moving_average, savitzky_golay_filter

import numpy as np
from scipy.stats import poisson


def call_peaks(chromosome, wsize=100, alpha=1e-5,
               apply_filter=False):
    """Call peaks from a given chromosome.

    Parameters:
        chromosome (:obj:´chipseq.Chromosome`):
            Input chromosome.
        wsize (int): Number of bases separating
            two tags.
        alpha (float): Significance threshold.
        apply_filter (bool): Whether to apply a
            filter on the signal beforehand.

    Returns:
        :obj:`chipseq.CalledPeaks`: Called peaks.
    """
    signal = chromosome.signal
    if apply_filter:
        # Apply prior Savitzky-Golay filter for
        # smoothing the signal
        filter_wsize = int(np.round(5000. / wsize))
        if filter_wsize % 2 == 0:
            # Filter's window size has to be an odd number
            filter_wsize += 1
        if len(signal) > wsize:
            signal = savitzky_golay_filter(
                    signal, filter_wsize)

    # Compute average of current window (size of 1 kb)
    lambda_1k = moving_average(
            signal, int(np.round(1000. / wsize)))

    # Beta is the relative size of current window w.r.t.
    # the whole chromosome
    beta = 1000. / len(chromosome)

    # Estimate background noise
    lambda_mean = np.mean(signal)
    lambda_bg = (lambda_mean - beta * lambda_1k) / (1. - beta)

    # Compute dynamic parameter lambda_local as in MACS peak caller
    lambda_5k = moving_average(
            signal, int(np.round(5000. / wsize)))
    lambda_10k = moving_average(
            signal, int(np.round(10000. / wsize)))
    lambda_local = np.maximum(
            np.maximum(lambda_bg, lambda_5k), lambda_10k)

    # Compute survival function
    log_sf = np.nan_to_num(poisson.logsf(signal, lambda_bg))

    # Identify segments of adjacent called peaks.
    # Segment starts correspond to 1 in `borders` and
    # segment ends correspond to 0 in `borders`.
    is_significant = (log_sf < np.log(alpha)).astype(np.int8)
    borders = np.empty(len(is_significant), dtype=np.int8)
    borders[1:] = is_significant[1:] - is_significant[:-1]
    borders[0] = int(is_significant[0])

    # Retrieve actual segment start and end positions
    # positions are measured in tags (not bases).
    start_positions = np.where(borders == 1)[0]
    end_positions = np.where(borders == -1)[0]
    if len(end_positions) == len(start_positions) - 1:
        end_positions = np.concatenate(
                [end_positions, [len(borders)-1]])
    
    # Compute scores
    scores = list()
    for start, end in zip(start_positions, end_positions):
        if end - start > 0:
            scores.append(-10. * np.min(log_sf[start:end]))
        else:
            # Peak with -inf score are removed afterwards
            scores.append(-np.inf)

    # Store peak calling information
    # and convert start-end positions from tags to bases
    return CalledPeaks(
        chromosome.name,
        chromosome.start_positions[start_positions], # conversion
        chromosome.end_positions[end_positions], # conversion
        scores, log_sf, lambda_bg, lambda_local, lambda_1k)
