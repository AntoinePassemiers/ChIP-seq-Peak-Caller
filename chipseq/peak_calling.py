# -*- coding: utf-8 -*-
# peak_calling.py
# author : Antoine Passemiers

from chipseq.chromosome import CalledPeaks
from chipseq.signal import moving_average

import numpy as np
from scipy.stats import poisson


def call_peaks(chromosome, wsize=10, alpha=1e-5):
    """Call peaks from a given chromosome.

    Parameters:
        chromosome (:obj:´chipseq.Chromosome`): Input chromosome.
        wsize (int): Size of the sliding window.
        alpha (float): Confidence threshold for
            rejecting null hypothesis.

    Returns:
        TODO
    """
    means_1k = moving_average(chromosome.signal, wsize)
    means_5k = moving_average(chromosome.signal, 50)
    means_10k = moving_average(chromosome.signal, 100)

    lambda_bg = np.mean(chromosome.signal)
    lambda_local = np.maximum(np.maximum(lambda_bg, means_5k), means_10k)

    log_sf = np.nan_to_num(poisson.logsf(means_1k, lambda_local))

    is_significant = (log_sf >= np.log(1. - alpha)).astype(np.int8)
    borders = np.empty(len(is_significant), dtype=np.int8)
    borders[1:] = is_significant[1:] - is_significant[:-1]
    borders[0] = int(is_significant[0])

    start_positions = np.where(borders == 1)[0]
    end_positions = np.where(borders == -1)[0]
    if len(end_positions) == len(start_positions) - 1:
        end_positions = np.concatenate([end_positions, [len(borders)-1]])
    
    scores = list()
    for start, end in zip(start_positions, end_positions):
        if end - start > 0:
            scores.append(10. * np.mean(log_sf[start:end]))
        else:
            scores.append(0.)

    return CalledPeaks(
        chromosome.name,
        chromosome.start_positions[start_positions],
        chromosome.end_positions[end_positions],
        scores)
