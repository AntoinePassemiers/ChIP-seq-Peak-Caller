# -*- coding: utf-8 -*-
# chromosome.py: Chromosome ADT
# author : Antoine Passemiers

from chipseq.signal import moving_average

import numpy as np


class Chromosome:
    """Chromosome represented by its measurements
    of local region enrichment.

    Attributes `start`, `end` and `signal` are assumed
    to be of the same length.

    Attributes:
        _name (str): Chromosome name.
        _start (:obj:`np.ndarray`): Array of start positions.
        _end (:obj:`np.ndarray`): Array of end positions.
        _signal (:obj:`np.ndarray`): Array of local coverages.
    """

    def __init__(self, name, start, end, signal):
        assert(len(start) == len(end))
        assert(len(end) == len(signal))
        self._name = str(name)
        self._start = np.asarray(start, dtype=np.int)
        self._end = np.asarray(end, dtype=np.int)
        self._signal = np.asarray(signal, dtype=np.float)

    def __len__(self):
        """Returns chromosome size, expressed in bases."""
        return self._end[-1] - self._start[0]

    @property
    def name(self):
        return self._name
    
    @property
    def start_positions(self):
        return self._start

    @property
    def end_positions(self):
        return self._end

    @property
    def signal(self):
        return self._signal
    

class CalledPeaks:
    """Called peaks of a given chromosome.

    Attributes `start`, `end` and `score` are assumed
    to be of the same length: the number of merged detected peaks.
    Attributes `log_sf`, `lambda_bg`, `lambda_local`, `lambda_1k`
    are assumed to be of the same length: the number of tags
    in the input signal.

    Attributes:
        _name (str): Chromosome name.
        _start (:obj:`np.ndarray`): Array of start positions.
        _end (:obj:`np.ndarray`): Array of end positions.
        _signal (:obj:`np.ndarray`): Array of enrichment scores.
        _log_sf ()
    """

    def __init__(self, name, start, end, score, log_sf, lambda_bg,
                 lambda_local, lambda_1k):
        assert(len(start) == len(end))
        assert(len(end) == len(score))
        self._name = str(name)
        self._start = np.asarray(start, dtype=np.int)
        self._end = np.asarray(end, dtype=np.int)
        self._score = np.asarray(score, dtype=np.float)
        self._log_sf = np.asarray(log_sf, dtype=np.float)
        self._lambda_bg = np.asarray(lambda_bg, dtype=np.float)
        self._lambda_local = np.asarray(lambda_local, dtype=np.float)
        self._lambda_1k = np.asarray(lambda_1k, dtype=np.float)

    def __len__(self):
        return self._start.shape[0]

    @property
    def name(self):
        return self._name
    
    @property
    def start_positions(self):
        return self._start

    @property
    def end_positions(self):
        return self._end

    @property
    def score(self):
        return self._score

    @property
    def log_sf(self):
        return self._log_sf

    @property
    def lambda_bg(self):
        return self._lambda_bg
    
    @property
    def lambda_local(self):
        return self._lambda_local

    @property
    def lambda_1k(self):
        return self._lambda_1k
        
