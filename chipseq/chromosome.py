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
    to be of the same length.

    Attributes:
        _name (str): Chromosome name.
        _start (:obj:`np.ndarray`): Array of start positions.
        _end (:obj:`np.ndarray`): Array of end positions.
        _signal (:obj:`np.ndarray`): Array of enrichment scores.
    """

    def __init__(self, name, start, end, score):
        assert(len(start) == len(end))
        assert(len(end) == len(score))
        self._name = str(name)
        self._start = np.asarray(start, dtype=np.int)
        self._end = np.asarray(end, dtype=np.int)
        self._score = np.asarray(score, dtype=np.float)

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
