# -*- coding: utf-8 -*-
# main.py
# author : Antoine Passemiers

from chipseq import call_peaks
from chipseq.io import BedGraphIO

import os
import matplotlib.pyplot as plt


DATA_FOLDER = 'data'


filepath = os.path.join(DATA_FOLDER, 'IP.bedGraph')
chromosomes = BedGraphIO.read(filepath)

called_peaks = [call_peaks(c) for c in chromosomes]

BedGraphIO.write(called_peaks, 'chipseq.out')
