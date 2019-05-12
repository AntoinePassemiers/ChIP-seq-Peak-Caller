# -*- coding: utf-8 -*-
# main.py
# author : Antoine Passemiers

from chipseq import call_peaks
from chipseq.io import BedGraphIO

import os
import numpy as np
import matplotlib.pyplot as plt


DATA_FOLDER = 'data'

# Significance level
ALPHA = 1e-5

# Parse chromosomes from bedGraph file
filepath = os.path.join(DATA_FOLDER, 'IP.bedGraph')
chromosomes = { c.name: c for c in BedGraphIO.read(filepath) }

# Keep only chromosomes of interest
CHR_NAMES = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrM', 'chrX', 'chrY']
COLORS = ['#69B578', '#3A7D44', '#254D32', '#856084', '#907AD6', '#4F518C', 'orangered', '#7FDEFF']
chromosomes = { name: chromosomes[name] for name in CHR_NAMES }

# Run peak calling algorithm on each chromosome
called_peaks = { c.name: call_peaks(c, alpha=ALPHA, apply_filter=False) \
        for c in chromosomes.values() }
called_peaks_filtered = { c.name: call_peaks(c, alpha=ALPHA, apply_filter=True) \
        for c in chromosomes.values() }

BedGraphIO.write(called_peaks.values(), 'chipseq.out')
BedGraphIO.write(called_peaks_filtered.values(), 'chipseq_filtered.out')

# Plot of detailed peak calls of some region of chromosome 2L
def plot_chr2L(called_peaks, subtitle):
    start, end = 0, 3000
    chr2L = chromosomes['chr2L']
    idx = chr2L.start_positions[start:end]
    peaks = called_peaks['chr2L']
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    fig.suptitle('Peaks called in region %i-%i of chromosome 2L (%s)' % \
            (idx[0], idx[-1], subtitle))
    ax1.plot(idx, chr2L.signal[start:end],
            label='signal', color='lightgreen')
    ax1.plot(idx, peaks.lambda_bg[start:end],
            label='lambda_BG', linestyle='--', color='#BBBBBB')
    ax1.plot(idx, peaks.lambda_1k[start:end],
            label='lambda_1k', color='darkgreen')
    ax1.plot(idx, peaks.lambda_local[start:end],
            label='lambda_local', color='blue')
    ax1.set_ylabel('Average coverage')
    ax1.legend()
    ax2.semilogy(idx, np.exp(peaks.log_sf[start:end]),
            label='p-values', color='orange')
    ax2.semilogy(idx, np.full(end-start, ALPHA),
            label='Significance threshold', linestyle='--', color='orangered')
    ax2.set_xlabel('Position (expressed in bases)')
    ax2.set_ylabel('Probability')
    ax2.legend()
    plt.show()
plot_chr2L(called_peaks, 'no prior filtering')
plot_chr2L(called_peaks_filtered, 'Savitzky-Golay filtering')

# Num called peaks w/wo prior filtering
xs = np.arange(len(CHR_NAMES))
ys = [len(called_peaks[name]) for name in CHR_NAMES]
zs = [len(called_peaks_filtered[name]) for name in CHR_NAMES]
fig, ax = plt.subplots()
width, opacity = .35, .9
ax.bar(xs, ys, width,
    alpha=opacity, color=COLORS[0], label='No prior filtering')
ax.bar(xs+width, zs, width,
    alpha=opacity, color=COLORS[1], label='S-G filtering')
ax.set_xlabel('Chromosome names')
ax.set_ylabel('Number of called peaks')
ax.set_title('Number of called peaks w/wo prior filtering')
ax.set_xticks(xs + width / 2)
ax.set_xticklabels(CHR_NAMES)
ax.legend()
plt.show()

# Plot sensitivity of num called peaks to significance threshold
thresholds = np.exp(np.arange(-20, -1))
print(thresholds)
n_called_peaks = { name: list() for name in CHR_NAMES }
for alpha in thresholds:
    # Run peak calling algorithm on each chromosome
    for name in chromosomes.keys():
        called_peaks = call_peaks(chromosomes[name], alpha=alpha)
        n_called_peaks[name].append(len(called_peaks))
fig = plt.figure()
ax = plt.gca()
for i, name in enumerate(n_called_peaks.keys()):
    ax.loglog(thresholds, n_called_peaks[name],
        marker='x', c=COLORS[i], label=name)
ax.set_xlabel('Significance threshold')
ax.set_ylabel('Number of called peaks')
plt.title('Sensitivity of called peaks to significance threshold')
plt.legend()
plt.show()
