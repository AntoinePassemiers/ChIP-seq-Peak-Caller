# -*- coding: utf-8 -*-
# io.py: Parsers and write operations
# author : Antoine Passemiers

from chipseq.chromosome import Chromosome

import numpy as np


class BedGraphIO:

    @staticmethod
    def read(filepath, sort=False):
        """Parses a bedGraph-formatted text file.

        Parameters:
            filepath (str): Location of bedGraph file.
            sort (bool, optional): Whether to sort entries.
                This parameter should be set to True if input
                file has unordered entries.

        Returns:
            list: List of Chromosome objects.
        """

        # Read bedGraph file
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Read entries from lines
        chrom = dict()
        for line in lines:
            if len(line) > 1:
                elements = line.split()
                chrom_name = elements[0]
                start = int(elements[1])
                end = int(elements[2])
                coverage = float(elements[3])
                entry = [start, end, coverage]
                if chrom_name not in chrom.keys():
                    chrom[chrom_name] = [entry]
                else:
                    chrom[chrom_name].append(entry)

        # Sort entries if asked
        if sort:
            for name in chrom.keys():
                chrom[name] = sorted(chrom[name], key=lambda x: x[0])

        # Instantiate chromosomes
        chromosomes = list()
        for name in chrom:
            start = np.array([x[0] for x in chrom[name]], dtype=np.int)
            end = np.array([x[1] for x in chrom[name]], dtype=np.int)
            coverage = np.array([x[2] for x in chrom[name]], dtype=np.float)
            chromosomes.append(Chromosome(name, start, end, coverage))
        return chromosomes

    @staticmethod
    def write(called_peaks_list, filepath):
        """Writes called peaks to a bedGraph-formatted text file.

        Parameters:
            called_peaks_list (list): List of `CalledPeaks` objects.
            filepath (str): Location of bedGraph file.
        """
        with open(filepath, 'w') as f:
            for called_peaks in called_peaks_list:
                name = called_peaks.name
                for i in range(len(called_peaks)):
                    start = called_peaks.start_positions[i]
                    end = called_peaks.end_positions[i]
                    score = called_peaks.score[i]
                    f.write('%s\t%i\t%i\t%f\n' % (name, start, end, score))
