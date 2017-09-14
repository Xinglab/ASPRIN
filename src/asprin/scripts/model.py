# Copyright (C) 2017 University of California, Los Angeles (UCLA)
# Emad Bahrami-Samani, Yi Xing
#
# Authors: Emad Bahrami-Samani, Yi Xing
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.

import fisher
import itertools as it
import operator
from scipy import stats
import numpy
import math
from collections import defaultdict

import mne
from mne import io
from mne.datasets import sample
from mne.stats import bonferroni_correction, fdr_correction

class Model :
  def __init__(self, minimum_coverage) :
    self.minimum_coverage = minimum_coverage

  def perform_test(self, genotype_info, clip_reads, rna_reads) :
    asprin_test = defaultdict(lambda: defaultdict(list))
    asprin_pvalues = []
    asprin_odds_ratio = []
    for chrom in genotype_info:
      for pos in genotype_info[chrom]:
        if (genotype_info[chrom][pos][2] != "none" and \
            clip_coverage[chrom][pos] >= self.minimum_coverage and \
            rna_coverage[chrom][pos] >= self.minimum_coverage) :
          asprin_test[chrom][pos] = \
            stats.fisher_exact([\
              [clip_reads[chrom][pos][genotype_info[chrom][pos][0]],\
               clip_reads[chrom][pos][genotype_info[chrom][pos][1]]],\
              [rna_reads[chrom][pos][genotype_info[chrom][pos][0]],\
               rna_reads[chrom][pos][genotype_info[chrom][pos][1]]]],\
              'two-sided')
          asprin_pvalues.append(asprin_test[chrom][pos][1])
          asprin_odds_ratio.append(asprin_test[chrom][pos][0])

    alpha = 0.1
    reject_fdr, asprin_qvalues = fdr_correction(asprin_pvalues, \
                                              alpha=alpha, method='indep')

    return asprin_qvalues, asprin_odds_ratio
