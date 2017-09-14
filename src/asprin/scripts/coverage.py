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

# standard python includes
import sys, os, argparse
import gzip,subprocess,re,logging,time,datetime,commands,argparse,random
from collections import defaultdict
import pysam
import progressbar
import thread
from multiprocessing.dummy import Pool as ThreadPool

class Coverage :
  def __init__(self, clipseq_fn, rnaseq_fn, \
               minimum_coverage, nthreads, all_variants):
               
    self.clipseq_fn = clipseq_fn
    self.rnaseq_fn = rnaseq_fn
    self.minimum_coverage = minimum_coverage
    self.nthreads = nthreads
    self.all_variants = all_variants

    self.genotype_info = defaultdict(lambda: defaultdict(list))
    self.snp_counter = 0

    self.bar = progressbar.ProgressBar()
    self.progress_counter = 0

    self.clip_seq_file = pysam.AlignmentFile(clipseq_fn)
    self.clip_reads = defaultdict(lambda: defaultdict(list))
    self.clip_coverage = defaultdict(lambda: defaultdict())

    self.rna_seq_file = pysam.AlignmentFile(rnaseq_fn)
    self.rna_reads = defaultdict(lambda: defaultdict(list))
    self.rna_coverage = defaultdict(lambda: defaultdict())

  def clipseq_reads_coverage_chrom(self, chrom):
    for pos in self.genotype_info[chrom]:
      if (self.genotype_info[chrom][pos][2] != "none") :
        for pileupcolumn in self.clip_seq_file.pileup(chrom, pos-1, pos, \
                                             truncate=True, stepper="nofilter"):
          for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
              self.clip_reads[chrom][pos][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
        self.clip_coverage[chrom][pos] = \
         self.clip_reads[chrom][pos][self.genotype_info[chrom][pos][0]] + \
         self.clip_reads[chrom][pos][self.genotype_info[chrom][pos][1]]
      self.progress_counter += 1
      self.bar.update(self.progress_counter)

  def rnaseq_reads_coverage_chrom(self, chrom):
    for pos in self.genotype_info[chrom]:
      if (self.genotype_info[chrom][pos][2] != "none" and \
          self.clip_coverage[chrom][pos] >= self.minimum_coverage) :
        for pileupcolumn in self.rna_seq_file.pileup(chrom, pos-1, pos, \
                                             truncate=True, stepper="nofilter"):
          for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
              self.rna_reads[chrom][pos][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
        self.rna_coverage[chrom][pos] = \
          self.rna_reads[chrom][pos][self.genotype_info[chrom][pos][0]] + \
          self.rna_reads[chrom][pos][self.genotype_info[chrom][pos][1]]
      self.progress_counter += 1
      self.bar.update(self.progress_counter)

  def clipseq_reads_coverage(self):
    sys.stderr.write('Reading CLIP-seq mapped reads file from: ' + \
                      self.clipseq_fn + '\n')
    self.bar = progressbar.ProgressBar(maxval=self.snp_counter, \
      widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
    self.progress_counter = 0
    self.bar.start()
    pool = ThreadPool(self.nthreads)
    pool.map(self.clipseq_reads_coverage_chrom, self.genotype_info)
    self.bar.finish()

  def rnaseq_reads_coverage(self):
    sys.stderr.write('Reading RNA-seq mapped reads file from: ' + \
                      self.rnaseq_fn + '\n')
    self.bar = progressbar.ProgressBar(maxval=self.snp_counter, \
      widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
    self.progress_counter = 0
    self.bar.start()
    pool = ThreadPool(self.nthreads)
    pool.map(self.rnaseq_reads_coverage_chrom, self.genotype_info)
    self.bar.finish()

  def initialize_tables(self):
    for chrom in self.genotype_info:
      for pos in self.genotype_info[chrom]:
        self.clip_reads[chrom][pos] = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
        self.rna_reads[chrom][pos] = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
