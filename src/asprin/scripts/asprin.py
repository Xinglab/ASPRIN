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

#asprin includes
from genotype import Genotype
from coverage import Coverage
from model import Model

# standard python includes
import sys, os, argparse
import gzip,subprocess,re,logging,time,datetime,commands,argparse,random
from collections import defaultdict

def write_output(genotype_info, clip_reads, clip_coverage, \
                 rna_reads, rna_coverage, asprin_qvalues, \
                 asprin_odds_ratio, minimum_coverage) :
  i = 0
  for chrom in genotype_info:
    for pos in genotype_info[chrom]:
      if (genotype_info[chrom][pos][2] != "none" and \
          clip_coverage[chrom][pos] >= minimum_coverage and \
          rna_coverage[chrom][pos] >= minimum_coverage) :

        print str(chrom) + "\t" + str(pos) + "\t" + \
              str(genotype_info[chrom][pos][2]) + "\t" + \
              str(genotype_info[chrom][pos][0]) + "\t" + \
              str(genotype_info[chrom][pos][1]) + "\t[" + \
              str(clip_reads[chrom][pos]['A']) + " , " + \
              str(clip_reads[chrom][pos]['C']) + " , " + \
              str(clip_reads[chrom][pos]['G']) + " , " + \
              str(clip_reads[chrom][pos]['T']) + "]\t[" + \
              str(rna_reads[chrom][pos]['A']) + " , " + \
              str(rna_reads[chrom][pos]['C']) + " , " + \
              str(rna_reads[chrom][pos]['G']) + " , " + \
              str(rna_reads[chrom][pos]['T']) + "]\t" + \
              str(asprin_qvalues[i]) + "\t" + str(asprin_odds_ratio[i])
        i += 1


def main() :
  parser = argparse.ArgumentParser(\
                  description='ASPRIN: Allele Specific Protein-RNA Interaction')
  group = parser.add_argument_group('required arguments')
  group.add_argument("-g", metavar="Genotype file", dest="genotype")
  group.add_argument("-c", metavar="CLIP-seq mapped reads file", dest="clipseq")
  group.add_argument("-r", metavar="RNA-seq mapped reads file", dest="rnaseq")
  parser.add_argument("-s", metavar="dbSNP VCF file", dest="dbsnp")
  parser.add_argument("-e", metavar="RADAR RNA editing database file", \
                      dest="radar")
  parser.add_argument("-t", type=int, default=25, dest="nthreads", \
                      metavar="Number of threads (default: 20)")
  parser.add_argument("-a", default="False", help="Use all the variants", \
                      action="store_true", dest="allvariants")

  args = parser.parse_args()

  if len(sys.argv) == 1 :
    parser.print_help()
  else :
    if (not args.genotype):
      sys.stderr.write("ASPRIN: -g genotype file is required\n")
      sys.exit()
    elif (not args.clipseq):
      sys.stderr.write("ASPRIN: -c CLIP-seq mapped reads file is required\n")
      sys.exit()
    elif (not args.rnaseq):
      sys.stderr.write("ASPRIN: -r RNA-seq mapped reads file is required\n")
      sys.exit()
    else :
      if (not args.dbsnp):
        dbsnp_fn = ""
      else :
        dbsnp_fn = args.dbsnp
      if (not args.radar):
        radar_fn = ""
      else :
        radar_fn = args.radar

      minimum_coverage = 10
      x = Coverage(args.clipseq, args.rnaseq, minimum_coverage, \
                   args.nthreads, args.allvariants)

      sys.stderr.write('Loading genotype information: ' + args.genotype +' \n')

      x.genotype_info, \
      x.snp_counter = Genotype(args.genotype, dbsnp_fn, radar_fn, \
                               args.allvariants).read_genotype_info()

      x.initialize_tables()
      x.clipseq_reads_coverage()
      x.rnaseq_reads_coverage()

      asprin_qvalues, asprin_odds_ratio = \
        Model(minimum_coverage).perform_test(x.genotype_info, \
                                             x.clip_reads, \
                                             x.rna_reads)

      write_output(x.genotype_info, x.clip_reads, x.clip_coverage, \
                   x.rna_reads, x.rna_coverage, asprin_qvalues, \
                   asprin_odds_ratio, minimum_coverage)

      sys.stderr.write('\n')
