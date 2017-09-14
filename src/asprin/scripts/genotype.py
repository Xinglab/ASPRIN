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
import progressbar

class Genotype :
  def __init__(self, genotype_fn, dbsnp_fn, radar_fn, all_variants) :
    self.genotype_fn = genotype_fn
    self.dbsnp_fn = dbsnp_fn
    self.radar_fn = radar_fn
    self.dbsnp_counter = 0
    self.radar_counter = 0
    self.all_variants = all_variants
    self.dbsnp_size = 323138224
    self.radar_size = 2576459

  def read_genotype_info(self):

    VCF_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL',\
                  'FILTER','INFO','FORMAT','SAMPLE']
    genotype_info = defaultdict(lambda: defaultdict(list))
    snp_counter = 0

    variant_id = "none"
    is_variant_id = False
    if ((self.radar_fn == "" and self.dbsnp_fn == "") or self.all_variants):
      is_variant_id = True

    if self.genotype_fn.split('.')[-1]=="gz":
      vcf_in = gzip.open(self.genotype_fn)
    else:
      vcf_in = open(self.genotype_fn)

    for line in vcf_in:
      if line.startswith('#'):
        continue
      result = {}
      fields = line.rstrip().split()
      for i,col in enumerate(VCF_HEADER):
        result[col] = fields[i]
      infos = [x for x in result['INFO'].split(';') if x.strip()]
      for i in infos:
        if '=' in i:
          key,value = i.split('=')
          result[key] = value
      geno = result['SAMPLE'].split(':')[0]
      if ('|' in geno or '/' in geno):
        if ('|' in geno):
          gt1 = int(geno.split('|')[0])
          gt2 = int(geno.split('|')[1])
        else :
          gt1 = int(geno.split('/')[0])
          gt2 = int(geno.split('/')[1])
        alt = [result['ALT']]
        alleles = [result['REF']] + alt
        if all([re.match(r'[ACGT]',i) for i in alleles]) and \
           all([len(i)==1 for i in alleles]):
          if ((gt1 == 0 and gt2 == 1) or (gt1 == 1 and gt2 == 0)):
            snp_counter += 1
            if (is_variant_id) : variant_id = "var_" + str(snp_counter)
            genotype_info[result['CHROM']][int(result['POS'])] = \
                                              [alleles[0],alleles[1],variant_id]
    vcf_in.close()

    DBSNP_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    RADAR_HEADER = ['chromosome','position','gene','strand','annot1','annot2',\
                    'alu?','non_alu_repetitive?','conservation_chimp',\
                    'conservation_rhesus','conservation_mouse']

    if (self.radar_fn != ""):
      sys.stderr.write('Readig RADAR editing events: ' + self.radar_fn + '\n')
      radar_in = open(self.radar_fn)

      bar = progressbar.ProgressBar(maxval=self.radar_size, \
        widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
      bar.start()
      for line in radar_in:
        fields = line.rstrip().split()
        if fields[0] == "chromosome":
          continue
        result = {}
        for i,col in enumerate(RADAR_HEADER):
          result[col] = fields[i]
        if (result['chromosome'] in genotype_info and \
            int(result['position']) in genotype_info[result['chromosome']]) :
          genotype_info[result['chromosome']][int(result['position'])][2] = \
                                                  "re" + str(self.radar_counter)
        self.radar_counter += 1
        if (self.radar_counter <= self.radar_size) :
          bar.update(self.radar_counter)
      radar_in.close()
      bar.finish()
      sys.stderr.write('\n')

    if (self.dbsnp_fn != ""):
      sys.stderr.write('Reading dbSNP SNPs: ' + self.dbsnp_fn +' \n')
      dbsnp_in = open(self.dbsnp_fn)

      bar = progressbar.ProgressBar(maxval=self.dbsnp_size, \
        widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
      bar.start()
      for line in dbsnp_in:
        if line.startswith('#'):
          continue
        result = {}
        fields = line.rstrip().split()
        for i,col in enumerate(DBSNP_HEADER):
          result[col] = fields[i]
        chrom = "chr" + result['CHROM']
        if (chrom in genotype_info and \
            int(result['POS']) in genotype_info[chrom]) :
          genotype_info[chrom][int(result['POS'])][2] = result['ID']
        self.dbsnp_counter += 1
        if (self.dbsnp_counter <= self.dbsnp_size) :
          bar.update(self.dbsnp_counter)
      dbsnp_in.close()
      bar.finish()
      sys.stderr.write('\n')

    return genotype_info, snp_counter
