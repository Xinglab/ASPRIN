#!/usr/bin/env python

import os,sys,argparse
import progressbar

def find_regions(exon_start_loci, exon_end_loci, cds_start, cds_end):

  for i in range(len(exon_start_loci)):
    if (exon_start_loci[i] <= cds_start < exon_end_loci[i]):
       m = i
    elif (cds_start == exon_end_loci[i]):
       m = i + 1
    if (exon_start_loci[i] <= cds_end <= exon_end_loci[i]):
       n = i
       break
  new_start = [cds_start] + exon_start_loci[m + 1:n + 1]
  new_end = exon_end_loci[m:n] + [cds_end]
  utr5_start = exon_start_loci[:m + 1]
  utr5_end = exon_end_loci[:m] + [cds_start]
  utr3_start = [cds_end] + exon_start_loci[n + 1:]
  utr3_end = exon_end_loci[n:]
  cds_length = 0
  utr5_length = 0
  utr3_length = 0
  for i in range(len(new_start)):
    cds_length += new_end[i] - new_start[i]
  for i in range(len(utr5_start)):
    utr5_length += utr5_end[i] - utr5_start[i]
  for i in range(len(utr3_start)):
    utr3_length += utr3_end[i] - utr3_start[i]
  return  [new_start, \
           new_end, \
           cds_length, \
           utr5_start, \
           utr5_end, \
           utr5_length, \
           utr3_start, \
           utr3_end, \
           utr3_length]


parser = argparse.ArgumentParser(\
                description='Finding the genomic regions that ASPRIN SNPs occur')
group = parser.add_argument_group('required arguments')
group.add_argument("-g", metavar="Gene annotation file from https://www.gencodegenes.org/", dest="gene_annot")
group.add_argument("-i", metavar="ASPRIN results", dest="asprin_file")
parser.add_argument("-o", metavar="output file (default: std)", dest="output_file")

args = parser.parse_args()

if len(sys.argv) == 1 :
  parser.print_help()
else :
  if (not args.gene_annot):
    sys.stderr.write("Gene annotation file is required from: https://www.gencodegenes.org/\n")
    sys.exit()
  elif (not args.asprin_file):
    sys.stderr.write("ASPRIN result file is required\n")
    sys.exit()
  else :
    if (args.output_file):
      sys.stdout = open(args.output_file, 'w')

    gene = {}
    gene_strand = {}

    for line in open(args.gene_annot):
      line = line.rstrip('\n\r')
      ref_record = line.split('\t')
      if (ref_record[0] == "#bin"): continue
      transcript = ref_record[1] 
      chrom = ref_record[2]
      strand = ref_record[3]
      cds_start = int(ref_record[6])
      cds_end = int(ref_record[7])
      transcript_start = ref_record[9].split(",")[:-1]
      transcript_end = ref_record[10].split(",")[:-1]
      gene_name = ref_record[12]
      for i in range(len(transcript_start)):
        transcript_start[i] = int(transcript_start[i])
        transcript_end[i] = int(transcript_end[i])
      updated_regions = find_regions(transcript_start, transcript_end, cds_start, cds_end)
      gene[transcript] = [chrom] + updated_regions[:] + [cds_start] + [cds_end] + [strand] + [gene_name]
      gene_strand[gene_name] = strand

    asprin_in = open(args.asprin_file)
    asprin_size = sum(1 for line in asprin_in)
    asprin_in.close()
    bar = progressbar.ProgressBar(maxval=asprin_size, \
      widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
    bar.start()
    asprin_counter = 0
    asprin_in = open(args.asprin_file)

    for line in asprin_in:
      line = line.rstrip('\n\r')
      asprin_result = line.split('\t')
      if (asprin_result[0] == "chromosome"): continue
      chrom = str(asprin_result[0])
      site = int(asprin_result[1])
      transcript = []
      region_type = []
      gene_name = []

      for key in gene.keys():
        transcript_chrom = gene[key][0]
        if (chrom != transcript_chrom): continue
        cds_start = []
        cds_end = []
        region = ""
        utr5_start = []
        utr5_end = []
        utr3_start = []
        utr3_end = []
      
        strand = gene[key][12]
        transcript_gene_name = gene[key][13]

        if(gene[key][4][0] <= site < gene[key][8][-1]):
          if (not gene[key][13] in gene_name):
            gene_name.append(gene[key][13])

        if (gene[key][4][0] <= site < gene[key][5][-1]):
          utr5_start = gene[key][4]
          utr5_end = gene[key][5]
          transcript.append(key)
          for j in range(len(utr5_start)):
            if (utr5_start[j] <= site < utr5_end[j]):
              if (gene[key][1] == gene[key][2]): region = "Noncoding"
              elif (strand=="+"): region = "5UTR"
              else: region = "3UTR"
              break
            elif (site < utr5_start[j]):
              if (((utr5_start[j] - site) > 500) and ((site - utr5_end[j - 1]) > 500)) :
                region = "Distal_Intron"
              else :
                first = (utr5_start[j] - site)
                last = (site - utr5_end[j - 1])
                if (((first <= last) and (strand=="+")) or ((last <= first) and (strand=="-"))) :
                  region = "Upstream_Proximal"
                elif (((last <= first) and (strand=="+")) or ((first <= last) and (strand=="-"))) :
                  region = "Downstream_Proximal"
              break
          region_type.append(region)
        elif (gene[key][7][0] <= site < gene[key][8][-1]):
          utr3_start = gene[key][7]
          utr3_end = gene[key][8]
          transcript.append(key)
          for j in range(len(utr3_start)):
            if (utr3_start[j] <= site < utr3_end[j]):
              if (gene[key][1] == gene[key][2]): region = "Noncoding"
              elif (strand == "+"): region = "3UTR"
              else: region = "5UTR"
              break
            elif (site < utr3_start[j]):
              if (((utr3_start[j] - site) > 500) and  ((site - utr3_end[j - 1]) > 500)) :
                region = "Distal_Intron"
              else :
                first = (utr3_start[j] - site)
                last = (site - utr3_end[j - 1])
                if (((first <= last) and (strand=="+")) or ((last <= first) and (strand=="-"))) :
                  region = "Upstream_Proximal"
                elif (((last <= first) and (strand=="+")) or ((first <= last) and (strand=="-"))) :
                  region = "Downstream_Proximal"
              break
          region_type.append(region)
        elif (gene[key][10] <= site < gene[key][11]):
          cds_start = gene[key][1]
          cds_end = gene[key][2]
          transcript.append(key)
          for j in range(len(cds_start)):
            if (cds_start[j] <= site < cds_end[j]):
              region = "CDS"
              break
            elif (site < cds_start[j]):
              if (((cds_start[j] - site) > 500) and ((site - cds_end[j - 1]) > 500)) :
                region = "Distal_Intron"
              else :
                first = (cds_start[j] - site)
                last = (site - cds_end[j - 1])
                if (((first <= last) and (strand=="+")) or ((last <= first) and (strand=="-"))) :
                  region = "Upstream_Proximal"
                elif (((last < first) and (strand=="+")) or ((first < last) and (strand=="-"))) :
                  region = "Downstream_Proximal"
              break
          region_type.append(region)

      region_consensus = ""
      if ("CDS" in region_type): region_consensus = "CDS"
      elif ("5UTR" in region_type): region_consensus = "5UTR"
      elif ("3UTR" in region_type): region_consensus = "3UTR"
      elif ("Noncoding" in region_type): region_consensus = "Noncoding"
      elif ("Upstream_Proximal" in region_type): region_consensus = "Upstream_Proximal"
      elif ("Downstream_Proximal" in region_type): region_consensus = "Downstream_Proximal"
      elif ("Distal_Intron" in region_type): region_consensus = "Distal_Intron"

      if (len(gene_name) > 0) :
        print('\t'.join([line , ','.join(gene_name), region_consensus, \
                                gene_strand[gene_name[0]]]))
      else :
        print('\t'.join([line , "NA", "Intergenic", "NA"]))

      asprin_counter += 1
      if (asprin_counter <= asprin_size) :
        bar.update(asprin_counter)

    asprin_in.close()
    bar.finish()

