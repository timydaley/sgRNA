import argparse
import sys


def main():
  parser = argeparse.ArgumentParser(description = 'A program to count sgRNAs in a sequenced read file from a reference')
  parser.add_argument('-o', '--output_filename', dest = 'output_filename',
                      type = argparse.FileType('w'),
                      help = 'output file name')
  parser.add_argument('-r', '--reference_filename', dest = 'ref_sgRNAs',
                      type = argparse.FileType('r'),
                      help = 'file name of reference sgRNA sequences')
  parser.add_argument('-i', '--input_filename', dest = 'input_filename',
                      type = argparse.FileType('r'),
                      help = 'input fastq file name')
  parser.add_argument('-m', '--mismatches', dest = 'n_mismatches',
                      default = 0, type = int, 
                      help = 'number of mismatches to allow')
  parser.add_argument('-V', '--VERBOSE', action = "store_true",
                      dest = 'VERBOSE', default = False,
                      help = "run in verbose mode")
  args = parser.parse_args()

  # read reference sgRNAs, save them in a hash table
  sgRNAs_set = set(line.strip() for line in args.ref_sgRNAs.readline())
  sgRNAs_list = list(sgRNAs_set)

  # read fastq file, iteratively check for existance in sgRNA_set
  indx = 0
  

if __name__ == "__main__":
  main()
