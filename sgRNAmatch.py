from __future__ import print_function
import argparse
import sys
import re


def main():
  parser = argparse.ArgumentParser(description = 'A program to count sgRNAs in a sequenced read file from a reference')
  parser.add_argument('-o', '--output_filename', dest = 'output_filename',
                      type = argparse.FileType('w'),
                      help = "output file name")
  parser.add_argument('-r', '--reference_filename', dest = 'ref_sgRNAs',
                      type = argparse.FileType('r'),
                      help = "file name of reference sgRNA sequences")
  parser.add_argument('-s', '--matching_seq', dest = 'matching_seq',
                      help = "sequence following the sgRNA to use for matching")
  parser.add_argument('-i', '--input_filename', dest = 'input_filename',
                      type = argparse.FileType('r'),
                      help = "input fastq file name")
  parser.add_argument('-m', '--mismatches', dest = 'n_mismatches',
                      default = 0, type = int, 
                      help = "number of mismatches to allow")
  parser.add_argument('-l', '--length', dest = 'length',
                      default = 20, type = int,
                      help = "length of sgRNA")
  parser.add_argument('-V', '--VERBOSE', action = "store_true",
                      dest = 'VERBOSE', default = False,
                      help = "run in verbose mode")
  args = parser.parse_args()

  if args.VERBOSE:
    print("sequence to match = ", args.matching_seq, file = sys.stderr)

  # read reference sgRNAs, save them in a hash table
  #sgRNAs_set = set(line.strip() for line in args.ref_sgRNAs.readline())
  #sgRNAs_list = list(sgRNAs_set)
  #if args.VERBOSE: 
  #  print("size of reference sequence set = ", len(sgRNAs_set), file = sys.stderr)
  
  # read fastq file, iteratively check for existance in sgRNA_set
  indx = 0
  for line in args.input_filename.readlines() :
    indx += 1
    # 2nd of overy 4 line is the sequence in a fastq file
    if indx % 4 == 2 :
      print(line, file = sys.stderr)
      seq_search = re.compile(line)
      seq_match = seq_search.search(args.matching_seq)
      # match
      if seq_match:
        match_loc = seq_match.start()
        if args.VERBOSE:
          print("match location found at position ", match_loc, file = sys.stderr)
      # no match
      else:
        print("No match", file = sys.stderr)
  

if __name__ == "__main__":
  main()
