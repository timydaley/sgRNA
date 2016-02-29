from __future__ import print_function
import argparse
import sys
import re
import math
from collections import defaultdict


# concatenate test_seq first
def exact_match(test_seq, reference_seqs, length) :
  IS_MATCH = False
  # check
  if len(test_seq) >= length :
    test_seq = test_seq[len(test_seq) - length : len(test_seq)]
    IS_MATCH = test_seq in reference_seqs
  return(IS_MATCH)

# use pigeonhole principle to test for edit distance one
def edit_distance_one_match(test_seq, reference_seqs, length) :
  IS_MATCH = False
  # check length
  if len(test_seq) >= length :
    # first half of sequence


                                






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
  parser.add_argument('-V', '--VERBOSE', action = "store_true",
                      dest = 'VERBOSE', default = False,
                      help = "run in verbose mode")
  args = parser.parse_args()

  args.matching_seq = args.matching_seq.replace(' ', '')
  if args.VERBOSE:
    print("sequence to match = ", args.matching_seq, file = sys.stderr)

  # read reference sgRNAs, save them in a hash table
  sgRNAs_set = set(line.strip() for line in args.ref_sgRNAs.readline())
  sgRNAs_list = list(sgRNAs_set)
  # length equals length of sgRNA
  length = len(sgRNAs_list[0])
  # make sure all reference sgRNAs are the same length
  for sgRNA in sgRNAs_list :
    assert len(sgRNA) == length
  if args.VERBOSE:
    print("size of reference sequence set = ", len(sgRNAs_set), file = sys.stderr)
  
  args.output_filename.write("obs_sequence\tref_sequence\tmismatches\n")
                             
  # read fastq file, iteratively check for existance in sgRNA_set
  indx = 0
  for line in args.input_filename.readlines() :
    line = line.rstrip('\n')
    indx += 1
    # 2nd of overy 4 line is the sequence in a fastq file
    if indx % 4 == 2 :
      if args.VERBOSE :
        print(line, file = sys.stderr)
      seq_search = re.compile(args.matching_seq)
      seq_match = seq_search.search(line)
      # match, process preceeding sequence
      if seq_match:
        match_loc = seq_match.start()
        if args.VERBOSE:
          print("match location found at position ", match_loc, file = sys.stderr)
        # test to see if position is correct to find a match
        if match_loc >= length + 1
          test_seq = line[(match_loc - length - 1):(match_loc - 1)]
          if exact_match(test_seq, reference_seqs, length) :
            args.output_filename.write("%s\t%s\t0\n" % (test_seq, test_seq))
          # test for edit distance 1 away
          elif args.n_mismatches > 0 :
            ref_seqs_first_half = []
            for seq in reference_seqs :
              ref_seqs_first_half.append(seq[1:floor(length/2)])
            # second half of sequence
            ref_seqs_second_half = []
            for seq in reference_seqs :
              ref_seqs_second_half.append(seq[ceil(length/2):len(seq)])
            # construct hash tables (dictionary in python) with
            ref_seqs_halves_dict = defaultdict(zip(ref_seqs_first_half + ref_seqs_second_half,
                                                   reference_seqs + reference_seqs))
            # testing
            if args.VERBOSE :
              for key, val in ref_seqs_halves_dict.iteritems() :
                print(str(key), " ", str(val), file = sys.stderr)
  
      # no match
      else:
        if args.VERBOSE :
          print("No match", file = sys.stderr)

  

if __name__ == "__main__":
  main()
