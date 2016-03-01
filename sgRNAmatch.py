from __future__ import print_function
import argparse
import sys
import re
import math
from collections import defaultdict
# why write edit distance function when I can import it?
import Levenshtein


# concatenate test_seq first
def exact_match(test_seq, reference_seqs, length) :
  IS_MATCH = False
  # check
  if len(test_seq) >= length :
    test_seq = test_seq[len(test_seq) - length : len(test_seq)]
    IS_MATCH = test_seq in reference_seqs
  return(IS_MATCH)







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
                      help = "number of mismatches to allow, maximum of 2")
  parser.add_argument('-V', '--VERBOSE', action = "store_true",
                      dest = 'VERBOSE', default = False,
                      help = "run in verbose mode")
  args = parser.parse_args()

  # ensure n_mismatches <= 2
  assert (args.n_mismatches < 3)  & (args.n_mismatches >= 0)
  
  # strip whitespace from matching sequence
  args.matching_seq = args.matching_seq.replace(' ', '')
  # testing check
  if args.VERBOSE:
    print("sequence to match = ", args.matching_seq, file = sys.stderr)

  # read reference sgRNAs, save them in a hash table
  sgRNAs_list = list(line.strip() for line in args.ref_sgRNAs.readlines())
  sgRNAs_set = set(sgRNAs_list)
  # testing check
  if args.VERBOSE:
    print("size of reference sequence set = ", len(sgRNAs_list), file = sys.stderr)
  # length equals length of sgRNA
  length = len(sgRNAs_list[0])
  # make sure all reference sgRNAs are the same length
  for sgRNA in sgRNAs_list :
    if len(sgRNA) != length :
      print("sgRNA = ", sgRNA, ", expected length = ", length, ", true length = ", len(sgRNA), file = sys.stderr)
    assert len(sgRNA) == length

  if args.n_mismatches > 0 :
    # constuct a hash table containing the two halfs of the reference sgRNAs
    # reduces the number of comparisons
    ref_seqs_first_half = []
    for seq in sgRNAs_list :
      ref_seqs_first_half.append(seq[0:(length/2)])
    # second half of sequence
    ref_seqs_second_half = []
    for seq in sgRNAs_list :
      ref_seqs_second_half.append(seq[(length/2):len(seq)])
    # construct hash tables (dictionary in python) with
    ref_seqs_halves_dict = defaultdict(list)
    for ref_seq, sgRNA in zip(ref_seqs_first_half + ref_seqs_second_half, sgRNAs_list + sgRNAs_list) :
      ref_seqs_halves_dict[ref_seq].append(sgRNA)

    # testing
    #if args.VERBOSE :
    #  for key, val in ref_seqs_halves_dict.iteritems() :
    #    print(str(key), " ", str(val), file = sys.stderr)


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
      seq_match = seq_search.finditer(line)
      print(seq_match, file = sys.stderr)
      match_locs = [[m.start()] for m in seq_match]
      # match, process preceeding sequence
      if len(match_locs) > 0:
        # first match
        match_loc = match_locs[0]
        # if match_loc is too small, move to next match
        if match_loc < length & len(match_locs) > 1:
          match_loc = match_locs[1]
        match_loc = match_loc[0]
        # afraid to go further
        # testing check
        if args.VERBOSE:
          print("match location found at position ", match_loc, file = sys.stderr)
        # test to see if position is correct to find a match
        if match_loc >= length :
          test_seq = line[(match_loc - length):match_loc]
          if exact_match(test_seq, sgRNAs_set, length) :
            args.output_filename.write("%s\t%s\t0\n" % (test_seq, test_seq))
          # test for edit distance 1 away
          elif args.n_mismatches > 0 :
            first_half = test_seq[0:(length/2)]
            second_half = test_seq[(length/2):len(seq)]
            matches = ref_seqs_halves_dict[first_half] + ref_seqs_halves_dict[second_half]
            for i in matches :
              d = Levenshtein.distance(i, test_seq)
              if d <= 1 :
                args.output_filename.write("%s\t%s\t%s\n" % (test_seq, i, d))
  
      # no match
      else:
        if args.VERBOSE :
          print("No match", file = sys.stderr)

  

if __name__ == "__main__":
  main()
