import os
import sys
import argparse
import pysam
import math
from tqdm import tqdm
from sortedcontainers import SortedSet
import numpy as np

def main():

    parser = argparse.ArgumentParser(description="bampolish - A tool to normalise coverage in long read sequencing pipelines")

    parser.add_argument("-i", "--input",
                        help="Input sam or bam filename")
    parser.add_argument("-o", "--output",
                        help="Output sam, bam or bed filename", default="output.bam")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Verbose output", default=False)
    parser.add_argument("-w", "--window",
                        help="Desired window size to calculate average coverage for", default=1000, type=int)
    parser.add_argument("-c", "--coverage",
                        help="Desired coverage limit for sequence for read", default=10, type=int)

    args = parser.parse_args()

    # Reading in of sam/bam file
    file_ext = args.input.split('.')[-1]

    if file_ext != 'sam' and file_ext != 'bam':
        # Early exit if wrong file type provided
        sys.stderr.write("ERROR: Provided file is not of sam or bam type\n")
        exit()
    else:
        # Read in files
        inFile = pysam.AlignmentFile(args.input, "r") # Autodetect SAM/BAM
        outFile = pysam.AlignmentFile(args.output, "wb", template=inFile) # Default out as BAM (need to change later)

        # Perform operations
        # coverage(inFile, outFile, args.coverage)
        # calculate_average_coverage_windows(inFile, args.window)
        CIGARCoverage(inFile)

#This method is greedy but reads in low coverage areas may be missed if the current set is full
#In future, windowing with a method to approach some average rather than a strict cutoff should be used
def coverage(inFile, ouFile, maxCoverage):


    curr = SortedSet()
    mapped = 0
    filtered = 0
    for r in inFile.fetch(until_eof=True):
        if(r.is_unmapped): continue
        mapped += 1

        #Attempt to find read that ends before this one
        itr = curr.irange(maximum=r.reference_start)
        try:
            ending = itr.__next__()

            #Some read is ending, replace it in the current set
            curr.discard(ending)
            curr.add(r.reference_end)
            outFile.write(r)
            filtered += 1
        except StopIteration:
            if(len(curr) < maxCoverage):
                #There is still room to include this read
                curr.add(r.reference_end)
                outFile.write(r)
                filtered += 1

    outFile.close()
    inFile.close()
    print("Reduced BAM from " + str(mapped) + " to " + str(filtered) + " reads")

def calculate_average_coverage_windows(inFile, windowSize):

    # Build dictionary of references
    ref_dic = {}
    for r, l in zip(*[inFile.references, inFile.lengths]):
        ref_dic[r] = l

    for reference, length in tqdm(ref_dic.items()):
        start = 0
        end = windowSize

        print(reference + " " + str(length))

        total = 0
        for window in range(math.floor(length/windowSize)):
            count = 0
            num_windows = 0

            for pileupcolumn in inFile.pileup(reference, start, end):
                count += pileupcolumn.n

            # print("The average coverage in window " + str(start) + " - " + str(end) + " is " + str(count))
            total += count/windowSize
            start = end + 1
            end = start + windowSize
            num_windows += 1

        average = total/num_windows
        print("The average coverage for the file is " + str(average))

def CIGARCoverage(inFile):

    # Create array to hold coverage
    coverage_array = np.zeros(250000000, dtype=int)

    for ref in inFile.references:
        for read in tqdm(inFile.fetch(ref, until_eof=True), total=inFile.mapped):
            if not read.is_unmapped:
                # Only if read is mapped, perform these operations
                counter = read.reference_start
                for cigPair in read.cigartuples:
                    # Reading each CIGAR alignment which is of the form (operation, length)
                    if (cigPair[0] == 0) or (cigPair[0] == 4) or (cigPair[0] == 7) or (cigPair[0] == 8):
                        # Match, soft-clip, == and mismatch all increase coverage
                        coverage_array[counter : counter + cigPair[1]] += 1

if __name__ == '__main__':
    main()
