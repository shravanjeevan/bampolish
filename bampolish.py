import os
import sys
import argparse
import pysam
from sortedcontainers import SortedSet

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
        #coverage(inFile, outFile, args.coverage)
        calculate_average_coverage_windows(inFile, args.window)

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

    for references, lengths in ref_dic.items():
        start = 0
        end = windowSize

        for window in range(lengths):
            count = 0
            total = 0
            windowCount = 0
            for pileupcolumn in inFile.pileup(references, start, end):
                count += pileupcolumn.n

            print("The average coverage in window " + str(start) + " - " + str(end) + " is " + str(count))
            total += count
            start = end + 1
            end = start + 1000
            windowCount += 1
            average = total/windowCount
            print("The average coverage for the file is " + str(average))


if __name__ == '__main__':
    main()
