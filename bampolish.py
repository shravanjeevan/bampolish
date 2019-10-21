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
                        help="Output sam, bam or bed filename")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Verbose output", default=False)
    parser.add_argument("-w", "--window",
                        help="Preffered window size to calculate average coverage for", default=1000)
    parser.add_argument("-c", "--coverage",
                        help="Desired coverage limit for sequence for read", default=10)

    args = parser.parse_args()
    coverage(args.input, args.output, args.coverage)

#TODO obtain overall window average rather than strictly cut
def coverage(inputBam, outputBam, maxCoverage):
    #Maintains a set of up to maxCoverage read end points
    curr = SortedSet() 

    #Iterate over mapped reads
    bf = pysam.AlignmentFile(inputBam, "rb")
    
    mapped = 0
    filtered = 0
    for r in bf.fetch(until_eof=True):
        if(r.is_unmapped): continue
        mapped += 1

        #Attempt to find a current read that ends before this read starts
        itr = curr.irange(maximum=r.reference_start)

        try:
            ending = itr.__next__()
            curr.discard(ending)
            curr.add(r.reference_end)
            filtered += 1
            
            #TODO Write this read to memory

        except StopIteration:
            #Add read to attempt to meet coverage
            if(len(curr) < maxCoverage):
                curr.add(r.reference_end)
                filtered += 1

    #TODO Write all reads still in curr to memory

    print("Reduced BAM from " + str(mapped) + " to " + str(filtered) + " reads")

if __name__ == '__main__':
    main()

