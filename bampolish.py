import os
import sys
import argparse
import pysam
from sortedcontainers import SortedSet

def main(argv):

    parser = argparse.ArgumentParser(description="bampolish - A tool to normalise coverage in long read sequencing pipelines")

    parser.add_argument("-i", "--input",
                        help="Input sam or bam filename")
    parser.add_argument("-o", "--output",
                        help="Output sam, bam or bed filename")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Verbose output")
    parser.add_argument("-w", "--window",
                        help="Preffered window size to calculate average coverage for")
    parser.add_argument("-c", "--coverage",
                        help="Desired coverage limit for sequence for read")

    args = parser.parse_args()
    coverage(args.i, args.o, args.c)

#TODO obtain overall window average rather than strictly cut
def coverage(inputBam, outputBam, maxCoverage = 10) {
    #Maintains a set of up to maxCoverage read end points
    curr = SortedSet() 

    #Iterate over mapped reads
    bf = pysam.AlignemFile(inputBam, "rb")
    for r in bf.fetch():
        #Attempt to find a current read that ends before this read starts
        itr = curr.irange(maximum=r.reference_start))
        if(itr):
            #Replace finished read with newly started read
            curr.discard(itr.next())
            curr.add(r.reference_end)
            count += 1
            
            #TODO Write this read to memory

        else if len(curr) < maxCoverage:
            #Add read to attempt to meet coverage
            curr.add(r.reference_end)
            count += 1

    #TODO Write all reads still in curr to memory

    print("Reduced BAM from " + bf.mapped + " to " + count + " reads")
}

if __name__ == '__main__':
    main()

