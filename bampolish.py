import os
import sys
import argparse
import pysam
import math
from tqdm import tqdm
from sortedcontainers import SortedSet
import numpy as np
from coveragetree import CoverageTree
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="bampolish - A tool to normalise coverage in long read sequencing pipelines")

    parser.add_argument("-i", "--input",
                        help="Input sam or bam filename")
    parser.add_argument("-o", "--output",
                        help="Output sam or bam filename")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Verbose output", default=False)
    parser.add_argument("-w", "--windowPower",
                        help="Desired window power for coverage calculation", default=5, type=int)
    parser.add_argument("-c", "--coverage",
                        help="Desired coverage limit for sequence for read", default=10, type=int)
    parser.add_argument("-p", "--progressbar", action='store_true',
                        help="Display progress bar in terminal for all operations", default=False)
    parser.add_argument("-g", "--graphOutput", action='store_true',
                        help="Produce a column graph using matplotlib on output coverage", default=False)
    parser.add_argument("-n", "--ndeviation",
                        help="Number of standard deviations that output should be filtered to", default=4, type=int)
    parser.add_argument("-s", "--stddeviationflag", action='store_true',
                        help="Use standard deviations as the method for filtering reads", default=False)
    args = parser.parse_args()

    # Reading in of sam/bam file
    file_ext = args.input.split('.')[-1]

    if file_ext != 'sam' and file_ext != 'bam':
        # Early exit if wrong file type provided
        sys.stderr.write("ERROR: Provided file is not of sam or bam type\n")
        exit()
    else:
        # Base filename
        basename = args.input.split('.')[0]

        # Read in files
        inFile = pysam.AlignmentFile(args.input, "r") # Autodetect SAM/BAM
        outName = str(basename + "_out.bam")
        outFile = pysam.AlignmentFile(outName, "wb", template=inFile) # Default out as BAM (need to change later)

        # Perform operations
        # coverage(inFile, outFile, args.coverage, args.verbose)
        # calculate_average_coverage_windows(inFile, args.window)

        # Build coverage tree

        for ref in inFile.references:
            total = inFile.count(reference=ref, until_eof=True)
            if total == 0: continue

            # Currently we remake tree each time, but memory should be reused in future
            ct = CoverageTree(inFile, ref, total, args.windowPower, args.verbose, args.progressbar, args.ndeviation, args.stddeviationflag)

            # Write output
            ct._outputFiltered(outFile)
            outFile.close()
            pysam.index(outName)

            ct._outputSpikes(filename=str(basename + "_spikes.bed"))

            # CIGARCoverage(inFile, args.verbose, args.progressbar)

            if args.graphOutput:
                plotOutput(array, basename)
                array = ct._getLeaves()
        # Close all files
        inFile.close()
        outFile.close()


def greedyCoverage(inFile, outFile, maxCoverage, verboseFlag):
    '''
    This method is greedy but reads in low coverage areas may be missed if the
    current set is full.
    '''

    if verboseFlag:
        print("Reducing coverage using the Greedy method")

    curr = SortedSet()
    mapped = 0
    filtered = 0
    for r in inFile.fetch(until_eof=True):
        if(r.is_unmapped): continue
        mapped += 1

        # Attempt to find read that ends before this one
        itr = curr.irange(maximum=r.reference_start)
        try:
            ending = itr.__next__()
            # Some read is ending, replace it in the current set
            curr.discard(ending)
            curr.add(r.reference_end)
            outFile.write(r)
            filtered += 1
        except StopIteration:
            if(len(curr) < maxCoverage):
                # There is still room to include this read
                curr.add(r.reference_end)
                outFile.write(r)
                filtered += 1

    if verboseFlag:
        print("Reduced BAM from " + str(mapped) + " to " + str(filtered) + " reads")

def calculate_average_coverage_windows(inFile, windowSize, verboseFlag):
    '''
    Calculates average using windows of specified size. Slow but produces exact result
    '''

    # Build dictionary of references
    ref_dic = {}
    for r, l in zip(*[inFile.references, inFile.lengths]):
        ref_dic[r] = l

    for reference, length in tqdm(ref_dic.items()):
        start = 0
        end = windowSize

        if verboseFlag:
            print(reference + " " + str(length))

        total = 0
        for window in range(math.floor(length/windowSize)):
            count = 0
            num_windows = 0

            for pileupcolumn in inFile.pileup(reference, start, end):
                count += pileupcolumn.n

            if verboseFlag:
                print("The average coverage in window " + str(start) + " - " + str(end) + " is " + str(count))
            total += count/windowSize
            start = end + 1
            end = start + windowSize
            num_windows += 1

        average = total/num_windows

        if verbose:
            print("The average coverage for the file is " + str(average))

def CIGARCoverage(inFile, verboseFlag, progressbarFlag):
    '''
    Calculates coverage by analysing CIGAR strings. The method is documented here:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030888/#sup1
    https://github.com/brentp/mosdepth#howitworks
    '''

    # Create array to hold coverage
    coverage_array = np.zeros(250000000, dtype=int)

    for ref in inFile.references:
        if verboseFlag:
            print("===== Processing " + str(ref) + " =====")
        for read in tqdm(inFile.fetch(ref, until_eof=True), total=inFile.count(reference=ref, until_eof=True), disable=not progressbarFlag):
            if not read.is_unmapped:
                # Only if read is mapped, perform these operations
                counter = read.reference_start
                for cigPair in read.cigartuples:
                    # Reading each CIGAR alignment which is of the form (operation, length)
                    if (cigPair[0] == 0) or (cigPair[0] == 4) or (cigPair[0] == 7) or (cigPair[0] == 8):
                        # Match, soft-clip, == and mismatch all increase coverage
                        coverage_array[counter : counter + cigPair[1]] += 1
                    # Increment counter
                    counter += cigPair[1]

    # Find cumulative sum to calculate the total coverage array
    coverage_array_final = coverage_array.cumsum()

    return coverage_array_final

def plotOutput(array, basename):
    '''
    Plots a column graph of coverage using matplotlib
    '''

    window = len(array)/576
    avgArray = []
    elementSum = 0
    avg = 0
    i = 0
    xVal = []

    for x in array:
        elementSum += array[x]
        i += 1
        if i == round(window):
            avg = elementSum/round(window)
            avgArray.append(round(avg))
            i = 0
            elementSum = 0

    for z in range(575):
        xVal.append(z)

    ax = plt.bar(xVal, avgArray)
    plt.ylabel('Coverage')
    plt.xlabel('Bases (Windowed)')
    plt.title('Normalised Coverage of ' + basename)

    plt.savefig(str(basename + "_plot.png"))


if __name__ == '__main__':
    main()
