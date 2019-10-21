import os
import sys
import argparse


def main():

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


if __name__ == '__main__':
    main()
