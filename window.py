import os
import sys
import argparse
import pysam
from sortedcontainers import SortedSet


inFile = pysam.AlignmentFile(inName, "rb")

start = 0
end = 100

for window in range(pysam.AlignmentFile.reference_lengths)

	for pileupcolumn in samfile.pileup("chr1", star, end):
		count+=pileupcolumn.n	

	print("The average coverage in window" + start + "-"+ end + "is"+count)

	total+=count
	start = end+1
	end = start+100
	windowCount+=1


average = total/windowCount

print("The average coverage for the file is " + average)








