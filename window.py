import os
import sys
import argparse
import pysam
from sortedcontainers import SortedSet


inFile = pysam.AlignmentFile(inName, "rb")

ref_dic = {}
for r, l in zip(*[inFile.references, inFile.lengths]):
    ref_dic[r] = l

for references, lengths in ref_dic.items:

	start = 0
	end = 1000

	for window in range(lengths):

		for pileupcolumn in samfile.pileup("references", star, end):
			count+=pileupcolumn.n	

		print("The average coverage in window" + start + "-"+ end + "is"+count)

		total+=count
		start = end+1
		end = start+100
		windowCount+=1


	average = total/windowCount

	print("The average coverage for the file is " + average)








