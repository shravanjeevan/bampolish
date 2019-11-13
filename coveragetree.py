import sys
import pysam
import numpy as np
import math
from tqdm import tqdm
import pdb
import random

class CoverageTree:
    #Number of times the chromosome length and read index resolutions are halved
    WINDOW_POWER = 5 #2^WINDOW_POWER is the window size but most of saving is memory rather than time

    #Maximum length of chromosome (~250 million bases in human)
    MAX_CHR_LENGTH = 250000000 >> WINDOW_POWER

    #Size allocation of the binary tree (padded to the next smallest power of two for balance)
    MAX_N = 2*pow(2, math.ceil(math.log(MAX_CHR_LENGTH, 2)))

    def __init__(self, filename):
        print("Constructing tree, chromosome scaled to " + str(self.MAX_CHR_LENGTH))

        self.file = pysam.AlignmentFile(filename, "r")

        #Stores sum of all subtree coverages (faster but more memory than maintaining average at each node)
        #TODO: removing this doubles python performance - if we don't need interval averages, remove
        self.sumarray = np.zeros(self.MAX_N, dtype=np.int64)

        #Stores the minimum of all subtree coverages
        self.minarray = np.zeros(self.MAX_N, dtype=np.int32)

        #Used to prevent inefficient propagation during queries
        self.lazyadd = np.zeros(self.MAX_N, dtype=np.int32)

        i = 10000
        #Iterate over all mapped reads
        for r in tqdm(self.file.fetch("chr1"), total=self.file.mapped):
            if(r.is_unmapped): continue

            #Scale reads to reduce resolution
            readStart = r.reference_start >> self.WINDOW_POWER
            readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            #Ensure interval is inclusive-exclusive
            if(readStart == readEnd): readEnd += 1

            #Incriment all bases in range by one
            self._update(readStart, readEnd, 1)

            #i -= 1
            #if i == 0: break

        #Calculate median of mapped bases only
        #leaves = self._getLeaves()
        #leavesSorted = np.sort(leaves[leaves > 0])
        #self.median = 0 if len(leavesSorted) == 0 else leavesSorted[math.floor(len(leavesSorted)/2)]

        #sumAns, minAns = self._query(0, self.MAX_CHR_LENGTH)  
        #self.mean = sumAns / self.MAX_CHR_LENGTH
        #print("Average coverage  = "  + str(self.mean))

        #self._testCoverage()
        
    #Reports coverage of the interval of the first mapped read to ensure non-zero
    def _testCoverage(self):
        firstRead = next(r for r in self.file.fetch("chr1") if not r.is_unmapped)
        firstReadStart = firstRead.reference_start >> self.WINDOW_POWER
        firstReadEnd = (firstRead.reference_end + 1) >> self.WINDOW_POWER
        if(firstReadStart == firstReadEnd): firstReadEnd += 1
        sumAns, minAns = self._query(firstReadStart, firstReadEnd)        
        self._updateAll()
        leaves = self._getLeaves(firstReadStart, firstReadEnd)
        for leaf in leaves:
            print(str(leaf) + ' ')
        leavesSorted = np.sort(leaves[leaves > 0])
        self.median = 0 if len(leavesSorted) == 0 else leavesSorted[math.floor(len(leavesSorted)/2)]
        print("Average coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " + str(sumAns / (firstReadEnd - firstReadStart)))
        print("Min coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " +  str(minAns))
        print("Median (mapped) coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " + str(self.median))

    #Updates the sum/min arrays according to the lazy counter
    def _recalculate(self, id, l, r):
        #self.sumarray[id] = self.lazyadd[id] * (r - l)
        self.minarray[id] = self.lazyadd[id]
        if(r-l != 1):
            #self.sumarray[id] += self.sumarray[id * 2] + self.sumarray[id * 2 + 1]
            self.minarray[id] += min(self.minarray[id * 2], self.minarray[id * 2 + 1])

    #Updates lazy counters
    def _update_lazy(self, id, v, l, r):
        self.lazyadd[id] += v
        self._recalculate(id, l, r)

    #Pushes the lazy counter down the tree, recalculating as it goes
    def _propagate(self, id, l, r):
        mid = (l + r)/2
        self._update_lazy(id * 2, self.lazyadd[id], l, mid)
        self._update_lazy(id * 2 + 1, self.lazyadd[id], mid, r)
        self.lazyadd[id] = 0

    #Incriments all base coverages within the given interval
    def _update(self, uL, uR, v, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        if(uL == cLeft and uR == cRight):
            self._update_lazy(i, v, cLeft, cRight)
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        if(uL < mid): self._update(uL, min(uR, mid), v, i * 2, cLeft, mid)
        if(uR > mid): self._update(max(uL, mid), uR, v, i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)
    
    #Non lazy update of tree, so that leaves are correct, super slow though so we need CIGAR method
    def _updateAll(self, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        #pdb.set_trace()
        if(cRight - cLeft == 1):
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        self._updateAll(i * 2, cLeft, mid)
        self._updateAll(i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)

    #Returns sum and minimum coverages within the given interval
    def _query(self, qL, qR, i=1, cLeft=0,  cRight=MAX_CHR_LENGTH):
        if(qL == cLeft and qR == cRight):
            return self.sumarray[i], self.minarray[i]
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        sumAns = 0
        minAns = math.inf
        if(qL < mid): 
            sumL, minL = self._query(qL, min(qR, mid), i * 2, cLeft, mid)
            sumAns += sumL
            minAns = min(minAns, minL)

        if(qR > mid): 
            sumR, minR = self._query(max(qL, mid), qR, i * 2 + 1, mid, cRight)
            sumAns += sumR
            minAns = min(minAns, minR)
        return sumAns, minAns

    #TODO: Fix changes not getting pushed to leaves (which is good but can't extract this way then)
    #Returns a modifiable numpy subarray of the leaves, which correspond to individual base position coverages
    def _getLeaves(self, cLeft=0, cRight=MAX_CHR_LENGTH):
        return self.sumarray[self.MAX_N-self.MAX_CHR_LENGTH+cLeft : self.MAX_N-self.MAX_CHR_LENGTH+cRight]

    #Reconstructs tree in linear time from modified leaves, bottom up
    def _resetFromLeaves(self, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        self.lazyadd[i] = 0

        if(cRight - cLeft == 1):
            self.minarray[i] = self.sumarray[i]
            return
        if(uL < mid): self._resetFromLeaves(i * 2, cLeft, mid)
        if(uR > mid): self._resetFromLeaves(i * 2 + 1, mid, cRight)

        self.minarray[i] = min(self.minarray[i*2], self.minarray[i*2+1])
        self.sumarray[i] = self.sumarray[i*2] + self.sumarray[i*2+1]

    #Writes high coverage intervals to a bed file
    def _outputSpikes(self, filename="spikes.bed"):
        threshold = self.median
        i=0
        wasSpike = false
        isSpike = false
        with open(filename,'wb') as outFile:
            for c in self._getLeaves():
                isSpike = (c > threshold)
                if(isSpike and not wasSpike):
                    #Start spike
                    outFile.write(i + '\t')
                elif(wasSpike and not isSpike):
                    #End spike
                    outFile.write(i + '\n')
                wasSpike = isSpike
                i += 1

        if(isSpike): 
            outFile.write(i + '\n')

    #Writes reads to a file without going below the median where possible
    def _outputFiltered(self, filename="filtered.bam"):
        outFile = pysam.AlignmentFile(filename, "wb", template=self.file)
        coverages = self._getLeaves()
        #threshold = self.mean/2 #TODO use n med-MADs instead

        for r in tqdm(self.file.fetch("chr1"), total=self.file.mapped):
            if r.is_unmapped: continue

            #Scale reads to reduce resolution
            readStart = r.reference_start >> self.WINDOW_POWER
            readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            #Ensure interval is inclusive-exclusive
            if(readStart == readEnd): readEnd += 1

            #Write if the read is part of minimum
            sumAns, minAns = self._query(readStart, readEnd)

            #Linear falloff
            desiredCov = 20

            chanceToKeep = min(1, desiredCov/minAns)
            if(chanceToKeep >= 1 or (chanceToKeep > 0 and random.uniform(0, 1) < chanceToKeep)):
                outFile.write(r)
                #self._update(readStart, readEnd, -1)
            outFile.close()
            pysam.index(filename)
        
if __name__ == '__main__':
    ct = CoverageTree(sys.argv[1])

    #array = ct._getLeaves() #TODO Shravan and Varsha write the base coverages into this array (size 250 million but you don't have to fill it) but will need to remove constructor range updates
    #ct._resetFromLeaves()

    ct._outputFiltered("output.bam")
