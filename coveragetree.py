import sys
import pysam
import numpy as np
import math
from tqdm import tqdm

class CoverageTree:
    #Number of times the chromosome length and read index resolutions are halved
    WINDOW_POWER = 0

    #Maximum length of chromosome (~250 million bases in human)
    MAX_CHR_LENGTH = 250000000 >> WINDOW_POWER

    #Size allocation of the tree padded to the next smallest power of two
    MAX_N = 2*pow(2, math.ceil(math.log(MAX_CHR_LENGTH, 2)))

    def __init__(self, filename):
        print("Constructing tree, chromosome scaled to " + str(self.MAX_CHR_LENGTH))

        self.file = pysam.AlignmentFile(filename, "r")

        #Stores sum of all subtree coverages (faster but more memory than maintaining average at each node)
        self.sumarray = np.zeros(self.MAX_N, dtype=np.int64)

        #Stores the minimum of all subtree coverages
        self.minarray = np.zeros(self.MAX_N, dtype=np.int32)

        #Used to prevent inefficient propagation during queries
        self.lazyadd = np.zeros(self.MAX_N, dtype=np.int32)

        #Iterate over all mapped reads
        #for r in tqdm(self.file.fetch("chr1"), total=self.file.mapped):
        #    if(r.is_unmapped): continue

            #Scale reads to reduce resolution
         #   readStart = r.reference_start >> self.WINDOW_POWER
         #   readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            #Ensure interval is inclusive-exclusive
         #   if(readStart == readEnd): readEnd += 1

            #Incriment all bases in range by one
            #self._update(readStart, readEnd, 1)

        self._testCoverage()

        
    #Reports coverage of the interval of the first mapped read to ensure non-zero
    def _testCoverage(self):
        firstRead = next(r for r in self.file.fetch("chr1") if not r.is_unmapped)
        firstReadStart = firstRead.reference_start >> self.WINDOW_POWER
        firstReadEnd = ((firstRead.reference_end + 1) >> self.WINDOW_POWER)
        if(firstReadStart == firstReadEnd): firstReadEnd += 1
        sumAns, minAns = self._query(firstReadStart, firstReadEnd)
        print("Average coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " + str(sumAns / (firstReadEnd - firstReadStart)))
        print("Min coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " +  str(minAns))
        print("Median coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " + str(np.sort(self._getBaseCoverages())[math.floor(self.MAX_CHR_LENGTH/2)]))

    #Updates the sum/min arrays according to the lazy counter
    def _recalculate(self, id, l, r):
        self.sumarray[id] = self.lazyadd[id] * (r - l)
        self.minarray[id] = self.lazyadd[id]
        if(r-l != 1):
            self.sumarray[id] += self.sumarray[id * 2] + self.sumarray[id * 2 + 1]
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
    
    #Returns a numpy subarray of the leaves, which correspond to individual base position coverages
    def _getBaseCoverages(self):
        return self.sumarray[self.MAX_N - self.MAX_CHR_LENGTH:]
        
if __name__ == '__main__':
    ct = CoverageTree(sys.argv[1])
