import sys
import pysam
import numpy as np
import math
from tqdm import tqdm
import pdb
import random
import matplotlib.pyplot as plt


class CoverageTree:
    #Number of times the chromosome length and read index resolutions are halved
    WINDOW_POWER = 8 #2^WINDOW_POWER is the window size but most of saving is memory rather than time

    #Maximum length of chromosome (~250 million bases in human)
    MAX_CHR_LENGTH = 47000000 >> WINDOW_POWER

    #Size allocation of the binary tree (padded to the next smallest power of two for balance)
    MAX_N = 2*pow(2, math.ceil(math.log(MAX_CHR_LENGTH, 2)))

    def __init__(self, filename):
        print("Constructing tree")
        print("(Chromosome scaled to " + str(self.MAX_CHR_LENGTH) + ")")

        self.file = pysam.AlignmentFile(filename, "r")

        #Stores the minimum of all subtree coverages
        self.minarray = np.zeros(self.MAX_N, dtype=np.int32)

        #Used to prevent inefficient propagation during queries
        self.lazyadd = np.zeros(self.MAX_N, dtype=np.int32)

        #Iterate over all mapped reads
        for r in tqdm(self.file.fetch("chr21"), total=self.file.mapped):
            if(r.is_unmapped): continue

            #Scale reads to reduce resolution
            readStart = r.reference_start >> self.WINDOW_POWER
            readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            #Ensure interval is inclusive-exclusive
            if(readStart == readEnd): readEnd += 1

            #Incriment all bases in range by one
            self._update(readStart, readEnd, 1)

        print("Pushing range updates to leaves")
        self.t = tqdm(total=self.MAX_CHR_LENGTH)
        self._updateAll()
        del self.t

        #Calculate median of mapped bases only
        leaves = self._getLeaves()
        leavesSorted = np.sort(leaves[leaves > 0])
        self.median = 0 if len(leavesSorted) == 0 else leavesSorted[math.floor(len(leavesSorted)/2)]
        print("(Median = " + str(self.median) + ")")

    #Reports coverage of the interval of the first mapped read to ensure non-zero
    def _testCoverage(self):
        firstRead = next(r for r in self.file.fetch("chr21") if not r.is_unmapped)
        firstReadStart = firstRead.reference_start >> self.WINDOW_POWER
        firstReadEnd = (firstRead.reference_end + 1) >> self.WINDOW_POWER
        if(firstReadStart == firstReadEnd): firstReadEnd += 1
        minAns = self._query(firstReadStart, firstReadEnd)        
        self._updateAll()
        leaves = self._getLeaves(firstReadStart, firstReadEnd)
        for leaf in leaves:
            print(str(leaf) + ' ')
        leavesSorted = np.sort(leaves[leaves > 0])
        self.median = 0 if len(leavesSorted) == 0 else leavesSorted[math.floor(len(leavesSorted)/2)]
        print("Min coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " +  str(minAns))
        print("Median (mapped) coverage [" + str(firstReadStart) + ", " + str(firstReadEnd) + ") = " + str(self.median))

    #Updates the sum/min arrays according to the lazy counter
    def _recalculate(self, id, l, r):
        self.minarray[id] = self.lazyadd[id]
        if(r-l != 1):
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
    
    #Non lazy update of tree, so that leaves are correct
    def _updateAll(self, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        if(cRight - cLeft == 1):
            self.t.update(1)
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        self._updateAll(i * 2, cLeft, mid)
        self._updateAll(i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)

    #Returns sum and minimum coverages within the given interval
    def _query(self, qL, qR, i=1, cLeft=0,  cRight=MAX_CHR_LENGTH):
        if(qL == cLeft and qR == cRight):
            return self.minarray[i]
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        minAns = math.inf
        if(qL < mid): minAns = min(minAns, self._query(qL, min(qR, mid), i * 2, cLeft, mid))
        if(qR > mid): minAns = min(minAns, self._query(max(qL, mid), qR, i * 2 + 1, mid, cRight))
        return minAns

    #TODO: Fix changes not getting pushed to leaves (which is good but can't extract this way then)
    #Returns a modifiable numpy subarray of the leaves, which correspond to individual base position coverages
    def _getLeaves(self, cLeft=0, cRight=MAX_CHR_LENGTH):
        return self.lazyadd[self.MAX_N-self.MAX_CHR_LENGTH+cLeft : self.MAX_N-self.MAX_CHR_LENGTH+cRight]

    #Reconstructs tree in linear time from modified leaves, bottom up
    def _resetFromLeaves(self, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        self.lazyadd[i] = 0
        if(cRight - cLeft == 1):
            return
        if(uL < mid): self._resetFromLeaves(i * 2, cLeft, mid)
        if(uR > mid): self._resetFromLeaves(i * 2 + 1, mid, cRight)

        self.minarray[i] = min(self.minarray[i*2], self.minarray[i*2+1])

    #Calculates the median absolute deviation of base coverages
    def _getMad():
        coverages = self._getLeaves()
        deviations = np.empty(len(coverages), dtype=np.int32)

        for i in range(len(coverages)):
            deviations[i] = abs(self.median - coverages[i])

        deviations.sort()
        return deviations[math.floor(len(deviations)/2)]

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
    def _outputFiltered(self, filename):
        print("Writing filtered reads to " + filename)
        outFile = pysam.AlignmentFile(filename, "wb", template=self.file)

        for r in tqdm(self.file.fetch("chr21"), total=self.file.mapped):
            if r.is_unmapped: continue

            #Scale reads to reduce resolution
            readStart = r.reference_start >> self.WINDOW_POWER
            readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            #Ensure interval is inclusive-exclusive
            if(readStart == readEnd): readEnd += 1

            #Write if the read is part of minimum
            minAns = self._query(readStart, readEnd)

            if(minAns <= 20): #TODO use n med-MADs instead
                outFile.write(r)
            else:
                self._update(readStart, readEnd, -1)

        outFile.close()
        pysam.index(filename)

#    def _printCoverage(array):
#        avgArray = []
#        count = 0
#
#        for x in array: 
#            elementSum += array[x]
#            i+=1
#            if i == 10000 :
#                avg = elementSum/10000
#                avgArray[count] = avg
#                count += 1
#                i = 0
#
#        plt.plot(avgArray)
#        plt.ylabel('Coverage')
#        plt.show()

   
if __name__ == '__main__':
    ct = CoverageTree(sys.argv[1])

    #TODO Shravan and Varsha write the base coverages into this array (size 250 million but you don't have to fill it) but will need to remove constructor range updates
    array = ct._getLeaves()
    
    window = len(array)/576
    avgArray = []
    elementSum = 0
    avg = 0
    i= 0
    xVal = []

    #print(len(array))
    #for y in range(len(array)): 
    	#print (array[y])
    #print(round(window))

    for x in array:
    	elementSum += array[x]
    	i+=1
    	if i == round(window) :
    		avg = elementSum/round(window)
    		avgArray.append(round(avg))
    		i = 0

  
    #for y in range(len(avgArray)): 
    	#print (avgArray[y])

    for z in range(575):
    	xVal.append(z)

    #print(len(xVal))
    #print(len(avgArray))
  
    ax = plt.bar(xVal, avgArray)
    plt.ylabel('Coverage')
    plt.xlabel('Base Window')
    plt.title('Coverage of Chromosome')
    plt.savefig("coverage.png")
    #graph = _printCoverage(array)
    #ct._resetFromLeaves()
    ct._outputFiltered("output2.bam")












