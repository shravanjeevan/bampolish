import sys
import pysam
import numpy as np
import math
from tqdm import tqdm
import pdb
import random
import matplotlib.pyplot as plt

class CoverageTree:
    CUTOFF_MEDMAD = 1
    CUTOFF_SD = 2

    def __init__(self, inFile, ref, total, windowPower, verboseFlag, progressbarFlag, ndeviation, stddeviationflag):
        self.file = inFile
        self.ref = ref
        self.ndeviation = ndeviation
        self.stddeviationflag = stddeviationflag

        # Number of times the chromosome length and read index resolutions are halved
        self.WINDOW_POWER = windowPower

        self.WINDOW_SIZE = math.floor(math.pow(2, self.WINDOW_POWER))
        self.verboseFlag = verboseFlag
        self.progressbarFlag = progressbarFlag

        # Maximum length of chromosome (~250 million bases in human)
        # Should pass in length too, this calculation is slow
        self.CHR_LENGTH = self.file.lengths[self.file.references.index(ref)] >> self.WINDOW_POWER
        self.MAX_CHR_LENGTH = pow(2, math.ceil(math.log(self.CHR_LENGTH, 2)))

        # Size allocation of the binary tree (padded to the next smallest power of two for balance)
        self.MAX_N = 2*pow(2, math.ceil(math.log(self.MAX_CHR_LENGTH, 2)))

        if self.verboseFlag:
            print("Constructing tree")
            print("(Chromosome scaled to " + str(self.MAX_CHR_LENGTH) + ")")

        # Stores the minimum of all subtree coverages
        self.minarray = np.zeros(self.MAX_N, dtype=np.int32)

        # Used to prevent inefficient propagation during queries
        self.lazyadd = np.zeros(self.MAX_N, dtype=np.int32)

        # Iterate over all mapped reads in ref
        if self.verboseFlag:
            print("===== Processing " + str(ref) + " =====")
        for r in tqdm(self.file.fetch(ref), total=total, disable=not self.progressbarFlag):
            if(r.is_unmapped): continue

            # Scale reads to reduce resolution
            readStart = r.reference_start >> self.WINDOW_POWER
            readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

            # Ensure interval is inclusive-exclusive
            if(readStart == readEnd): readEnd += 1

            # Increment all bases in range by one
            self._update(readStart, readEnd, 1)

        if self.verboseFlag:
            print("Pushing range updates to leaves")
        self.t = tqdm(total=self.MAX_CHR_LENGTH)

        self.leaves = np.empty(self.MAX_CHR_LENGTH, dtype=np.int32)
        self._updateAll()
        self.t.close()
        del self.t

        # Calculate median of mapped bases only
        self.median = self._getMedian()
        if self.verboseFlag:
            print("(Median = " + str(self.median) + ")")

        self.mad = self._getMad()
        if self.verboseFlag:
            print("(MAD = " + str(self.mad) + ")")

        self.mean = self._getMean()
        self.sd = self._getSD()


    def _initRefFromFile():
        return

    def _initRefFromCigar():
        return

    # Updates the sum/min arrays according to the lazy counter
    def _recalculate(self, id, l, r):
        self.minarray[id] = self.lazyadd[id]
        if(r-l != 1):
            self.minarray[id] += min(self.minarray[id * 2], self.minarray[id * 2 + 1])

    # Updates lazy counters
    def _update_lazy(self, id, v, l, r):
        self.lazyadd[id] += v
        self._recalculate(id, l, r)

    # Pushes the lazy counter down the tree, recalculating as it goes
    def _propagate(self, id, l, r):
        mid = (l + r)/2
        self._update_lazy(id * 2, self.lazyadd[id], l, mid)
        self._update_lazy(id * 2 + 1, self.lazyadd[id], mid, r)
        self.lazyadd[id] = 0

    # Increments all base coverages within the given interval
    def _update(self, uL, uR, v, i=1, cLeft=0, cRight=None):
        cRight = cRight or self.MAX_CHR_LENGTH
        if(uL == cLeft and uR == cRight):
            self._update_lazy(i, v, cLeft, cRight)
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        if(uL < mid): self._update(uL, min(uR, mid), v, i * 2, cLeft, mid)
        if(uR > mid): self._update(max(uL, mid), uR, v, i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)

    # Non lazy update of tree, so that leaves are correct
    def _updateAll(self, i=1, cLeft=0, cRight=None):
        cRight = cRight or self.MAX_CHR_LENGTH
        if(cRight - cLeft == 1):
            self._recalculate(i, cLeft, cRight)
            self.t.update(1)
            #self.leaves[cLeft] = self.minarray[i]
            #pdb.set_trace()
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        self._updateAll(i * 2, cLeft, mid)
        self._updateAll(i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)

    # Returns sum and minimum coverages within the given interval
    def _query(self, qL, qR, i=1, cLeft=0,  cRight=None):
        cRight = cRight or self.MAX_CHR_LENGTH
        if(qL == cLeft and qR == cRight):
            return self.minarray[i]
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        minAns = math.inf
        if(qL < mid): minAns = min(minAns, self._query(qL, min(qR, mid), i * 2, cLeft, mid))
        if(qR > mid): minAns = min(minAns, self._query(max(qL, mid), qR, i * 2 + 1, mid, cRight))
        return minAns

    # TODO: Fix changes not getting pushed to leaves (which is good but can't extract this way then)
    # Returns a modifiable numpy subarray of the leaves, which correspond to individual base position coverages
    def _getLeaves(self):
        return self.minarray[self.MAX_N-self.MAX_CHR_LENGTH : self.MAX_N-self.MAX_CHR_LENGTH+self.CHR_LENGTH]
        #return self.leaves

    # Reconstructs tree in linear time from modified leaves, bottom up
    def _resetFromLeaves(self, i=1, cLeft=0, cRight=None):
        cRight = cRight or self.MAX_CHR_LENGTH
        self.lazyadd[i] = 0
        if(cRight - cLeft == 1):
            return
        if(uL < mid): self._resetFromLeaves(i * 2, cLeft, mid)
        if(uR > mid): self._resetFromLeaves(i * 2 + 1, mid, cRight)

        self.minarray[i] = min(self.minarray[i*2], self.minarray[i*2+1])

    def _getMean(self):
        leaves = self._getLeaves()
        return np.average(leaves[leaves > 0])

    def _getSD(self):
        coverages = self._getLeaves()[self._getLeaves() > 0]
        squareDiffSum = 0
        for i in range(len(coverages)):
            squareDiffSum += math.pow(coverages[i] - self.mean, 2)
        return math.sqrt(squareDiffSum/len(coverages))

    def _getMedian(self):
        leaves = self._getLeaves()
        leavesSorted = np.sort(leaves[leaves > 0])
        return 0 if len(leavesSorted) == 0 else leavesSorted[math.floor(len(leavesSorted)/2)]

    # Calculates the median absolute deviation of base coverages
    def _getMad(self):
        coverages = self._getLeaves()[self._getLeaves() > 0]
        deviations = np.empty(len(coverages), dtype=np.int32)

        for i in range(len(coverages)):
            deviations[i] = abs(self.median - coverages[i])

        deviations.sort()
        return deviations[math.floor(len(deviations)/2)]

    # Writes high coverage intervals to a bed file
    def _outputSpikes(self, filename):
        baseline = self.mean if self.stddeviationflag else self.median
        deviation = self.sd if self.stddeviationflag else self.mad
        threshHi = baseline + self.ndeviation * deviation
        threshLo = baseline - self.ndeviation * deviation

        spikeStart = 0
        spikeCurr = 0
        cumsum = 0
        wasHiSpike = False
        isHiSpike = False
        wasLoSpike = False
        isLoSpike = False
        spike_id = 0
        leaves = self._getLeaves()
        with open(filename,'w') as outFile:
            for c in leaves:
                isHiSpike = (c > threshHi)
                isLoSpike = (c < threshLo)
                if((isHiSpike and not wasHiSpike)):# or (isLoSpike and not wasLoSpike)):
                    #Start spike
                    spikeStart = spikeCurr
                    cumsum = 0
                elif((isHiSpike and spikeCurr == len(self._getLeaves())-1) or (wasHiSpike and not isHiSpike)):# or (wasLoSpike and not isLoSpike)):
                    #End spike
                    outFile.write(self.ref + '\t' + str(spikeStart << self.WINDOW_POWER) + '\t' + str(spikeCurr << self.WINDOW_POWER) + '\t' + "spike_" + str(spike_id) + '\t' + str(cumsum/(spikeCurr - spikeStart)) + '\n')
                    spike_id += 1
                wasHiSpike = isHiSpike
                wasLoSpike = isLoSpike
                spikeCurr += 1
                cumsum += c

    # Writes reads to a file without going below the median where possible
    # Decay (0, 1] increases the reads included above the otherwise strict cutoff
    def _outputFiltered(self, outFile, decay=0):
        baseline = self.mean if self.stddeviationflag else self.median
        deviation = self.sd if self.stddeviationflag else self.mad
        threshHi = baseline + self.ndeviation * deviation

        if self.verboseFlag:
            print("Writing filtered reads to output file")

        for ref in self.file.references:
            total = self.file.count(reference=ref, until_eof=True)
            if total == 0: continue
            if self.verboseFlag:
                print("===== Processing " + str(ref) + " =====")
            for r in tqdm(self.file.fetch(ref), total=total, disable=not self.progressbarFlag):
                if r.is_unmapped: continue

                # Scale reads to reduce resolution
                readStart = r.reference_start >> self.WINDOW_POWER
                readEnd = (r.reference_end + 1) >> self.WINDOW_POWER

                # Ensure interval is inclusive-exclusive
                if(readStart == readEnd): readEnd += 1

                # Write if the read is part of minimum
                minAns = self._query(readStart, readEnd)

                if(minAns <= (threshHi if decay <= 0 else threshHi + math.log(threshHi - minAns)/decay)):
                    outFile.write(r)
                else:
                    self._update(readStart, readEnd, -1)
