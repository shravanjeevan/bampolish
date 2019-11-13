import sys
import pysam
import numpy as np
import math
from tqdm import tqdm

#TODO: change sumarray to averagearray and minarry (queries from [x, x+1) should still return correct coverage in O(logn), can direct access array for O(1))
class CoverageTree:
    WINDOW_POWER = 0 #How many times to divide resolution by 2
    MAX_CHR_LENGTH = 250000000 >> WINDOW_POWER
    MAX_N = 2*pow(2, math.ceil(math.log(MAX_CHR_LENGTH, 2)))
    def __init__(self, filename):
        self.file = pysam.AlignmentFile(filename, "r")
        
        #Initiallize lazy nodes
        self.sumarray = np.zeros(self.MAX_N)
        self.lazyadd = np.zeros(self.MAX_N)
        print("Constructing tree, chromosome scaled to " + str(self.MAX_CHR_LENGTH))

        #Range updates
        for r in tqdm(self.file.fetch("chr1")):
            if(r.is_unmapped): continue
            if(r.reference_end >= self.MAX_CHR_LENGTH): continue
            self._update(r.reference_start >> self.WINDOW_POWER, (r.reference_end + 1) >> self.WINDOW_POWER, 1)

    def _recalculate(self, id, l, r):
        self.sumarray[id] = self.lazyadd[id] * (r-l)
        if(r-l != 1):
            self.sumarray[id] += self.sumarray[id * 2]
            self.sumarray[id] += self.sumarray[id * 2 + 1]

    def _update_lazy(self, id, v, l, r):
        self.lazyadd[id] += v
        self._recalculate(id, l, r)

    def _propagate(self, id, l, r):
        mid = (l + r)/2
        self._update_lazy(id * 2, self.lazyadd[id], l, mid)
        self._update_lazy(id * 2 + 1, self.lazyadd[id], mid, r)
        self.lazyadd[id] = 0

    def _update(self, uL, uR, v, i=1, cLeft=0, cRight=MAX_CHR_LENGTH):
        if(uL == cLeft and uR == cRight):
            self._update_lazy(i, v, cLeft, cRight)
            return
        self._propagate(i, cLeft, cRight)
        mid = math.floor((cLeft + cRight) / 2)
        if(uL < mid): self._update(uL, min(uR, mid), v, i * 2, cLeft, mid)
        if(uR > mid): self._update(max(uL, mid), uR, v, i * 2 + 1, mid, cRight)
        self._recalculate(i, cLeft, cRight)
        
    def _query(self, qL, qR, i=1, cLeft=0,  cRight=MAX_CHR_LENGTH):
        if(qL == cLeft and qR == cRight):
            return sumarray[i]
        self._propagate(i, cLeft, cRight)
        mid = (cLeft + cRight) / 2
        ans = 0
        if(qL < mid): ans += self._query(qL, min(qR, mid), i * 2, cLeft, mid)
        if(qR > mid): ans += self._query(max(qL, mid), qR, i * 2 + 1, mid, cRight)
        return ans
        
if __name__ == '__main__':
    ct = CoverageTree(sys.argv[1])
