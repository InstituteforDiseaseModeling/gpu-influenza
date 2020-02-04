#!/usr/bin/env python3

import numpy as np

data = np.fromfile('transmissions.bin', dtype=np.uint32)
rows = data.shape[0] // 3
data = data.reshape((rows,3))
np.savetxt('transmissions.csv', data, fmt='%u', delimiter=',', header='infectedById,id,timeInfected', comments='')

data = np.fromfile('recoveryCounts.bin', dtype=np.uint32)
rows = data.shape[0]
np.savetxt('recoveryCounts.csv', np.dstack((np.arange(0,rows), data))[0], fmt='%u', delimiter=',', header='day,count', comments='')
