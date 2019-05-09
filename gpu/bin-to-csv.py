#!/usr/bin/env python3

import numpy as np

data = np.fromfile('transmissions.bin', dtype=np.uint32)
rows = data.shape[0] // 3
data = data.reshape((rows,3))
np.savetxt('transmissions.csv', data, fmt='%u', delimiter=',', header='infectedById,id,timeInfected', comments='')