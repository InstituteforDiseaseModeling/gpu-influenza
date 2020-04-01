#! /usr/bin/env python3

import numpy as np

data = np.fromfile('transmissions.bin', dtype=np.int32)
rows = data.shape[0] // 4
data = data.reshape((rows,4))
np.savetxt('transmissions.csv', np.hstack((np.arange(0,rows).reshape((rows,1)),data)), fmt='%u', delimiter=',', header='id,infectedById,timeExposed,timeInfectious,timeRecovered', comments='')
