""" gpu_influenza_test

"""

# Standard packages
import numpy
import pandas
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend("Agg")
import unittest

# IDM packages
import gpuflu
from gpuflu.interface import gpu_influenza as gpu_influenza
from gpuflu.interface import gpu_influenza_sampleI as gpu_influenza_sampleI


# Parameters
modelParams = pandas.DataFrame( { "population"  : [100, 200],
                                  "connections" : [  5, 4],
                                  "exposures"   : [  0, 1],
                                  "infections"  : [  3, 5],
                                  "time"        : [ 40, 30],
                                  #"p"           : [1.0, 1.0]
                                 }
                                )
seed = 101
expected_result = pandas.read_csv("gpu_influenza_output.csv")
sampleI_expected_result = pandas.read_csv("gpu_influenza_sampleI_output.csv")


class  TestGpuInfluenza(unittest.TestCase):
    def test_gpu_influenza(self):
        numpy.random.seed(seed)
        # Run model
        y = gpu_influenza(modelParams)

        # Draw some plots
        plt.figure()  # Full model output
        ax1 = plt.subplot(3, 1, 1)
        y.filter(regex="^E", axis=1).transpose().plot(title="exposed", ax=ax1)
        expected_result.filter(regex="^E", axis=1).transpose().plot(title="exposed", ax=ax1)

        ax2 = plt.subplot(3, 1, 2)
        y.filter(regex="^I", axis=1).transpose().plot(title="infectious", ax=ax2)
        expected_result.filter(regex="^I", axis=1).transpose().plot(title="infectious", ax=ax2)

        ax3 = plt.subplot(3, 1, 3)
        y.filter(regex="^R", axis=1).transpose().plot(title="recovered", ax=ax3)
        expected_result.filter(regex="^R", axis=1).transpose().plot(title="recovered", ax=ax3)

        plt.tight_layout(pad=3.0)
        plt.savefig("gpu_influenza.png")
        plt.close()

        self.assertTrue(y.equals(expected_result))

    def test_gpu_influenza_sampleI(self):
        numpy.random.seed(seed)
        # Run model
        ySampled = gpu_influenza_sampleI(modelParams)

        plt.figure()  # Sampled output
        ax1 = plt.subplot(2, 1, 1)
        expected_result.filter(regex="^I", axis=1).transpose().plot(title="infectious", ax=ax1)

        ax2 = plt.subplot(2, 1, 2)
        ySampled.transpose().plot(title="infectious (sampled)", ax=ax2)
        sampleI_expected_result.transpose().plot(title="infectious (sampled)", ax=ax2)

        plt.tight_layout(pad=3.0)
        plt.savefig("gpu_influenza_sampleI.png")
        plt.close()

        assert ySampled.equals(sampleI_expected_result)


if __name__ == '__main__':
    unittest.main()

