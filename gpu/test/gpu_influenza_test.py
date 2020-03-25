""" gpu_influenza_test

"""

# Standard packages
import sys
import os
import numpy
import pandas
import unittest
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend("Agg")
from pathlib import Path

# IDM packages
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(str(Path(os.path.join(test_dir, "./../")).resolve().absolute()))
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
modelParams_1_connection = pandas.DataFrame( { "population"  : [500],
                                               "connections" : [  1],
                                               "exposures"   : [  0],
                                               "infections"  : [  1],
                                               "time"        : [500],
                                               #"p"           : [1.0, 1.0]
                                 }
                                )

modelParams_n_connections = pandas.DataFrame( { "population"  : [500],
                                                "connections" : [500],
                                                "exposures"   : [  0],
                                                "infections"  : [  1],
                                                "time"        : [ 70],
                                                #"p"           : [1.0, 1.0]
                                 }
                                )
seed = 101

output_folder = "expected_output"
expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_output.csv"))
sampleI_expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_sampleI_output.csv"))

one_connection_expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_1_connection_output.csv"))
sampleI_one_connection_expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_sampleI_1_connection_output.csv"))

n_connections_expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_n_connections_output.csv"))
sampleI_n_connections_expected_result = pandas.read_csv(os.path.join(output_folder, "gpu_influenza_sampleI_n_connections_output.csv"))


class  TestGpuInfluenza(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        # move intermediate result
        filenames = os.listdir(test_dir)
        test_name = str(self._testMethodName)
        if not os.path.exists(test_name):
            os.mkdir(test_name)
        for item in filenames:
            if item.endswith(".bin"):
                os.rename(os.path.join(test_dir, item), os.path.join(test_dir, test_name, item))

    def test_gpu_influenza(self):
        numpy.random.seed(seed)
        # Run model
        y = gpu_influenza(modelParams)

        plot_gpu_influenza(y, expected_result, plot_name="gpu_influenza.png")

        self.assertTrue(tolerance_test(y, expected_result, 1))

    def test_gpu_influenza_sampleI(self):
        numpy.random.seed(seed)
        # Run model
        ySampled = gpu_influenza_sampleI(modelParams)

        plot_gpu_influenza_sampleI(ySampled, expected_result, sampleI_expected_result, "gpu_influenza_sampleI.png")

        self.assertTrue(tolerance_test(ySampled, sampleI_expected_result, 1))

    def test_gpu_influenza_1_connection(self):
        numpy.random.seed(seed)
        # Run model
        y = gpu_influenza(modelParams_1_connection)

        plot_gpu_influenza(y, one_connection_expected_result, plot_name="gpu_influenza_1_connection.png")

        self.assertTrue(y.equals(one_connection_expected_result))

    def test_gpu_influenza_sampleI_1_connection(self):
        numpy.random.seed(seed)
        # Run model
        ySampled = gpu_influenza_sampleI(modelParams_1_connection)
        # with open("gpu_influenza_sampleI_1_connection_output.csv", "w") as file:
        #     ySampled.to_csv(file, index=False, header=True, line_terminator="\n")

        plot_gpu_influenza_sampleI(ySampled, one_connection_expected_result, sampleI_one_connection_expected_result,
                                   "gpu_influenza_1_connection_sampleI.png")

        self.assertTrue(ySampled.equals(sampleI_one_connection_expected_result))

    def test_gpu_influenza_n_connections(self):
        numpy.random.seed(seed)
        # Run model
        y = gpu_influenza(modelParams_n_connections)

        plot_gpu_influenza(y, n_connections_expected_result, plot_name="gpu_influenza_n_connection.png")

        self.assertTrue(tolerance_test(y, n_connections_expected_result, 20))

    def test_gpu_influenza_sampleI_n_connections(self):
        numpy.random.seed(seed)
        # Run model
        ySampled = gpu_influenza_sampleI(modelParams_n_connections)

        plot_gpu_influenza_sampleI(ySampled, n_connections_expected_result, sampleI_n_connections_expected_result,
                                   "gpu_influenza_n_connections_sampleI.png")

        self.assertTrue(tolerance_test(ySampled, sampleI_n_connections_expected_result, 20))


def plot_gpu_influenza(y, expected_result, plot_name):
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
    plt.savefig(plot_name)
    plt.close()


def plot_gpu_influenza_sampleI(ySampled, expected_result, sampleI_expected_result, plot_name):
    plt.figure()  # Sampled output
    ax1 = plt.subplot(2, 1, 1)
    expected_result.filter(regex="^I", axis=1).transpose().plot(title="infectious", ax=ax1)

    ax2 = plt.subplot(2, 1, 2)
    ySampled.transpose().plot(title="infectious (sampled)", ax=ax2)
    sampleI_expected_result.transpose().plot(title="infectious (sampled)", ax=ax2)

    plt.tight_layout(pad=3.0)
    plt.savefig(plot_name)
    plt.close()


def tolerance_test(y, expected_result, tolerance):
    diff = expected_result.subtract(y)
    # print(diff.abs().max(axis=1))
    res = diff.abs() <= tolerance
    return res.all(axis=None)


if __name__ == '__main__':
    unittest.main()

