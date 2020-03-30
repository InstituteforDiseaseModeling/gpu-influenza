from setuptools import setup
import os


# Todo: config block size, nvidia architecture
os.system( "make" )

setup( name        = "gpuflu",
       version     = "0.1",
       description = "Network-based compartmental model for NVIDIA GPUs",
       url         = "https://github.com/InstituteforDiseaseModeling/gpu-influenza/",
       packages    = ["gpuflu"],
       zip_safe    = False,
       install_requires = [ "numpy",
                            "scipy",
                            "pandas",
                            "matplotlib",
                            "seaborn",
                            "mkl",
                            "jupyter",
                            "jupyter-core",
                            "jupyter-console",
                            "jupyter-client"
                           ],
       include_package_data = True
     )
