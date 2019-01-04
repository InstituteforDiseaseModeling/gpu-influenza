#!/usr/bin/env python

import yaml
import argparse
import os

print('''
=============================
Matrix-based Individual Model
Flu Prototype
=============================
''')

logger = None

# @profile
def main(configfile):

    from flu.simulation import Simulation
    from flu.logs import setup_logging
    from flu.resource import check_resources
    from flu.report import Report

    #read configs
    configs = yaml.load(open(configfile))

    #setup logging levels
    setup_logging("flu", configs["logging"])

    #check cpu and memory resources
    check_resources(configs)

    sim = Simulation(configs)
    report = Report(sim)

    #check output directory and create if it doesn't exist
    try:
        os.makedirs(configs["reports"]["output_directory"])
    except OSError:
        if not os.path.isdir(configs["reports"]["output_directory"]):
            raise

    if configs["reports"]["plots"]["enabled"] == True:
        #create the output directory
        #check plot directory and create if it doesn't exist
        try:
            os.makedirs(configs["reports"]["plots"]["output_directory"])
        except OSError:
            if not os.path.isdir(configs["reports"]["plots"]["output_directory"]):
                raise

    #simulation starts here

    logger.info("simulation start.")

    sim.populate()

    seeding_time = configs["simulation"]["seeding_time"]

    for t in range(configs["simulation"]["total_time"]):

        if t in seeding_time:
            sim.seed()

        sim.update(t)
        report.update(t)

    logger.info("simulation end.")

    #reports and plots
    report.write_reports()

    if configs["reports"]["plots"]["enabled"] == True:
        report.write_plots()

    return


def configure_logging(logfile):

    import logging
    global logger

    logger = logging.getLogger("flu")
    logger.setLevel(logging.INFO)

    if logfile: # != "":
        ch = logging.FileHandler(logfile, mode='w')
    else:
        ch = logging.StreamHandler()

    formatter = logging.Formatter("%(asctime)s | %(name)-18s| %(levelname)-8s| %(message)s", "%H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog='flu_prototype',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="Flu Prototype")
    parser.add_argument('--config', '-C', metavar='CONFIG_FILE_NAME', help="config file name", default="config.yaml")
    parser.add_argument('--log', '-L', metavar='LOG_FILE_NAME', help="log file name", default="")
    input_args = parser.parse_args()

    configure_logging(input_args.log)

    main(input_args.config)
