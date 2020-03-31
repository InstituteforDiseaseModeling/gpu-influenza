from builtins import str
import logging
log = logging.getLogger(__name__)

import psutil

def check_resources(params):
    #check available CPU and memory
    number_cpus = psutil.cpu_count()
    log.info("CPU = "+str(number_cpus))
    meminfo = psutil.virtual_memory()
    log.info("MEM = "+str(meminfo.total >> 30)+" GB, available = "+str(meminfo.available >> 30)+" GB.")
