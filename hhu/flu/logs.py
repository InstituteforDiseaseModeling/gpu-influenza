import logging
log = logging.getLogger(__name__)

def setup_logging(default_name, logging_params):
    #setup logging

    levels = {}
    levels["debug"]=logging.DEBUG
    levels["info"]=logging.INFO
    levels["warning"]=logging.WARNING
    levels["error"]=logging.ERROR
    levels["critical"]=logging.CRITICAL

    for key in list(logging_params.keys()):
        if key == "default":
            level_name = "flu"
        else:
            level_name = "flu."+key
        log = logging.getLogger(level_name)
        log.setLevel(levels[logging_params[key]])
        log = logging.getLogger(__name__)
        log.info("set "+ key + " logLevel to "+logging_params[key]+".")

    if "default" not in list(logging_params.keys()):
       log = logging.getLogger(__name__)
       log.info("default log level not set, use info.")

    return
