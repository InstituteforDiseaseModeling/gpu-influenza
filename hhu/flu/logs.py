import logging
log = logging.getLogger(__name__)

def setup_logging(default_name, logging_params):

    levels = {}
    levels["debug"]=logging.DEBUG
    levels["info"]=logging.INFO
    levels["warning"]=logging.WARNING
    levels["error"]=logging.ERROR
    levels["critical"]=logging.CRITICAL

    for key in logging_params.keys():

        level_name = 'flu' if key == 'default' else 'flu.'+key

        log = logging.getLogger(level_name)
        log.setLevel(levels[logging_params[key]])
        log = logging.getLogger(__name__)
        log.info( 'set {0} logLevel to {1}.'.format(key, logging_params[key]) )

    if "default" not in logging_params.keys():
       log = logging.getLogger(__name__)
       log.info("default log level not set, use info.")

    return
