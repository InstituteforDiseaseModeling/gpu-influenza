import logging
log = logging.getLogger(__name__)
import numpy as np
import numpy.random as ran

def get_value_from_distribution(params):

    params_distribution = params["distribution"]
    result = 0
    if params_distribution == "normal":
        loc = params["loc"]
        scale = params["scale"]
        result = ran.normal(loc, scale)
    elif params_distribution == "exponential":
        scale = params["scale"]
        result = ran.exponential(scale)
    elif params_distribution == "fixed":
        result = params["value"]
    elif params_distribution == "poisson":
        result = ran.poisson(params["mean"])
    elif params_distribution == "lognormal":
        result = ran.lognormal(params["mean"], params["sigma"])
    elif params_distribution == "binned_uniform":
        #binned uniform distribution
        #find bin
        ran_num = [ ran.random() ]
        bin_value = np.digitize(ran_num, params["cumulative_p"])[0]
        result = ran.uniform(params["bins"][bin_value], params["bins"][bin_value+1])
    elif params_distribution == "choice":
        result = ran.choice(params["selection"], p=params["p"])

    else:
        raise RuntimeError('error: unknown distribution.')

    return result
