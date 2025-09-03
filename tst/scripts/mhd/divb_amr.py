import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True

logger = logging.getLogger('athena' + __name__[7:])

def run(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.run('mhd/ffc_divb.athinput')

def analyze():
    logger.debug('Analyzing test ' + __name__)
    data = athena_read.athdf('divb_amr.out1.00001.athdf')
    bx = data['B1c']
    by = data['B2c']
    x1v = data['x1v']
    x2v = data['x2v']
    dx1 = x1v[1] - x1v[0]
    dx2 = x2v[1] - x2v[0]
    divb = np.gradient(bx, dx1, axis=2) + np.gradient(by, dx2, axis=1)
    max_divb = np.max(np.abs(divb))
    logger.debug('Max divB %e', max_divb)
    return max_divb < 1e-10
