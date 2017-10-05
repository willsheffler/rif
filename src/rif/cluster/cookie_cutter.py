from _rif._cluster._cookie_cutter import *
import numpy as np

metrics = dict()
metrics['L2'] = lambda x, y: np.sqrt(np.sum((x - y)**2))
metrics['hamming'] = lambda x, y: np.sum(x != y) / x.shape[0]
metrics['ndiff'] = lambda x, y: np.sum(x != y)


def cookie_cutter(x, thresh, metric='L2'):
    "use the cpp version, it's 20x faster"
    if not isinstance(x, np.ndarray):
        x = x.as_matrix()
    metric = metrics[metric]
    centers = list()
    for i, r in enumerate(x):
        for j, c in centers:
            if metric(r, c) <= thresh:
                break
        else:
            centers.append((i, r))
    return np.array([i for i, r in centers])
