import os


def homedir():
    return os.environ['HOME'] if 'HOME' in os.environ else '/'


def resource_path(path):

    candidate_dirs = [
        os.path.join(homedir(), 'rif_data'),
        os.path.join(homedir(), 'data'),
        os.path.join(homedir(), 'data/rif'),
    ]
    for d in candidate_dirs:
        tryme = os.path.join(d, path)
        print('trial:', tryme)
        if os.path.exists(tryme):
            return tryme
    raise IOError(path + ' not found in any known data location')
