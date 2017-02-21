import os


def locate_resource_file(fname):
    return os.environ['HOME'] + '/rif_data/' + fname
