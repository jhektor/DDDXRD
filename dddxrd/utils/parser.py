import json
import yaml
import os

def parse_fio(file):
    """ Parse a .fio file from the fastsweep macro at P21.2"""


def parse_parameters(file):
    """ Parse parameters stored in a json file. Returns a dictionary"""
    if file.endswith('.json'):
        parameters = _parse_json(file)
    elif file.endswith(('.yml','.yaml')):
        parameters = _parse_yaml(file)
    return parameters

def _parse_json(file):
    with open(file) as f:
        pars = json.load(f)
    return pars

def _parse_yaml(file):
    with open(file) as f:
        pars = yaml.safe_load(f)
    return pars

def find_padding(pars):
    """ Finds the padding of image filenames, im_0001 return 4"""
    ndigits = None
    for i in range(7):
        pth = os.path.join(pars['image_path'],pars['image_stem'])
        first_im_name = '{}{:0>{pad}}{}'.format(pth,pars['first_image'],pars['filetype'],pad=i)
        if os.path.exists(first_im_name):
            ndigits = i
            break
    if ndigits is None:
        raise RuntimeError('format string problem, revisit')
    return ndigits

    