import json
import yaml
import os
import numpy as np

def parse_fio(file,nchannels=1):
    """ Parse a .fio file from the fastsweep macro at P21.2. 
    Not tested on multiple channel sweeps"""
    scaninfo = {} #dictonary to store scan info in
    colnames=[]
    colformats=[]
    with open(file) as f:
        for i,line in enumerate(f):
            if line.startswith('fastsweep'):
                scaninfo['nbr_images'] = line.split(':')[1].split('/')[0]
            elif 'acquisition started' in line:
                scaninfo['starttime'] = line.split('started at')[-1]
            elif ' = ' in line:
                kv = line.split(' = ')
                try:
                    scaninfo[kv[0]]=float(kv[1])
                except ValueError:
                    scaninfo[kv[0]]=kv[1].rstrip()
            elif line.startswith(' Col '):
                colnames.append(line.split(' ')[3])
            elif line.startswith('  '): #first line of data
                break
        headerlines=i
        #read the data
        scaninfo['data'] = np.genfromtxt(file,dtype=None, names=colnames,skip_header=headerlines,encoding=None)     

    return scaninfo

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

def main():
    fiofile = '/Users/al8720/Box/projects/Inconel/eh3scan1_00055.fio'
    info = parse_fio(fiofile)
    for keys in info:
        print(keys,info[keys])

if __name__ == "__main__":
    main()