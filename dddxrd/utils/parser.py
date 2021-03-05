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
    scaninfo['fiofile'] = file
    with open(file) as f:
        for i,line in enumerate(f):
            if line.startswith('fastsweep'):
                scaninfo['nbr_images'] = int(line.split(':')[1].split('/')[0])
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

        #Compute some additional stuff here
        scaninfo['beamheight'] = float(scaninfo['eh3st']+scaninfo['eh3sb']) #only if eh3 slits define the beam
        scaninfo['beamwidth'] = float(scaninfo['eh3si']+scaninfo['eh3so'])
        idx = np.nonzero(scaninfo['data']['type'] == 'exposure')[0][0] # index of first exposure
        scaninfo['omegastep']=float(np.diff(scaninfo['data']['end'])[0])
        scaninfo['startomega'] = float(scaninfo['data']['end'][idx]-0.5*scaninfo['omegastep'])
        scaninfo['filetype'] = '.{}'.format(str(scaninfo['data']['filename']).split('.')[-1][:-2])
        scaninfo['first_image'] = int(scaninfo['data']['filename'][idx].split(scaninfo['filetype'])[0].split('_')[-1]) #this is not general
        # Not general
        stemlist = scaninfo['data']['filename'][idx].split(scaninfo['filetype'])[0].split('_')[:-1] #this is not general
        stem = ''
        for s in stemlist:
            stem += '{}_'.format(s)
        scaninfo['stem'] = stem
    return scaninfo

def print_scan_info(scaninfo):
    """ Print info from a fio file"""
    for keys in scaninfo:
        if (keys == 'data') or keys.startswith('hasep212'):
            continue
        else:
            print(keys,scaninfo[keys])

def init_yaml_from_fio(fiofile,yamlfile,save_file=True,return_dict=False):
    """ Initializes a yaml file with parameters from a fio file"""
    scaninfo = parse_fio(fiofile)
    yaml_dict={}
    yaml_dict['image_path'] = scaninfo['channel1_FileDir1']
    yaml_dict['first_image'] = scaninfo['first_image']
    yaml_dict['image_stem'] = scaninfo['stem']
    yaml_dict['nbr_images'] = scaninfo['nbr_images']
    yaml_dict['filetype'] = scaninfo['filetype']
    yaml_dict['omegastep'] = scaninfo['omegastep']
    yaml_dict['startomega'] = scaninfo['startomega']
    if save_file:
        dump_yaml(yamlfile,yaml_dict)
    if return_dict:
        return yaml_dict
    return
def dump_yaml(yamlfile,yaml_dict):
    with open(yamlfile, 'w') as file:
        yaml.dump(yaml_dict, file)
    return
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
    print_scan_info(info)
    yamlfile = '/Users/al8720/Box/projects/Inconel/peaksearch_fio.yaml'
    init_yaml_from_fio(fiofile,yamlfile)

if __name__ == "__main__":
    main()