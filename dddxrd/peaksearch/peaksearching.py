import os
import subprocess
import dddxrd.utils.parser as parser

def run_peaksearch(par_file):
    """ Wrapper for the ImageD11 peaksearch.py script"""
    pars = parser.parse_parameters(par_file)
    first_im = pars['first_image']
    last_im = pars['first_image'] + pars['nbr_images'] - 1
    ndigits = parser.find_padding(pars)
    impath = os.path.join(pars['image_path'],pars['image_stem'])
    outpath = os.path.join(pars['output_dir'],'peaks')
    os.makedirs(outpath, exist_ok=True)
    # construct the command for peaksearch.py
    command = ('peaksearch.py -n {} -F {} -f {:d} -l {:d} -o {} -d {} -p Y --ndigits {:d} -S {:.3f} -T {:.3f} '.format(
        impath,pars['filetype'],first_im,last_im,outpath,
        pars['dark_image'],ndigits,pars['omegastep'], pars['startomega']
        ))
    # Adds threshold values to command
    for t in pars['thresholds']:
        command += '-t {:d} '.format(t)
    # Adds keyword args
    if 'kwargs' in pars:
        command += '{} '.format(pars['kwargs'])
    print('Running peaksearch with the following command:')
    print(command)
    try:
        subprocess.call(command, shell=True)
    except AttributeError as a:
        print('peaksearch.py ended with error. It seems to work nonetheless.', a)
    return

def merge_peaks(par_file):
    ''' Wrapper for ImageD11 merge_flt.py'''
    pars = parser.parse_parameters(par_file)
    if 'merged_name' in pars:
        outfile = os.path.join(pars['output_dir'],pars['merged_name'])
    else:
        outfile = os.path.join(pars['output_dir'],'merged.flt')
    inp = os.path.join(pars['output_dir'],'peaks')
    print('Merging flt files matching {}'.format(inp))
    if 'par_file' in pars:
        par_file = pars['par_file']
    else:
        par_file = 'junk'
    command = 'merge_flt.py {} {} {} {:d} '.format(par_file,inp,outfile,pars['pixel_tol']) + ('{:d} '*len(pars['thresholds'])).format(*pars['thresholds'])
    print(command)
    subprocess.call(command, shell=True)

def main(par_file):
    run_peaksearch(par_file)
    merge_peaks(par_file)
    return

if __name__ == "__main__":
    parfile = 'dddxrd/tests/peaksearch.yaml'
    main(parfile)
