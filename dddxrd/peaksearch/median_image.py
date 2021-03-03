import numpy as np
from PIL import Image
import fabio
from dddxrd.utils.parser import parse_parameters
import glob
import os

def median_image(par_file=None, pars=None):
    """ Computes the median image using parameters in par_file"""
    if (pars is None) and (par_file is not None):
        pars = parser.parse_parameters(par_file)
    elif (pars is None) and (par_file is None):
        print('Must supply either par_file or pars to run_peaksearcher')
        raise ValueError
    print("Looking for images in: {}".format(pars['image_path']))
    glob_pattern = os.path.join(pars['image_path'],pars['stem'])
    img_list = glob.glob(glob_pattern)
    print("Found {:d} images".format(len(img_list)))
    imgs = []
    for im in img_list:
        apath=os.path.abspath(im)
        imgs.append(fabio.open(apath).data)
    "compute median image"
    mim = np.median(imgs,axis=0)
    tiff = Image.fromarray(mim.astype(np.uint16))
    outfile = os.path.abspath(os.path.join(pars['output_path'],pars['output_stem']))
    tiff.save(outfile)
    print('Saved median image in {}'.format(outfile))
    return tiff

def main():
    par_file = 'dddxrd/tests/median_image.yaml'
    median_image(par_file)

if __name__ == "__main__":
    main()