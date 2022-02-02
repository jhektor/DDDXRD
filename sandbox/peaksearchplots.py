import matplotlib.pyplot as plt
import numpy as np
import fabio
import dddxrd.utils.parser as parser
from ImageD11 import columnfile
import os
import glob

def max_projection(ims):
    """Computes the maximum projection of all images in a FileSeries"""
    max_im = np.zeros(ims.get_frame(0).data.shape)
    for i in range(0,ims.nframes,1):
        im = ims.get_frame(i).data
        max_im = np.where(im>max_im,im,max_im) 
    return max_im

def make_np_array(ims):
    imlist = []
    for iml in ims:
        print(iml)
        for i in range(899):
            imlist.append(iml.get_frame(i).data)
    return np.array(imlist)
    

yaml = '/Users/al8720/projects/fisk/P21_2021/processed/scripts/3dxrd/peaksearch.yaml'
pars = parser.parse_parameters(yaml)

# Corrections for peak coordinates (stich the detectors)
pixs = 0.15 #pixelsize
npix = 2880
ty = [-180.6/pixs,199.1/pixs,-270.8/pixs,291.83/pixs]
tz = [276.0/pixs,-283.4/pixs,-197.8/pixs,189.2/pixs]

omegastep = 0.2
omegas = np.arange(omegastep/2,180,omegastep)
print(omegas)
#detector = 1
#load data
images = []
peaks = []
for detector in range(1,5):
    if detector ==3:
        digs = 6
    else:
        digs = 5
    #find first image
    impath = '{}{:d}'.format(pars['image_path'],detector)
    imgs = sorted(glob.glob('{}/*.cbf'.format(impath)))
    img0 = imgs[1].split('_')[-1].split('.cbf')[0]  #first frame is clearing
    pars['first_image'] = int(img0)
    first_im = '{path}{detector}/{stem}{imnum:0{digs}}.cbf'.format(path=pars['image_path'],detector=detector,stem=pars['image_stem'],imnum=pars['first_image'],digs=digs) 
    print(os.path.isfile(first_im))
    images.append(fabio.open_series(first_filename=first_im,single_frame=True))
    peakpath = impath.replace('raw','processed')
    peaklist = os.path.join(peakpath,'merged_{:d}.flt'.format(detector))
    c = columnfile.columnfile(peaklist)
    c.sc -= ty[detector-1]
    c.fc -= tz[detector-1]
    peaks.append(c)



for j in [99,100,101,102,103]:
    fig,axs = plt.subplots(2,2)
    axes = axs.ravel()
    fig.suptitle('Omega: {:.3f} - {:.3f}'.format(omegas[j-1],omegas[j]))
    for i in range(4):
        ims = images[i]
        pks = peaks[i]
        #bin omegas
        idx = np.digitize(pks.omega,omegas)
        idxj = np.where(idx==j)[0]
        axes[i].imshow(ims.get_frame(j).data,origin='lower')
        axes[i].scatter(pks.fc[idxj],pks.sc[idxj],s=80,edgecolors='red',facecolors='none')
        axes[i].set_title('Detector {:d}'.format(i))

plt.show()