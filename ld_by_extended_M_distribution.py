import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fftengine

import astropy.io.fits as pyfits


class deflector(object):
    """
    initialize the deflector using a surface density (convergence) map
    the boolean variable pad indicates whether zero-padding is used
    or not
    """
    def __init__(self,filekappa,pad=False):
        kappa,header=pyfits.getdata(filekappa,header=True)
        self.kappa=kappa
        self.nx=kappa.shape[0]
        self.ny=kappa.shape[1]
        self.pad=pad
        if (pad):
            self.kpad()
        self.kx,self.ky=self.kernel()

    '''
    implement the kernel function Kimplement the kernel function K
    '''
    def kernel(self):
        x = np.linspace(-0.5,0.5,self.kappa.shape[0])
        y = np.linspace(-0.5,0.5,self.kappa.shape[1])
        kx,ky = np.meshgrid(x,y)
        norm = (kx**2 + ky**2 +1e-12)
        kx = kx / norm
        ky = ky / norm
        return (kx,ky)

    '''
    compute the deflection angle maps by convolving
    the surface density with the kernel function
    Note that the returned values will be in pixel units
    '''
    def angles(self):
        # FFT of the surface density and of the two components of the kernel
        kappa_ft = fftengine.rfftn(self.kappa,axes=(0,1))
        kernelx_ft = fftengine.rfftn(self.kx,axes=(0,1),s = self.kappa.shape)
        kernely_ft = fftengine.rfftn(self.ky,axes=(0,1),s = self.kappa.shape)
        # perform the convolution in Fourier space and transform the result
        # back in real space. Note that a shift needs to be applied using fftshift
        alphax = 1.0/np.pi * fftengine.fftshift(fftengine.irfftn(kappa_ft*kernelx_ft))
        alphay = 1.0/np.pi * fftengine.fftshift(fftengine.irfftn(kappa_ft*kernely_ft))
        return (alphax,alphay)
    
    '''
    returns the surface-density (convergence) of the deflector
    '''
    def kmap(self):
        return (self.kappa)

    '''
    performs zero-padding
    '''
    def kpad(self):
        # add zeros around the original array
        def padwithzeros(vector,pad_width,iaxis,kwargs):
            vector[:pad_width[0]] = 0
            vector[-pad_width[1]:] = 0
            return vector
        # use the pad method from numpy.lib to add zeros (padwithzeros)
        # in a frame with thicksness self.kappa.shape[0]
        self.kappa=np.lib.pad(self.kappa, self.kappa.shape[0],padwithzeros)

    '''
    crop the maps to remove zero-padded areas and get back to the original region.
    '''
    def mapCrop(self,mappa):
        xmin = np.int(0.5*(self.kappa.shape[0]-self.nx))
        ymin = np.int(0.5*(self.kappa.shape[1]-self.ny))
        xmax=xmin+self.nx
        ymax=ymin+self.ny
        mappa=mappa[xmin:xmax,ymin:ymax]
        return(mappa)

df=deflector('data/kappa_2.fits')
angx_nopad,angy_nopad=df.angles()
kappa=df.kmap()

fig,ax = plt.subplots(1,3,figsize=(16,8))
ax[0].imshow(kappa,origin="lower",cmap = 'jet')
ax[0].set_title('convergence')
ax[1].imshow(angx_nopad,origin="lower",cmap = 'jet')
ax[1].set_title('angle 1')
ax[2].imshow(angy_nopad,origin="lower",cmap = 'jet')
ax[2].set_title('angle 2')
plt.show()

'''
ddf=deflector('data/kappa_2.fits',pad=True)
angx,angy=df.angles()
kappa=df.kmap()

angx=df.mapCrop(angx)
angy=df.mapCrop(angy)

fig,ax = plt.subplots(1,3,figsize=(16,8))
ax[0].imshow(kappa,origin="lower",cmap = 'jet')
ax[0].set_title('convergence')
ax[1].imshow(angx,origin="lower",cmap = 'jet')
ax[1].set_title('angle 1')
ax[2].imshow(angy,origin="lower",cmap = 'jet')
ax[2].set_title('angle 2')
plt.show()

'''

