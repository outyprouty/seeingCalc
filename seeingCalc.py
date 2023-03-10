from matplotlib import pyplot as plt
from glob import glob
from astropy.io import fits
from photutils import centroids
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button,TextBox
from numpy import log, exp, pi
import numpy
from scipy import optimize as opt
import logging

from sys import argv

numpy.seterr(divide = 'ignore')
rad2arcsec = 3600*180/pi


def gauss(x, amp, mu, sig):
       return amp*exp(-(x-mu)**2/(2.*sig**2))

class ImgHolder:

    def __init__(self, directory, log):

        self.log = log
        self.fig, self.imgAx = plt.subplots(figsize=(12,8))
        self.fig.number = "main"
        self.directory = directory
        
        logger.info('Set working directory as {}'.format(directory))
        self.file = glob("{}/*.fit*".format(self.directory))[0]
        self.header = fits.open(self.file)[0].header
        self.naxis = int(self.header['NAXIS'])
        if self.naxis == 2:
            self.data = fits.open(self.file)[0].data
        elif self.naxis == 3:
            self.data = fits.open(self.file)[0].data[0]
        else:
            #some unknown option, fix this -- is sloopy
            print('error')
        self.fwhm = []

        self.dx = 100

        self.f = 6500 #mm for DFM Scope
        self.pixScale = float(self.header['XPIXSZ'])
        self.instrument = self.header['INSTRUME']

        self.ASperPix = self.pixScale/(self.f*1000.0) * rad2arcsec #Using small angle approx and conv to arcsec
        
        img = self.imgAx.imshow(log(self.data+1), cmap='gray', origin='lower')
        self.imgAx.set_title(self.instrument)
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        cid = self.fig.canvas.mpl_connect('close_event', self.onclose)
        
        axFocLTB = self.fig.add_axes([0.1,0.05,0.1,0.025])
        focLTB   = TextBox(axFocLTB, "focal length [mm]", initial=self.f)
        focLTB.on_submit(self.updateFL)

        axPSclTB = self.fig.add_axes([0.1,0.095,0.1,0.025])
        pSclTB   = TextBox(axPSclTB, "Pixel Scale [um]", initial=self.pixScale)
        pSclTB.on_submit(self.updatePixScale)

        axSFTB = self.fig.add_axes([0.3,0.05,0.1,0.025])
        sFTB   = TextBox(axSFTB, "sub frame dim [p]", initial=self.dx)
        sFTB.on_submit(self.updateSF)
        
        axButt = self.fig.add_axes([0.5,0.05,0.1,0.065])
        reCentButt = Button(axButt, "Re-Center")
        reCentButt.on_clicked(self.recenter)
        
        axFind = self.fig.add_axes([0.6,0.05,0.1,0.065])
        findButt = Button(axFind, "Find FWHM")
        findButt.on_clicked(self.findfwhm)
        
        axClear = self.fig.add_axes([0.8,0.05,0.1,0.065])
        clearButt = Button(axClear, "Clear")
        clearButt.on_clicked(self.clear)
        
        plt.show()

    def onclose(self, event):
        
        #grab all patches
        for i,p in enumerate(self.imgAx.patches):
            #extract subframe
            x = int(p.get_x() + self.dx)
            y = int(p.get_y() + self.dx)
            
            subFrame = self.data[y-self.dx:y+self.dx,x-self.dx:x+self.dx]
            
            #get vert data
            vData = subFrame[:,self.dx]
            coeffv, var  = opt.curve_fit(gauss,range(self.dx*2),vData, p0=[5,self.dx,1])

            #get hori data
            hData = subFrame[self.dx,:]
            coeffh, var  = opt.curve_fit(gauss,range(self.dx*2),hData, p0=[5,self.dx,1])
            
            fwhmPix = (abs(coeffh[2])+abs(coeffv[2]))*0.5 * 2.355
            fwhmArcSec = fwhmPix * self.ASperPix * rad2arcsec

            fig1, axArr= plt.subplots(figsize=(8,6),nrows=1, ncols=2, sharey=True)
            axArr[0].step(range(x-self.dx,x+self.dx), hData, color='b', alpha=0.5)
            axArr[1].step(range(y-self.dx,y+self.dx), vData, color='r', alpha=0.5)
            axArr[0].plot(range(x-self.dx,x+self.dx), gauss(range(self.dx*2),*coeffh), 'b.')
            axArr[1].plot(range(y-self.dx,y+self.dx), gauss(range(self.dx*2),*coeffv), 'r.')
            axArr[0].set_title("Source: Horizontal Data & Gaussian Fit")
            axArr[1].set_title("Source: Vertical Data & Gaussian Fit")
            axArr[0].set_xlabel("Pixel Location (Horizontal) [p]")
            axArr[1].set_xlabel("Pixel Location (Vertical) [p]")
            axArr[0].set_ylabel("Profile")
            axArr[0].grid(1)
            axArr[1].grid(1)
            plt.savefig("{}/subFrameGaussFit_{:03d}.png".format(self.directory, i))
        
        fig, ax = plt.subplots(figsize=(16,10))
        ax.imshow(log(self.data+1), cmap='gray', origin='lower')
        ax.set_title("Camera: {}; Pixel Scale: {:0.2f}[um/p] Scope focalLen: {}[mm]".format(\
            self.instrument,self.pixScale, self.f))
        for i,p in enumerate(self.imgAx.patches):
            ax.add_patch(Rectangle((p.get_x(), p.get_y()), p.get_width(), p.get_height(), edgecolor='r', fill=False))
            ax.text(p.get_x()+2*self.dx+5, p.get_y()+2*self.dx+5, "{:0.3f} \'\'".format(self.fwhm[i]), color='r')
        plt.savefig("{}/main-figure.png".format(self.directory))
        plt.close('all')

    def updateSF(self, text):
        self.dx = int(text)
        

    def updatePixScale(self, text):
        self.pixScale = float(text)
        self.ASperPix = self.pixScale/(self.f*1000.0) * rad2arcsec
        print(self.ASperPix)

    def updateFL(self, text):
        self.f = float(text)
        self.ASperPix = self.pixScale/(self.f*1000.0) * rad2arcsec
        print(self.ASperPix)

    def clear(self, event):
        self.imgAx.clear()
        self.imgAx.imshow(log(self.data+1), cmap='gray', origin='lower')
        self.imgAx.set_title(self.instrument)

        plt.gcf().canvas.draw_idle()

    def onclick(self, event):
        x, y = int(event.xdata), int(event.ydata)
        ax = event.inaxes
        ax.add_patch(Rectangle((x-self.dx,y-self.dx), self.dx*2, self.dx*2, edgecolor='r', fill=False))
        self.fwhm.append(-999)
        plt.gcf().canvas.draw_idle()

    def recenter(self, event):
        for p in self.imgAx.patches:
            #Extract subframe
            x = int(p.get_x() + p.get_width()/2)
            y = int(p.get_y() + p.get_height()/2)
            
            subFrame = self.data[y-self.dx:y+self.dx,x-self.dx:x+self.dx]
            
            #Calculate Centroid
            center = centroids.centroid_com(subFrame)
            
            #Update rectangle
            p.set(x=x+center[0]-2*self.dx, y=y+center[1]-2*self.dx)

            #Redraw rectangle        
            plt.gcf().canvas.draw_idle()

    def findfwhm(self, event):
        for i, p in enumerate(self.imgAx.patches):
            #extract subframe
            x = int(p.get_x() + p.get_width()/2)
            y = int(p.get_y() + p.get_height()/2)

            subFrame = self.data[y-self.dx:y+self.dx,x-self.dx:x+self.dx]
            
            #get vert data
            vData = subFrame[:,self.dx]
            coeffv, var  = opt.curve_fit(gauss,range(self.dx*2),vData, p0=[5,self.dx,1])

            #get hori data
            hData = subFrame[self.dx,:]
            coeffh, var  = opt.curve_fit(gauss,range(self.dx*2),hData, p0=[5,self.dx,1])
            
            fwhmPix = (abs(coeffh[2])+abs(coeffv[2]))*0.5 * 2.355
            fwhmArcSec = fwhmPix * self.ASperPix

            self.fwhm[i] = fwhmArcSec
            
            self.imgAx.text(x+self.dx+5, y+self.dx+5, "{:0.3f} \'\'".format(fwhmArcSec), color='r')
            
        #Redraw rectangle        
        plt.gcf().canvas.draw_idle()
 
from logging.config import fileConfig

logging.basicConfig(filename='seeingCalc.log', level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()
logger.info('Logger instance created')
#DEBUG INFO WARNING ERROR CRITICAL

directory = argv[1]
logger.info('Received working directory CLI input as {}'.format(directory))

iH = ImgHolder(directory, log)
logger.info('Attempting to start seeingCalc session; ImgHolder Obj')

