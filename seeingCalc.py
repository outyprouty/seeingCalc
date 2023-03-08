from matplotlib import pyplot as plt
from glob import glob
from astropy.io import fits
from photutils import centroids
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button,TextBox
from numpy import log, exp, pi
from scipy import optimize as opt
import logging

from sys import argv

verb = '-v' in argv
dx, dy = 50,50

rad2arcsec = 3600*180/pi

def gauss(x, amp, mu, sig):
       return amp*exp(-(x-mu)**2/(2.*sig**2))

class ImgHolder:

    def __init__(self, directory):
        self.fig, self.imgAx = plt.subplots(figsize=(12,8))
        self.fig.number = "main"
        self.directory = directory
        self.file = glob("{}/*.fit*".format(self.directory))[0]
        self.data = fits.open(self.file)[0].data[0]
        self.fwhm = []

        self.f = 6500 #mm for DFM Scope
        self.pixScale = float(fits.open(self.file)[0].header['XPIXSZ'])
        self.instrument = fits.open(self.file)[0].header['INSTRUME']

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
        
        axButt = self.fig.add_axes([0.2,0.05,0.1,0.065])
        reCentButt = Button(axButt, "Re-Center")
        reCentButt.on_clicked(self.recenter)
        
        axFind = self.fig.add_axes([0.3,0.05,0.1,0.065])
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
            x = int(p.get_x() + dx)
            y = int(p.get_y() + dy)
            
            subFrame = self.data[y-dy:y+dy,x-dx:x+dx]
            
            #get vert data
            vData = subFrame[:,dx]
            coeffv, var  = opt.curve_fit(gauss,range(dy*2),vData, p0=[5,dy,1])

            #get hori data
            hData = subFrame[dy,:]
            coeffh, var  = opt.curve_fit(gauss,range(dx*2),hData, p0=[5,dx,1])
            
            fwhmPix = (abs(coeffh[2])+abs(coeffv[2]))*0.5 * 2.355
            fwhmArcSec = fwhmPix * self.ASperPix * rad2arcsec

            fig1, ax1= plt.subplots()
            ax1.step(range(x-dx,x+dx), hData, color='b', label='hor; {:0.2f}p'.format(coeffh[2]))
            ax1.step(range(y-dy,y+dy), vData, color='r', label='ver; {:0.2f}p'.format(coeffv[2]))
            ax1.plot(range(x-dx,x+dx), gauss(range(dx*2),*coeffh), 'b.')
            ax1.plot(range(y-dy,y+dy), gauss(range(dy*2),*coeffv), 'r.')
            ax1.legend(loc=0)
            ax1.set_title("SubFrame Gaussian Fit")
            ax1.set_xlabel("Pixel Location")
            ax1.set_ylabel("Profile")
            ax1.grid(1)
            plt.savefig("{}/subFrameGaussFit_{:03d}.png".format(self.directory, i))
        
        fig, ax = plt.subplots(figsize=(16,10))
        ax.imshow(log(self.data+1), cmap='gray', origin='lower')
        ax.set_title("Camera: {}; Pixel Scale: {:0.2f}[um/p] Scope focalLen: {}[mm]".format(\
            self.instrument,self.pixScale, self.f))
        for i,p in enumerate(self.imgAx.patches):
            ax.add_patch(Rectangle((p.get_x(), p.get_y()), p.get_width(), p.get_height(), edgecolor='r', fill=False))
            ax.text(p.get_x()+2*dx+5, p.get_y()+2*dy+5, "{:0.3f} \'\'".format(self.fwhm[i]), color='r')
        plt.savefig("{}/main-figure.png".format(self.directory))
        plt.close('all')
            

    def updatePixScale(self,text):
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
        imgAx.set_title(self.instrument)

        plt.gcf().canvas.draw_idle()

    def onclick(self, event):
        x, y = int(event.xdata), int(event.ydata)
        if verb: print("{:d}, {:d}".format(x, y))
        ax = event.inaxes
        ax.add_patch(Rectangle((x-dx,y-dy), dx*2, dy*2, edgecolor='r', fill=False))
        self.fwhm.append(-999)
        plt.gcf().canvas.draw_idle()

    def recenter(self, event):
        for p in self.imgAx.patches:
            #Extract subframe
            x = int(p.get_x() + dx)
            y = int(p.get_y() + dy)
            if verb: print("Current center: {},{}".format(x, y))
            
            subFrame = self.data[y-dy:y+dy,x-dx:x+dx]
            
            #Calculate Centroid
            center = centroids.centroid_com(subFrame)
            if verb: print(center[0], center[1])
            
            #Update rectangle
            p.set(x=x+center[0]-2*dx, y=y+center[1]-2*dy)

            #Redraw rectangle        
            plt.gcf().canvas.draw_idle()

    def findfwhm(self, event):
        for i, p in enumerate(self.imgAx.patches):
            #extract subframe
            x = int(p.get_x() + dx)
            y = int(p.get_y() + dy)

            subFrame = self.data[y-dy:y+dy,x-dx:x+dx]
            
            #get vert data
            vData = subFrame[:,dx]
            coeffv, var  = opt.curve_fit(gauss,range(dy*2),vData, p0=[5,dy,1])

            #get hori data
            hData = subFrame[dy,:]
            coeffh, var  = opt.curve_fit(gauss,range(dx*2),hData, p0=[5,dx,1])
            
            fwhmPix = (abs(coeffh[2])+abs(coeffv[2]))*0.5 * 2.355
            fwhmArcSec = fwhmPix * self.ASperPix

            self.fwhm[i] = fwhmArcSec
            
            self.imgAx.text(x+dx+5, y+dy+5, "{:0.3f} \'\'".format(fwhmArcSec), color='r')
            
        #Redraw rectangle        
        plt.gcf().canvas.draw_idle()
 


directory = argv[1]

iH = ImgHolder(directory)


