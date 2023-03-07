from matplotlib import pyplot as plt
from astropy.io import fits
from photutils import centroids
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button
from numpy import log, exp
from scipy import optimize as opt

from sys import argv

verb = '-v' in argv
dx, dy = 50,50

def gauss(x, amp, mu, sig):
       return amp*exp(-(x-mu)**2/(2.*sig**2))

class ImgHolder:

    def __init__(self, imgAx, data):
        self.imgAx = imgAx
        self.data = data
        
        img = self.imgAx.imshow(log(data+1), cmap='gray', origin='lower')
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        
        axButt = fig.add_axes([0.7,0.05,0.1,0.075])
        reCentButt = Button(axButt, "Re-Center")
        reCentButt.on_clicked(self.recenter)
        
        axFind = fig.add_axes([0.2,0.05,0.1,0.075])
        findButt = Button(axFind, "Find FWHM")
        findButt.on_clicked(self.findfwhm)
        
        
        plt.show()

    def onclick(self, event):
        x, y = int(event.xdata), int(event.ydata)
        if verb: print("{:d}, {:d}".format(x, y))
        ax = event.inaxes
        ax.add_patch(Rectangle((x-dx,y-dy), dx*2, dy*2, edgecolor='r', fill=False))

        plt.gcf().canvas.draw_idle()

    def recenter(self, event):
        for p in self.imgAx.patches:
            #Extract subframe
            x = int(p.get_x() + dx)
            y = int(p.get_y() + dy)
            if verb: print("Current center: {},{}".format(x, y))
            
            subFrame = data[y-dy:y+dy,x-dx:x+dx]
            
            #Calculate Centroid
            center = centroids.centroid_com(subFrame)
            if verb: print(center[0], center[1])
            
            #Update rectangle
            p.set(x=x+center[0]-2*dx, y=y+center[1]-2*dy)

            #Redraw rectangle        
            plt.gcf().canvas.draw_idle()

    def findfwhm(self, event):
        for p in self.imgAx.patches:
            #extract subframe
            x = int(p.get_x() + dx)
            y = int(p.get_y() + dy)

            subFrame = data[y-dy:y+dy,x-dx:x+dx]
            
            #get vert data
            vData = subFrame[:,dx]
            coeffv, var  = opt.curve_fit(gauss,range(dy*2),vData, p0=[5,dy,1])

            #get hori data
            hData = subFrame[dy,:]
            coeffh, var  = opt.curve_fit(gauss,range(dx*2),hData, p0=[5,dx,1])
            #fig1 = plt.figure()
            #plt.step(range(dx*2), hData, color='b', label='hor')
            #plt.step(range(dy*2), vData, color='r', label='ver')
            #plt.plot(range(dx*2), gauss(range(dx*2),*coeffh), 'b.')
            #plt.plot(range(dy*2), gauss(range(dy*2),*coeffv), 'r.')
            #plt.legend(loc=0)
            #plt.grid(1)
            #plt.savefig("subFrameGaussFit.png")
        
            self.imgAx.text(x+dx+5, y+dy+5, "{:0.3f} pix".format((abs(coeffh[2])+abs(coeffv[2])*0.5*2.355)), color='r')
            
        #Redraw rectangle        
        plt.gcf().canvas.draw_idle()
 

fig, ax = plt.subplots(figsize=(16,10))
data = fits.open("M42-trap-color_00001.fits")[0].data[0]

iH = ImgHolder(ax, data)


#fig.canvas.mpl_disconnect(cid)
