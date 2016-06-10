#!/usr/bin/env python
'''
Created on February 10, 2016
The program to create cut out stamp around object coordinates provided in
input as RA,Dec in degrees and stamp size in degrees
 The program calculates pixel coordinates of the object and cuts out
 a stamp in such way that the object will be in center of the stamp. 

@author: kuropat
'''
import os
import sys
import fitsio
import exceptions
import string
import getopt
#import pyfits
import numpy as np
from despyastro import wcsutil
#from pyfits._release import PYFITS_HOMEPAGE_BASE_URL

class MakeCutOutStamp(exceptions.Exception):
    
    def __init__(self,ext, infile,outfile,stampS):
        self.degPerPix = 0.2636/3600.
        self.stampsize = string.atof(stampS) # stamp size in degrees
        self.stamppix = (int)(self.stampsize/self.degPerPix)
        if self.stamppix % 2 != 0: self.stamppix += 1
        self.infile = infile
        self.outfile = outfile
        self.ext = ext # extention to cutout 0 image 1 weight
        self.resh = {}
        self.wcsInfo = []
        self.objectName = ''
        self.band = ''
        self.im_type = ''
       
    """ Set info about the stamp object name, band, type - image or weight """    
    def setStampAtt(self,objectName,band,Im_type):
        self.objectName = objectName
        self.bane = band
        self.im_type = Im_type

    """ Create stamp image """    
    def makeStamp(self,image,xobj,yobj):
        s = (self.stamppix,self.stamppix)
        stamp = np.zeros(s, dtype='f4')
        halfstamp = (int)(self.stamppix/2);
        " galaxy should be in the center of stamp "
        galX = self.stamppix/2;
        galY = self.stamppix/2;
        deltaX = 0;
        deltaY = 0;
        self.tileX = (int)(xobj - halfstamp)
        self.tileY = (int)(yobj - halfstamp)

        ''' if stamp goes outside of image shift it back 
          in this case shift is positive '''
        if (self.tileX + self.stamppix) > self.ncol:
            deltaX = self.tileX+self.stamppix - self.ncol
            self.tileX = self.ncol - self.stamppix
        elif  self.tileX < 0:
            deltaX = self.tileX;
            self.tileX = 0;            
        if (self.tileY + self.stamppix) > self.nrow:
            deltaY = self.tileY+self.stamppix - self.nrow;
            self.tileY = self.nrow - self.stamppix;
        elif self.tileY < 0:
            deltaY = self.tileY;
            self.tileY = 0;
        galX -= deltaX;
        galY -= deltaY;

        sy = (int)(self.tileY)
        ey = (int)(self.tileY + self.stamppix)
        sx = (int)(self.tileX)
        ex = (int)(self.tileX + self.stamppix)
        subStamp = np.copy(image[sy:ey:1,sx:ex:1])
        xs = 0; ys = 0; x1s = 0; y1s = 0
        xe = (int)(self.stamppix); x1e = (int)(self.stamppix); ye = (int)(self.stamppix); y1e = (int)(self.stamppix)
        if deltaX == 0:
            xs = 0
            xe = (int)( self.stamppix   )
            x1s = 0
            x1e = (int)(self.stamppix  ) 
        elif deltaX < 0:
            xs = 0
            x1s =(int)(-deltaX)
            x1e = (int)(self.stamppix  )
            xe = x1e - x1s            
        else:
            xs = (int)(deltaX)
            xe = (int)(self.stamppix )
            x1s = 0
            x1e = (int)(self.stamppix - deltaX  ) 
 
        if deltaY == 0:
            ys = 0
            ye = (int)(self.stamppix   )
            y1s = 0
            y1e = (int)(self.stamppix )
        elif deltaY < 0:
            ys = 0
            y1s = (int)(-deltaY)
            y1e = (int)(self.stamppix )
            ye = y1e - y1s
        else:
            ys = (int)(deltaY)
            ye = (int)(self.stamppix )
            y1s =  0
            y1e = (int)(self.stamppix -deltaY  )
        stamp[y1s:y1e:1,x1s:x1e:1] = subStamp[ys:ye:1,xs:xe:1]

        return stamp

    
    """ create header of the new stamp """
    def copyHeader(self, master):
        keys = master.keys()
        res = self.resh
        for key in keys:
            if key == "NAXIS1":
                res.update({key:self.stamppix})
                continue
            elif ( key == "NAXIS2"):
                res.update({key:self.stamppix})
                continue
            elif (key == "CHECKSUM"):
                continue
            elif (key == "DATASUM"):
                continue    
            elif (key == "PCOUNT"):
                continue
            elif (key == "GCOUNT"):
                continue        
            elif (key == "END"):
                continue
            elif ( key == "XTENSION"):
                res.update({"SIMPLE":True})
                continue
            elif (key == "CRPIX1"):
                crpix1 =  master.get(key) - (self.posX - self.stamppix / 2.)
                res.update({"CRPIX1":crpix1})
            elif ( key == "CRPIX2"):
                crpix2 = master.get(key) - (self.posY - self.stamppix / 2.)
#                print 'crpix2=%d \n' % crpix2
                res.update({"CRPIX2":crpix2})
            else:
                res.update({key:master.get(key)})
        """ Now add stamp description comments """
        res.update({'OBJECT':self.objectName})
        res.update({'BAND':self.band})
        res.update({'TYPE':self.im_type})
        
        return res

        
        '''
        Create RGB image
        '''
        
    def produce(self,RAS,DecS):
        outpath = os.path.normpath(self.outfile)
        outdir = os.path.dirname(self.outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if os.path.exists(outpath):
            os.remove(outpath)
        self.RA = string.atof(RAS)
        self.Dec = string.atof(DecS)
        fitsfile = os.path.normpath(self.infile)
        fits=fitsio.FITS(fitsfile)
        if fitsfile.find('.fz') > 0 :
            prihdr = fits[1].read_header()
            self.ncol = prihdr['ZNAXIS1']
            self.nrow = prihdr['ZNAXIS2']
            data1 = fits[self.ext+1].read()
        else:
            prihdr = fits[0].read_header()
            self.ncol = prihdr['NAXIS1']
            self.nrow = prihdr['NAXIS2']
            data1 = fits[self.ext].read()

        w = wcsutil.WCS(prihdr)
        objx,objy = w.sky2image(self.RA,self.Dec)
        self.posX = objx
        self.posY = objy

        data = self.makeStamp(data1,objx,objy)

        self.writeStamp(data,prihdr)

    
    def writeStamp(self,data,prihdr):
        fits1 = fitsio.FITS(self.outfile,'rw')
        
        hdict = self.copyHeader(prihdr)
 
        fits1.write(data, header=hdict)
 
 
    def clean(self):
        print 'Remove orig. file %s \n' % self.infile
        try:
            os.remove(self.infile)
        except:
            print "remove failed on file %s \n" % infile
            
    def getPos(self):
        return (self.posX,self.posY)
           
if __name__ == "__main__":
    print sys.argv
    outFile=""
    infile=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:e:r:d:",
                                  ["infile=","outfile=","stampsize=","extension","RA=","Dec="])
    except getopt.GetoptError:
        print "MakeCutOutStamp.py -i <infile> -o <outfile> -s <stampsize> -e <extension> -r <RA> -d <Dec>" 
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "MakeCutOutStamp.py -i <infile> -o <outfile> -s <stampsize> -e <extension> -r <RA> -d <Dec>" 
            sys.exit()
        elif opt in ("-i","--infile"):
            print "got -i arg= %s"%arg
            infile=arg
        elif opt in ("-o","--outile"):
            print "got -o arg= %s"%arg
            outfile=arg            
        elif opt in ("-s","--stampsize"):
            print "got -s arg= %s"%arg
            stampS=arg  
        elif opt in ("-e","--extension"):
            print "got -e arg= %s"%arg
            ext=arg  
        elif opt in ("-r","--RA"):
            print "got -r arg= %s"%arg
            RA=arg 
        elif opt in ("-d","--Dec"):
            print "got -d arg= %s"%arg
            Dec=arg  
    mrgbs = MakeCutOutStamp(ext,infile,outfile,stampS)
    mrgbs.produce(RA,Dec)
    mrgbs.clean()

    sys.exit(0)
