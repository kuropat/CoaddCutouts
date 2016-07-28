#!/usr/bin/env python
'''
Created on February 23, 2016
The program to create psf image file from psfexcat file in given RA,Dec 
Input: RA,Dec in degrees.
 The program calculates pixel coordinates of the object and construct
 fits image for given point 
 

@author: kuropat

'''
import os
import sys
import fitsio
import exceptions
import string
import getopt
import numpy as np
#from despyastro import wcsutil

class CreatePSFFile(exceptions.Exception):
    
    def __init__(self, infile,outfile):
        self.degPerPix = 0.27/3600.
        self.infile = infile
        self.outfile = outfile
        self.resh = {}
        self.wcsInfo = []
        self.objectName = ''
        self.band = ''
        self.posx = 0.
        self.posy = 0.


    def setPSFAtt(self,objectName,band):
        self.objectName = objectName
        self.band = band
        
        '''
        Create PSF image
        '''
        
    def produce(self,posX,posY):
        outpath = os.path.normpath(self.outfile)
        if os.path.exists(outpath):
            os.remove(outpath)
        self.posX = posX
        self.posY = posY
#
        if not os.path.exists(self.infile):
            self.posX = 0.
            self.posY = 0.
            return
        fitsfile = os.path.normpath(self.infile)
        fits=fitsio.FITS(fitsfile)
        prihdr = fits[1].read_header()

        x_psf = self.posX
        y_psf = self.posY
#
        x0 = prihdr['POLZERO1']
        xscale = prihdr['POLSCAL1']
        y0 = prihdr['POLZERO2']
        yscale = prihdr['POLSCAL2']
        psf_fwhm = prihdr['PSF_FWHM']
        psf_samp = prihdr['PSF_SAMP']
        psf_axis = prihdr['PSFNAXIS']
        psf_ny = prihdr['PSFAXIS1']
        psf_nx = prihdr['PSFAXIS2']
        psf_nc = prihdr['PSFAXIS3']
#
        pim = fits[1]['PSF_MASK'][0]  # this is array 6,15,15
        psf_im = np.zeros((psf_ny,psf_nx))
        x1 = (x_psf - x0)/xscale;
        y1 = (y_psf - y0)/yscale;

        for  i in range(0,psf_ny):
            for j in range(0,psf_nx):
                psf_im[i,j] = pim[0,i,j] + pim[1,i,j]*x1 + pim[2,i,j]*x1*x1 + pim[3,i,j]*y1 + pim[4,i,j]*y1*y1 + pim[5,i,j]*x1*y1

        self.writePSF(psf_im)

            
    def writePSF(self,psf_im):
        fits1 = fitsio.FITS(self.outfile,'rw')
        fits1.create_image_hdu(psf_im)
        header = fits1[0].read_header()
#        header['OBJECT']=self.objectName
#        header['BAND']=self.band
        fits1.write(psf_im,header=header)
        fits1[0].write_keys({'OBJECT':self.objectName,'BAND':self.band})
        fits1.close()

           
if __name__ == "__main__":
    print sys.argv
    outFile=""
    infile=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:x:y:",
                                  ["infile=","outfile=","X=","Y="])
    except getopt.GetoptError:
        print "CreatPSFFile.py -i <infile> -o <outfile> -x <x> -y <y>" 
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "CreatePSFFile.py -i <infile> -o <outfile>  -x <x> -y <y>" 
            sys.exit()
        elif opt in ("-i","--infile"):
            print "got -i arg= %s"%arg
            infile=arg
        elif opt in ("-o","--outile"):
            print "got -o arg= %s"%arg
            outfile=arg
        elif opt in ("-x","--X"):
            print "got -x arg= %s"%arg
            x=string.atof(arg) 
        elif opt in ("-y","--Y"):
            print "got -y arg= %s"%arg
            y=string.atof(arg)  
    mkpsf = CreatePSFFile(infile,outfile)
    mkpsf.produce(x,y)

    sys.exit(0)
