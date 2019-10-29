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
import psfex
import getopt
import numpy as np
#

class CreatePSFFile(exceptions.Exception):
    
    def __init__(self, infile,outfile):
        self.degPerPix = 0.263/3600.
        self.infile = infile
        self.outfile = outfile
        self.resh = {}
        self.wcsInfo = []
        self.objectName = ''
        self.band = ''



    def setPSFAtt(self,objectName,band):
        self.objectName = objectName
        self.band = band
        
        '''
        Create PSF image
        '''
        
    def produce(self,posX,posY):
        print "CreatePSFFile input %s \n" % self.infile
        if not os.path.isfile(self.infile):
            print " CreatePSFFile input file %s do not exists \n" % self.infile
            return

#
        outpath = os.path.normpath(self.outfile)
        if os.path.exists(outpath):
            os.remove(outpath)
        self.posX = posX
        self.posY = posY
#        fitsfile = os.path.normpath(self.posX,self.posY)

        pex=psfex.PSFEx(self.infile) 
        image=pex.get_rec(self.posX,self.posY)
#
        print "writing PSF image to file %s \n" % self.outfile
        self.writePSF(image)
            
    def writePSF(self,psf_im):
        fits1 = fitsio.FITS(self.outfile,'rw')
        fits1.create_image_hdu(psf_im)
        header = fits1[0].read_header()

        fits1.write(psf_im,header=header)
        fits1[0].write_keys({'OBJECT':self.objectName,'BAND':self.band,'TYPE':'PSF'})
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
