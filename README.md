# Coadd Cutouts
 This repository contains programs to create fits cutouts from coadded DES tiles
and  DES SE images.

## CoaddCutouts

The program MakeCoaddCutoutsY6.py is designed to be used
by a user that have have DES stack installed and able to 
use easyaccess program. The program is using some modules
from easyaccess to communicate with the DES database.

The program optimizes  IO operations. For this the list of objects is analysed and
list of objects that belong to a tile is created. So input files for given tile and 
 neighboring tiles if necessary will be copied to data/ sub-directory and be used 
 for making cut outs for all the objects.
 For each object cutout stamps are created in workdir/ sub-directory for each band.
These files will be then moved to a sub-directory named with the object name.
Required PSF files will be created in the last one.
This intermediate directory is used when we need to coadd stamps that belongs to neighboring
tiles. The coadded stamps will be moved to resulting sub-directory.

The size of the cut-out stamp is given in arcsec and will be corrected to provide exact number
of pixels.  

There are two options for presenting results:
 res_flag = 0 Created cut-out images, weight masks and PSF images for each band will be stored in 
 a sub-directory for each object.
* res_flag = 1 a multi-extension  fits file will be created for each tile.
 The file will contain 3 image extensions for each object (IMAGE,WEIGHT,PSF)
 Each image in the file contains following header cards:
* OBJECT - name of the object
* BAND - name of the band (g,r,i,z,Y)
* TYPE - (IMAGE,WEIGHT,PSF)
 
 The DESCoaddCutL.sh gives a simple example of the program usage.

## FitsCutouts

The program MakeFitsCutoutsY6.py creates cutouts from DES Y6 SE images. The
program selects best SE images and produced fits files similar to one produced
by MakeCoaddCutoutsY6.py

The DESFitsCutL.sh script gives a simple example of the program usage.


## PSF image

The CreatePSFFile.py is used to create PSF image at given coordinate. It uses
Erin Sheldon's program to create PSF from psfex file. The program is used for
both coadd and single epoch cutouts.

## MakeCutOutStamp.py
 
 It is a simple program to create subimage from given extension in fits file.
 
