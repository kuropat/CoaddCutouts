#!/usr/bin/env python

"""
 The program to create a set of fits cutouts for provided RA,Dec coordinates
 Usage:
    MakeCoaddCutoutsY6.py   -l <objlist> -s <stampsize> -r <res_flag>
    Here:
          <objlist> list of objects in the form objName,RA,Dec (RA,Dec in degrees)
          <stampsize> is size of the cutout stamp in arcsec
          <res_flag> if 0 each object cutouts are stored in subdirectory with the object name,
                     if 1 results are stored in multi-extension fits file named by tile 
    By N. P. Kuropatkin  May 11 2016
    Updated By N. Kuropatkin April 30 2018 to use Erin Sheldons CreatePSFFile 
    

"""

import os
import shutil
import sys
import string
import getopt
import subprocess


import fitsio
from time import sleep
from MakeCutOutStamp import MakeCutOutStamp
from CreatePSFFile import CreatePSFFile
import urllib2, ssl
import config_ea as config_mod
import cx_Oracle
#from scipy.stats.mstats_extras import hdquantiles
try:
    from termcolor import colored
except:
    def colored(line, color):
        return line





class MakeCoaddCutoutsY6():
    
    def __init__(self, objlist, ssize, resType ):
        '''
        Constructor
        '''
        urlbase = "https://desar2.cosmology.illinois.edu/DESFiles/"
        self.objlist = objlist
        self.dbname = 'dessci'
#
        self.resType = string.atoi(resType)
        self.pixscale=0.2636
        self.spix = int(string.atof(ssize)/self.pixscale)  # stamp size in pixels
        if self.spix % 2 != 0:
            self.spix +=1
        self.ssize = self.spix*self.pixscale
        self.ssize /= 3600. # stamp size in degrees
        print "Pix size=%d arcsec size = %f \n" % (self.spix, self.ssize*3600.)
        self.autocommit = True
        self.quiet = False
        self.bands = ['g','r','i','z','Y']
        self.curdir = os.getcwd()
        desfile = os.getenv("DES_SERVICES")
        if not desfile: desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
        self.desconfig = config_mod.get_desconfig(desfile, self.dbname)
        self.connectDB()
        self.opener = self.get_connected(urlbase,self.user,self.password)
        self.tiles = {}
        self.objects = {}
        self.tilesForObj = {}
        self.objForTile = {}


    """ establish connection to the database """
    def connectDB(self):        
        self.user = self.desconfig.get('db-' + self.dbname, 'user')
        self.dbhost = self.desconfig.get('db-' + self.dbname, 'server')
        self.port = self.desconfig.get('db-' + self.dbname, 'port')
        self.password = self.desconfig.get('db-' + self.dbname, 'passwd')
        kwargs = {'host': self.dbhost, 'port': self.port, 'service_name': self.dbname}
        self.dsn = cx_Oracle.makedsn(**kwargs)
        if not self.quiet: print('Connecting to DB ** %s ** ...' % self.dbname)
        connected = False
        
        for tries in range(3):
            try:
                self.con = cx_Oracle.connect(self.user, self.password, dsn=self.dsn)
                if self.autocommit: self.con.autocommit = True
                connected = True
                break
            except Exception as e:
                lasterr = str(e).strip()
                print(colored("Error when trying to connect to database: %s" % lasterr, "red"))
                print("\n   Retrying...\n")
                sleep( 8 )
        if not connected:
            print('\n ** Could not successfully connect to DB. Try again later. Aborting. ** \n')
            os._exit(0)
        self.cur = self.con.cursor()




    """ make file path for file copy """
    def makeFilePath(self,pathname):
        pathT = """https://desar2.cosmology.illinois.edu/DESFiles/desarchive/"""
        fpath = pathT+pathname
        return fpath


    
    def get_objects(self):
        f = open(self.objlist, 'rb')
        for line in f:
            tokens = string.split(line,',')
            self.objects.update({tokens[0]:[string.atof(tokens[1]),string.atof(tokens[2])]})
        f.close()

    """ Creates dictionary of object name - list of involved tiles   """
    def get_tilegeom(self):
        stampHS = self.ssize/2.
#
        obj_rec = []
        for okey in self.objects:
            ObjName = okey
            obj_rec = self.objects[okey]
            SRA = obj_rec[0]
            SDEC = obj_rec[1]
            SRAmin = SRA - stampHS
            SRAmax = SRA + stampHS
            SDecmin = SDEC - stampHS
            SDecmax = SDEC + stampHS
            

            tileinfo = self.get_tile(SRA,SDEC)
            tile_list = tileinfo.keys()
            print tile_list[0]
            RACMIN = tileinfo[tile_list[0]][0]["RACMIN"]
            DECCMIN = tileinfo[tile_list[0]][0]["DECCMIN"]
            RACMAX = tileinfo[tile_list[0]][0]["RACMAX"]
            DECCMAX = tileinfo[tile_list[0]][0]["DECCMAX"]

            for tilename in tile_list:
#
                if  tilename not in self.tiles:
                    self.tiles.update({tilename:tileinfo[tilename]})
                else:
                    self.tiles[tilename]= self.tiles[tilename] + (tileinfo[tilename]) 
 
            # Check if we need neighbor tiles 
            if SRAmin < RACMIN and SDecmin > DECCMIN : # West
                tileinfo1 = self.get_tile(SRAmin,SDEC) # dict of list of dict.
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename]= self.tiles[tilename] + (tileinfo1[tilename])   

                    
            elif  SRAmin < RACMIN and SDecmin < DECCMIN: # Left bottom corner
                tileinfo1 = self.get_tile(SRAmin,SDEC) # West
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename] =  self.tiles[tilename]+(tileinfo1[tilename])

                tileinfo2 = self.get_tile(SRAmin,SDecmin) # South - West
                tileL = tileinfo2.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo2[tilename]})
                    else:
                        self.tiles[tilename]= self.tiles[tilename] + (tileinfo2[tilename])
                        
                tileinfo3 = self.get_tile(SRA,SDecmin) # South 
                tileL = tileinfo3.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo3[tilename]})
                    else:
                        self.tiles[tilename] =self.tiles[tilename] +(tileinfo3[tilename])

                    
            elif  SRAmin < RACMIN and SDecmax > DECCMAX: # Left upper corner
                tileinfo1 = self.get_tile(SRAmin,SDEC) # west
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo1[tilename])
                    
                tileinfo2 = self.get_tile(SRAmin,SDecmax) # West - North
                tileL = tileinfo2.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo2[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo2[tilename])

                    
                tileinfo3 = self.get_tile(SRA,SDecmax) # North    
                tileL = tileinfo3.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo3[tilename]})
                    else:
                        self.tiles[tilename]= self.tiles[tilename] + (tileinfo3[tilename])      
                tilename = tileinfo3.keys()[0]

                    
            elif SRAmin > RACMIN and SRAmin <RACMAX and  SDecmin < DECCMIN: # bottom
                tileinfo1 = self.get_tile(SRA,SDecmin) # South
                tileL = tileinfo1.keys()
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename]= self.tiles[tilename]+(tileinfo1[tilename])

                
            elif  SRAmax > RACMAX and SDecmax > DECCMAX: # right upper corner
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename]= self.tiles[tilename] + (tileinfo1[tilename])

            
                tileinfo2 = self.get_tile(SRAmax,SDecmax) # North - East
                tileL = tileinfo2.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo2[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo2[tilename])

                    
                tileinfo3 = self.get_tile(SRA,SDecmax) # North
                tileL = tileinfo3.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo3[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo3[tilename])

                
            elif  SRAmax > RACMAX and SDecmin < DECCMIN: # Right lower corner
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo1[tilename])

                
                tileinfo2 = self.get_tile(SRAmax,SDecmin) # South - East
                tileL = tileinfo2.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo2[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo2[tilename])

                tileinfo3 = self.get_tile(SRA,SDecmin) # South
                tileL = tileinfo3.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo3[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo3[tilename])

            elif SRAmax > RACMAX and SDecmax < DECCMAX and  SDecmin > DECCMIN: # Right
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tileL = tileinfo1.keys()
#
                tile_list= list(set(tile_list)|set(tileL))
                for tilename in tileL:
                    if tilename not in self.tiles :
                        self.tiles.update({tilename:tileinfo1[tilename]})
                    else:
                        self.tiles[tilename] = self.tiles[tilename] + (tileinfo1[tilename])
#

                if tilename not in self.tiles:
                    self.tiles.update({tilename:tileinfo1[tilename]})
                if tilename not in tile_list: tile_list.append(tilename)
                    
            self.tilesForObj.update({ObjName:tile_list}) 

#        keys_ob = self.tilesForObj.keys()
#        return self.tiles
    
    """ query DB for tile containing given point """
    def get_tile(self,RA,DEC):
        tile_inf = {}
        file_tile = {}
        tiledat = [] # list of dict.
        tilelist=[]
        filelist=[]
        
#
        query_geom = """ select c.FILENAME,c.TILENAME,c.BAND,c.MAG_ZERO,c.RA_CENT,c.DEC_CENT,c.NAXIS1,c.NAXIS2,c.PFW_ATTEMPT_ID,c.RACMIN,c.RACMAX,c.DECCMIN,c.DECCMAX
            from Y6A1_COADD c where c.RACMIN<{SRA} and c.DECCMIN<{SDEC} and c.RACMAX>{SRA} and c.DECCMAX>{SDEC}
            and c.BAND is not NULL and c.BAND not like 'det' """        
        try:
            query = query_geom.format(SRA = str(RA),SDEC = str(DEC))
            self.cur.execute(query)            
            desc = [d[0] for d in self.cur.description]
            lines = self.cur.fetchall()
            for line in lines:      # number of lines correspond to number of files
                if line != None:
                    tiledat.append(dict(zip(desc,line)))
                else:
                    print "No tile found \n"

        except Exception as e:
            query = query_geom.format(SRA = str(RA+360.),SDEC = str(DEC))
            self.cur.execute(query)
            desc = [d[0] for d in self.cur.description]
            lines = self.cur.fetchall()
            for line in lines:
                if line != None:
                    tiledat.append(dict(zip(desc,line)))
                else:
                    print "no tile found \n"

        for tileI in tiledat:
            tilename = tileI['TILENAME']
            filename = tileI['FILENAME']
            if filename not in filelist:
                filelist.append(filename)
            if tilename not in tilelist:
                tilelist.append(tilename)
                file_tile.update({tilename:[tileI]})
                tile_inf.update({tilename: [tileI]})
            else:
                file_tile[tilename].append(tileI)
                tile_inf[tilename].append(tileI)
#        if len(tilelist) > 0:
#            tile_inf.update({tilelist[0]: tiledat})
#        else:
#            print "tilelist is empty \n"
#        print " file_tile for Ra,Dec %f,%f \n" % (RA,DEC)
#        print tile_inf
#        print " tiles for the coord \n"
#        print file_tile.keys()
#        print " files for coord \n"
#        print filelist
        return tile_inf    
    
    """ Creates a list of objects in a tile """
    def get_objForTile(self):
        obj_dict = self.tilesForObj.copy()
        for okey in self.objects:  # loop on objects
            objID = okey
            obj_list = []
            try:
                tileListForObj =  self.tilesForObj[objID]         
                tileID0 =  tileListForObj[0] # get first tile
                keys = obj_dict.keys()  # list of tiles in object
                if objID in keys:
                    for key in keys:  # loop on objects
                        tileID = (obj_dict[key])[0]
                        if tileID == tileID0:
                            obj_list.append(key)
                self.objForTile.update({tileID0:obj_list})
            except Exception as e:
                lasterr = str(e).strip()
                print(colored("Error when searching for tile info: %s" % lasterr, "red"))

        
    
    """ get data for all objects in the tile """
    def getObjData(self,key):
        objlist = self.objForTile[key]
        tile_list = []
        for obj in objlist:
            tiles = self.tilesForObj[obj]
            for otile in tiles:
                if otile not in tile_list:
                    tile_list.append(otile)
        " Now we have a list of tiles lets get data" 
        for tileN in tile_list:
                self.getTileData(tileN)
        print "All data for tile are collected \n" 
               
    """ get data for given tile: find list of images and list of psfcat files """        
    def getTileData(self,tilename):
        print "Get tile data for tile %s \n" % tilename
        if tilename=='None': return
        if tilename== None: return
        queryT = """select filename,archive_name,path,compression from Y6A1_FILE_ARCHIVE_INFO where filename='{FILENAME}'"""
        tileinfo = self.tiles[tilename]

        for tileI in tileinfo:
            filename = tileI['FILENAME']
            psfcatF = string.split(filename,'.fits')[0]+'_psfcat.psf'

            query = queryT.format(FILENAME = filename)
            dirpath = self.curdir+'/data/'
            filepaths = []
            catpaths = []
            try:
                self.cur.execute(query)
                res = self.cur.fetchone()
                fpath = res[2]
                comp = res[3]
                filepath = self.makeFilePath( fpath);
                filepaths.append(filepath+'/'+filename+comp)
            except:
                print(' Failed on file %s \n' %filename)
        
            query1 = queryT.format(FILENAME = psfcatF)
            try:
                self.cur.execute(query1)
                res = self.cur.fetchone()
               
                cpath = res[2]
                

                catpath = self.makeFilePath(cpath);
                catpaths.append(catpath+'/'+psfcatF)  
            except:
                print(' Failed on file %s \n' %filename)

            for fpath in filepaths:
                fname = string.split(fpath,'/')[-1]
                lname = dirpath+fname
                if self.file_exists(fpath): 
                    self.get_file(self.opener,fpath,lname)
                    
                    command = ['funpack','-D',('%s'%lname)]
#                    print command
                    if lname.find('.fz') > 0 and os.path.isfile(lname) :
                        try:
                            subprocess.check_output(command)
                            if os.path.exists(lname):
                                os.remove(lname)    
                        except subprocess.CalledProcessError as e:
                            print "error %s" % e 
            for cpath in catpaths:
                cname = string.split(cpath,'/')[-1]
                lcname = dirpath+cname
                if self.file_exists(cpath): 
                    self.get_file(self.opener,cpath,lcname)

    """ create a CSV list of files for swarp """
    def makeList(self,im_list):
        first = True
        res = ''
        for imName in im_list:
            if first:
                res +=imName
                first = False
            else:
                res += (','+imName)
        return res

    """ call swarp to create a stamp if more than one tile """
    def SWARPcaller(self,images,weights,imName,weightName,RA,Dec):
 
        command = ['swarp',"%s"%images]
        command +=['-NTHREADS','1','-c','default.swarp']

        command +=["-PIXELSCALE_TYPE","MANUAL","-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER_TYPE","MANUAL","-CENTER","%s,%s"%(RA,Dec)]
        command +=['-IMAGE_SIZE',"%f,%f"%(self.spix,self.spix)]

        command +=["-RESAMPLING_TYPE","LANCZOS3","-RESAMPLE","Y"]
        command +=["-SUBTRACT_BACK","Y","-BLANK_BADPIXELS","N"]
        command +=["-DELETE_TMPFILES","Y"]
        command +=["-FSCALASTRO_TYPE","VARIABLE"]
        command +=["-COMBINE","Y","-COMBINE_TYPE","WEIGHTED"] 
        command +=["-IMAGEOUT_NAME",imName]
        command +=["-WEIGHTOUT_NAME",weightName]
        command +=["-RESAMPLE_SUFFIX","_im.fits"]
        command +=["-WEIGHT_TYPE","MAP_WEIGHT","-WEIGHT_IMAGE","%s"%weights]
        command +=["-HEADER_ONLY","N","-VERBOSE_TYPE","QUIET"]
        command +=["-WRITE_XML","N"]
#        print command
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e


        
    
    def writeResFile(self,resOutFile,objName,band,tile):
        """ Read created stamps and PSF and append them to
         results file with modifying headers if need """
        resOutPath = self.curdir+'/'+resOutFile
        first = False
        if  not os.path.exists(resOutPath): first = True
        respath = self.curdir+'/'+objName+'/'
        imfile = respath+objName+'_'+band+'_image.fits'
        im_type = 'IMAGE'
        if os.path.exists(imfile):
            fits1 = fitsio.FITS(imfile,'rw')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,objName,band,im_type,first)
            first = False
        imfile = respath+objName+'_'+band+'_weight.fits'
        im_type = 'WEIGHT'
        if os.path.exists(imfile):
            fits1 = fitsio.FITS(imfile,'rw')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,objName,band,im_type,first)
        psffile = respath+objName+'_'+band+'_psf.fits'
        im_type = 'PSF'
        if os.path.exists(psffile):
            fits1 = fitsio.FITS(psffile,'rw')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,objName,band,im_type,first)
        
    def writeFits(self,outPath,im_data,header,objName,band,im_type,first):
        if not first: fits2 = fitsio.FITS(outPath,'rw')
        header['BAND'] = band
        header['OBJECT'] = objName
        header['TYPE'] = im_type
        if first:
            fitsio.write(outPath,im_data,header=header)
        else:
            fits2.write(im_data,header=header)

            
    """ create cutouts in the workdir """
    def makeTileCutouts(self,tilename):
#        resOutFile = tilename+'_cutouts.fits'
        tileM = tilename
        objs = self.objForTile[tilename]
        for ob in objs:      # Loop on objects in the tile
            print "Start working on object %s \n" % ob
            objname = ob
            resOutFile = objname+'_cutouts.fits'
            RA = self.objects[ob][0]
            Dec = self.objects[ob][1]            
            imTab = {'g':[],'r':[],'i':[],'z':[],'Y':[]}
            weightTab = {'g':[],'r':[],'i':[],'z':[],'Y':[]}
            tiles = self.tilesForObj[objname]
            
            respath = self.curdir+'/'+objname+'/'
            if not os.path.exists(respath):
                os.mkdir(respath)
            " files for object "
            files = []

            for tileN in tiles:
                if tileN == 'None': continue
                tile_dat = self.tiles[tileN]

                for tileinfo in tile_dat:
                    band = tileinfo['BAND']
                    if band == None or band == 'det':
                        continue
                    fileN = tileinfo['FILENAME']
                    if fileN not in files:
                        files.append(fileN)                    

                    tile = tileinfo['TILENAME']
                    stampname = self.curdir+'/workdir/'+tile+'_'+ob+'_'+band+'_image.fit'
                    infile = self.curdir+'/data/'+ fileN
                    mkfs = MakeCutOutStamp(0,infile,stampname,self.ssize)
                    mkfs.setStampAtt(ob,band,'IMAGE')
                    mkfs.produce(RA,Dec)
                    (posX,posY) = mkfs.getPos()
                    bandfiles = imTab[band]
                    if stampname not in bandfiles:
                        (imTab[band]).append(stampname)
                    
                    stampname = self.curdir+'/workdir/'+tile+'_'+ob+'_'+band+'_weight.fit'
                    mkfs = MakeCutOutStamp(2,infile,stampname,self.ssize)
                    mkfs.setStampAtt(ob,band,'WEIGHT')
                    mkfs.produce(RA,Dec)
                    bandweights = weightTab[band]
                    if stampname not in bandweights:
                        (weightTab[band]).append(stampname)
                    if tile == tileM:  # create PSF image for the object
                        psfcatF = string.split(fileN,'.fits')[0] + '_psfcat.psf'
                        psfFile = self.curdir+'/data/'+psfcatF
#
                        resFile = respath+objname+'_'+band+'_psf.fits'
                        mkpsf = CreatePSFFile(psfFile,resFile)
                        mkpsf.setPSFAtt(ob, band)
                        mkpsf.produce(posX,posY)
                        print " created PSF file %s \n" % resFile

            bands = imTab.keys()

            for band in bands:

                if len(imTab[band])>=1:
                                          
                    images = self.makeList(imTab[band])
                    weights = self.makeList(weightTab[band])
                    resWeight = respath+objname+'_'+band+'_weight.fits'
                    resIm = respath+objname+'_'+band+'_image.fits'
                    self.SWARPcaller(images,weights,resIm,resWeight,RA,Dec)
                if self.resType == 1:
                    """ Now add created images to multi-extension fits file """
                    self.writeResFile(resOutFile,ob,band,tile)
        """ remove object data if type == 1 """
        if self.resType == 1:
            for ob in objs:                 
                obj_path = self.curdir + '/'+ob + '/'
                shutil.rmtree(obj_path)

        """ now clean up workdir """            
        dirpath = self.curdir+'/workdir/'
        for filename in os.listdir(dirpath):
            filepath = os.path.join(dirpath, filename)
            try:
                shutil.rmtree(filepath)
            except OSError:
                os.remove(filepath)
        
        
    def run(self):
        self.get_tilegeom()
        self.get_objForTile()
#        print " self.tiles \n"
#        print self.tiles
        print " tilesForObj  \n"
        print self.tilesForObj 
        print " objForTile \n"
        print self.objForTile     
        dirpath = self.curdir+'/data/'
        " Now try to make cutouts based on tile info"

        for key in self.objForTile:          # loop on tiles
            self.getObjData(key)            
            self.makeTileCutouts(key)
            
            " Now remove data files for given tile"       
            for filename in os.listdir(dirpath):
                filepath = os.path.join(dirpath, filename)
                try:
                    shutil.rmtree(filepath)
                except OSError:
                    os.remove(filepath)

        
         
    """ The method to create connection with url server.
        The user set password and create OpenerDirector to use the authentication
        for all consecutive transfers """
                
    def get_connected(self,urlbase,username,password):
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        opener = urllib2.OpenerDirector()
        try:
            password_mgr.add_password(None, urlbase, username, password)

            handler = urllib2.HTTPBasicAuthHandler(password_mgr)
            sslhandler =urllib2.HTTPSHandler(context=ctx)
        # create "opener" (OpenerDirector instance)
            opener = urllib2.build_opener(sslhandler,handler)
            print " created opener"
            urllib2.install_opener(opener)
        except urllib2.HTTPError as e:
            print e.code
            print e.read()
        else:
            print "Everething is OK opener is created"
                   
            return opener

    """ Simple method to check if the url file exists """
    def file_exists(self,infile):
        print " file_exists %s \n" % infile
        request = urllib2.Request(infile)
        request.get_method = lambda : 'HEAD'
        try:
            response = urllib2.urlopen(request)
            return True
        except:
            print "file exists failed \n"
            return False
        
    """ This method copies a single file from a url using provided opener """ 
    def get_file(self,opener,infile, outfile):
        success = False
        comp_file = outfile
        " Check if the file already exists "
        if outfile.find('.fz') > 0:
            uncompout = string.split(outfile,'.fz')[0]
        else:
            uncompout = outfile
        if os.path.isfile(uncompout):
            print "file %s already exists \n" % uncompout
            return True
        "     "
        file_name = outfile.split('/')[-1]
        if not os.path.exists(outfile):
                
            dirname   = os.path.dirname(outfile)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
                    
        f = open(comp_file, 'wb')
        file_size_dl = 0
        
        if opener :
            try:
                u = opener.open(infile)
            except urllib2.URLError as e:
                print "exception in get_file"
                if hasattr(e, 'reason'):
                    print 'We failed to reach a server.'
                    print 'Reason: ', e.reason
                elif hasattr(e, 'code'):
                    print 'The server couldn\'t fulfill the request.'
                    print 'Error code: ', e.code
            else:
                print "authentication is OK loading file"   
                meta = u.info()
                file_size = int(meta.getheaders("Content-Length")[0])
                print "Downloading: %s Bytes: %s" % (file_name, file_size)
                for tries in range(3):
                    copied = False
                    try:
                        file_size_dl = 0
                        block_sz = 8192
#                        block_sz = 1024
                        while True:
                            Fbuffer = u.read(block_sz)
                            if not Fbuffer:
                                copied = True
                                success = True
                                break

                            file_size_dl += len(Fbuffer)
                            f.write(Fbuffer)
                            status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
                            status = status + chr(8)*(len(status)+1)
                            print status,
                    except Exception as e:
                        lasterr = str(e).strip()
                        print(colored("Error when trying to copy file: %s" % lasterr, "red"))
                        print("\n   Retrying...\n")
                    if copied :
                        success = True
                        break
                    print 'Try again '
                    sleep( 8 )
                if not copied:
                    print('\n ** Failed to copy file  %s \n' % file_name)
                    return success
        f.close()
        return success
    


 
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    print nbpar
    if nbpar < 5:
        "Usage: MakeCoaddCutoutsY6.py <required inputs>"
        print "  Required inputs:"
        print "  -l <objlist> - list of objects to process like sl1_list.csv"
        print "  -s <stampsize> - size of stamp in arcsec like 30"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hl:s:r:",
                    ["objlist=","stampsize=","res_type="])
    except getopt.GetoptError:
        print "Usage: MakeCoaddCutoutsY6.py <required inputs>"
        print "  Required inputs:"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage: MakeCoaddCutoutsY6.py <required inputs>"
            print "  Required inputs:"
            print "  -l <objlist> - list of tiles to process"
            print "  -s <stampsize> - size of stamp in arcsec"
            print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
            sys.exit(2)
        elif opt in ("-l","--objlist"):
            list_flag = 1
            objlist = arg
        elif opt in ("-s","--stampsize"):
            stamp_flag = 1
            ssize = arg  
        elif opt in ("-r","--res_type"):
            res_flag = 1
            resType = arg  
    inpsum = list_flag + stamp_flag + res_flag
    print " list_flag=%d stamp_flag=%d res_flag=%d\n" %(list_flag,stamp_flag,res_flag)
    if inpsum != 3:
        print "Usage: MakeCoaddCutoutsY6.py <required inputs>"
        print "  Required inputs:"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(-2)
    cutM = MakeCoaddCutoutsY6(objlist,ssize,resType)
    cutM.get_objects()
    cutM.run()            