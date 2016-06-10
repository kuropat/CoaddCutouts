#!/usr/bin/env python

"""
 The program to create a set of fits cutouts for provided RA,Dec coordinates
 Usage:
    MakeCoaddCutouts.py   -t <tagname> -d <db archive>  -l <objlist> -s <stampsize> -r <res_flag>
    Here:
          <tabname> is prefix of  coadd table name like Y1A1  
          <db_archive> is the database name where we will search for files like destest or desoper
          <objlist> list of objects in the form objName,RA,Dec (RA,Dec in degrees)
          <stampsize> is size of the cutout stamp in arcsec
          <res_flag> if 0 each object cutouts are stored in subdirectory with the object name,
                     if 1 results are stored in multi-extension fits file named by tile 
    By N. P. Kuropatkin  May 11 2016

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
import urllib2
import config_ea as config_mod
import cx_Oracle
try:
    from termcolor import colored
except:
    def colored(line, color):
        return line





class MakeCoaddCutouts():
    
    def __init__(self,dbtab, dbname,  objlist, ssize, resType):
        '''
        Constructor
        '''
        urlbase = "https://desar2.cosmology.illinois.edu/DESFiles/"
        self.objlist = objlist
        self.dbname = dbname
        self.dbtab = dbtab
        self.resType = string.atoi(resType)
        self.pixscale=0.2636
        self.spix = int(string.atof(ssize)/self.pixscale)  # stamp size in pixels
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
        self.tileinfo = {}
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
        pathT = """https://desar2.cosmology.illinois.edu/DESFiles/desardata/"""
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
        stampPix = self.spix 

        tiles = {}
        obj_rec = []
        objkeys = self.objects.keys()
        for okey in objkeys:
            tile_list = []
            ObjName = okey
            obj_rec = self.objects.get(okey)
            SRA = obj_rec[0]
            SDEC = obj_rec[1]
            SRAmin = SRA - stampHS
            SRAmax = SRA + stampHS
            SDecmin = SDEC - stampHS
            SDecmax = SDEC + stampHS
            
            tileinfo = {}
            tileinfo = self.get_tile(SRA,SDEC)
            tilename = tileinfo.get("TILENAME")
            tile_list.append(tilename)
            RALL = tileinfo.get("RALL")
            DECLL = tileinfo.get("DECLL")
            RAUL = tileinfo.get("RAUL")
            DECUL = tileinfo.get("DECUL")
            RAUR = tileinfo.get("RAUR")
            DECUR = tileinfo.get("DECUR")
            RALR = tileinfo.get("RALR")
            DECLR = tileinfo.get("DECLR")
            URALL = tileinfo.get("URALL")
            UDECLL = tileinfo.get("UDECLL")
            URAUR = tileinfo.get("URAUR")
            UDECUR = tileinfo.get("UDECUR")
            keys = tiles.keys()
            if tilename not in keys:
                tiles.update({tilename:tileinfo})
            # Check if we need neighbor tiles 
            if SRAmin < RALL and SDecmin > DECLL: # West
                tileinfo1 = self.get_tile(SRAmin,SDEC)
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
            elif  SRAmin < RALL and SDecmin < DECLL: # Left bottom corner
                tileinfo1 = self.get_tile(SRAmin,SDEC) # West
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
                tileinfo2 = self.get_tile(SRAmin,SDecmin) # South - West
                tilename = tileinfo2.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo2})
                tileinfo3 = self.get_tile(SRA,SDecmin) # South 
                tilename = tileinfo3.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo3})
            elif  SRAmin < RAUL and SDecmax > DECUL: # Left upper corner
                tileinfo1 = self.get_tile(SRAmin,SDEC) # west
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
                tileinfo2 = self.get_tile(SRAmin,SDecmax) # West - North
                tilename = tileinfo2.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo2})
                tileinfo3 = self.get_tile(SRA,SDecmax) # North
                tilename = tileinfo3.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo3})
            elif SRAmin > RALL and SRAmin <RALR and  SDecmin < DECLR: # bottom
                tileinfo1 = self.get_tile(SRA,SDecmin) # South
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
            if SRAmax < RAUR and SRAmax > RAUL and SDecmax > DECUR:    # North
                tileinfo1 = self.get_tile(SRA,SDecmax)
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
            elif  SRAmax > RAUR and SDecmax > DECUR: # right upper corner
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
                tileinfo2 = self.get_tile(SRAmax,SDecmax) # North - East
                tilename = tileinfo2.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo2})
                tileinfo3 = self.get_tile(SRA,SDecmax) # North
                tilename = tileinfo3.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo3})
            elif  SRAmax > RALR and SDecmin < DECLR: # Right lower corner
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
                tileinfo2 = self.get_tile(SRAmax,SDecmin) # South - East
                tilename = tileinfo2.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo2})
                tileinfo3 = self.get_tile(SRA,SDecmin) # South
                tilename = tileinfo3.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo3})                    
            elif SRAmax > RAUR and SDecmax < DECUR and  SDecmin > DECLR: # Right
                tileinfo1 = self.get_tile(SRAmax,SDEC) # East
                tilename = tileinfo1.get("TILENAME")
                tile_list.append(tilename)
                keys = tiles.keys()
                if tilename not in keys:
                    tiles.update({tilename:tileinfo1})
                        
            self.tilesForObj.update({ObjName:tile_list})            
        keys_ob = self.tilesForObj.keys()
        return tiles
    
    """ query DB for tile containing given point """
    def get_tile(self,RA,DEC):
        query_geom = """ select TILENAME,PROJECT,RA,DEC,PIXELSIZE,NPIX_RA,NPIX_DEC,URALL,UDECLL,URAUR,UDECUR,
        RALL,DECLL,RAUL,DECUL,RAUR,DECUR,RALR,DECLR
                     from coaddtile where URALL<{SRA} and UDECLL<{SDEC} and URAUR>{SRA} and UDECUR>{SDEC}"""
#        tile_inf = {}
        query = query_geom.format(SRA = str(RA),SDEC = str(DEC))
        self.cur.execute(query)
        desc = [d[0] for d in self.cur.description]
        line = self.cur.fetchone()
        tile_inf = dict(zip(desc,line))
        return tile_inf
    
    """ Creates a list of objects in a tile """
    def get_objForTile(self):
        inobj = self.objects
        objkeys = inobj.keys()
        obj_dict = self.tilesForObj.copy()
        for okey in objkeys:  # loop ob objects
            objID = okey
            obj_list = []
            tileListForObj =  self.tilesForObj.get(objID)         
            tileID0 =  tileListForObj[0] # get first object
            
            keys = obj_dict.keys()  # list of tiles in object
            if objID in keys:
                for key in keys:  # loop on objects
                    tileID = (obj_dict.get(key))[0]
                    if tileID == tileID0:
                        obj_list.append(key)
            self.objForTile.update({tileID0:obj_list})

        
    
    """ get data for all objects in the tile """
    def getObjData(self,key):
        objlist = self.objForTile.get(key)
        tile_list = []
        for obj in objlist:
            tiles = self.tilesForObj.get(obj)
            for otile in tiles:
                if otile not in tile_list:
                    tile_list.append(otile)
        " Now we have a list of tiles lets get data" 
        for tile in tile_list:
                self.getTileData(tile)
        print "All data for tile are collected \n" 
               
    """ get data for given tile: find list of images and list of psfcat files """        
    def getTileData(self,tile):
        print "Get tile data for tile %s \n" % tile
        queryT = """SELECT IMAGENAME,PATH from {TABLE} where TILENAME='{TILENAME}' and
                 IMAGENAME not like '%det.fits'"""
        queryC = """SELECT PATH,CATALOGNAME from {TABLE} where TILENAME='{TILENAME}' and
                 CATALOGNAME like '%_psfcat.fits'"""

        query = queryT.format(TABLE = self.dbtab+'_COADD',TILENAME = tile)
        dirpath = self.curdir+'/data/'
        filepaths = []
        catpaths = []
        
        try:
            self.cur.execute(query)
                
            res = self.cur.fetchall()
            for row in res:
                fpath = row[1]
                filepath = self.makeFilePath( fpath);
                filepaths.append(filepath)  
        except:
            print(' Failed on tile %s \n' %tile)
        
        query1 = queryC.format(TABLE=self.dbtab+'_CATALOG',TILENAME = tile)
        try:
            self.cur.execute(query1)
            res = self.cur.fetchall()
            for row in res:
                cpath = row[0]
                " Have to substitute .fits with .psf due to database error "
                if string.find(cpath,'.fits') > 0 :
                    cpath = string.split(cpath,'.fits')[0]+'.psf'
                catpath = self.makeFilePath(cpath);
                catpaths.append(catpath)  
        except:
            print(' Failed on tile %s \n' %tile)
        for fpath in filepaths:
            fname = string.split(fpath,'/')[-1]
            lname = dirpath+fname

            if self.file_exists(fpath): 
                self.get_file(self.opener,fpath,lname)
                
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
        fits1 = fitsio.FITS(imfile,'rw')
        imhdr = fits1[0].read_header()
        data =  fits1[0].read()
        self.writeFits(resOutPath,data,imhdr,objName,band,im_type,first)
        first = False
        imfile = respath+objName+'_'+band+'_weight.fits'
        im_type = 'WEIGHT'
        fits1 = fitsio.FITS(imfile,'rw')
        imhdr = fits1[0].read_header()
        data =  fits1[0].read()
        self.writeFits(resOutPath,data,imhdr,objName,band,im_type,first)
        psffile = respath+objName+'_'+band+'_psf.fits'
        im_type = 'PSF'
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

            
    """ create cutouts in the wordir """
    def makeTileCutouts(self,tile):
        resOutFile = tile+'_cutouts.fits'
        objs = self.objForTile.get(tile)
        for ob in objs:      # Loop on objects in the tile
            objname = ob
            RA = self.objects.get(ob)[0]
            Dec = self.objects.get(ob)[1]            
            imTab = {'g':[],'r':[],'i':[],'z':[],'Y':[]}
            weightTab = {'g':[],'r':[],'i':[],'z':[],'Y':[]}
            tiles = self.tilesForObj.get(objname)
            tileM = tiles[0] # master tile
            respath = self.curdir+'/'+objname+'/'
            if not os.path.exists(respath):
                os.mkdir(respath)
            for tile in tiles:    # loop on all tiles participating in object
                for band in self.bands:
                    stampname = self.curdir+'/workdir/'+tile+'_'+ob+'_'+band+'_image.fit'
                    infile = self.curdir+'/data/'+tile+'_'+band+'.fits.fz'
                    mkfs = MakeCutOutStamp(0,infile,stampname,self.ssize)
                    mkfs.setStampAtt(ob,band,'IMAGE')
                    mkfs.produce(RA,Dec)
                    (posX,posY) = mkfs.getPos()
                    (imTab.get(band)).append(stampname) 
                    
                    stampname = self.curdir+'/workdir/'+tile+'_'+ob+'_'+band+'_weight.fit'
                    mkfs = MakeCutOutStamp(1,infile,stampname,self.ssize)
                    mkfs.setStampAtt(ob,band,'WEIGHT')
                    mkfs.produce(RA,Dec)
                    (weightTab.get(band)).append(stampname)
                    if tile == tileM:  # create PSF image for the object
                        psfFile = self.curdir+'/data/'+tile+'_'+band+'_psfcat.psf'
#                        outfileF = self.curdir+'/workdir/'+objname+'_'+band+'_psf.fits'
                        resFile = respath+objname+'_'+band+'_psf.fits'
                        mkpsf = CreatePSFFile(psfFile,resFile)
                        mkpsf.setPSFAtt(ob, band)
                        mkpsf.produce(posX,posY)

            bands = imTab.keys()
            for band in bands:
                if len(imTab.get(band)) == 1: # Just copy stamp in proper place
                    imageName = imTab.get(band)
                    weightName = weightTab.get(band)
                    resWeight = respath+objname+'_'+band+'_weight.fits'
                    resIm = respath+objname+'_'+band+'_image.fits'
                    shutil.move(imageName[0], resIm)
                    shutil.move(weightName[0], resWeight)
                else:
                    images = self.makeList(imTab.get(band))
                    weights = self.makeList(weightTab.get(band))
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
        dirpath = self.curdir+'/data/'
        " Now try to make cutouts based on tile info"
 
        keys = self.objForTile.keys()
        for key in keys:          # loop on tiles
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
       
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        opener = urllib2.OpenerDirector()
        try:
            password_mgr.add_password(None, urlbase, username, password)

            handler = urllib2.HTTPBasicAuthHandler(password_mgr)

        # create "opener" (OpenerDirector instance)
            opener = urllib2.build_opener(handler)
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
        "Usage: MakeCoaddCutouts.py <required inputs>"
        print "  Required inputs:"
        print "  -t <dbtab> -coadd table name like Y1A1"
        print "  -d <db archive> - archive name like desoper"
        print "  -l <objlist> - list of objects to process like sl1_list.csv"
        print "  -s <stampsize> - size of stamp in arcsec like 30"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"ht:d:l:s:r:",
                    ["tabname=","db_archive=","objlist=","stampsize=","res_type="])
    except getopt.GetoptError:
        print "Usage: MakeCoaddCutouts.py <required inputs>"
        print "  Required inputs:"
        print "  -t <dbtab> -coadd table name like Y1A1"
        print "  -d <db archive> like desoper"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage: MakeCoaddCutouts.py <required inputs>"
            print "  Required inputs:"
            print "  -t <dbtab> -coadd table name like Y1A1"
            print "  -d <db archive> like desoper"
            print "  -l <objlist> - list of tiles to process"
            print "  -s <stampsize> - size of stamp in arcsec"
            print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
            sys.exit(2)
        elif opt in ("-t","--dbtab"):
            t_flag = 1
            dbtab = arg 
        elif opt in ("-d","--db_archive"):
            db_flag = 1
            db_name = arg
        elif opt in ("-l","--objlist"):
            list_flag = 1
            objlist = arg
        elif opt in ("-s","--stampsize"):
            stamp_flag = 1
            ssize = arg  
        elif opt in ("-r","--res_type"):
            res_flag = 1
            resType = arg  
    inpsum = t_flag + db_flag + list_flag + stamp_flag + res_flag
    print " t_flag=%d db_flag=%d list_flag=%d stamp_flag=%d res_flag=%d\n" %(t_flag,db_flag,list_flag,stamp_flag,res_flag)
    if inpsum != 5:
        print "Usage: MakeCoaddCutouts.py <required inputs>"
        print "  Required inputs:"
        print "  -t <dbtab> -coadd table name like Y1A1"
        print "  -d <db archive> - database name like desoper"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <resType> - if 0 - each object data in separate directory 1- multi-extension fits file"
        sys.exit(-2)
    cutM = MakeCoaddCutouts(dbtab,db_name,objlist,ssize,resType)
    cutM.get_objects()
    cutM.run()            