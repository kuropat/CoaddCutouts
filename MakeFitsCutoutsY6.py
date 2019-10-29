#!/usr/bin/env python

"""
 The program to create a set of unrecfiyied fits cutouts for provided RA,Dec coordinates using dessci Y6 database
 Created fits cutouts will contain original nonlinear WCS coefficients.
 Usage:
    MakeFitsCutoutsY6.py   -w <workdir> -d <db archive>  -l <objlist> -s <stampsize>  -r <restype>
    Here: <workdir> is current working directory. It could be ./
          <db_archive> is the database name where we will search for files like desoper
          <objlist> list of objects in the form objName,RA,Dec (RA,Dec in degrees)
          <stampsize> is size of the cutout stamp in arcsec 
          <restype> is type of result 0 - cutout stamps are stored in subdirectory with the
                     object name, 1- cutouts are stored in multiextension fits file

    By N. P. Kuropatkin   08/28/2018
   

"""

import os
import sys
import math
import string
import shutil
import getopt
import subprocess
import fitsio
from despyastro import wcsutil

from CreatePSFFile import CreatePSFFile
#
import timeit
from time import sleep
import urllib2, ssl
import easyaccess
#
import cx_Oracle
from tornado.netutil import ssl_match_hostname
try:
    from termcolor import colored
except:
    def colored(line, color):
        return line





class MakeFitsCutoutsY6():
    
    def __init__(self, dbname, workdir, objlist,ssize, restype):
        '''
        Constructor
        '''
        urlbase = "https://desar2.cosmology.illinois.edu/DESFiles/"
        self.workdir = workdir
        self.objlist = objlist
        self.dbname = dbname
        self.pixscale=0.2636
        self.stamp_sec = ssize
        self.spix = int(float(ssize)/self.pixscale)  # stamp size in pixels
        self.ssize = self.spix*self.pixscale
        self.ssize /= 3600. # stamp size in degrees
        self.fedCut = self.ssize/4.
        self.restype = restype
        self.curdir = os.getcwd()
        self.autocommit = True
        self.quiet = False
        self.RA = ''
        self.Dec = ''
        desfile = os.getenv("DES_SERVICES")
        if not desfile: desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
        self.desconfig = easyaccess.config_ea.get_desconfig(desfile, self.dbname)
       
        self.connectDB()
        self.opener = self.get_connected(urlbase,self.user,self.password)
         

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
        
    """ The explist is not used in file selection, it is used just to reject
     coordinates with no exposures """
    def makeExpList(self,RA,Dec):
        explist = []
        ra = float(RA)
        dec = float(Dec)
        factor = math.cos(math.radians(dec))
        ramin = ra - 1.5/factor
        ramax = ra + 1.5/factor
        ramax1 = ramax
        ramax2 = ramax
        ramin1 = ramin
        ramin2 = ramin
        nexp = 0
        case = 0
        if ramin < 0. and ramax <= 0.:
            ramin1 =360. + ramin
            ramax1 = 360. + ramax
        elif ramin<0. and ramax >=0.:
            case = 1
            ramin1 = 0.
            ramax1 = ramax
            ramin2 = 360. + ramin
            ramax2 = 360.
        elif ramax > 360.0 and ramin >= 360.:
            ramin1 = ramin -360.
            ramax1 = ramax - 360.
        elif ramax > 360. and ramin <= 360.:
            case = 1
            ramax1 = ramax - 360.
            ramin1 = 0.
            ramax2 = 360.
            ramin2 = ramin
        decmin = dec -1.5
        decmax = dec + 1.5
#        print "RAMAx=%f RAMin=%f DECmax=%f DECmin=%f \n" % (ramax,ramin,decmax,decmin)
        queryT="""select expnum from Y6A1_exposure where radeg between
         {RAMIN} and {RAMAX} and decdeg between {DECMIN} and {DECMAX} and exptime < 400."""

        query = queryT.format(RAMAX=ramax1,RAMIN=ramin1, DECMAX=decmax,DECMIN=decmin)
        query2 = queryT.format(RAMAX=ramax2,RAMIN=ramin2, DECMAX=decmax,DECMIN=decmin)
#        print query
#        print query2
        try:
            if case == 0:
                self.cur.execute(query)             
                res = self.cur.fetchall()
#            print res

                for row in res:
                    explist.append(str(row[0]))
                    nexp+=1
                    if (nexp > 1000): break
            else:
                self.cur.execute(query)             
                res = self.cur.fetchall()
#            print res

                for row in res:
                    explist.append(str(row[0]))
                    nexp+=1
                    if (nexp > 1000): break
                self.cur.execute(query2)             
                res = self.cur.fetchall()
#            print res
                for row in res:
                    explist.append(str(row[0]))
                    nexp+=1
                    if (nexp > 1000): break
           
            if len(explist) == 0:
                print " No exposures found for this coordinates  \n" 
                return ''
#
            explist.sort()
            return explist
        except:
           
            print("Unexpected error:", sys.exc_info()[0])
            return []
        
    """ Create query string. In this query we select all images for all exposures available for
    given coordinates. There is no EXPLIST in this query """
    def makeQueryString(self,RA,Dec,explist):

        queryT = """  select i.filename,i.band,i.ccdnum,i.nite,f.t_eff,i.expnum,rac2,rac4,decc2,decc4,naxis1,naxis2,ra_cent,dec_cent 
        from Y6A1_image i,Y6A1_qa_summary f,Y6A1_proctag t
         where  t.tag like '%_FINALCUT'  and f.pfw_attempt_id = t.pfw_attempt_id and 
         i.pfw_attempt_id = t.pfw_attempt_id  and i.filetype='red_immask'  and f.t_eff > 0.3 and 
          {RA}< (CASE WHEN i.crossra0='Y'  THEN i.racmax+360.0  ELSE i.racmax  END) 
        and i.racmin < {RA} and i.deccmin < {DEC} and i.deccmax > {DEC}"""

        query = queryT.format(EXPLIST=explist,RA=RA,DEC=Dec)

#        print query
        return query

    """ make file path for file copy """
    def makeFilePath(self, expnum, rev, proc, Year,ccd):
        ccdnum = ccd[1:]
        prepath="https://desar2.cosmology.illinois.edu/DESFiles/desarchive/"
        pathquery="""select a.filename, a.path, a.compression
        from Y6A1_file_archive_info a, Y6A1_image i, Y6A1_proctag t
        where t.tag like '%_FINALCUT'
        and t.pfw_attempt_id=i.pfw_attempt_id
        and i.filetype='red_immask'
        and i.expnum={EXPNUM}
        and i.ccdnum={CCD}
        and i.filename=a.filename
        and a.archive_name='desar2home'"""


        

        query = pathquery.format(EXPNUM=expnum, CCD=ccdnum, REV=rev)
        print query
        comp = ''
        try:
            self.cur.execute(query)
            rows = self.cur.fetchall()
#            print rows
            if len(rows) > 0:
                subpath = rows[0][1]
                comp = rows[0][2]
#                print subpath
                fpath = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/"
                fpath +=subpath
            else:
                comp=''
                fpath=''
        except Exception as e1:
                lasterr = str(e1).strip()
                print(colored("Error when trying to connect to database: %s" % lasterr, "red"))
        print "File path=%s \n" % fpath
        return (fpath,comp)


    
    "query DB to get proper zeropoint "
    def getZeropoint(self,ccdname,expnum):
        ccdnum = string.atoi(ccdname[1:])
        queryT = """ select mag_zero,sigma_mag_zero,source,version,flag from Y6A1_zeropoint where expnum={EXPNUM} and ccdnum={CCDNUM}"""
        query = queryT.format(EXPNUM=expnum, CCDNUM=ccdnum)
#        print query
        mag_zero = 0.
        sigma_mag_zero = 0.
        try:
            self.cur.execute(query)
            rows = self.cur.fetchall()
#            print rows
            if len(rows) > 0:
                flagmin = 1000
                mag_zero = 0.
                sigma_mag_zero = 0.
                for row in rows:                     
                    source = row[2]
                    flag = int(row[4])
#                    print "source %s flag %d \n" % (source,flag)
                    if flag < flagmin and source == 'FGCM':
                        flagmin = flag
                        mag_zero = float(row[0])
                        sigma_mag_zero = float(row[1])
            else:
                mag_zero = 0.
                sigma_mag_zero = 0.
        except Exception as e1:
                lasterr = str(e1).strip()
                print(colored("Error when trying to connect to database: %s" % lasterr, "red"))
        print " At the end zeropoint mag_zero=%s \n" % mag_zero
        return (mag_zero,sigma_mag_zero)
  
    "query DB to get t_eff "
    def getTeff(self,expnum):      
        queryT = """select f.t_eff from Y6A1_qa_summary f, Y6A1_proctag t where  t.pfw_attempt_id=f.pfw_attempt_id and t.tag like '%_FINALCUT' and expnum={EXPNUM} """

        query = queryT.format(EXPNUM=expnum)
        teff = 0.
        try:
            self.cur.execute(query)
            row = self.cur.fetchone()
            teff = float(row[0])
        except Exception as e1:
                lasterr = str(e1).strip()
                print(colored("Error when trying to connect to database: %s" % lasterr, "red"))
        return teff 
    
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
    

    def  sortIm(self,imperobj):
        halfstamp = self.ssize/2
        raobj = float(self.RA)
        decobj = float(self.Dec)
        sramin = raobj - halfstamp
        sramax = raobj + halfstamp
        sdmin = decobj - halfstamp
        sdmax = decobj + halfstamp

        bands = ['g','r','i','z','Y','u']
        imforband = {}
        bandfiles = {'g':[],'r':[],'i':[],'z':[],'Y':[],'u':[]}
        for fileN in imperobj:
            try:
                select = False
                imb = imperobj[fileN][0]
                imt = float(imperobj[fileN][1])
                ramin= float(imperobj[fileN][11])
                ramax= float(imperobj[fileN][12])
            
                if ramin > ramax:
                    ramin -= 360.
                    raobj = float(self.RA) - 360.
                    sramin = raobj - halfstamp
                    sramax = raobj + halfstamp
                decmin = float(imperobj[fileN][13])
                decmax = float(imperobj[fileN][14])


                if sramax<= ramax and sramin >= ramin and sdmin >= decmin and sdmax <= decmax: # completely in CCD
                    select = True
                elif sramax <= (ramax + self.fedCut) and sramin >= ramin : # slightly on right 
                    
                    if sdmax <= decmax + self.fedCut and sdmin >= decmin: # in or slightly up
                        select = True
                    elif sdmin >= (decmin - self.fedCut) and sdmax <= decmax: # in or slightly down
                        select = True
                elif sramin >= (ramin - self.fedCut) and sramax <= ramax:
                   
                    if sdmax <= (decmax + self.fedCut) and sdmin >= decmin:
                        select = True
                    elif sdmin >= (decmin - self.fedCut) and sdmax <= decmax:
                        select = True
                else:
                    select = False
#                print select
                if select:
#                    print "appending file %s for band %s \n" % (fileN,imb)
                    (bandfiles[str(imb)]).append(fileN)
            except:
                print " problem with file %s \n" % fileN
                print imperobj[fileN]
                print(" error:", sys.exc_info()[0])
                continue
        " Now select them by best t_eff"
        for band in bands:
            fileList = bandfiles[band]
            tmax = -99.
            selList = []
            for fileN in fileList:
                imt = float(imperobj[fileN][1])
#                selList = []
#                if band == 'z':
#                print "fileN=%s teff=%f \n"  % (fileN,imt)
                if imt >= tmax:
                    tmax = imt
                    selList = [fileN] + imperobj[fileN]
            imforband[str(band)] = selList
        " at this point we mast select best image for each band "
#        print imforband
        return imforband  
     
    """ call swarp to create a stamp if more than one ccd """
    def SWARPcaller(self,flist,ra_cent,dec_cent,naxis1,naxis2,stampname):
        imName = stampname+'_rect.fits'
        weightName = stampname+'_rect_weight.fits'
#        print "imName=%s weightName=%s \n" % (imName,weightName)
        images = ''
        weights = ''
        first = True
        for fileN in flist:
            imext = "[0]"
            wext = "[2]"
            if first:
                images +=str(fileN+imext)
                weights += str(fileN+wext)
                first = False
            else:
                images += str(','+fileN+imext)
                weights += str(','+fileN+wext )
                
        command = ['swarp',"%s"%images]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','swarp.config')]       
        command +=["-PIXELSCALE_TYPE","MANUAL","-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER_TYPE","MANUAL","-CENTER","%f,%f"%(ra_cent,dec_cent)]
        command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
        command +=["-RESAMPLING_TYPE","LANCZOS3","-RESAMPLE","N"]
        command +=["-SUBTRACT_BACK","N","-BLANK_BADPIXELS","N"]
        command +=["-DELETE_TMPFILES","Y"]
        command +=["-FSCALASTRO_TYPE","VARIABLE"]
        command +=["-COMBINE","Y","-COMBINE_TYPE","AVERAGE"] 
        command +=["-IMAGEOUT_NAME",imName]
        command +=["-WEIGHTOUT_NAME",weightName]
        command +=["-RESAMPLE_SUFFIX","_im.fits"]
        command +=["-WEIGHT_TYPE","MAP_WEIGHT","-WEIGHT_IMAGE","%s"%weights]
        command +=["-HEADER_ONLY","N","-VERBOSE_TYPE","NORMAL"]
        command +=["-COPY_KEYWORDS","OBJECT,SKYSIGMA,SKYBRITE,FWHM"]
        command +=["-WRITE_XML","N"]
        
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
            
    def SEXcaller(self,stampname):
        image = stampname+'.fits'
        weight = stampname+'_weight.fits'
        catName = stampname+'_psfex.fits'
        print "Sex caller catName=%s \n" % catName
        command=['sex',"%s"%image,'-c',os.path.join('./','etc','sexforpsfex.config')]
        command+=['-WEIGHT_IMAGE',"%s"%weight]
        command+=['-CATALOG_NAME',"%s"%catName]
        command+=['-CATALOG_TYPE',"FITS_LDAC",'-WEIGHT_TYPE','MAP_WEIGHT']
        command+=['-PARAMETERS_NAME', os.path.join('./','etc','sex.param_psfex')]
        command+=['-FILTER_NAME', os.path.join('./','etc','gauss_3.0_7x7.conv')]
        command+=['-STARNNW_NAME', os.path.join('./','etc','sex.nnw')]
        command+=['-SATUR_LEVEL', '65000','-VERBOSE_TYPE','NORMAL','-DETECT_MINAREA' , '9']
        print "command= %s \n" % command
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
            

       
        
    """ Use swarp to make cutout """  
    def MakeStamp(self,flist,stampname,band,objname):

        
        images = ''
        weights = ''
        first = True
        for fileN in flist:
            imext = "[0]"
            wext = "[2]"
            if first:
                images +=str(fileN+imext)
                weights += str(fileN+wext)
                first = False
            else:
                images += str(','+fileN+imext)
                weights += str(','+fileN+wext )
#        images = ''
#        weights = ''
        print('IN MAKE STAMP band=%s \n') % band
        print(" Center stamp on RA=%s Dec=%s \n") % (self.RA,self.Dec)
        print images
        print weights

        image = stampname+'.fits'
        weight = stampname+'_weight.fits'
#        print "Size=%d \n" % self.spix
        command = ['swarp',"%s"%images]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','swarp.config')]
        command +=["-PIXELSCALE_TYPE","MANUAL","-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER_TYPE","MANUAL","-CENTER","%s,%s"%(self.RA,self.Dec)]
        command +=['-IMAGE_SIZE',"%d,%d"%(self.spix,self.spix)]
        command +=["-RESAMPLING_TYPE","LANCZOS3","-RESAMPLE","N"]
        command +=["-SUBTRACT_BACK","N","-BLANK_BADPIXELS","N"]
        command +=["-DELETE_TMPFILES","Y"]
        command +=["-FSCALASTRO_TYPE","VARIABLE"]
        command +=["-COMBINE","Y","-COMBINE_TYPE","AVERAGE"] 
        command +=["-IMAGEOUT_NAME",image]
        command +=["-WEIGHTOUT_NAME",weight]
        command +=["-RESAMPLE_SUFFIX","_im.fits"]
        command +=["-WEIGHT_TYPE","MAP_WEIGHT","-WEIGHT_IMAGE","%s"%weights]
        command +=["-HEADER_ONLY","N","-VERBOSE_TYPE","NORMAL"]
        command +=["-COPY_KEYWORDS","OBJECT,SKYSIGMA,SKYBRITE,FWHM"]
        command +=["-WRITE_XML","N"]
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        fits1 = fitsio.FITS(image,'rw')
        imhdr = fits1[0].read_header()
        imhdr['BAND'] = band
        imhdr['OBJECT'] = objname
        imhdr['TYPE'] = 'image'
        data = fits1[0].read()
        os.remove(image)
        fitsio.write(image,data,header=imhdr)
        fits2 = fitsio.FITS(weight,'rw')
        imhdr1 = fits2[0].read_header()
        imhdr1['BAND'] = band
        imhdr1['OBJECT'] = objname
        imhdr1['TYPE'] = 'weight'
        data2 = fits2[0].read()
        os.remove(weight)
        fitsio.write(weight,data2,header=imhdr1)       
        
    def PSFEXcaller(self,stampname):
        " Will create _psfex.psf"
        inputCat = stampname+'_psfex.fits'
        outcatPSF = stampname+'_psfex_out.fits'
        xmlFile = stampname+'_psfex.xml'
#         -PSF_DIR psf  
        command=['psfex','%s'%inputCat,'-c',os.path.join('./','etc','default.psfex')]
#        command+=['-WRITE_XML','Y', '-XML_NAME', '%s'% xmlFile , '-OUTCAT_TYPE', 'FITS_LDAC' ]
#        command+=[ '-OUTCAT_NAME','%s'% outcatPSF]
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        
            
    def clean(self,flist,stampname):
        for fileN in flist:
            try:
                os.remove(fileN)
            except:
                print "remove failed on file %s \n" % fileN               
        " clean all other files we do jnot need"
        try:
            os.remove(stampname+'_rect.fits')
            os.remove(stampname+'_rect_weight.fits')
            os.remove(stampname+'_psfex.fits')
            os.remove(stampname+'_psfex.psf')
        except:
                print "remove failed on file %s \n" % stampname+'_rect.fits'
		
    """ convert exposure number in the observing year"""
    def expToYear(self,expnum):
        Year = 'Y3'
        if expnum <= 200000: Year = 'SV'
        if expnum > 200000 and expnum <= 340000: Year = 'Y1'
        if expnum > 340000 and expnum <= 450000: Year = 'Y2'
	if expnum > 450000 and expnum <= 519543: Year = 'Y3'
	if expnum > 519543 and expnum <= 666742: Year = 'Y4'
	if expnum > 666742 and expnum <= 724364: Year = 'Y5'
	if expnum > 724364 : Year = 'Y6'
        return Year                
                
                
    def run(self):
        runstat = 0
        list_file = self.objlist
        outfile = ''
        for line in open(list_file):
#            line = line.rstrip('\n')
            line = line.rstrip()
            print(' line=%s \n' % line)
            tokens = line.split(',')
            objname = tokens[0]
            self.RA = tokens[1]
            self.Dec = tokens[2]
            byband = {}
            byband['g'] = []
            byband['r'] = []
            byband['i'] = []
            byband['z'] = []
            byband['Y'] = []
            byband['u'] = []
            t0 = timeit.default_timer()
            if not os.path.exists('./'+objname):
                os.mkdir(objname)
            explist = self.makeExpList(self.RA, self.Dec)
#            print explist
            if len(explist) >= 1:
                query = self.makeQueryString(self.RA, self.Dec,explist)
#                print query
                try:
                    t1 = timeit.default_timer()
                    self.cur.execute(query)             
                    res = self.cur.fetchall()
                    print res
                    t2 = timeit.default_timer()
                    print "Query time=%f \n" % (t2 - t1)
                    if len(res) == 0:
                        print " No records found for object %s \n" % line
                        runstat = -999
                        continue
                    " Create directory if not exists "
                    imperobj = {}
                    for row in res:
                        filename = row[0]
                        band = row[1]
                        ccdnum = row[2]
                        if ccdnum == 'c31': continue
                        nite = row[3]
                        teff = float(row[4])
                        expnum = row[5]
                        rac2 = float(row[6])
                        rac4 = float(row[7])
                        decc2= float(row[8])
                        decc4 = float(row[9])
                        naxis1 = int(row[10])
                        naxis2 = int(row[11])
                        ra_cent = float(row[12])
                        dec_cent = float(row[13])
                        print "naxis1=%d naxis2=%d ra_cent=%f dec_cent=%f rac2=%f rac4=%f decc2=%f decc4=%f \n" % (naxis1,naxis2,ra_cent,dec_cent,rac2,rac4,decc2,decc4)
                        # unpack filename
                        tokens2 = filename.split('_')
                        Dpart = tokens2[0]   # D00233091
                        band = tokens2[1]
                        ccdnum = tokens2[2]   # c01
                        expnum = string.atoi(Dpart[1:])
                        Year = self.expToYear(expnum)
                        revp = tokens2[3];   # r1984p01
                        tokens3 = revp.split('p');
                        rev = (tokens3[0])[1:];
                        proc = "p"+tokens3[1];
                        imrec = [band,teff,ccdnum,nite,expnum,naxis1,naxis2,ra_cent,dec_cent,rev,proc,rac2,rac4,decc2,decc4]
                        imperobj[filename] = imrec
                    " We have all images in the dictionary lets select best"
                    print imperobj
                    selfiles = self.sortIm(imperobj)
                    print 'Selected files \n'
                    print selfiles
                    
                    for band in selfiles:
                        if len(selfiles[band]) == 0: continue
                        filename = selfiles[band][0]
                        teff = float(selfiles[band][2])
                        ccdnum = selfiles[band][3]
                        nite= selfiles[band][4]
                        expnum = selfiles[band][5]
                        rev =  selfiles[band][10]
                        proc = selfiles[band][11]
                        revp = 'r'+rev+'p'+proc
                        Year = self.expToYear(expnum)
                        ccdnum = selfiles[band][3]
                        
                        ra_cent = float( selfiles[band][8])
                        dec_cent = float(selfiles[band][9])

                        (filepath,comp) = self.makeFilePath( expnum, rev, proc,Year,ccdnum);
                        print "band=%s filepath=%s comp=%s \n" % (band,filepath,comp)
                        infile = filepath+'/'+filename
                        if comp == '.fz': infile +=comp
                        outfileZ = './'+objname+'/'+filename+'.fz'
                        outfileF = './'+objname+'/'+filename
                        outfile = outfileF
#
                        psfFile = Dpart+'_'+band+'_'+ccdnum+'_'+revp+'_psfexcat.psf'
                        outfileP = './'+ objname +'/' + psfFile
                        copys = False

                        if self.file_exists(infile):
 
                            (mag_zero,sigma_mag_zero) = self.getZeropoint(ccdnum,expnum)
                            " throw away bad exposures "
                            if mag_zero == 0. or teff <= 0.3:
                                print "file rejected by mag_zero \n"
                                continue
                            tag = objname+'_'+str(expnum)+'_'+ccdnum+'_'+str(nite)+'_'+band+'_'+self.RA+'_'+self.Dec
                            psfFile= filename.split('immask')[0] + 'psfexcat.psf'
                            outfileP = './'+ objname +'/' + psfFile
                            psfpath = filepath.split('red')[0]+'psf'
                            psfinF = psfpath+'/'+psfFile
                            print "start coping psf file %s \n" % psfinF
                            print " psf path=%s \n" % psfpath
                            stat = self.get_file(self.opener,psfinF,outfileP)
                            print " got the file \n"
                            print stat
                            if comp == '.fz':
                                outfile = outfileZ
                                copys = self.get_file(self.opener,infile,outfileZ)
                                if not copys: continue
                                " Now uncompress the file as swarp do not work with .fz "
#
                                command = ['funpack',"%s"%outfileZ]
                                try:
                                    subprocess.check_output(command)
                                except subprocess.CalledProcessError as e:
                                    print "error %s"% e
                        
                                outfile = outfileF

                                try:
                                    os.remove(outfileZ)
                                except:
                                    print "remove failed on file %s \n" % outfileZ   
                            else:
                                outfile = outfileF
                                copys = self.get_file(infile,outfileF)
                            
                            if not copys: continue
                            stampname = './'+objname+'/'+tag
                            outPSF = './'+objname+'/'+tag+'_psf.fits'
                            filelist = [outfile]
                            self.MakeStamp(filelist,stampname,band,objname)
#                            print "After MakeStamp \n"                           
                            """ Now create PSF """                           
                            RAV = float(self.RA)
                            DecV = float(self.Dec)
                            print "Start conversion with file %s \n" % stampname+'.fits'
                            fitsfile = os.path.normpath(stampname+'.fits')
                            fits=fitsio.FITS(fitsfile)
                            prihdr = fits[0].read_header()
                            w = wcsutil.WCS(prihdr)
#                        print "RA=%f Dec=%f \n" % (RAV,DecV)
                            (objx,objy) = w.sky2image(RAV,DecV)
                            print "objx=%f objy=%f \n" % (objx,objy)

                            try:
                                psfcatF = stampname+'_psfex.psf'
                                print " psfcat.psf=%s outPSF=%s \n" % (psfcatF,outPSF)
                                mkpsf = CreatePSFFile(outfileP,outPSF)
                                mkpsf.setPSFAtt(objname, band)
                                mkpsf.produce(objx,objy)
                                try:
#                                    os.remove(outfileF)
                                    print "Now removing file %s \n" % outfileF
                                except:
                                    print "remove failed on file %s \n" % outfileF
                            except:
                                print " Failed on psf file %s \n" % outPSF
#                            self.clean(filelist,stampname)
                            byband[band].append((tag,mag_zero,sigma_mag_zero,teff))
                        else:
                            continue   

                except:
                    print(' Failed on object %s \n' %objname)
                    print("Unexpected error:", sys.exc_info()[0])
                    runstat +=1 
                    continue
            else:
                print " No exposures found for Ra=%s Dec=%s \n" % (self.RA,self.Dec)
            
            " Now we have cutouts in ./objname subdir lets make fits file "
            if self.restype == '1':
                keys=byband.keys()
                fitsout = objname+'_cutouts.fits'
                if  os.path.exists(fitsout):           
                    os.remove(fitsout)
                for key in keys:
                    tags = byband.get(key)
                    
                    for t in tags:
                        self.writeResFile(fitsout,objname,t,key)
                shutil.rmtree('./'+objname)
            t3 = timeit.default_timer()

            exectime = float((t3 - t0)/60.)
            print " end work with line %s  exectime=%.2f min \n" % (line,exectime)
        
        SystemExit(runstat)           

    def writeResFile(self,resOutFile,objName,t,band):
        """ Read created stamps and PSF and append them to
         results file with modifying headers if need """
        tag = t[0]
        resOutPath = self.curdir+'/'+resOutFile
        first = False
        if  not os.path.exists(resOutPath): first = True
        
        respath = self.curdir+'/'+objName+'/'
        imfile = respath+tag+'.fits'
        im_type = 'IMAGE'
        if os.path.exists(imfile):
            fits1 = fitsio.FITS(imfile,'r')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,t,band,im_type,first)
            first = False
        else:
            print "The file %s does not exists \n" % imfile
        imfile = respath+tag+'_weight.fits'
        im_type = 'WEIGHT'
        if os.path.exists(imfile):
            fits1 = fitsio.FITS(imfile,'rw')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,t,band,im_type,first)
        else:
            print "The file %s does not exists \n" % imfile
        psffile = respath+tag+'_psf.fits'
        im_type = 'PSF'
        if os.path.exists(psffile):
            fits1 = fitsio.FITS(psffile,'rw')
            imhdr = fits1[0].read_header()
            data =  fits1[0].read()
            self.writeFits(resOutPath,data,imhdr,t,band,im_type,first)
        else:
            print " The file %s does not exists \n" % psffile
            
            
    def writeFits(self,outPath,im_data,header,t,band,im_type,first):
        if not first: fits2 = fitsio.FITS(outPath,'rw')
        header['BAND'] = band
        header['OBJECT'] = t[0]
        header['MAG_ZERO'] = str(t[1])
        header['SIGMA_MZ'] = str(t[2])
        header['T_EFF'] = str(t[3])
        header['TYPE'] = im_type
        if first:
            fitsio.write(outPath,im_data,header=header)
        else:
            fits2.write(im_data,header=header)
        


         
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
            print "Everything is OK opener is created"             
            return opener

    """ Simple method to check if the url file exists """
    def file_exists(self,infile):
        request = urllib2.Request(infile)
        request.get_method = lambda : 'HEAD'
        try:
            response = urllib2.urlopen(request)
            return True
        except Exception as e1:
            lasterr = str(e1).strip()
            print(colored("Error when trying to copy file: %s" % lasterr, "red"))
            print "file exists failed for file %s \n" % infile
            return False
        
    """ This method copies a single file from a url using provided opener """ 
    def get_file(self,opener,infile, outfile):
        success = False
        comp_file = outfile
        file_name = outfile.split('/')[-1]
        dir_name = outfile.split(file_name)[0]
        outpath = os.path.realpath(dir_name)
        if not os.path.exists(outpath):
            print " path do not exists %s \n" % outpath    
            dirname   = os.path.dirname(outpath)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            else:
                print " directory %s exists \n" % dirname        

        fileres = os.path.normpath(outpath+'/'+file_name)
        f = open(fileres, 'wb')
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
                                f.close()
                                break

                            file_size_dl += len(Fbuffer)
                            f.write(Fbuffer)
                            status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
                            status = status + chr(8)*(len(status)+1)
#                            print status,
#                            copied = True
                    except Exception as e:
                        copied = False
                        lasterr = str(e).strip()
                        print(colored("Error when trying to copy file: %s" % lasterr, "red"))
                        print("\n   Retrying...\n")
                    if copied :
                        csize = os.path.getsize(fileres)
                        print " Copyed size = %d orig_size=%d \n" % (csize, file_size)
                        if csize == file_size:
                            success = True
                            f.close()
                            return success
                        else:
                            copied = False
                            success = False
                    print 'Try again '
                    sleep( 8 )
                if not copied:
                    print('\n ** Failed to copy file  %s \n' % file_name)
                    return success
        f.close()
#        print "end of getfile "
        return success
    


 
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 5:
        "Usage: MakeFitsCutoutsY6.py <required inputs>"
        print "  Required inputs:"
        print "  -w <workdir> -workdir name"
        print "  -d <db archive> - archive name like desoper"
        print "  -l <objlist> - list of objects to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <restype> - type of result 0 - subdir. 1. fits file"
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hw:d:l:s:r:",
                    ["workdir=","db_archive=","objlist=","stampsize=","restype="])
    except getopt.GetoptError:
        print "Usage: MakeFitsCutoutsY6.py <required inputs>"
        print "  Required inputs:"
        print "  -w <workdir> -workdir name"
        print "  -d <db archive> like desoper"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <restype> - type of result 0 - subdir. 1. fits file"
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage:MakeFitsCutoutsY6.py  <required inputs>"
            print "  Required inputs:"
            print "  -w <workdir> -workdir name"
            print "  -d <db archive> like desoper"
            print "  -l <objlist> - list of tiles to process"
            print "  -s <stampsize> - size of stamp in arcsec"
            print "  -r <restype> - type of result 0 - subdir. 1. fits file"
            sys.exit(2)
           
        elif opt in ("-w","--workdir"):
            work_flag = 1
            workdir = arg 
        elif opt in ("-d","--db_archive"):
            db_flag = 1
            db_name = arg
        elif opt in ("-l","--objlist"):
            list_flag = 1
            objlist = arg
        elif opt in ("-s","--stampsize"):
            stamp_flag = 1
            ssize = arg  
        elif opt in ("-r","--restype"):
            res_flag = 1
            rtype = arg  
    inpsum = work_flag + db_flag + list_flag + stamp_flag + res_flag

    if inpsum != 5:
        print "Usage:MakeFitsCutoutsY6.py  <required inputs>"
        print "  Required inputs:"
        print "  -w <workdir> -workdir name"
        print "  -d <db archive> - database name like desoper"
        print "  -l <objlist> - list of tiles to process"
        print "  -s <stampsize> - size of stamp in arcsec"
        print "  -r <restype> - type of result 0 - subdir. 1. fits file"       
        sys.exit(-2)
    print " Start with MakeFitsCutouts \n"
    cutM = MakeFitsCutoutsY6(db_name,workdir,objlist,ssize,rtype)
    cutM.run()            
