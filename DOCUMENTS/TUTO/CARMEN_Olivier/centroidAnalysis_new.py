
import os
import numpy
import aotools
from astropy.io import fits
import matplotlib.pylab as plt
import h5py
plt.ion()


def calcR0(cents,pxlscale,subD,lamda=500E-9):
    cents =cents.astype(numpy.float16)
    r0 = ((numpy.std(cents*pxlscale)*numpy.pi/(180.*3600.))**2/(0.162*(lamda**2)*(subD**(-1./3.))))**(-3./5.)
    return r0

def calcWFE(cents,pxlscale,d,D=4.2,lamda=500E-9):
    r0 = calcR0(cents,pxlscale,subD)
    WFE = numpy.sqrt(1.03*(D/r0)**(5./3.))*lamda/(1E-9*2*numpy.pi)
    return WFE

def variance_pupil(img,pupil=None):
    if pupil is None:
        pupil=numpy.where(img!=0.,1.,0.)   # cauition with strong scintillation, phase patters go -ve
    temp=numpy.where(pupil!=0,(img*pupil-(img*pupil).sum()/pupil.sum())**2,0.)
    return temp.sum()/pupil.sum()



if __name__ == '__main__':
    # data directory
    dataDir = 'Data_tomo/'
    filename = 'datatomo_2013-05-21_12h37m50s_val00.fits'

    # dataDir = '/Users/jo/source/data_mining/canary/bench/datatomo_july2013/'
    # filename = 'datatomo_2013-07-21_21h20m47s_.fits'

    fnlist = os.listdir(dataDir)
    fnlist.sort()
    #fnlist = [filename]

    LA_altitude = []
    NN_altitude = []
    altitude=[]
    for ifile in fnlist:
        if 's_val' in ifile:

            print dataDir+ifile
            data = fits.getdata(dataDir+ifile)
            header = fits.getheader(dataDir+ifile)
            print ifile
            print data.shape

            # index of wavefront sensor (truth sensor)
            IWFS = header['ITS']
            # centroid data
            cents = data[:,72*(IWFS-1):72*IWFS]
            # pixel scale
            pxlscale = header['PIXARC%i'%(IWFS)]
            # subaperture size
            subD = 4.2/float(header['WFSSUBX'].split()[IWFS-1])

            # Calculate Learn and apply slopes
            zmat=fits.getdata('GLOB_mrz_7x7.fits')
            la_mt=fits.getdata('mt_2013-05-25_21h13m02s_James.fits')
            #LAcents = la_mt

            # L&A cents from off-axis WFS
            LAcents = numpy.dot(data[:,:-72],la_mt)
            #file = h5py.File('Compare.h5', 'r')
            #LAcents = file["sim"][:,:]

            ########################
            # NN centroids
            # Read in NN centroids here (set to copy L&A centroids for now)
            #NNcents = fits.getdata('CARMEN_datatomo_2013-05-21_12h37m50s_val00.fits')
            # NNcents = fits.getdata('your_fits_file_here.fits') # size 10000x72
            file = h5py.File('Compare.h5', 'r')
            NNcents = file["red"][:,:]
            ########################

            # Calculate wavefronts
            zerns = aotools.zernikeArray(zmat.shape[1],100)
            phsTruthTotal = numpy.zeros((100,100),numpy.float)
            phsLATotal = numpy.zeros((100,100),numpy.float)
            phsLAresidual = numpy.zeros((100,100),numpy.float)
            phsNNTotal = numpy.zeros((100,100),numpy.float)
            phsNNresidual = numpy.zeros((100,100),numpy.float)

            TruthWFETotal = 0.
            LAWFETotal = 0.
            LAresidualWFETotal = 0.
            NNWFETotal = 0.
            NNresidualWFETotal = 0.

            # plot phase or not?
            # 1 for plotting
            # 0 for not plotting
            verbose = 0

            for i in range(data.shape[0]):
                # reconstruct zernike coefficients
                zCoeffTruth = numpy.dot(cents[i],zmat)* (2.1E9 * pxlscale / (3600.*180./numpy.pi)) *(2*numpy.pi/500.)
                zCoeffLA = numpy.dot(LAcents[i],zmat)* (2.1E9 * pxlscale / (3600.*180./numpy.pi)) * (2*numpy.pi/500.) 
                zCoeffNN = numpy.dot(NNcents[i],zmat)* (2.1E9 * pxlscale / (3600.*180./numpy.pi)) * (2*numpy.pi/500.) 


                # calculate phase
                phsTruth = (zerns*zCoeffTruth[:,None,None]).sum(0)
                phsLA = (zerns*zCoeffLA[:,None,None]).sum(0)
                phsNN = (zerns*zCoeffNN[:,None,None]).sum(0)

                phsTruthTotal += phsTruth
                phsLATotal += phsLA
                phsNNTotal += phsNN

                phsLAresidual += (phsTruth - phsLA)
                phsNNresidual += (phsTruth - phsNN)

                # calaculate WFE
                TruthWFE = numpy.sqrt(variance_pupil(phsTruth,zerns[0]))*500/(2*numpy.pi)
                TruthWFETotal += TruthWFE

                LAWFE = numpy.sqrt(variance_pupil(phsLA,zerns[0]))*500/(2*numpy.pi)
                LAresidualWFE = numpy.sqrt(variance_pupil(phsTruth-phsLA,zerns[0]))*500/(2*numpy.pi)
                LAWFETotal += LAWFE
                LAresidualWFETotal += LAresidualWFE

                NNWFE = numpy.sqrt(variance_pupil(phsNN,zerns[0]))*500/(2*numpy.pi)
                NNresidualWFE = numpy.sqrt(variance_pupil(phsTruth-phsNN,zerns[0]))*500/(2*numpy.pi)
                NNWFETotal += NNWFE
                NNresidualWFETotal += NNresidualWFE


                if verbose:
                    print 'Truth WFE: %.3f, L&A WFE: %.3f, L&A residual WFE: %.3f, NN WFE: %.3f, NN residual WFE: %.3f'%(TruthWFE,LAWFE,LAresidualWFE, NNWFE, NNresidualWFE)
                    plt.figure('phs truth')
                    plt.clf()
                    plt.imshow(phsTruth)
                    plt.colorbar()
                    plt.figure('phs')
                    plt.clf()
                    plt.subplot(221)
                    plt.imshow(phsLA)
                    plt.colorbar()
                    plt.title('L&A')
                    plt.subplot(222)
                    plt.imshow(phsTruth - phsLA)
                    plt.colorbar()
                    plt.title('L&A residual')
                    plt.subplot(223)
                    plt.imshow(phsNN)
                    plt.colorbar()
                    plt.title('NN')
                    plt.subplot(224)
                    plt.imshow(phsTruth - phsNN)
                    plt.colorbar()
                    plt.title('NN residual')

                    plt.draw()
                    raw_input()

            
            TruthWFETotal/=float(data.shape[0])
            LAWFETotal/=float(data.shape[0])
            LAresidualWFETotal/=float(data.shape[0])
            NNWFETotal/=float(data.shape[0])
            NNresidualWFETotal/=float(data.shape[0])

            print 'Truth WFE: %.3f, L&A WFE: %.3f, L&A residual WFE: %.3f, NN WFE: %.3f, NN residual WFE: %.3f'%(TruthWFETotal,LAWFETotal,LAresidualWFETotal,NNWFETotal,NNresidualWFETotal)

            LA_altitude.append(LAresidualWFETotal)
            NN_altitude.append(NNresidualWFETotal)
            altitude.append(int(ifile[-7:-5]))

    plt.figure('WFE Altitude')
    plt.plot(altitude,LA_altitude,'rx')
    plt.plot(altitude,NN_altitude,'kx')
    plt.xlabel('altitude (mm on bench)')
    plt.ylabel('WFE (nm)')
    #plt.savefig('~/Dropbox/ANN_altitude.pdf')
    plt.draw()
    raw_input()
    
    

    








