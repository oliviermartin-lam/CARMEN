## James Osborn, University of Durham
## Easy AO simulation
## Makes several assumptions so be careful 


import matplotlib.pyplot as plt
plt.ion()
import numpy 
import copy
from astropy.io import fits
import scipy
from scipy import ndimage,interpolate


import aotools
from aotools.turbulence import phasescreen
from aotools.turbulence import infinitephasescreen

# Measure WF slopes directly
def MagicShackHartmann(phs,nsubx,subap_n,pupil_map,magicSlopes):
    tmp1=phs*magicSlopes
    tmp2=phs*magicSlopes.transpose()
    cent_temp=numpy.zeros((2,nsubx,nsubx),numpy.float)
    for i in range(nsubx):
        for j in range(nsubx):
            if pupil_map[i,j] > 0:
                cent_temp[0,i,j]=tmp1[i*subap_n:(i+1)*subap_n,j*subap_n:(j+1)*subap_n].sum()
                cent_temp[1,i,j]=tmp2[i*subap_n:(i+1)*subap_n,j*subap_n:(j+1)*subap_n].sum()
    return cent_temp


# Generate the slope array for the magic WFS
def MagicSHslopes(npup,nsubx,pupil,pupil_map):
    slop=numpy.ones((npup,npup),numpy.float)
    slop[:]=(numpy.arange(npup)+1.)
    slop*=pupil
    n=npup/nsubx
    for i in range(nsubx):
        for j in range(nsubx):
            if pupil_map[i,j]>0:
                p=slop[i*n:(i+1)*n,j*n:(j+1)*n].sum()/float(pupil[i*n:(i+1)*n,j*n:(j+1)*n].sum())
                slop[i*n:(i+1)*n,j*n:(j+1)*n]-=p*pupil[i*n:(i+1)*n,j*n:(j+1)*n]
                s2=(slop[i*n:(i+1)*n,j*n:(j+1)*n]**2).sum()
                slop[i*n:(i+1)*n,j*n:(j+1)*n]/=numpy.sqrt(s2)
            else:
                slop[i*n:(i+1)*n,j*n:(j+1)*n]*=0.
    return slop

def zernikeSubapSlopes(zerns,gamx,gamy,nsubx):   # Calculate mean gradients of the Zernike functions over out WFS subapertures 
    nz=zerns.shape[0]
    npup=zerns.shape[1]
    n=npup/nsubx
    area=numpy.zeros((nsubx,nsubx),numpy.float)
    for i in range(nsubx):
        for j in range(nsubx):
            area[i,j]=zerns[0,i*n:(i+1)*n,j*n:(j+1)*n].sum() #*pupil[i,j]

    dzdx=numpy.dot(zerns.transpose(1,2,0),gamx.transpose())
    dzdy=numpy.dot(zerns.transpose(1,2,0),gamy.transpose())

    sub_grad=numpy.zeros((nz,2,nsubx,nsubx),numpy.float)		# Find mean zernike grads 

    for i in range(nsubx):
        for j in range(nsubx):
            if area[i,j]>0.:
                sub_grad[:,0,i,j]=dzdy[i*n:(i+1)*n,j*n:(j+1)*n,:].sum(0).sum(0)/area[i,j]
                sub_grad[:,1,i,j]=dzdx[i*n:(i+1)*n,j*n:(j+1)*n,:].sum(0).sum(0)/area[i,j]
    return sub_grad 

def makeImage(phs,pupil,nsamp):
    n=len(phs)
    amp=numpy.zeros((n,n),numpy.complex)                     # make complex amp in pupil
    amp.real=numpy.cos(phs)
    amp.imag=numpy.sin(phs)
    nfield=len(amp)
    img= numpy.absolute(numpy.fft.fftshift(numpy.fft.fft2(amp*pupil,s=(nsamp*nfield,nsamp*nfield))))**2 /(float(nsamp*nfield)**2)
    return img

# system parameters
D = float(4)                    # Telescope Diameter (m)
nsubx = int(16)                 # Number of SH subapertures across pupil 
npix = nsubx*10                 # Pixels across pupil
nscrn = int(npix*1.5)           # Number of pixels across seed phase screen
lamda = 500E-9                  # Wavelength (m)
fRate = float(100)              # Frame Rate (Hz)
nsteps = int(100)               # Number of iterations in simulation
nZerns = int(20)                # Number of zernike modes in correction
nLayers = 2                     # number of layers
r0Layers = [0.1,0.2]            # list of r0 for each layer (m)
wSpeedLayers = [6.5,10.]        # list of turbulence speed for each layer (m/s)
wDirLayers = [0.,45.]           # list of turbulence direction for each layer (degrees)
AOlatency = 1                   # integer frames
AOgain = 1./(AOlatency+1)       # needs to be chosen carefully for stability of loop
phsSeed = 10 #int(numpy.random.uniform()*100)   # manually set seed to a number to repeat the same phase screens - not working??

# npix/nsubx must be integer
assert (float(npix)/float(nsubx))%1 == 0

# select debug modes (useful for plotting)
debugAtmos = 0
debugAO = 0
debugFunction = 0

# Time varying function
windSpeedFunction = numpy.ones((nLayers,nsteps))
windDirectionFunction = numpy.ones((nLayers,nsteps))
r0Function = numpy.ones((nLayers,nsteps))

# Change turbulence speed, direction and strength during sequence
# Must set function for each layer manually
# wind speed and r0 multiply by original values
windSpeedFunction[0] = 0.5*numpy.sin((numpy.arange(nsteps)/float(nsteps))*2*numpy.pi)+1
r0Function[0] = 0.2*numpy.sin((numpy.arange(nsteps)/float(nsteps))*2*numpy.pi)+1
# wind direction is added to original value
windDirectionFunction[0] = 15*numpy.sin((numpy.arange(nsteps)/float(nsteps))*2*numpy.pi)

if debugFunction:
    plt.figure('functions')
    for ilayer in range(nLayers):
        plt.plot(windSpeedFunction[ilayer])
        plt.plot(windDirectionFunction[ilayer])
        plt.plot(r0Function[ilayer])
    plt.draw()
    raw_input()

# telescope
zerns = aotools.zernikeArray(nZerns,npix)
pupil = zerns[0]
subapertureMap = aotools.pupil.circle(nsubx/2,nsubx)
subapertureMap -= aotools.pupil.circle(2,nsubx)

# magic WFS
magicSlopes=MagicSHslopes(npix,nsubx,pupil,subapertureMap)
gammas = aotools.makegammas(aotools.zernIndex(nZerns)[0])
subGrads = zernikeSubapSlopes(zerns,gammas[0][:nZerns,:nZerns],gammas[1][:nZerns,:nZerns],nsubx)
z=subGrads.reshape(nZerns,2*nsubx*nsubx)
zmat=numpy.linalg.pinv(z).transpose()

# initiate arrays
phsLayers = numpy.zeros((nLayers,npix,npix),numpy.float)
coeffs = numpy.zeros((zerns.shape[0]),numpy.float)
centroidArray = numpy.zeros((nsteps,2*int(subapertureMap.sum())),numpy.float)
totalMoved = numpy.zeros((nLayers))
demandMoved = numpy.zeros((nLayers))
DMphs = numpy.zeros((npix,npix))
longImg = numpy.zeros((2*npix,2*npix))
longResImg = numpy.zeros((2*npix,2*npix))

# initiate lists
coeffRegister = []
phsScrnList = []
stepSizeList = [] 

# initiate variables
totalPhsWfe = 0.
totalResWfe = 0.
spScale = D/npix                # m/pixel

# Initiate phase screens
for ilayer in range(nLayers):
    tpScale = wSpeedLayers[ilayer]/fRate
    stepSizeList.append(tpScale/spScale)
    phsScrnList.append(infinitephasescreen.PhaseScreenKolmogorov(nscrn,spScale,r0Layers[ilayer],100.,random_seed=(phsSeed+ilayer)))

# calibrate WFS
zScale=[]
for izern in range(len(zerns)):
    calibCents = MagicShackHartmann(zerns[izern],nsubx,npix/nsubx,subapertureMap,magicSlopes)[::-1]
    coeffs = numpy.dot(zmat,calibCents.reshape(2*nsubx*nsubx))
    zScale.append(1./coeffs[izern])
zScale[0] = 0.                      # piston

# run loop
for istep in range(nsteps):

    # atmosphere
    phsLayers*=0.
    for ilayer in range(nLayers):

        # interpolate for wind speed
        stepSize = stepSizeList[ilayer] * windSpeedFunction[ilayer][istep]
        if stepSize<0:
            stepSize = 0
            print 'warning: requested step size is negative - setting to zero'
        demandMoved[ilayer] += stepSize
        totalOffset = demandMoved[ilayer]
        nRows = int(totalOffset - totalMoved[ilayer])
        interpRows = (totalOffset-totalMoved[ilayer])%1

        # print istep,stepSize,'total',totalOffset,'nRows',nRows,'interp',interpRows,'totalMoved',totalMoved[ilayer],'total interp',totalMoved[ilayer]+interpRows

        # add rows
        for iRows in range(nRows):
            phsScrnList[ilayer].add_row()
            totalMoved[ilayer]+=1
        
        # get phs
        phs = copy.copy(phsScrnList[ilayer].scrn)

        # interpolate for sub-pixel shift
        f = interpolate.interp2d(numpy.arange(phs.shape[0]),numpy.arange(phs.shape[1]),phs)
        phs2 = f(numpy.arange(phs.shape[0]),numpy.arange(phs.shape[1])-interpRows)

        # rotate with wind direction
        if wDirLayers[ilayer]+windDirectionFunction[ilayer][istep]!=0:
            phs2 = ndimage.rotate(phs2,wDirLayers[ilayer]+windDirectionFunction[ilayer][istep],reshape=False)

        # scale r0
        phs2 = phs2*r0Function[ilayer][istep]**(-5./6.)

        # cut and add to integrated phase
        phs2 = phs2[nscrn/2-npix/2:nscrn/2+npix/2,nscrn/2-npix/2:nscrn/2+npix/2]
        phsLayers[ilayer] = phs2

    # Add phase from all layers
    phsTemp = phsLayers.sum(0)

    if debugAtmos:
        plt.figure('scrn',figsize=(10,5))
        # plt.clf()
        for ilayer in range(nLayers):
            ax = plt.subplot(1,nLayers,ilayer+1)
            ax.cla()
            plt.imshow(phsLayers[ilayer])
            # plt.subplot(122)
            # plt.plot(numpy.arange(0,10),phs[nscrn/2-npix/2:nscrn/2+npix/2,nscrn/2-npix/2:nscrn/2+npix/2][:10,0],'rx')
            # plt.plot(numpy.arange(0,10)-interpRows,phs2[:10,0],'b.')
            plt.pause(0.05)
            # plt.draw()
            # raw_input()
    
    # remove piston
    phsTemp = (phsTemp - (phsTemp*pupil).sum()/pupil.sum())*pupil

    # remove DM phase
    partialCorrectedPhs = copy.copy(phsTemp) - DMphs

    # magic WFS
    slopes = MagicShackHartmann(partialCorrectedPhs,nsubx,npix/nsubx,subapertureMap,magicSlopes)[::-1]

    # populate centroid array
    centroidArray[istep,:int(subapertureMap.sum())] = numpy.reshape(slopes[0][subapertureMap==1],(int(subapertureMap.sum())))
    centroidArray[istep,int(subapertureMap.sum()):] = numpy.reshape(slopes[1][subapertureMap==1],int(subapertureMap.sum()))

    # predicted slopes
    predSlopes = slopes

    # calculate zernike coefficients
    coeffs = numpy.dot(zmat,predSlopes.reshape(2*nsubx*nsubx))*zScale
    coeffRegister.append(coeffs)

    # apply phase to DM
    if len(coeffRegister)>AOlatency:
        # apply correction
        DMphs = DMphs + AOgain*(coeffRegister[0][:,None,None]*zerns).sum(0)
        coeffRegister.pop(0)
     
     # calculate WFE - only after correction has started, 3 is a bit arbitrary?
    if istep>AOlatency*3: 
        phsWFE = numpy.sqrt((phsTemp**2).sum()/pupil.sum())*lamda/(1E-9*2*numpy.pi)
        residualWFE = numpy.sqrt((partialCorrectedPhs**2).sum()/pupil.sum())*lamda/(1E-9*2*numpy.pi)
        totalPhsWfe += phsWFE
        totalResWfe += residualWFE
    else:
        phsWFE = numpy.sqrt((phsTemp**2).sum()/pupil.sum())*lamda/(1E-9*2*numpy.pi)
        residualWFE = phsWFE
    
    # real-time print outs
    if istep%50==0 and istep>0:
        print istep
        print 'WFE (uncorrected): %.3f, WFE (corrected): %.3f'%(totalPhsWfe/float(istep-AOlatency*3),totalResWfe/float(istep-AOlatency*3))

    if debugAO == 1:
        if istep>AOlatency*3: 
            img = makeImage(phsTemp,pupil,2)
            resImg = makeImage(partialCorrectedPhs,pupil,2)
            longImg += img
            longResImg +=resImg
        plt.figure('phs',figsize=(10,4))
        plt.clf()
        plt.subplot(131)
        plt.imshow(phsTemp)
        plt.colorbar()
        plt.title('uncorrected phase')
        plt.subplot(132)
        plt.imshow(DMphs)
        plt.colorbar()
        plt.title('reconstructed phase')
        plt.subplot(133)
        plt.imshow(partialCorrectedPhs)
        plt.colorbar()
        plt.title('corrected phase')
        plt.figure('image',figsize=(10,5))
        plt.clf()
        plt.subplot(121)
        plt.imshow(longImg/(istep*pupil.sum()))
        plt.title('long image')
        plt.colorbar()
        plt.subplot(122)
        plt.title('long corrected image')
        plt.imshow(longResImg/(istep*pupil.sum()))
        plt.colorbar()

        plt.pause(0.05)
        # plt.draw()
        # raw_input()

totalPhsWfe /= float(nsteps-AOlatency*3)
totalResWfe /= float(nsteps-AOlatency*3)
print istep
print 'WFE (uncorrected): %.3f, WFE (corrected): %.3f'%(phsWFE,residualWFE)

plt.figure('centroids')
plt.imshow(centroidArray)
plt.draw()

fits.writeto('centroids_'+str(phsSeed)+'.fits',centroidArray,overwrite=True)

print 'done'
raw_input()
