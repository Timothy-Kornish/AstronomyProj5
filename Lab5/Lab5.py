import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle
from photutils import daofind, aperture_photometry, CircularAperture

############ Look at individual images ##############

dir_t = '/Volumes/TIM/Transits/'
dir_trans='ABand_VIC_x_o_20141120_0'

dir_d = 'Darks_VIC_x_d_20141120_072037.835094_0010.fits'
dir_f ='flatR_VIC_x_o_20141120_010910.893597_0007.fits' 

dark= fits.getdata(dir_t + dir_d)
darkavg = np.mean(dark)
darkmed = np.median(dark)

flat = fits.getdata(dir_t+dir_f)
flatmed = np.median(flat)
mediansignal = np.median(flatmed - darkmed)
flats =(flatmed-darkmed)/mediansignal
normflat = np.median(flats)

temporary = np.zeros((2048, 2048,2))
temporary[:,:,0] = fits.getdata(dir_t + dir_trans + '13549.575626_0000.fits')  #'13644.460464_0008.fits'
temporary[:,:,1] = fits.getdata(dir_t + dir_trans + '34744.580663_0160.fits')
sky = temporary[:,:,0][830:930,1270:1370]
sky_med= np.median(sky)

temporary[:,:,0] = (temporary[:,:,0]-darkavg)/normflat
temporary[:,:,1] = (temporary[:,:,1]-darkavg)/normflat

plt.figure(1)
plt.imshow(temporary[:,:,0],vmin = 0, vmax = 300, origin = 'lower', interpolation ='nearest')
plt.show()

plt.figure(2)
plt.imshow(temporary[:,:,1],vmin = 00, vmax =300, origin = 'lower', interpolation ='nearest')
plt.show()
FWHM = 3
aperture_r= 6
aperture_r1 =8
aperture_r2 = 7
thresh = 25
exposure_time = 5.0

########### before ~ 180 degree flip ##############
target_b = np.loadtxt(dir_t+'transit_b', dtype = str)

edit = 0
edit_up = 430
hour_b = np.zeros(target_b.size - edit )
minute_b = np.zeros(target_b.size - edit )
second_b = np.zeros(target_b.size - edit )
time_b = np.zeros(target_b.size - edit )

mag_b= np.zeros(target_b.size - edit )
mag_b1= np.zeros(target_b.size - edit )
mag_b2= np.zeros(target_b.size - edit )

phot_table_b = np.zeros( target_b.size - edit )
phot_table_b1 = np.zeros( target_b.size - edit )
phot_table_b2 = np.zeros( target_b.size - edit )

for i in range(0,target_b.size - edit):
    image_b, hd1 = fits.getdata(dir_t + target_b[i], header=True)
    hd=hd1['UTC']
    hour_b[i] = float(hd[9:11])
    minute_b[i] = float((hd[11:13]))/60.0
    second_b[i] = float((hd[13:22]))/3600.0
    time_b[i] = hour_b[i]+minute_b[i]+second_b[i]
    
    image_b =( image_b -darkavg)/normflat
    slice_comp1 = image_b[1237:1337,760:860]   #[795:825,1272:1302]
    slice_target  = image_b[650:750,620:720]  #[660:690,690:720]
    slice_comp2  = image_b[1315:1415,1330:1430]  #[1365:1395,1350:1380]
    
    sources_b = daofind(slice_target, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_b = zip(sources_b['xcentroid'], sources_b['ycentroid'])
    apertures_b = CircularAperture(positions_b, r= aperture_r)
    phot_table_b = aperture_photometry(slice_target, apertures_b)

    mag_b[i] =np.max( phot_table_b['aperture_sum'])/exposure_time
    
    sources_b1 = daofind(slice_comp1, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_b1 = zip(sources_b1['xcentroid'], sources_b1['ycentroid'])
    apertures_b1 = CircularAperture(positions_b1, r= aperture_r1)
    phot_table_b1 = aperture_photometry(slice_comp1, apertures_b1)

    mag_b1[i] =np.max( phot_table_b1['aperture_sum'])/exposure_time

    sources_b2 = daofind(slice_comp2, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_b2 = zip(sources_b2['xcentroid'], sources_b2['ycentroid'])
    apertures_b2 = CircularAperture(positions_b2, r= aperture_r2)
    phot_table_b2 = aperture_photometry(slice_comp2, apertures_b2)

    mag_b2[i] =np.max( phot_table_b2['aperture_sum'])/exposure_time


pickle.dump(mag_b ,open('/Volumes/TIM/Lab5/trans/mag_b.txt', 'wb'))
pickle.dump(mag_b1,open('/Volumes/TIM/Lab5/trans/mag_b1.txt', 'wb'))
pickle.dump(mag_b2,open('/Volumes/TIM/Lab5/trans/mag_b2.txt', 'wb'))
pickle.dump(time_b,open('/Volumes/TIM/Lab5/trans/time_b.txt', 'wb'))

mags_b = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b.txt','rb'))
mags_b1 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b1.txt','rb'))
mags_b2 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b2.txt','rb'))

time_b = pickle.load(open('/Volumes/TIM/Lab5/trans/time_b.txt', 'rb'))

plt.figure(1)
plt.plot( time_b, mags_b, 'bo')
plt.xlabel('UT (hours)')
plt.ylabel('Flux (Counts/Second)')
plt.xlim(1,8)
plt.show(1)


########### after ~ 180 degree flip ##############
target_a = np.loadtxt(dir_t+'transit_a', dtype = str)

edit = 0
hour_a = np.zeros(target_a.size - edit )
minute_a = np.zeros(target_a.size - edit )
second_a = np.zeros(target_a.size - edit )
time_a = np.zeros(target_a.size - edit )

phot_table_a = np.zeros( target_a.size-edit )
phot_table_a1 = np.zeros( target_a.size-edit )
phot_table_a2 = np.zeros( target_a.size-edit )
mag_a= np.zeros(target_a.size - edit )
mag_a1= np.zeros(target_a.size - edit )
mag_a2= np.zeros(target_a.size - edit )


for i in range(0,target_a.size - edit):
    image_a, hd1 = fits.getdata(dir_t + target_a[i], header=True)
    hd=hd1['UTC']
    hour_a[i] = float(hd[9:11])
    minute_a[i] = float((hd[11:13]))/60.0
    second_a[i] = float((hd[13:22]))/3600.0
    time_a[i] = hour_a[i]+minute_a[i]+second_a[i]

    
    image_a =( image_a-darkavg)/normflat
    slice_comp1 = image_a[342:442,1335:1435]   #[795:825,1272:1302]
    slice_target  = image_a[900:1000,1550:1650]  #[660:690,690:720]
    slice_comp2  = image_a[343:443,760:860]  #[1365:1395,1350:1380]
    
    sources_a = daofind(slice_target, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_a = zip(sources_a['xcentroid'], sources_a['ycentroid'])
    apertures_a = CircularAperture(positions_a, r= aperture_r)
    phot_table_a = aperture_photometry(slice_target, apertures_a)

    mag_a[i] =np.max( phot_table_a['aperture_sum'])/exposure_time
    
    sources_a1 = daofind(slice_comp1, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_a1 = zip(sources_a1['xcentroid'], sources_a1['ycentroid'])
    apertures_a1 = CircularAperture(positions_a1, r= aperture_r1)
    phot_table_a1 = aperture_photometry(slice_comp1, apertures_a1)

    mag_a1[i] =np.max( phot_table_a1['aperture_sum'])/exposure_time
    
    sources_a2 = daofind(slice_comp2, fwhm = FWHM, threshold = thresh, exclude_border=True)
    positions_a2 = zip(sources_a2['xcentroid'], sources_a2['ycentroid'])
    apertures_a2 = CircularAperture(positions_a2, r= aperture_r2)
    phot_table_a2 = aperture_photometry(slice_comp2, apertures_a2)

    mag_a2[i] =np.max( phot_table_a2['aperture_sum'])/exposure_time

pickle.dump(mag_a,open('/Volumes/TIM/Lab5/trans/mag_a.txt', 'wb'))
pickle.dump(mag_a1,open('/Volumes/TIM/Lab5/trans/mag_a1.txt', 'wb'))
pickle.dump(mag_a2,open('/Volumes/TIM/Lab5/trans/mag_a2.txt', 'wb'))
pickle.dump(time_a,open('/Volumes/TIM/Lab5/trans/time_a.txt', 'wb'))

mags_a1 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_a1.txt','rb'))
mags_a2 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_a2.txt','rb'))

mags_b1 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b1.txt','rb'))
mags_b2 = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b2.txt','rb'))
mags_b = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b.txt','rb'))
time_b = pickle.load(open('/Volumes/TIM/Lab5/trans/time_b.txt', 'rb'))

mags_a = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_a.txt','rb'))
time_a = pickle.load(open('/Volumes/TIM/Lab5/trans/time_a.txt', 'rb'))

mag_Acomp =mags_a1+ mags_a2
mag_Bcomp =mags_b1+ mags_b2

mag_Amed = mag_Acomp/np.median(mag_Acomp)
mag_Bmed = mag_Bcomp/np.median(mag_Bcomp)

mag_Adetrend =mags_a/mag_Amed
mag_Bdetrend = mags_b/mag_Bmed 
"""
plt.plot(time_b,mag_Bdetrend, 'ro')
plt.plot(time_a,mag_Adetrend, 'ro')

plt.plot(time_b,mags_b,'bo')
plt.plot( time_a,mags_a,'bo')

plt.plot(time_b,mags_b1, 'go')
plt.plot( time_a,mags_a1, 'go')
plt.plot(time_b,mags_b2,'go')
plt.plot( time_a,mags_a2,'go')
"""
plt.figure(1)
plt.plot(time_b,mag_Bdetrend, 'ro')
plt.plot(time_a,mag_Adetrend, 'ro')
plt.xlabel('UT (hours)')   
plt.ylabel('Flux (Counts/Second)')
plt.grid()
plt.xlim(1.4,7.4)
plt.ylim(18000,40000)
plt.show()

plt.figure(2)
plt.plot(time_b,mag_Bdetrend/37700, 'ro')
plt.plot(time_a,mag_Adetrend/37700, 'ro')
plt.xlabel('UT (hours)')   
plt.ylabel('Flux (% drop)')
plt.grid()
plt.xlim(1.4,7.4)
plt.ylim(0.45,1.05)
plt.show(0)


min_1 = np.min(time_a) # 5.871743425
min_2 = np.min(time_b) # 1.810587467

period = 2*(5.871743425- 1.810587467)



# code below will build a model of an eclipsing binary system
# PHysics Of Eclipsing BinariEs (PHOEBE)

import phoebe

# read in the inscrutable "active" file required by the code
pset,lcset,rvset = phoebe.wd.lcin_to_ps('/Volumes/TIM/Lab5/test01lcin.active',version='wd2003')
lcset['phinc'] = 0.01

# set some model parameters:
period = 0.331893 # days   (this is the period of AB And according to
#7.965432 hours

pset['teff1'] = 0.5670      # effective temperature of star 1 (multiply by 10000 K)
pset['teff2'] = 0.5370      # effective temperature of star 2 (multiply by 10000 K)
pset['period'] = period    # period (days)
pset['incl'] = 87.3      # orbital inclination
pset['ecc'] = 0.0        # orbital eccentricity
pset['rm'] = 0.75       # mass ratio, secondary/primary

# to print all parameters used in the model:
print(pset)

# build the model
curve,params = phoebe.wd.lc(pset,request='curve',light_curve=lcset)
image,params = phoebe.wd.lc(pset,request='image',light_curve=lcset)

# plot it up
p = plt.figure(5)
p = plt.subplot(211)
p = plt.plot(curve['indeps'],curve['lc'],'r-')
p = plt.subplot(212,aspect='equal')
p = plt.plot(image['y'],image['z'],'ko',ms=1)

#------------------------------#

dark= fits.getdata(dir_t + dir_d)
darkavg = np.mean(dark)
darkmed = np.median(dark)

flat = fits.getdata(dir_t+dir_f)
flatmed = np.median(flat)
mediansignal = np.median(flatmed - darkmed)
flats =(flatmed-darkmed)/mediansignal
normflat = np.median(flats)

temporary = np.zeros((2048, 2048,2))

temporary[:,:,0],hd1 = fits.getdata(dir_t + dir_trans + '13549.575626_0000.fits', header =True)  #'13644.460464_0008.fits'
temporary[:,:,1],hd2 = fits.getdata(dir_t + dir_trans + '34744.580663_0160.fits',header = True)
sky = temporary[:,:,0][830:930,1270:1370]
sky_med= np.median(sky)

temporary[:,:,0] = (temporary[:,:,0]-darkavg)/normflat
temporary[:,:,1] = (temporary[:,:,1]-darkavg)/normflat

fits.writeto('/Volumes/TIM/Lab5/trans/image_b.fits',temporary[:,:,0], header =hd1)
fits.writeto('/Volumes/TIM/Lab5/trans/image_a.fits',temporary[:,:,1], header =hd2)

import aplpy


astro1= fits.getdata('/Volumes/TIM/Lab5/trans/astro/image_a1.fits')
astro = aplpy.FITSFigure('/Volumes/TIM/Lab5/trans/astro/image_a1.fits')
fig = astro
plt.imshow(astro1, vmin=0, vmax = 500,cmap ='gray')
plt.show()


plt.imshow(astro1, vmin=0, vmax = 500,cmap ='gray')
plt.show()

fig.add_label(0.85,0.47, 'Target Star',color = 'red', relative = True) # target star
fig.add_label(0.785,0.2, 'Comparison Star 1',color = 'red', relative =True) 
fig.add_label(0.50, 0.205, 'Comparison Star 2' ,color = 'red', relative = True)
fig.add_grid()
fig.grid.set_linestyle('dotted')

astro2= fits.getdata('/Volumes/TIM/Lab5/trans/astro/image_b1.fits')
astro = aplpy.FITSFigure('/Volumes/TIM/Lab5/trans/astro/image_b1.fits')
fig = astro
plt.figure(2)
plt.imshow(astro2, vmin=0, vmax = 600,cmap ='gray')
plt.show()

fig.add_label(0.39,.35, 'Target Star',color = 'red', relative = True) # target star
fig.add_label(0.50,0.64, 'Comparison Star 1',color = 'red', relative =True) 
fig.add_label(0.78,0.68, 'Comparison Star 2' ,color = 'red',  relative = True)
fig.add_grid()
fig.grid.set_linestyle('dotted')
