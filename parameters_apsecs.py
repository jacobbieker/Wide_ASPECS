'''User inputs'''
#files & directories
input_catalog = '/net/amstelmeer/data3/MAGPHYS/example_catalog.fits'

#cosmology
H0 = 70.0 #km/s/mpc
omega_matter = 0.3
omega_lambda = 0.7

#catalog stuff
redshift = 'z' #columnname for column containing source redshifts
rows_to_fit = range(0,63499) #iterable containing the row numbers of the sources to fit (e.g. range(0,99999999))

#catalog photometry
flux_conversion = 1.0e-6 #all fluxes are multiplied by this number (e.g. to convert from uJy in the catalog to Jy used in the fitting by multiplying with 1.0e-6)
filter_fluxes = [] #columnname for column containing flux in jansky
filter_errors = [] #columnname for column containing error in jansky
filter_refnrs = [] #reference number in the magphys filters.log file
filter_wavelengths = [] #filter wavelength in micrometers
filter_bandwidths = [] #filter width in micrometers (optinal; used for plotting)
#filter
filter_fluxes.append('fU38')
filter_errors.append('efU38')
filter_refnrs.append(249)
filter_wavelengths.append(0.3676)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fU')
filter_errors.append('efU')
filter_refnrs.append(363)
filter_wavelengths.append(0.37651)
filter_bandwidths.append(0.0)

#filter
# VIMOS R, since Rc is from WFI
filter_fluxes.append('fR')
filter_errors.append('efR')
filter_refnrs.append(364)
filter_wavelengths.append(0.63755)
filter_bandwidths.append(0.0)

#filter
# For sure from WFI, no Rc filter, so R filter from Hildebrandt et al. (2006)
filter_fluxes.append('fRc')
filter_errors.append('efRc')
filter_refnrs.append(252)
filter_wavelengths.append(0.64284)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fB')
filter_errors.append('efB')
filter_refnrs.append(250)
filter_wavelengths.append(0.44562)
filter_bandwidths.append(0.0)

#filter
# Also WFI
filter_fluxes.append('fV')
filter_errors.append('efV')
filter_refnrs.append(251)
filter_wavelengths.append(0.53401)
filter_bandwidths.append(0.0)

#filter
# Also WFI
filter_fluxes.append('fI')
filter_errors.append('efI')
filter_refnrs.append(253)
filter_wavelengths.append(0.81454)
filter_bandwidths.append(0.0)

#filter
# ISAAC J band b/c in GOOD-S
filter_fluxes.append('fJ')
filter_errors.append('efJ')
filter_refnrs.append(211)
filter_wavelengths.append(1.24279)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fH')
filter_errors.append('efH')
filter_refnrs.append(212)
filter_wavelengths.append(1.63661)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fKs')
filter_errors.append('efKs')
filter_refnrs.append(213)
filter_wavelengths.append(2.15217)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA427')
filter_errors.append('efIA427')
filter_refnrs.append(265)
filter_wavelengths.append(0.4263)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA505')
filter_errors.append('efIA505')
filter_refnrs.append(268)
filter_wavelengths.append(0.5063)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA527')
filter_errors.append('efIA527')
filter_refnrs.append(269)
filter_wavelengths.append(0.5261)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA574')
filter_errors.append('efIA574')
filter_refnrs.append(270)
filter_wavelengths.append(0.5764)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA624')
filter_errors.append('efIA624')
filter_refnrs.append(271)
filter_wavelengths.append(0.6233)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA679')
filter_errors.append('efIA679')
filter_refnrs.append(272)
filter_wavelengths.append(0.6781)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA738')
filter_errors.append('efIA738')
filter_refnrs.append(274)
filter_wavelengths.append(0.7361)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIA767')
filter_errors.append('efIA767')
filter_refnrs.append(275)
filter_wavelengths.append(0.7684)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f125w')
filter_errors.append('ef125w')
filter_refnrs.append(328)
filter_wavelengths.append(1.2486)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f140w')
filter_errors.append('ef140w')
filter_refnrs.append(329)
filter_wavelengths.append(1.3923)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f606w')
filter_errors.append('ef606w')
filter_refnrs.append(94)
filter_wavelengths.append(0.5959)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f814w')
filter_errors.append('ef814w')
filter_refnrs.append(95)
filter_wavelengths.append(0.818644)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f435w')
filter_errors.append('ef435w')
filter_refnrs.append(214)
filter_wavelengths.append(0.4328)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f775w')
filter_errors.append('ef775w')
filter_refnrs.append(216)
filter_wavelengths.append(0.7705)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f850lp')
filter_errors.append('ef850lp')
filter_refnrs.append(217)
filter_wavelengths.append(0.9048)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f160w')
filter_errors.append('ef160w')
filter_refnrs.append(298)
filter_wavelengths.append(1.5419)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIRAC1')
filter_errors.append('efIRAC1')
filter_refnrs.append(153)
filter_wavelengths.append(3.5634)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIRAC2')
filter_errors.append('efIRAC2')
filter_refnrs.append(154)
filter_wavelengths.append(4.5110)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIRAC3')
filter_errors.append('efIRAC3')
filter_refnrs.append(155)
filter_wavelengths.append(5.7593)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fIRAC4')
filter_errors.append('efIRAC4')
filter_refnrs.append(156)
filter_wavelengths.append(7.9594)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f24um')
filter_errors.append('ef24um')
filter_refnrs.append(157)
filter_wavelengths.append(24.0)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fnu_1mm')
filter_errors.append('efnu_1mm')
filter_refnrs.append(350) #possibly incorrect
filter_wavelengths.append(1000)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('f125w')
filter_errors.append('ef125w')
filter_refnrs.append(328)
filter_wavelengths.append(1.23049)
filter_bandwidths.append(0.0)

#filter
# This is Wircam
filter_fluxes.append('ftJ')
filter_errors.append('eftJ')
filter_refnrs.append(277)
filter_wavelengths.append(1.24815)
filter_bandwidths.append(0.0)

#filter
# This is Wircam, TENIS Ks band
filter_fluxes.append('ftK')
filter_errors.append('eftK')
filter_refnrs.append(279)
filter_wavelengths.append(2.13378)
filter_bandwidths.append(0.0)

##filter
#filter_fluxes.append('')
#filter_errors.append('')
#filter_refnrs.append()
#filter_wavelengths.append()
#filter_bandwidths.append()

#magphys
deltaz = 0.005 #MUST BE MULTIPLE OF 0.0025!!!!!!!!!!!   0.005 is good enough for most applications
overwrite_colors = True
overwrite_fits = True
overwrite_catalog = True
first_time_magphys = False #change this only if you are running a fresh download of magphys for the very first time
save_catalog_every_n_sources = 10

'''End user inputs'''


'''Default Settins'''
import os
NAME = os.path.basename(__file__) #not necessary to change
NAME = NAME[11:-3]
savedir = '/net/amstelmeer/data3/MAGPHYS/'+'%s/'%(NAME) #can be left like this
magphysdir = '/net/amstelmeer/data3/MAGPHYS/magphys_highz/' #does not have to be changed

'''Fixed inputs'''
colorsdir = savedir+'colors/'
plotsdir = savedir+'plots/'
datadir = savedir+'data/'

'''Sorting filters'''
filter_wavelengths,filter_fluxes,filter_errors,filter_refnrs,filter_bandwidths = (list(t) for t in zip(*sorted(zip(filter_wavelengths,filter_fluxes,filter_errors,filter_refnrs,filter_bandwidths)))) #sorting by effective wavelenth is required











