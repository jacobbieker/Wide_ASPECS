from astropy.io import fits
import numpy as np
from catalog_builder import build_catalog
from astropy.table import Table

hdu_list = fits.open("data/magphys_in.fits")
#print(hdu_list.info())
#print(hdu_list[1].header)
#print(hdu_list[1].columns.info())
#print(hdu_list[1].data)

data_field = hdu_list[1].data

list_of_used_ones = []
for element in data_field:
    list_of_used_ones.append(element[0])

full_catalog = build_catalog()

filter_fluxes = []
filter_errors = []
filter_refnrs = []
filter_wavelengths = []
filter_bandwidths = []

#filter
filter_fluxes.append('fU38')
filter_errors.append('efU38')
filter_refnrs.append(249)
filter_wavelengths.append(0.3676)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fU')
filter_errors.append('efU')
filter_refnrs.append(222)
filter_wavelengths.append(0.36402)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fR')
filter_errors.append('efR')
filter_refnrs.append(225)
filter_wavelengths.append(0.6442)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fB')
filter_errors.append('efB')
filter_refnrs.append(223)
filter_wavelengths.append(0.43363)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fV')
filter_errors.append('efV')
filter_refnrs.append(224)
filter_wavelengths.append(0.5354)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fI')
filter_errors.append('efI')
filter_refnrs.append(226)
filter_wavelengths.append(0.7936)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fJ')
filter_errors.append('efJ')
filter_refnrs.append(450)
filter_wavelengths.append(1.2369)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fH')
filter_errors.append('efH')
filter_refnrs.append(451)
filter_wavelengths.append(1.6464)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fK')
filter_errors.append('efK')
filter_refnrs.append(452)
filter_wavelengths.append(2.2109)
filter_bandwidths.append(0.0)

#filter
filter_fluxes.append('fKs')
filter_errors.append('efKs')
filter_refnrs.append(108)
filter_wavelengths.append(2.13708)
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

print(full_catalog.columns)
size_f = 0
print("Missing from Filter set")
for name in full_catalog.columns:
    if "flag_" not in name and "f" in name and name not in filter_fluxes:
        print(name)
        size_f += 1
print(size_f)

print("In Filter Set but not in Full Catalog")
size_f = 0
for name in filter_fluxes:
    if  name not in full_catalog.columns:
        print(name)
        size_f += 1
print(size_f)

print("In Filter Set")
size_f = 0
for name in filter_fluxes:
    if "ef" not in name and "flag_" not in name and "f" in name and name in full_catalog.columns \
            and "sfr" not in name:
        print(name)
        size_f += 1
print(size_f)

print("Length of Filter Set: " + str(len(filter_fluxes)))
for name in filter_fluxes:
    print(name)

# TODO: B,V,I,J,H,1.1mm,tJ,Rc,R,U,tK,K