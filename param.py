#!/usr/bin/python
import math
import numpy as np
import string
import urllib
import urllib.request
import os
import argparse
import pandas as pd


def _output(header=False, parameters=None, switch='on'):
    """Create the output file 'stellar_characterization.dat''
    Input
    -----
    header    - Only use True if this is for the file to be created
    """

    if header and (switch == 'on'):
        hdr = ['star_name star_teff star_erteff star_logg star_erlogg star_metal star_ermetal star_mass_padova star_ermass_padova star_radius_padova star_erradius_padova star_age star_erage star_logg_padova star_erlogg_padova logg_hip erlogg_hip Rt dRt Mcal dMcal']
        if not os.path.isfile('stellar_characterization.dat'):
            with open('stellar_characterization.dat', 'w') as output:
                output.write('\t'.join(hdr)+'\n')

    elif header and (switch == 'off'):
        hdr = 'star_name star_teff star_erteff star_logg star_erlogg star_metal star_ermetal Rt dRt Mcal dMcal'
        if not os.path.isfile('stellar_characterization.dat'):
            with open('stellar_characterization.dat', 'w') as output:
                output.write('\t'.join(hdr)+'\n')

    else:
        with open('stellar_characterization.dat', 'a') as output:
            output.write('\t'.join(map(str, parameters))+'\n')


def bcflow(teff):
    """Table 1 from Torres 2010: Bolometric Corrections by Flower (1996) as a Function of Temperature:
    BC_V = a + b(log T_eff) + c(log T_eff)^2 + sdot sdot sdot
    Coefficient	log T_eff < 3.70	3.70 < log T_eff < 3.90	   log T_eff>3.90
    a	       -0.190537291496456E+05	-0.370510203809015E+05	-0.118115450538963E+06
    b	        0.155144866764412E+05	 0.385672629965804E+05	 0.137145973583929E+06
    c	       -0.421278819301717E+04	-0.150651486316025E+05	-0.636233812100225E+05
    d	        0.381476328422343E+03	 0.261724637119416E+04	 0.147412923562646E+05
    e	                    ... 	    -0.170623810323864E+03	-0.170587278406872E+04
    f	                    ... 	            ... 	         0.788731721804990E+02
    """

    a = [-0.190537291496456e+05, -0.370510203809015e+05, -0.118115450538963e+06]
    b = [0.155144866764412e+05,   0.385672629965804e+05,  0.137145973583929e+06]
    c = [-0.421278819301717e+04, -0.150651486316025e+05, -0.636233812100225e+05]
    d = [0.381476328422343e+03,   0.261724637119416e+04,  0.147412923562646e+05]
    e = [-0.170623810323864e+03, -0.170587278406872e+04]
    f = [0.788731721804990e+02]

    lteff= np.log10(teff)
    if lteff < 3.7:
        bc = a[0] + (b[0]*lteff) + (c[0]*(lteff**2)) + (d[0]*(lteff**3))
    elif (lteff >= 3.7) and (lteff < 3.9):
        bc = a[1] + (b[1]*lteff) + (c[1]*(lteff**2)) + (d[1]*(lteff**3)) + (e[0]*(lteff**4))
    elif lteff >= 3.9:
        bc = a[2] + (b[2]*lteff) + (c[2]*(lteff**2)) + (d[2]*(lteff**3)) + (e[1]*(lteff**4)) + (f[0]*(lteff)**5)
    return bc


def logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass):
    """Calculate the trigonometric logg and error"""
    #np.geterr()
    if mass == 'nan':
        logg, dlogg = 'nan', 'nan'

    else:
        e = 2.718281828
        logg  = 4.44 + np.log10(mass) + (4.0*np.log10(teff/5777.)) + (0.4*(v + bc)) + (2.0*np.log10(par/1000.0)) + 0.108
        logg  = np.round(logg, 2)
        dlogg = np.sqrt(((dmass*np.log10(e))/mass)**2 + ((4.*dteff*np.log10(e))/teff)**2 + ((2.*0.05*np.log10(e))/par)**2)
        dlogg = np.round(dlogg, 2)
    return logg, dlogg


def get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, temp, er_temp, metal, er_metal, minage, maxage):
    """Returns mass, radius, and age from padova interface.
       Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
       Parallax is in mas.
       Be careful with the constrains in age!
    """

    url = 'http://stev.oapd.inaf.it/cgi-bin/param_1.3'
    #These are the parameters in the webpage to tune
    form_data = {'param_version': '1.3', 'star_name': star, 'star_teff': int(temp), 'star_sigteff': int(er_temp), 'star_feh': round(metal,2), 'star_sigfeh': round(er_metal,2),  'star_vmag': round(vmag,3), 'star_sigvmag': round(er_vmag,3), 
                 'star_parallax': round(parallax,3), 'star_sigparallax': round(er_parallax,3),  'isoc_kind': 'parsec_CAF09_v1.1', 'kind_interp': '1', 'kind_tpagb': '0', 'kind_pulsecycle': '0', 'kind_postagb': '-1', 
                 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat', 'sfr_file': 'tab_sfr/sfr_const_z008.dat', 'sfr_minage': minage, 'sfr_maxage': maxage, 'flag_sismo': '0', 'photsys_file': 'tab_mag_odfnew/tab_mag_ubvrijhk.dat', 
                 'submit_form': 'Submit' }
    print('Parameters for %s from Padova interface.' % star)

    for key, value in form_data.items():
        print('%s -> %s' % (key, value))
    urllib.request.urlretrieve(url, 'parameters.html', lambda x,y,z:0, urllib.parse.urlencode(form_data).encode("utf-8"))

    #write results
    with open('parameters.html') as f:
        
        line = f.readlines()[15]

    line = line.replace(' .<p>Results for ', '')
    line = line.replace('Mass=', '')
    line = line.replace('&plusmn;', ' , ')
    line = line.replace(' ','')
    line = line.replace('<i>R</i>=', '')
    line = line.replace(':Age=', ',')
    line = line.replace('<i>M</i>&#9737', '')
    line = line.replace('<i>R</i>&#9737', '')
    line = line.replace('log<i>g</i>=', '')
    line = line.replace('(cgs)', '')
    line = line.replace('Gyr', '')
    line = line.split(',')

    age =      float(line[1])
    erage =    float(line[2])
    mass =     float(line[3])
    ermass =   float(line[4])
    logg_p =   float(line[5])
    erlogg_p = float(line[6])
    radius =   float(line[7])
    erradius = float(line[8])
    print('---------------------%s done!----------------------------' % star)
    return mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p


def mass_torres(hd, T, logg, fe, dT, dlogg, dfe):
    """Calculate masses and radii from the calibration of Torres et al. 2010"""

    # coefficients
    a  = [1.5689, 1.3787, 0.4243, 1.139, -0.1425,  0.01969, 0.1010]
    da = [0.058,  0.029,  0.029,  0.240, 0.011,    0.0019,  0.014]
    b  = [2.4427, 0.6679, 0.1771, 0.705, -0.21415, 0.02306, 0.04173]
    db = [0.038,  0.016,  0.027,  0.13,  0.0075,   0.0013,  0.0082]

    X = np.log10(T) - 4.1
    dX = dT/T
    log_M = a[0] + (a[1]*X) + (a[2]*(X**2)) + (a[3]*(X**3)) + (a[4]*(logg**2)) + (a[5]*(logg**3)) + (a[6]*fe)
    log_R = b[0] + (b[1]*X) + (b[2]*(X**2)) + (b[3]*(X**3)) + (b[4]*(logg**2)) + (b[5]*(logg**3)) + (b[6]*fe)
    # must check the errors
    dlog_M = np.sqrt((da[0]**2) + ((a[1]*dX)**2) + ((da[1]*X)**2) + ((da[2]*(X**2))**2) + ((a[2]*2*X*dX)**2) + ((da[3]*(X**3))**2) + ((a[3]*3*X*X*dX)**2) + ((da[4]*(logg**2))**2) + ((a[4]*2*logg*dlogg)**2) + ((da[5]*(logg**3))**2) + ((a[5]*3*logg*logg*dlogg)**2) + ((da[6]*fe)**2) + ((a[6]*dfe)**2))
    dlog_R = np.sqrt((db[0]**2) + ((b[1]*dX)**2) + ((db[1]*X)**2) + ((db[2]*(X**2))**2) + ((b[2]*2*X*dX)**2) + ((db[3]*(X**3))**2) + ((b[3]*3*X*X*dX)**2) + ((db[4]*(logg**2))**2) + ((b[4]*2*logg*dlogg)**2) + ((db[5]*(logg**3))**2) + ((b[5]*3*logg*logg*dlogg)**2) + ((db[6]*fe)**2) + ((b[6]*dfe)**2))
    Mt = np.power(10,log_M)
    Rt = np.power(10,log_R)
    dMt = dlog_M*Mt*np.power(10,(log_M-1.0))
    dRt = dlog_R*Rt*np.power(10,(log_R-1.0))
    # Apply Santos et al. (2013) correction
    Mcal = (0.791*(Mt**2.0)) - (0.575*Mt) + 0.701
    dMcal = np.sqrt(((0.791*Mt*dMt)**2) + ((0.575*dMt)**2))

    Mcal  = np.round(Mcal,2)
    dMcal = np.round(dMcal,2)
    Rt    = np.round(Rt, 2)
    dRt   = np.round(dRt,2)
    return hd, Rt, dRt, Mcal, dMcal

if __name__ == '__main__':

    # params = pd.read_csv('star', delimiter='\t', comment='#', header=1, usecols=(range(11)), names = ['star', 'teff', 'erteff', 'logg', 'erlogg', 'feh', 'erfeh', 'vmag', 'ervmag', 'par', 'erpar'])
    params = pd.read_csv('star', delimiter='\t')
    print(params)

    for i, x in enumerate(params.star[:]):
            #Padova mass, radius, age, and logg from isochrones
            mass_padova, ermass_padova, radius_padova, erradius_padova, age, erage, logg_p, erlogg_p = get_mass_radius(params.star[i], params.vmag[i],  params.ervmag[i], params.par[i], params.erpar[i], params.teff[i], params.erteff[i], params.feh[i], params.erfeh[i], params.sfr_minage[i], params.sfr_maxage[i])
            #Trigonometric logg. Uses the get_mass_radius function
            bc = bcflow(params.teff[i])
            logg_hip, erlogg_hip = logg_trigomonetric(params.teff[i], mass_padova, params.vmag[i], bc, params.par[i], params.erpar[i], params.erteff[i], ermass_padova)
            #Mass and radius from Torres calibration
            name, Rt, dRt, Mcal, dMcal = mass_torres(params.star[i], params.teff[i], params.logg[i], params.feh[i], params.erteff[i], params.erlogg[i], params.erfeh[i])

            #Save results
            parameters = [name, params.teff[i], params.erteff[i], params.logg[i], params.erlogg[i], params.feh[i], params.erfeh[i], mass_padova, ermass_padova, radius_padova, erradius_padova, age, erage, logg_p, erlogg_p, logg_hip, erlogg_hip, Rt, dRt, Mcal, dMcal]
            _output(parameters=parameters)
            os.remove('parameters.html')

 
