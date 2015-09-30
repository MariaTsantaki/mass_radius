#!/usr/bin/python
import math
import numpy as np
import string
import urllib
import os, sys
import argparse

parser = argparse.ArgumentParser(description='set flag if there is information on the parallax and V mag.')
parser.add_argument('name', type=str, help ='File with the data')
parser.add_argument('plx_switch', type=str, help ='Parallaxes and Vmags are included by default')
args = parser.parse_args()

def bcflow(teff):
    """Table 1 from Torres 2010: Bolometric Corrections by Flower (1996) as a Function of Temperature: BC_V = a + b(log T_eff) + c(log T_eff)^2 + sdot sdot sdot

    Coefficient	log T_eff < 3.70	3.70 < log T_eff < 3.90	   log T_eff>3.90	
    a	       -0.190537291496456E+05	-0.370510203809015E+05	-0.118115450538963E+06	
    b	        0.155144866764412E+05	 0.385672629965804E+05	 0.137145973583929E+06	
    c	       -0.421278819301717E+04	-0.150651486316025E+05	-0.636233812100225E+05	
    d	        0.381476328422343E+03	 0.261724637119416E+04	 0.147412923562646E+05	
    e	                    ... 	-0.170623810323864E+03	-0.170587278406872E+04	
    f	                    ... 	            ... 	 0.788731721804990E+02	
    """

    lteff= np.log10(teff)
    new_bc = []
    for teff_value in lteff:
        if teff_value<3.7:
	    bc=-0.190537291496456e+05+(0.155144866764412e+05*teff_value)-(0.421278819301717e+04*(teff_value**2))+(0.381476328422343e+03*(teff_value**3))
            new_bc.append(bc)
	elif teff_value>=3.7 and teff_value<3.9:
	    bc=-0.370510203809015e+05+(0.385672629965804e+05*teff_value)-(0.150651486316025e+05*(teff_value**2))+(0.261724637119416e+04*(teff_value**3))-(0.170623810323864e+03*(teff_value**4)) 
            new_bc.append(bc)
	elif teff_value>=3.9:
	    bc=-0.118115450538963e+06+(0.137145973583929e+06*teff_value)-(0.636233812100225e+05*(teff_value**2))+(0.147412923562646e+05*(teff_value**3))-(0.170587278406872E+04*(teff_value**4))+(0.788731721804990e+02*(teff_value)**5)
	    new_bc.append(bc)
    return new_bc

def logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass):
    """Calculate the trigonometric logg and error"""
    #np.geterr()
    e = 2.718281828
    logg = 4.44 + np.log10(mass) + 4.*np.log10(teff/5777.) + 0.4*(v+bc) + 2.*np.log10(par/1000.)+0.108
    logg = np.round(logg,2)
    dlogg = np.sqrt(((dmass*np.log10(e))/mass)**2 + ((4.*dteff*np.log10(e))/teff)**2 + ((2.*0.05*np.log10(e))/par)**2)
    dlogg = np.round(dlogg,2)
    return logg, dlogg

def get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, temp, er_temp, metal, er_metal):
    """Returns mass, radius, and age from padova interface.
       Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
       Be careful with the constrains in age!
    """
    url = 'http://stev.oapd.inaf.it/cgi-bin/param'
    #These are the parameters in the webpage to tune
    form_data = {'param_version': '1.3', 'star_name': star, 'star_teff': temp, 'star_sigteff': er_temp, 'star_feh': metal, 'star_sigfeh': er_metal,  'star_vmag': vmag, 'star_sigvmag': er_vmag, 'star_parallax': parallax, 'star_sigparallax': er_parallax,  'isoc_kind': 'parsec_CAF09_v1.1', 'kind_interp': '1', 'kind_tpagb': '0', 'kind_pulsecycle': '0', 'kind_postagb': '-1', 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat', 'sfr_file': 'tab_sfr/sfr_const_z008.dat', 'sfr_minage': '0.1e9', 'sfr_maxage': '12.0e9', 'flag_sismo': '0', 'submit_form': 'Submit' }
    print form_data
    urllib.urlretrieve(url, 'parameters.html', lambda x,y,z:0, urllib.urlencode(form_data))

    #write results
    with open('parameters.html') as f:
        line = f.readlines()[19]

    line = line.replace(' .<p>Results for ', '')
    line = line.replace('Mass=', '')
    line = line.replace('&plusmn;', ' , ')
    line = line.replace(' ','')
    line = line.replace('<i>R</i>=','')
    line = line.replace(':Age=',',')
    line = line.replace('<i>M</i>&#9737','')
    line = line.replace('<i>R</i>&#9737','')
    line = line.replace('log<i>g</i>=','')
    line = line.replace('(cgs)','')
    line = line.replace('Gyr','')
    line = line.split(',')

    name = str(line[0])
    age = float(line[1])
    erage = float(line[2])
    mass = float(line[3])
    ermass = float(line[4])
    logg_p = float(line[5])
    erlogg_p = float(line[6])
    radius = float(line[7])
    erradius = float(line[8])
    print name, ' done'

    return mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p


def mass_torres(hd, T, logg, fe, dT, dlogg, dfe):
    """Calculate masses and radii from the calibration of Torres et al. 2010"""

    #coefficients   
    a1 = 1.5689; a2 = 1.3787; a3 = 0.4243; a4 = 1.139; a5 = -0.1425 ; a6 = 0.01969 ; a7 = 0.1010; b1 = 2.4427; b2 = 0.6679; b3 =  0.1771; b4 = 0.705; b5 = -0.21415; b6 = 0.02306; b7 = 0.04173; da1 = 0.058; da2 = 0.029; da3 = 0.029; da4 = 0.240; da5 = 0.011; da6 = 0.0019; da7 = 0.014; db1 = 0.038; db2 = 0.016; db3 = 0.027; db4 = 0.13; db5 = 0.0075; db6 = 0.0013; db7 = 0.0082
    X = np.log10(T)-4.1
    dX = dT/T
    log_M = a1 + a2*X + a3*(X**2) + a4*(X**3) + a5*(logg**2) + a6*(logg**3) +a7*fe
    #must check the errors
    dlog_M = np.sqrt((da1**2) + ((a2*dX)**2) + ((da2*X)**2) + ((da3*(X**2))**2) + ((a3*2*X*dX)**2) + ((da4*(X**3))**2) + ((a4*3*X*X*dX)**2) + ((da5*(logg**2))**2) + ((a5*2*logg*dlogg)**2) + ((da6*(logg**3))**2) + ((a6*3*logg*logg*dlogg)**2) + ((da7*fe)**2) + ((a7*dfe)**2))
    log_R = b1 + b2*X + b3*(X**2) + b4*(X**3) + b5*(logg**2) + b6*(logg**3) + b7*fe
    dlog_R = np.sqrt((db1**2) + ((b2*dX)**2) + ((db2*X)**2) + ((db3*(X**2))**2) + ((b3*2*X*dX)**2) + ((db4*(X**3))**2) + ((b4*3*X*X*dX)**2) + ((db5*(logg**2))**2) + ((b5*2*logg*dlogg)**2) + ((db6*(logg**3))**2) + ((b6*3*logg*logg*dlogg)**2) + ((db7*fe)**2) + ((b7*dfe)**2))
    Mt = np.power(10,log_M)
    Rt = np.power(10,log_R)
    dMt = dlog_M*Mt*np.power(10,(log_M-1.0))
    dRt = dlog_R*Rt*np.power(10,(log_R-1.0))
    Mcal = 0.791*Mt*Mt - 0.575*Mt +0.701
    dMcal = np.sqrt(((0.791*Mt*dMt)**2) + ((0.575*dMt)**2))

    Mcal = np.round(Mcal,2)
    dMcal = np.round(dMcal,2)
    Rt = np.round(Rt, 2)
    dRt = np.round(dRt,2)

    return hd, Rt, dRt, Mcal, dMcal

#main program
if args.plx_switch == 'on': 
    print 'Parallax info'

    star = np.genfromtxt(args.name, dtype=None, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9,10), names = ['star', 'teff', 'erteff', 'logg', 'erlogg', 'feh', 'erfeh', 'V', 'eV', 'Plx', 'e_Plx'])
    star_name = star['star']
    star_vmag = star['V']
    star_ervmag = star['eV']
    star_par = star['Plx']
    star_erpar = star['e_Plx']
    star_teff = star['teff'] 
    star_erteff = star['erteff']
    star_logg = star['logg']
    star_erlogg = star['erlogg']
    star_metal = star['feh']
    star_ermetal = star['erfeh']

    #Calculate parameters
    #Padova mass, radius, age, and logg from isochrones
    padova_params = []

    for i, x in enumerate(star_name[:]):
        padova_params.append(get_mass_radius(star_name[i], star_vmag[i], star_par[i], star_erpar[i], star_teff[i], star_erteff[i], star_metal[i], star_ermetal[i]))

    star_mass_padova = np.column_stack(padova_params)[0]
    star_ermass_padova = np.column_stack(padova_params)[1]
    star_radius_padova = np.column_stack(padova_params)[2]
    star_erradius_padova = np.column_stack(padova_params)[3]
    star_age = np.column_stack(padova_params)[4]
    star_erage = np.column_stack(padova_params)[5]
    star_logg_padova = np.column_stack(padova_params)[6]
    star_erlogg_padova = np.column_stack(padova_params)[7]

    #Trigonometric logg. Uses the get_mass_radius function
    bc = bcflow(star_teff)
    logg_hip, erlogg_hip = logg_trigomonetric(star_teff, star_mass_padova, star_vmag, bc, star_par, star_erpar, star_erteff, star_ermass_padova)

    #Mass and radius from Torres calibration
    name, Rt, dRt, Mcal, dMcal = mass_torres(star_name, star_teff, star_logg, star_metal, star_erteff, star_erlogg, star_ermetal)
    header = 'star_name, star_teff, star_erteff, star_logg, star_erlogg, star_metal, star_ermetal, star_mass_padova, star_ermass_padova, star_radius_padova, star_erradius_padova, star_age, star_erage, star_logg_padova, star_erlogg_padova, logg_hip, erlogg_hip, Rt, dRt, Mcal, dMcal'
    data = np.column_stack((star_name, star_teff, star_erteff, star_logg, star_erlogg, star_metal, star_ermetal, star_mass_padova, star_ermass_padova, star_radius_padova, star_erradius_padova, star_age, star_erage, star_logg_padova, star_erlogg_padova, logg_hip, erlogg_hip, Rt, dRt, Mcal, dMcal))
    np.savetxt('stellar_characterization.dat', data, delimiter='\t', fmt="%s", header=header) 

elif args.plx_switch == 'off': 
    print 'Mass and Radius from Torres calibration'
    star = np.genfromtxt(args.name, dtype=None, skiprows=2, usecols=(0,1,2,3,4,5,6), names = ['star', 'teff', 'erteff', 'logg', 'erlogg', 'feh', 'erfeh'])

    star_name = star['star']
    star_teff = star['teff'] 
    star_erteff = star['erteff']
    star_logg = star['logg']
    star_erlogg = star['erlogg']
    star_metal = star['feh']
    star_ermetal = star['erfeh']

    #Calculate parameters
    #Mass and radius from Torres calibration
    name, Rt, dRt, Mcal, dMcal = mass_torres(star_name, star_teff, star_logg, star_metal, star_erteff, star_erlogg, star_ermetal)
    header='star_name, star_teff, star_erteff, star_logg, star_erlogg, star_metal, star_ermetal, Rt, dRt, Mcal, dMcal'
    data = np.column_stack((star_name, star_teff, star_erteff, star_logg, star_erlogg, star_metal, star_ermetal, Rt, dRt, Mcal, dMcal))
    np.savetxt('stellar_characterization.dat', data, delimiter='\t', fmt="%s", header=header) 

else:
    print 'Type "on" in there are parallaxes, "off" in not'
