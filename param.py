#!/usr/bin/python
import numpy as np
import urllib
import urllib.request
import ssl
import os
import time
import pandas as pd
import re
from html import unescape

def _output(header=False, parameters=None):
    """Create the output file 'stellar_characterization.dat''
    Input
    -----
    header    - Only use True if this is for the file to be created
    """

    hdr = ['name,teff,erteff,metal,ermetal,mass_iso,ermass_iso,radius_iso,erradius_iso,age_iso,erage_iso,logg_iso,erlogg_iso,radius_cal,erradius_cal,mass_cal,ermass_cal,logg_trig,erlogg_trig']
    if not os.path.isfile('stellar_characterization.dat'):
            with open('stellar_characterization.dat', 'w') as output:
                output.write(','.join(hdr)+'\n')
    parameters = [item for sublist in parameters for item in (sublist if isinstance(sublist, tuple) else (sublist,))]         
    with open('stellar_characterization.dat', 'a') as output:
        output.write(','.join(map(str, parameters))+'\n')

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

def mass_torres(T, logg, fe, dT, dlogg, dfe):
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
    torres =  Rt, dRt, Mcal, dMcal
    return torres

def get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, temp, er_temp, metal, er_metal, minage, maxage):
    """Returns mass, radius, and age from padova interface.
       Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
       Parallax is in mas.
       Be careful with the constrains in age!
    """
    if os.path.exists('parameters.html'):
        os.remove('parameters.html')

    url = 'https://stev.oapd.inaf.it/cgi-bin/param_1.3'
    form_data = {
        'param_version': '1.3', 'star_name': star, 
        'star_teff': str(int(temp)), 
        'star_sigteff': str(int(er_temp)), 
        'star_feh': str(round(metal, 2)), 
        'star_sigfeh': str(round(er_metal, 2)),  
        'star_vmag': str(round(vmag, 2)), 
        'star_sigvmag': str(round(er_vmag, 2)), 
        'star_parallax': str(round(parallax, 3)), 
        'star_sigparallax': str(round(er_parallax, 3)), 
        'isoc_kind': 'parsec_CAF09_v1.1', 
        'flag_sismo': '0', 
        'kind_interp': '1', 
        'kind_tpagb': '0', 
        'kind_pulsecycle': '0', 
        'kind_postagb': '-1', 
        # 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat', 
        # 'sfr_file': 'tab_sfr/sfr_const_z008.dat', 
        # 'photsys_file': 'tab_mag_odfnew/tab_mag_ubvrijhk.dat',
        "sfr_maxage": "12.0e9", 
        "sfr_minage": "0.1e9",
        'submit_form': 'Submit' }

    print('Parameters for %s from Padova interface.' % star)
    for key, value in form_data.items():
        print('%s -> %s' % (key, value))
        
    # Create an insecure context
    context = ssl._create_unverified_context()

    # Use the insecure context when making a request
    response = urllib.request.urlopen(url, context=context)
    urllib.request.urlretrieve(url, 'parameters.html', lambda x, y, z:0, urllib.parse.urlencode(form_data).encode("utf-8"))
    return

def read_params():

    time.sleep(60)    
    # Read the HTML file
    with open("parameters.html", "r", encoding="utf-8") as f:
        line = f.readlines()[15]

    # Decode HTML entities
    text = unescape(line)

    # Extract the results
    pattern = (
    r"Results for\s+(.*?)\s*:"
    r".*?Age\s*=\s*(.*?)\s*±\s*(.*?)\s*Gyr"
    r".*?Mass\s*=\s*(.*?)\s*±\s*(.*?)\s*<i>"
    r".*?log<i>g</i>\s*=\s*(.*?)\s*±\s*(.*?)\s*\(cgs\)"
    r".*?<i>R</i>\s*=\s*(.*?)\s*±\s*(.*?)\s*<i>"
    )
    
    match = re.search(pattern, text, re.DOTALL)
    if match:
        ex = True
        star, age, erage, mass, ermass, logg_p, erlogg_p, radius, erradius = match.groups()
        star = star.strip()  
        age, erage, mass, ermass, logg_p, erlogg_p, radius, erradius = map(float, [age, erage, mass, ermass, logg_p, erlogg_p, radius, erradius])
        print('----------------Results for ' + str(star) + ' ----------------------------')
        print(f"Mass   = {mass:.3f} ± {ermass:.3f} M_sun")
        print(f"Radius = {radius:.3f} ± {erradius:.3f} R_sun")
        print(f"Age    = {age:.3f} ± {erage:.3f} Gyr")
        print(f"logg   = {logg_p:.3f} ± {erlogg_p:.3f} dex")
    else:
        ex = False
        print("Could not find the stellar parameters in the HTML file.")
        mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p = (np.nan,) * 8 
    print('--------------------------------------------------')
        
    results = mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p
    return results, ex
    
if __name__ == '__main__':

    ssl._create_default_https_context = ssl._create_unverified_context
    ffile = 'Ariel_stellar_catalogue_v0.1_20251001_param_input.csv'
    params = pd.read_csv(ffile, comment='#')

    for i, x in enumerate(params.name[:]):
        # Initiate
        ex = False        
        results = (np.nan,) * 8
        torres = (np.nan,) * 4
        logg_hip, erlogg_hip = (np.nan,) * 2
        name, vmag, evmag, parallax, eparallax, teff, eteff, feh, efeh = params.name[i], params.vmag[i], params.evmag[i], params.parallax[i], params.eparallax[i], params.teff[i], params.eteff[i], params.feh[i], params.efeh[i]

        # Try 3 times
        for attempt in range(3):            
            # param mass, radius, age, and logg from isochrone [add a slight push]
            get_mass_radius(name, vmag, evmag, parallax, eparallax, teff, eteff+attempt, feh, efeh, 0.1e9, 12.0e9)
            results, ex = read_params()
            if ex:
                break

        if not ex:
            print('Try again with increased errors')
            get_mass_radius(name, vmag, evmag, parallax, eparallax, teff, eteff+30, feh, efeh+0.03, 0.1e9, 12.0e9)
            results, ex = read_params()
   
        if ex:
            #Trigonometric logg. Uses the get_mass_radius function
            bc = bcflow(teff)
            logg_hip, erlogg_hip = logg_trigomonetric(teff, results[0], vmag, bc, parallax, eparallax, eteff, results[1])
            #Mass and radius from Torres calibration
            torres = mass_torres(teff, logg_hip, feh, eteff, erlogg_hip, efeh)
        else: 
            print('Unable to recover parameters')

        #Save results
        parameters = [x, teff, eteff, feh, efeh, results, torres, logg_hip, erlogg_hip]
        _output(parameters=parameters)

 
