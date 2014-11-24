#!/usr/bin/python

import urllib
import os, sys

#Create a file 'star' and fill with the parameters: star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
#This code will get the mass and radius of the stars from the Padova interface in the 'mass_radius.txt' file.

def get_mass_radius(star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal):
	"""Returns mass and radius from padova interface.
	   Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
	"""
	url = 'http://stev.oapd.inaf.it/cgi-bin/param'
	#These are the parameters in the webpage to tune
	form_data = {'param_version': '1.3', 'star_name': star, 'star_teff': temp, 'star_sigteff': er_temp, 'star_feh': metal, 'star_sigfeh': er_metal,  'star_vmag': vmag, 'star_sigvmag': '0.0', 'star_parallax': parallax, 'star_sigparallax': er_parallax,  'isoc_kind': 'parsec_CAF09_v1.1', 'kind_interp': '1', 'kind_tpagb': '0', 'kind_pulsecycle': '0', 'kind_postagb': '-1', 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat', 'sfr_file': 'tab_sfr/sfr_const_z008.dat', 'sfr_minage': '0.1e9', 'sfr_maxage': '12.0e9', 'flag_sismo': '0', 'submit_form': 'Submit' }
	print form_data
	urllib.urlretrieve(url, 'parameters.html', lambda x,y,z:0, urllib.urlencode(form_data))

	#write results
	with open('parameters.html') as f:
		line = f.readlines()[19]

	mass = []
	radius = []
	name = []
	ermass = []
	erradius = []

       	line = line.replace(' .<p>Results for ', '')
 	line = line.replace('Mass=', '')
 	line = line.replace('&plusmn;', ' , ')
	line = line.replace(' ','')
	line = line.replace('<i>R</i>=','')
	line = line.replace(':Age=',',')
	line = line.replace('<i>M</i>&#9737','')
	line = line.replace('<i>R</i>&#9737','')
 	line = line.split(',')
       	name.append(str(line[0]))
       	mass.append(str(line[3]))
       	ermass.append(str(line[4]))
       	radius.append(str(line[7]))
       	erradius.append(str(line[8]))

	return name, mass, ermass, radius, erradius
################################################################################
#'star' file is the input file with my parameters 
with open('star') as f:

		lines = f.readlines()[1:]

hd = [] 	#star name
v = [] 		#vmag
par = [] 	#parallax in mas
er_par = []	#error in parallax
teff = []	#temperature
er_teff = []	#error in temperature
fe = []		#metallicity
er_fe = []	#error in metallicity

for line in lines:
 	line = line.split()
       	hd.append(str(line[0]))
       	v.append(float(line[1]))
       	par.append(float(line[2]))
       	er_par.append(float(line[3]))
       	teff.append(float(line[4]))
       	er_teff.append(float(line[5]))
       	fe.append(float(line[6]))
       	er_fe.append(float(line[7]))

#Write results in the 'mass_radius.txt' file
out_file = open('mass_radius.txt','w')
out_file.writelines('star' + '\t' + 'mass' + '\t' + 'er_mass' + '\t' + 'radius' + '\t' + 'er_radius' + '\n')
for i, x in enumerate(hd[:]):
	name, mass, er_mass, radius, er_radius = get_mass_radius(hd[i], v[i], par[i], er_par[i], teff[i], er_teff[i], fe[i], er_fe[i])
	out_file.writelines(str(name[0]) + '\t' + str(mass[0]) + '\t' + str(er_mass[0]) + '\t' + str(radius[0]) + '\t' + str(er_radius[0]) + '\n' )
out_file.close()
print '--------------------------THE END------------------------'
os.system('rm parameters.html')
