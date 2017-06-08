mass_radius
===========

Get the stellar mass, radius and age from the Padova interface
This code also includes mass and radius calculations using the Torres et al. 2010 calibration.  

Create a file 'star' and fill with the parameters separated by tabs: star name, temperature, error in temperature, metallicity, error in metallicity, V magnitude, error in V, parallax, error in parallax.
The file 'star' is given as an example.
This code will return the mass and radius from the Padova interface (http://stev.oapd.inaf.it/cgi-bin/param) in the 'stellar_characterization.dat' file.

In case there is no information on the parallaxes, masses and radii will be calculated using only the calibration. Therefore, the user must specify with a flag

Run in case parallaxes are included:

    python mass_radius_age.py on 

Run in case parallaxes are not included:

    python mass_radius_age.py off
    
The user inputs the name of the file.
