mass_radius
===========

Get the stellar mass, radius and age from the Padova interface
This code also includes mass and radius calculations using the Torres et al. 2010 calibration.  

Create a file 'star' and fill with the parameters separated by tabs.

**Important: The 'star' must have a specific format. Tab separated values and 2 lines as header.**
The columns should be in order of the following: 

    star name, temperature, error in temperature, gravity, error in gravity, metallicity, error in metallicity, V magnitude, error in V, parallax, error in parallax 

In case gravity and error in gravity are not known fill the columns with zeros and ignore the logg from Torres et al. calibration. 

In case parallax is not known, mass and radius will be calculated by Torres et. al. 2010 and the file should include the following columns: 

    star name, temperature, error in temperature, gravity, error in gravity, metallicity, error in metallicity 
    
This will work if all values are filled.

The file 'star' is given as an example.
This code will return the mass and radius from the Padova interface (http://stev.oapd.inaf.it/cgi-bin/param) in the 'stellar_characterization.dat' file.

In case there is no information on the parallaxes, masses and radii will be calculated using only the calibration. Therefore, the user must specify with a flag

Run in case parallaxes are included:

    python mass_radius_age.py star on 

Run in case parallaxes are not included:

    python mass_radius_age.py star off
    
The user inputs the name of the file.
