mass_radius
===========

Get the stellar mass, radius and age from the Padova interface
This code also includes mass and radius calculations using the Torres et al. 2010 calibration.  

Create a file 'star' and fill with the parameters separated by tabs.

**Important: The 'star' must have a specific format (see example).**

In case gravity and error in gravity are not known fill the columns with zeros and ignore the logg from Torres et al. calibration. 
    
This will work if all values are filled.

The file 'star' is given as an example.
This code will return the mass and radius from the Padova interface (http://stev.oapd.inaf.it/cgi-bin/param) in the 'stellar_characterization.dat' file.

Run:

    python param.py 

