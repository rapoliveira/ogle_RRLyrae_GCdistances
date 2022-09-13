=======================================================================
     Over 78000 RR Lyrae Stars in the OGLE Galactic Disk Fields

 I. Soszynski, A. Udalski, M. K. Szymanski, P. Pietrukowicz, P. Mroz,
    J. Skowron, S. Kozlowski, R. Poleski, D. Skowron, K. Ulaczyk,
           K. Rybicki, P. Iwanek, M. Wrona, M. Gromadzki

=======================================================================
soszynsk@astrouw.edu.pl

This directory contains the OGLE catalog of RR Lyrae stars identified
in the fields toward the Galactic bulge.

The directory structure is as follows:

README            - this file

ident.dat         - identification of stars

RRab.dat          - parameters of RRab stars
RRc.dat           - parameters of RRc stars
RRd.dat           - parameters of RRd stars

phot/I/           - I-band photometry of individual objects
phot/V/           - V-band photometry of individual objects
phot.tar.gz       - gzipped phot/ directory

fcharts/          - finding charts of individual objects
                       
remarks.txt       - remarks on selected objects

pap.pdf           - PDF version of the paper Soszynski et al. (2019), 
                    Acta Astron. 69, 321, (arXiv:2001.00025) describing
                    the catalog.


Format of the file ident.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 19  A19   ---     Star's ID
 22- 25  A4    ---     Subtype of RR Lyr star (RRab, RRc, RRd, aRRd)
 28- 29  I2    h       Right ascension, equinox J2000.0 (hours)
 31- 32  I2    m       Right ascension, equinox J2000.0 (minutes)
 34- 38  F5.2  s       Right ascension, equinox J2000.0 (seconds)
     40  A1    ---     Declination, equinox J2000.0 (sign)
 41- 42  I2    deg     Declination, equinox J2000.0 (degrees)
 44- 45  I2    arcmin  Declination, equinox J2000.0 (arc minutes)
 47- 50  F4.1  arcsec  Declination, equinox J2000.0 (arc seconds)
 53- 68  A16   ---     OGLE-IV ID
 70- 84  A15   ---     OGLE-III ID
 86-100  A15   ---     OGLE-II ID
102-     A     ---     Other designation (from VSX)
--------------------------------------------------------------------------

Format of the files RRab.dat and RRc.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 19  A19   ---     Star's ID
 22- 27  F6.3  mag     Intensity mean I-band magnitude
 29- 34  F6.3  mag     Intensity mean V-band magnitude
 37- 46  F10.8 days    Period
 48- 57  F10.8 days    Uncertainty of period
 60- 69  F10.5 days    Time of maximum brightness (HJD-2450000)
 72- 76  F5.3  mag     I-band amplitude (maximum-minimum)
 79- 83  F5.3  ---     Fourier coefficient R_21
 85- 89  F5.3  ---     Fourier coefficient phi_21
 92- 96  F5.3  ---     Fourier coefficient R_31
 98-102  F5.3  ---     Fourier coefficient phi_31
--------------------------------------------------------------------------

Format of the files RRd.dat and aRRd.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 19  A19   ---     Star's ID
 22- 27  F6.3  mag     Intensity mean I-band magnitude
 29- 34  F6.3  mag     Intensity mean V-band magnitude
 37- 46  F10.8 days    First-overtone period
 48- 57  F10.8 days    Uncertainty of first-overtone period
 60- 69  F10.5 days    Time of maximum brightness (HJD-2450000)  
 72- 76  F5.3  mag     I-band amplitude (maximum-minimum)
 79- 83  F5.3  ---     Fourier coefficient R_21
 85- 89  F5.3  ---     Fourier coefficient phi_21
 92- 96  F5.3  ---     Fourier coefficient R_31
 98-102  F5.3  ---     Fourier coefficient phi_31
105-114  F10.8 days    Fundamental-mode period
116-125  F10.8 days    Uncertainty of fundamental-mode period
128-137  F10.5 days    Time of maximum brightness (HJD-2450000)
140-144  F5.3  mag     I-band amplitude (maximum-minimum)
147-151  F5.3  ---     Fourier coefficient R_21
153-157  F5.3  ---     Fourier coefficient phi_21
160-164  F5.3  ---     Fourier coefficient R_31
166-170  F5.3  ---     Fourier coefficient phi_31
--------------------------------------------------------------------------

Finding charts are gzipped Postscript images. The names of the files are in
the form: OGLE-BLG-RRLYR-NNNNN.ps.gz. The finding charts are 60"x60" subframes
of the I-band reference frames centered on the star. White cross marks the
star. North is up and East to the left.


Any presentation of the scientific analysis or usage of the data from the
catalog of RR Lyrae stars in the Galactic bulge should include the appropriate
reference(s) to OGLE paper(s).

