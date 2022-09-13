=======================================================================
     Over 78000 RR Lyrae Stars in the OGLE Galactic Bulge Fields

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

pap1.pdf          - PDF version of the paper Soszynski et al. (2014),
                    Acta Astron. 64, 177, (arXiv:1410.1542) describing the
                    first part of the catalog.
pap2.pdf          - PDF version of the paper Soszynski et al. (2019), 
                    Acta Astron. 69, 321, (arXiv:2001.00025) describing the
                    second part of the catalog.


Format of the file ident.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 20  A20   ---     Star's ID
 23- 26  A4    ---     Subtype of RR Lyr star (RRab, RRc, RRd, aRRd)
 29- 30  I2    h       Right ascension, equinox J2000.0 (hours)
 32- 33  I2    m       Right ascension, equinox J2000.0 (minutes)
 35- 39  F5.2  s       Right ascension, equinox J2000.0 (seconds)
     41  A1    ---     Declination, equinox J2000.0 (sign)
 42- 43  I2    deg     Declination, equinox J2000.0 (degrees)
 45- 46  I2    arcmin  Declination, equinox J2000.0 (arc minutes)
 48- 51  F4.1  arcsec  Declination, equinox J2000.0 (arc seconds)
 54- 69  A16   ---     OGLE-IV ID
 71- 85  A15   ---     OGLE-III ID
 87-101  A15   ---     OGLE-II ID
103-     A     ---     Other designation (from VSX)
--------------------------------------------------------------------------

Format of the files RRab.dat and RRc.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 20  A20   ---     Star's ID
 23- 28  F6.3  mag     Intensity mean I-band magnitude
 30- 35  F6.3  mag     Intensity mean V-band magnitude
 38- 47  F10.8 days    Period
 49- 58  F10.8 days    Uncertainty of period
 61- 70  F10.5 days    Time of maximum brightness (HJD-2450000)
 73- 77  F5.3  mag     I-band amplitude (maximum-minimum)
 80- 84  F5.3  ---     Fourier coefficient R_21
 86- 90  F5.3  ---     Fourier coefficient phi_21
 93- 97  F5.3  ---     Fourier coefficient R_31
 99-103  F5.3  ---     Fourier coefficient phi_31
--------------------------------------------------------------------------

Format of the files RRd.dat and aRRd.dat:
--------------------------------------------------------------------------
 Bytes  Format Units   Description
--------------------------------------------------------------------------
  1- 20  A20   ---     Star's ID
 23- 28  F6.3  mag     Intensity mean I-band magnitude
 30- 35  F6.3  mag     Intensity mean V-band magnitude
 38- 47  F10.8 days    First-overtone period
 49- 58  F10.8 days    Uncertainty of first-overtone period
 61- 70  F10.5 days    Time of maximum brightness (HJD-2450000)  
 73- 77  F5.3  mag     I-band amplitude (maximum-minimum)
 80- 84  F5.3  ---     Fourier coefficient R_21
 86- 90  F5.3  ---     Fourier coefficient phi_21
 93- 97  F5.3  ---     Fourier coefficient R_31
 99-103  F5.3  ---     Fourier coefficient phi_31
106-115  F10.8 days    Fundamental-mode period
117-126  F10.8 days    Uncertainty of fundamental-mode period
129-138  F10.5 days    Time of maximum brightness (HJD-2450000)
141-145  F5.3  mag     I-band amplitude (maximum-minimum)
148-152  F5.3  ---     Fourier coefficient R_21
154-158  F5.3  ---     Fourier coefficient phi_21
161-165  F5.3  ---     Fourier coefficient R_31
167-171  F5.3  ---     Fourier coefficient phi_31
--------------------------------------------------------------------------

Finding charts are gzipped Postscript images. The names of the files are in
the form: OGLE-BLG-RRLYR-NNNNN.ps.gz. The finding charts are 60"x60" subframes
of the I-band reference frames centered on the star. White cross marks the
star. North is up and East to the left.


Any presentation of the scientific analysis or usage of the data from the
catalog of RR Lyrae stars in the Galactic bulge should include the appropriate
reference(s) to OGLE paper(s).


Updates:

2014-10-16  Reclassified OGLE-BLG-RRLYR-22097: RRd -> RRc
            (thanks to Radoslaw Smolec)

2017-12-06  828 additional RR Lyrae stars (from OGLE-BLG-RRLYR-38290
            to OGLE-BLG-RRLYR-39122; Soszynski et al. 2017, Acta
            Astron. 67, 297).
            Light curves supplemented with new datapoints, up to
            August 2017.
            Several stars reclassified and removed from the list.

2020-01-02  29075 additional RR Lyrae stars (from OGLE-BLG-RRLYR-39123
            to OGLE-BLG-RRLYR-68197; Soszynski et al. 2019, Acta
            Astron. 69, 321).
            OGLE-BLG-RRLYR-16703 reclassified and removed from the list.
