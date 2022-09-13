from PyAstronomy import pyasl

# Create class instance
r = pyasl.Ramirez2005()

# Which color bands are available
print("Available color bands: ", r.availableBands())

# Convert B-V to effective temperature and back
bv = 0.307
feh = -1.10
teff = r.colorToTeff("B-V", bv, feh)
bv1 = r.teffToColor("B-V", teff, feh)
# Watch out for differences between input bv and the output bv1
print("B-V = ", bv, ", Teff = ", teff, ", bv1 = ", bv1, ", bv-bv1 = ", bv-bv1)