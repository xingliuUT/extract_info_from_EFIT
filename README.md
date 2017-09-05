# extract_info_from_EFIT

## Description
This project contains Python code for extracting useful info from [EFIT](https://fusion.gat.com/theory/Efit) files.

`read_EFIT.py` takes an EFIT file and construct a dictionary of:
* resolution `nw` (Rgrid) and `nh` (Zgrid)
* position of the magnetic axis: `rmaxis`, `zmaxis`
* poloidal magnetic flux on (R, Z) grid: `psirz`
* physical quantities: F ( = Btor * R), Pressure, FFprime, Pprime, q, jtor
* magnetic field strength: Btor, Bpol
* coordinates: `psipn`, `rhotn`

`calc_volume_from_EFIT.py` contains methods:
* `dVolume` returns the total volume between two given flux surfaces
* `totalV` scans from magnetic axis to the boundary to write out volume elements versus `psipn` coordinate in EFIT file
* `surfaceArea` returns the total surface area given the flux surface `psipn`
* `totalN` multiplies volume element with density profile to get the total stored particle numbers

## Usage

## Note
TODO: find total magnetic field at a given (R,Z), not limited in psipn < 1
