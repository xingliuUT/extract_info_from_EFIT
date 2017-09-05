# extract_info_from_EFIT

## Description
This project contains Python code for extracting useful info from [EFIT](https://fusion.gat.com/theory/Efit) files.

`read_EFIT.py` takes an EFIT file and construct a dictionary of:
* resolution `nw` (Rgrid) and `nh` (Zgrid)
* position of the magnetic axis: `rmaxis`, `zmaxis`
* poloidal magnetic flux on (R, Z) grid: `psirz`
* physical quantities: F ( = Btor * R), Pressure, FFprime, Pprime, q, jtor
* magnetic field strength: Btor, Bpol
* coordinates: `psipn`, 'rhotn`

## Usage

## Note
TODO: find total magnetic field at a given (R,Z), not limited in psipn < 1
