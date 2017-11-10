# SMIP
Satellite Model Integrated Product

process_smip.py is the top level code to collocate quantities from other satellite (IMERG) and reanalysis (ECMWF) products to the native SEVIRI geotemporal grid. The basis of the code works on daily intervals, and hence, the inputs are YEAR MONTH DAY and version number. The SEVIRI data required for this product requires a merged output from the ORAC postprocessor and derived products FORTRAN module and is run separately using an IDL code maintained by Matt Christensen U. Oxford.
