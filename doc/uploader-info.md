# exciting-parser

## Main file

The main file is INFO.OUT. The parsing starts from this one and it fails if it is not present.

## Other files coming from a ground-state calculation

 The following files must be in the same folder as INFO.OUT:

 * dos.xml: density of states.
 * bandstructure.xml: band structure.
 * FERMISURF.bxsf: Fermi surface.
 * EIGVAL.OUT: eigenvalues.
 * input.xml: 

## Other files from GW calculation

 GW output files must be in the same folder as the ground-state calculation.

 The GW output files that are parsed are:

 * GW\_INFO.OUT or GWINFO.OUT (it depends on the code version): it is the main GW output.
 * EVALQP.DAT or EVALQP.TXT (it depends on the code version): contains the eigenvalues.
 * TDOS\-QP.OUT: density of states.
 * bandstructure\-qp.dat or BAND\-QP.OUT (it depends on the code version): band structure
 * BANDLINES.OUT: start and end point of the segments in the band structure plot.
 * SELFC.DAT: self energy
 * input.xml or input\-gw.xml or input\_gw.xml: most of the parameters that characterize the GW calculation are retrieved from this file.
