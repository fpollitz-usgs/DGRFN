		Program Package DGRFN

		by Fred F. Pollitz

---------------------------------------------------------------------------------------

These programs are an implementation of the Direct Green's Function method
described by Friederich and Dalkolmo (1995) and Dalkolmo (1993).
They solve the seismic wave equation in a spherically layered isotropic medium
using a decomposition into spheroidal and toroidal motions.  For each spherical
harmonic degree l and azimuthal order number m, the (l,m) response function is
deternined subject to jumps in the displacement-stress vector at the source radius,
a zero-traction boundary condition at Earth's surface, and a homogeneous isotropic
elastic solid at the base of the specified Earth model.

DGRFN uses two subroutines (with modifications) provided in the GEMINI-2.2 release
authored by Joerg R. Dalkolmo.  The remainder is an adaptation of the static-deformation
code COGRAV, authored by Fred F. Pollitz and described by Pollitz (1996), to the
case of non-zero complex frequency to simulate seismic wave propagation, including
the effects of anelasticity.

DGRFN presently handles only an isotropic medium; it is intended to modify it
soon to handle a laterally isotropic medium.

There are two main programs:

1) WAVE0:	Computes response functions for spherical harmonic degrees
		from l=0 up to input l=[lmax], and azimuthal order numbers
		m=0, 1, and 2.  The source and observation-pt depths are fixed.

2) WAVE1:	Computes seismic velocity time series at specified observation points
		and source geometry, assuming a distribution of point sources
		along a horizontal line segment at the source depth specified
		in WAVE0 and the observation-pt depth specified in WAVE0.
		A step-function source is assumed.

---------------------------------------------------------------------------------------
	
		EXAMPLE 1

A homogeneous structure with kappa=50 GPa, mu=30 GPa, and rho=3 g-cm^{-3} 
is represented by 'earth.modelHOMO':

  14   300.000
[# layers, Earth radius in km]
   42.000   39.000    3.000    5.000    3.000   25.000   62.500
[layer 1: lower depth, upper depth (km), rho (g-cm^(-3)), kappa (10^10 Pa),
mu (10^10 Pa), Q_beta, Q_kappa]
   39.000   36.000    3.000    5.000    3.000   25.000   62.500
[layer 2: lower depth, upper depth (km), rho (g-cm^(-3)), kappa (10^10 Pa),
mu (10^10 Pa), Q_beta, Q_kappa]
   36.000   33.000    3.000    5.000    3.000   25.000   62.500
   33.000   30.000    3.000    5.000    3.000   25.000   62.500
   30.000   27.000    3.000    5.000    3.000   25.000   62.500
   27.000   24.000    3.000    5.000    3.000   25.000   62.500
   24.000   21.000    3.000    5.000    3.000   25.000   62.500
   21.000   18.000    3.000    5.000    3.000   25.000   62.500
   18.000   15.000    3.000    5.000    3.000   25.000   62.500
   15.000   12.000    3.000    5.000    3.000   25.000   62.500
   12.000    9.000    3.000    5.000    3.000   25.000   62.500
    9.000    6.000    3.000    5.000    3.000   25.000   62.500
    6.000    3.000    3.000    5.000    3.000   25.000   62.500
    3.000    0.000    3.000    5.000    3.000   25.000   62.500
[layer 14: lower depth, upper depth (km), rho (g-cm^(-3)), kappa (10^10 Pa),
mu (10^10 Pa), Q_beta, Q_kappa]

In general, the layer thickness (which may be variable and not necessarily constant
as in this example) should be no greater than the shortest wavelength
expected to be synthesized.  A rough estimate of this wavelength is
(2 * pi * R)/(lmax + 0.5), where R is Earth's model radius and lmax is
the maximum spherical harmonic degree

DGRFN (through WAVE0) will implement a homogeneous sphere below the deepest
layer listed in 'earth.model' -- assigning it the same elastic parameters and
density as this deepest layer.

The example consists of running two command files:  

1) wave0.xEX1 calculates
the response functions at spherical harmonic degrees from 0 to 2400 and a set 
of complex-valued angular frequencies just below the real frequency-axis. 
The command file is

cp earth.modelHOMO earth.model
[copy 'earth.modelHOMO' onto 'earth.model' which is read in by WAVE0

wave0 << ! > /dev/null 
2400 [maximum spherical harmonic degree]
19. [source depth (km)]
29. [observation pt. depth (km)]
0.07 [time interval dt, used to calculate the frequency spacing]  
256 [nlen=number of time samples, and synthetic records will go from t=0
to t=0.07*256=17.92 sec.  The shortest period computed for the Green's functions
is 3*dt=0.21 sec.]

In general, dt should
be <= one-seventh the corner period intended to be used in WAVE1.

In general, the maximum spherical harmonic degree lmax should be chosen to be 
somewhat greater than the degree of the fundamental mode associated with a period of
twice the corner period intended to be used in WAVE1.

2) wave1.xEX1 calculates three-component velocity seismograms.
The source depth and observation depth are fixed at the values input to WAVE0.
The command file has several independent runs of WAVE1.  One section reads:

wave1 << !
0.7 [corner period (sec)]
1 [# sources]
0. 0. 0.01 180. 60. 0. 1.e-6 [latitude,longitude,segment length (km),segment strike (deg.
CW from due North], dip (deg.), dislocation rake (deg.), scalar seismic moment (10^20 N m).
The latitude and longitude are the location of the segment endpoint closest to
the strike direction.]
1 [# observation points]
0. 0.1694 [latitude and longitude of observation point]
!

The output seismograms are in 
wave1.outx
wave1.outy
wave1.outz
where x, y, and z denote the local East, North, and Up directions, respectively.
The time series are given in (t,v) pairs, where t is in sec after the origin time
and v is in cm/yr.
Output time series are oversampled by a factor of 4 relative to the number
of frequency samples used in the inverse Fourier transform.
If N observation points are specified, then each output file contains
4*N*(nlen/2) lines, one set of 4*(nlen/2) lines for each observation point in the
order in which they were read in.

The above command file computes the response to a shear dislocation on a plane
defined to have unit normal pointing from the source to the receiver.  The non-trivial
displacement is in the y-direction (i.e. north component).  

wave1.xEX1 computes the response to a shear dislocation or an isotropic source
at a set of (dip, receiver location)
pairs with source-receiver distance equal to 10, 20, 30, and 40 km, the dip always
being chosen so that the dislocation plane has unit normal pointing from the source 
(at 19 km depth) to the 
receiver (at 29 km depth).

The sections of wave1.xEX1 that do an isotropic source have the form

wave1 << !
0.7
1
0. 0. 0.01 180. 60.00 0. -1.e-6
1.e-6 1.e-6 1.e-6 0. 0. 0.
1
0. 0.1694
!

In this case, the negative value of the input seismic moment (-1.e-6) means
that WAVE1 will read in an additional line that has the moment tensor components
Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r=Up, t=theta=South, p=Phi=South,
in units of 10^20 N m.

The resulting three-component velocity time series have been resolved into longitudinal
and transverse components and compared with the analytic solution (see validation.pdf).


		EXAMPLE 2

In this example we compute synthetic seismograms for the M5.2 April 18, 2008 
Mount Carmel, IL earthquake on a layered model appropriate for the western US.
The layered structure is represented by 'earth.model-Nuttli-1969':

  10  6371.000
   42.000   40.000    3.350   12.464    7.463  500.000 1250.000
   40.000   35.000    2.790    6.952    4.178  500.000 1250.000
   35.000   30.000    2.790    6.952    4.178  500.000 1250.000
   30.000   25.000    2.790    6.952    4.178  500.000 1250.000
   25.000   20.000    2.790    6.952    4.178  500.000 1250.000
   20.000   15.000    2.580    5.423    3.251  250.000  625.000
   15.000   10.000    2.580    5.423    3.251  250.000  625.000
   10.000    5.000    2.580    5.423    3.251  250.000  625.000
    5.000    1.000    2.580    5.423    3.251  250.000  625.000
    1.000    0.000    2.150    2.981    1.796  100.000  250.000

1) wave0.xEX2 calculates
the response functions at spherical harmonic degrees from 0 to 22000 and a set 
of complex-valued angular frequencies just below the real frequency-axis. 

wave0 << ! > /dev/null 
11400 [maximum spherical harmonic degree]
11. [source depth (km)]
0. [observation pt. depth (km)]
0.48 [time interval dt, used to calculate the frequency spacing]  
200 [nlen=number of time samples, and synthetic records will go from t=0
to t=0.48*200=96 sec.  The shortest period computed for the Green's functions
is 3*dt=1.44 sec.]

2) wave1.xEX2 calculates three-component velocity seismograms at 6 stations
of the Cooperative New Madrid Seismic Network.
The source depth and observation depth are fixed at the values input to WAVE0.

wave1 << ! > /dev/null
4.0 [corner period (sec)]
1 [# point sources]
38.472 -87.918 0.01 296. 84. 4. 9.3e-4 [latitude,longitude,segment length (km),segment strike (deg.
CW from due North], dip (deg.), dislocation rake (deg.), scalar seismic moment (10^20 N m).
The latitude and longitude are the location of the segment endpoint closest to
the strike direction.]
6 [# observation points]
37.96510 -87.66600 USIN [latitude and longitude of observation point; station name (optional)]
38.73380 -88.09910 OLIL
37.97160 -87.52970 EVIN
37.75270 -88.43730 HAIL
38.13000 -87.93600 NHIN
38.42980 -87.78170 WVIL
!

Resulting three-component seismograms at these 6 stations (and the result
of a source relocation) and corresponding observed seismograms,
courtesy of Chris Kramer, are shown in the figure 'EX2-figure.pdf'

---------------------------------------------------------------------------------------

		VALIDATION

Using the results of three-component velocity computed with wave1.xEX1 on a homogeneous
sphere, validation.pdf contains a comparison between these Direct Green's Function
solutions and the corresponding solutions on a homogeneous full space (Stokes solution)
as given by Aki and Richards (1980).

DGRFN has also been used to replicate Figure 13 of Kuhn (1985) (vertical- 
component displacement wavefield resulting from an explosive source
observed at depth), which is based on numerical finite difference.
These same calculations, including both vertical and horizontal components, 
also agree closely with AXITRA (Bouchon, 1981), which is based on a frequency-
wavenumber integration.

---------------------------------------------------------------------------------------

		REFERENCES

Bouchon, M., 1981. A simple method to calculate Green's functions for elastic layered 
media, Bull. Seismol. Soc. Am., 71, 959–971.

Dalkolmo, J., 1993. Synthetische Seimogramme fur eine spharisch
symmetrische, nichtrotierende Erde durch direkte Berechnung
der Greenschen Funktion, Diploma thesis, Institute of
Geophysics, Stuttgart University.

Friederich, W. and Dalkolmo, J., 1995.  Complete synthetic seismograms for a spherically
symmetric earth by a numerical computation of the Green's function in the frequency
domain, Geophys. J. Int., 122, 537-550.

Kuhn, M. J., A numerical study of Lamb's Problem, Geophysical Prospecting, 33, 
1103-1137, 1985.

Pollitz, F. F., 1996. Coseismic deformation from earthquake faulting on a
layered spherical earth. Geophys. J. Int., 125, 1-14.
