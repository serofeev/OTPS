# OTPS
Oregon State University Tidal Prediction Software
AUTHORS:
  Gary Egbert & Lana Erofeeva
  College of Atmospheric and Oceanic Sciences
  104 COAS Admin. Bldg.
  Oregon State University
  Corvallis, OR 97331-5503
  
  E-mail:  egbert@coas.oregonstate.edu                                      
  Fax:     (541) 737-2064
  Ph.:     (541) 737-2947                                        
  https://www.tpxo.net/

COPYRIGHT: OREGON STATE UNIVERSITY, 2010
(see the file COPYRIGHT for license agreement)
=====================================================================
OTPS UPDATES 2021
- 2q1 and m3 constituents added;
- added chek of minor constituents already in the model. If a minor
  constituent is in the model or setup.inp file, it is not inferred.
=====================================================================
OTPS UPDATES 2020
- new load file obtained from TPXO9 included;
- ability to directly access single bytes of TPXO(atlas) bin files
  added; this significantly reduces required RAM and speeds up output;
- matlab script tpxo_atlas2local added to OTPS/matlab. The script can be
  used as plug-in to TMD matlab toolbox: 
  https://www.esr.org/research/polar-tide-models/tmd-software/
  TMD works with multi-constituent OTPS bin format. tpxo_atlas2local 
  can be used to make a multi-constituent outcut from the single constituent
  atlas files i.e. provides the same functionality as extract_local_model.f90
******************************************************************************
1. INTRODUCTION

OTPS accomplishes 3 tasks with 3 executables:
- extracting harmonic constants from a barotropic tidal solution
  in OTPS binary format (see below) at given locations - extract_HC;

- predicting tides at given times and locations from any local,
  global or atlas solution given in OTPS binary format - predict_tide;

- extracting a local solution from global TPXO8/9-atlas in OTPS format-
  extract_local_model.

Global, regional and global-atlas barotropic inverse solutions in OTIS format
are available for download at:

  https://www.tpxo.net/

******************************************************************************
2. OTPS BINARY FORMAT

The files are big_endian binary, recorded on Linux. They
can be read sequentially or by direct access.

There is one header record, which gives the grid size (n grid cells in
longitude and m grid cells in latitude), number of constituents (nc),
limits of the area (theta_lim(2),phi_lim(2)), and constituent names
(c_id(nc)). The header is followed by nc records, each giving the
either elevation fields (m, ocean tide) OR transports fields (m^2/s)
for one constituent.

Elevations (h) and transports/currents are given as complex amplitudes,
so that the partial tide for a single constituent of frequency w is given by

     h(t,x) = Re [ h(x) exp { i [w (t - t0) + V0(t0)] } ]

where V0(t0) is the astronomical argument for the constituent at t0.
Note that with the usual conventions, amplitude and phase are given
by    amp = | h |    phase = atan (-Im(h)/Re(h)) .

To read with a simple FORTRAN program:
         .
         .
         .
      integer n,m,nc
      real theta_lim(2),phi_lim(2)
      character*4  c_id(21)
      complex, allocatable: h(:,:,:),uv(:,:,:)
      complex, allocatable: u(:,:,:),v(:,:,:)
         .
         .
         .
      open(unit=1,file='h_tpxo9.v1',form='unformatted',
     *     status='old')
ccc  (header)
      read(1) n,m,nc,theta_lim,phi_lim,c_id(1:nc) 
      allocate(h(nc,n,m))
      do ic=1,nc
       read(1)h(ic,:,:)
      enddo   
         .
         .
         .
      allocate(uv(2,n,m),u(nc,n,m),v(nc,n,m))
      open(unit=1,file='u_tpxo9.v1',form='unformatted',)
     *     status='old')
ccc  (skip header)
      read(1)
      do ic=1,nc 
       read(1) uv
       u(ic,:,:)=uv(1,:,:)
       v(ic,:,:)=uv(2,:,:) 
      enddo 
         .
         .
         .

  you now have elevations (m/s) in array h  and transports in
  arrays u and v for nc constituents on a  n x m grid, more precisely:

==> u(1,.,.) gives the complex amplitude of  zonal transport 
            (m**2/s) to the East ... call this U
==> u(2,.,.) gives meridional transport
            (m**2/s) to the North  ... call this V

See also matlab scripts h_in.m, h_out.m, u_in,m, uv_out.m in OTPS/matlab
on how to read/write OTPS binary files.

******************************************************************************
3. Arakawa C-GRID

There are "n" divisions in longitude, "m" in latitude. The elevations
and transports are given on a C-grid. Elevation nodes are at the grid
cell center. U nodes are at the center of left (east) side, and
V nodes are at the center of bottom (south) side of grid cells.
Then lats and lons (in degrees) for individual grid nodes are:

For h-nodes   lat(i,j) =  theta_lim(1) + (j-.5)*dy
              lon(i,j) =  phi_lim(1) + (i-.5)*dx

For U-nodes:  lat(i,j) =  theta_lim(1) + (j-.5)*dy
              lon(i,j) =  phi_lim(1) + (i- 1)*dx

For V-nodes:  lat(i,j) =  theta_lim(1) + (j- 1)*dy
              lon(i,j) =  phi_lim(1) + (i-.5)*dx

Here dx,dy is resolution in longitude and latitude correspondingly.
Most of our solutions are given on C grids uniform in lats,lons, but
some of them (i.e. Arctic solutions) are given on C grids uniform in
kilometers. The OTPS is applicable for these solutions also.
 
See matlab scripts grd_in.m, grd_out.m in OTPS/matlab on how
to read/write grid* files.

******************************************************************************
4. SETUP file "setup.inp"

This is an example of the input file for extract_HC and predict_tide,
included with OTPS. You can change the file name and edit for your
path, input file names, and parameter values.

DATA/Model_tpxo9.v1        ! 1. tidal model control file
lat_lon_time               ! 2. latitude/longitude/time file
z                          ! 3. z/U/V/u/v
m2,s2,n2,k2,k1,o1,p1,q1    ! 4. tidal constituents to include
AP                         ! 5. AP/RI
oce                        ! 6. oce/geo
1                          ! 7. 1/0 correct for minor constituents
tmp                        ! 8. output file (ASCII)


Comments on lines 1-7

1. Tidal model control file (ASCII, supplied with the model, starting
   from "Model_") contains AT LEAST 3 lines:

   elevation model file name
   transport model file name
   bathymetry grid file name
   {name of function to convert x,y to lon,lat}

   Please show FULL path to the model files in the "Model_*" file
   unless they are located in /DATA (download default).

   4th line in Model_* file is used ONLY for models calculated on
   uniform grid in km. Converting functions are provided with 
   current version of OTPS. Try upgrading OTPS if a function is
   missing.

2. latitude, longitude and time file is an ASCII file, consisting 
   eighter of 2 OR 8 columns:
   latitude (degrees) longitude (degrees) [yyyy mm dd hh mm ss]
   Common sign convention:
   lat>0 - degrees North, lon>0 - degrees East
   lat<0 - degrees South, lon<0 - degrees West
   [yyyy mm dd hh mm ss] - year month day hour seconds GMT
   These 6 columns are needed for tide predictions only.
   You may leave them empty if only extracting HC.
   If you want tidal predictions at the same times, but different locations,
   you may also provide only 2 columns lat_lon file. BUT then you have to
   provide also 6 column "time" file, consisting of times, given as
    yyyy mm dd hh mm ss
   In this case predict_tide usage will be:

   predict_tide -ttime_file<setup.inp

   This option is useful, if you want to obtain time series at open 
   boundary nodes.
   
3. Extract HC/predict tide for:
   z/U/V/u/v - elevation(m) /WE transport(m^2/s)/ SN transport(m^2/s)/
                             WE velocity (cm/s) / SN velocity (cm/s)
   NOTE: Changed, Nov 2004: now for any of u/v/U/V all components
         (that is U V u v) are calculated by predict_tide.

4. Constituent names should be in LOWER case and separated by comma.
   Leave the line blank, if ALL model constituents are included

5. AP/RI - output amplitude and phase(GMT) OR real/imaginary parts.
           Only used when extracting HC.
           May leave blank, if predicting tide.

6. oce/geo - extract ocean/geocentric HC for elevations only.
             Geocentric tide is appropriate for comparison to or
             correction of altimetry data.
             May leave blank if not z on line 3.

7. 1/0 - Do/Not correct for the 16 minor constituents:
   2q1, sigma1, rho1, m1, chi1, pi1, phi1, theta1, j1, oo1, 2n2,
   mu2, nu2, lambda2, l2, t2
 
   Inference only works if 8 major tidal constituents
   are included: m2,s2,n2,k2,k1,o1,p1,q1

   If a minor constituent is in the model and included, then it is
   not inferred. For example, TPXO9v5a includes 2n2,oo1, mu2, nu2, l2, j1, 2q1.
   If you leave line 4 blank (i.e. include all constituents), then only sigma1,
   rho1, m1, chi1, pi1, phi1, theta1,lambda2 and t2 will be inferred.

*****************************************************************************
5. SETUP file "setup.local"

This is an example of the input file for executable extract_local_model,
included with OTPS. You can change the file name and edit for your
path, input file names, and parameter values.

DATA/Model_tpxo9_atlas_v5 ! 
DATA/Model_Hawaii         ! Local area Model name
 16  28                   ! Lat limits (degrees N)
185 208                   ! Lon limits (degrees E, range -180 +360)

Here control file Model_tpxo9_atlas_v5 includes 3 lines with file names
(please edit for your own path!) for tpxo9-atlas elevation, transport and grid
files. Since tpxo9-atlas consists of single constituent files, wild card
(i.e. * ) usage is allowed in atlas names. Sample file DATA/Model_tpxo9_atlas_v5
is included with OTPS and looks like:

DATA/TPXO9_atlas_v5/h_*_tpxo9_atlas_30_v5
DATA/TPXO9_atlas_v5/u_*_tpxo9_atlas_30_v5
DATA/TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5

The Model_* control file on the second line of "setup.local" 
(i.e. DATA/Model_Hawaii in the example above) should contain 3 lines with
path and desired names of output elevaion, transports and grid files, where you 
intend to place the outcut from tpxo9-atlas. Output files will contain
all constituents included in tpxo9-atlas you have downloaded and placed
into ONE special directory, as for example OTPS/DATA/TPXO9_atlas_v5.

Instead of using executable extract_local_model you can use matlab script 
tpxo_atlas2local.m. See usage from matlab. Included script is for linux
or Windows, the script for MacOS is available for download
at www.tpxo.net

******************************************************************************
6. OTHER input files

   extract_HC - you need lat_lon file (example is provided);
   predict_tide - you need lat_lon_time file (example is provided) OR
                  lat_lon and time file; sample matlab script mk_lltime.m
                  creating lat_lon_time file, is provided in OTPS/matlab
   
7. COMPILING and RUNNING OTPS

Example of makefile for gfortran compiler is provided with OTPS.
To compile:

make extract_HC
make extract_local_model
make predict_tide

If you need to create executables for other compiler, please provide options
to deal with big endian bytes encoding similar to gfortran options 
"-fconvert=swap -frecord-marker=4". For example, for ifort these 
options are: "-convert big_endian -assume byterecl",
for pgf "-Mbyteswapio".

AFTER compiling successfully and editing Model_* files for your path,
editing setup.inp or setup.local, creating your own input files 
run executable as:

extract_HC<setup.inp                   # have lat_lon file ready
predict_tide<setup.inp                 # have lat_lon_time file ready
extract_local_model<setup.local        # have two Model_*
                                       # control files ready

Examples of setup.local and setup.inp as well as examples of DATA/Model_*
control files are provided with OTPS.

******************************************************************************
Contacts:  Gary Egbert:   Gary.Egbert@oregonstate.edu
           Lana Erofeeva: Svetlana.Erofeeva@oregonstate.edu
                          lana@tpxo.net
