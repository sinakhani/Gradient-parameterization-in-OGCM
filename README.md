# Gradient-parameterization-in-OGCM

In this project, we introduce a new subgrid-scale parameterization for representing unresolved volume/eddy transports and momentum fluxes based on a Taylor series expansion of resolved varibales in ocean general circulation models (OGCM). An apriori study using high-resolution MOM6 outputs with 1/16 and 1/8 deg resolutions in channel and Neverworld bathymetry, is performed. Results of different ocean bathymetry and spatial resolutions, from eddy-permitting to non-eddying, are shown.     

In FBC (flat bottom channel in 1/8 deg), CWR (channel with bottom topography in 1/8 deg) and NWR (Neverworld in 1/8 deg) directories, the coarse-grained degree is given as (NUM - 1) / 8, where NUM is the tail number at the end of each directory (for example, Filt_09 refers to a coarse-grained resolution of (9-1)/8 = 1 degree).

In NWR_HR (high resolution 1/16 deg  Neverworld) directory, the coarse-grained degree is (NUM - 1) / 16, where NUM is the tail number at the end of each directory (for example, Filt_25 refers to a coarse-grained resolution of (25-1)/16 = 1.5 degree).

***Please contact sina.khani@austin.utexas.edu for access to high-resolution Neverworld data at 1/16 deg***
