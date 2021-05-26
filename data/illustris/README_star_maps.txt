-------------------------------------------------------------------------------------
HDF5 ILLUSTRIS DATA FILE
-------------------------------------------------------------------------------------

The following table describes the content of the sfh_*.hdf5 files. A few comments:

- Unit system
  Comoving coordinates, referring to position in box     kpc / h (comoving)
  Halo masses                                            Msun / h
  All other coordinates (e.g. within galaxy)             kpc (physical)
  All other masses                                       Msun
  Times                                                  Gyr
  Star formation rate                                    Msun / yr
  Black hole accretion rate                              Msun / yr
  All velocities                                         km / s
  
- Cosmology
  All calculations are performed assuming the Illustris cosmology, with H0 = 70.4,
  Om = 0.2726, sigma8 = 0.8090.

- Times
  All times are measured since the big bang, i.e. the beginning of the universe 
  is at t = 0 and the time today is t = 13.75. 

- Halos / subhalos / groups
  In the Illustris catalogs, “groups” refer to friends-of-friends groups, and subhalos
  to their sub-groups. Here, each subhalo corresponds to a galaxy that belongs to a 
  group (which may have more than one subhalo in it). All fields denotet cat_sh_*
  below refer to properties of the subhalo, fields denoted cat_grp_* refer to the 
  group.

- Dimensions of fields
  In the table below, nh is the number of halos (or galaxies), and nt is the number
  of time bins. These bins are not tied to snapshots, because stellar particles remember
  their exact birth time independent of snapshot files. The properties extracted along
  the main branch of the tree, however, are tied to snapshots, and the number of 
  snapshots extracted (counting back from the last snapshot) is nsnaps. The ordering 
  in the tree arrays is time-forward though. Currently, the dimensions are set to:
  
  nt     = 100
  ntree  = 98  (going back to snapshot 38, z = 8)
  nenv   = 7

- Metallicity: the metallicity is given in absolute units, i.e. mass_metals / mass_total
  as opposed to solar metallicity.


--------------------------------------------------------------------------------------------------------------------------
FIELD                            DIMENSIONS          UNITS         EXPLANATION
--------------------------------------------------------------------------------------------------------------------------
cat_grp_id                       nh                  -             Group ID 
cat_grp_m200c                    nh                  Msun/h        Group mass within R200c
cat_grp_n                        nh                  -             The number of subhalos in this SH’s group
cat_grp_is_primary               nh                  -             Is this the primary subhalo?
cat_grp_is_primary_at_peak       nh                  -             Was this the primary subhalo at peak SFR?
--------------------------------------------------------------------------------------------------------------------------
cat_sh_id                        nh                  -             Subhalo ID
cat_sh_mdm                       nh                  Msun          Mass in dark matter (bound)
cat_sh_mgas                      nh                  Msun          Mass in gas (bound)
cat_sh_mstar                     nh                  Msun          Mass in stars according to catalog (bound)
cat_sh_sfr                       nh                  Msun/yr       Star formation rate according to catalog
cat_sh_mbh                       nh                  Msun          Mass in black holes
cat_sh_bh_mdot                   nh                  Msun/yr       Mass accretion rate of the black hole
cat_sh_halfmrad_stars            nh                  kpc           The stellar half-mass radius
cat_sh_metallicity               nh                  -             The metallicity 
cat_sh_vmax                      nh                  km/s          Maximum circular velocity
cat_sh_phot_U                    nh                  -             Magnitude in U-band       
cat_sh_phot_B                    nh                  -             Magnitude in B-band       
cat_sh_phot_V                    nh                  -             Magnitude in V-band       
cat_sh_phot_K                    nh                  -             Magnitude in K-band       
cat_sh_phot_g                    nh                  -             Magnitude in g-band       
cat_sh_phot_r                    nh                  -             Magnitude in r-band       
cat_sh_phot_i                    nh                  -             Magnitude in i-band       
cat_sh_phot_z                    nh                  -             Magnitude in z-band       
--------------------------------------------------------------------------------------------------------------------------
map_stars                        nh*3*3*npix*npix    Msun/kpc^2    Surf. dens. of M* for SH/group/fuzz, xy-xz-yz proj.
map_stars_nptl                   nh * 3              -             Number of stellar ptl in SH/group/fuzz samples
map_stars_npixel                 1                   -             The number of pixels in the stellar maps
map_stars_size                   1                   kpc           The physical side length of the stellar maps
--------------------------------------------------------------------------------------------------------------------------


