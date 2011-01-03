/* ------- file: -------------------------- writeinputdata_p.c ---------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Thu Dec 30 22:58:15 2010 --

       --------------------------                      -----------RH-- */

/* --- Writes input data (including input, atmos, geom) to output file */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

#define INPUTDATA_FILE "output_indata.ncdf"
#define ARR_STRLEN 20

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern InputData input;
extern char messageStr[];
extern NCDF_Atmos_file infile;
extern MPI_data mpi;
extern IO_data io; 

/* ------- begin --------------------------   init_ncdf_indata.c  --- */
void init_ncdf_indata(void)
/* Creates the netCDF file for the input data */
{
  const char routineName[] = "init_ncdf_indata";
  int    ierror, ncid, ncid_input, ncid_atmos, ncid_mpi, nx_id, ny_id, 
         nspace_id, nhydr_id, nelem_id, nrays_id, nproc_id, naa_id, atstr_id,
         temp_var, ne_var, vz_var, vturb_var, B_var, gB_var, chiB_var, nh_var,
         ew_var, ab_var, eid_var, mu_var, wmu_var, height_var,x_var, y_var,
         xnum_var, ynum_var, tm_var, it_var, conv_var, dm_var, ntsk_var, host_var,
         ft_var, dimids[4];
  bool_t PRD_angle_dep, XRD;
  char   startJ[MAX_LINE_SIZE], StokesMode[MAX_LINE_SIZE], angleSet[MAX_LINE_SIZE];
  FILE  *test;

  /* Check if we can open the file */
  if ((test = fopen(INPUTDATA_FILE, "w")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", INPUTDATA_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file  */
  if ((ierror = nc_create_par(INPUTDATA_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  
  /* Write atmos.ID as global attribute */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  // ideas for more global attributes? (root group)


  /* Create groups */
  if ((ierror = nc_def_grp(ncid, "input", &ncid_input))) ERR(ierror,routineName);
  if ((ierror = nc_def_grp(ncid, "atmos", &ncid_atmos))) ERR(ierror,routineName);
  if ((ierror = nc_def_grp(ncid, "mpi",   &ncid_mpi)))   ERR(ierror,routineName);

  /* --- Definitions for the INPUT group --- */
  /* dimensions */
  // None?

  /* variables*/
  // None?

  /* attributes */
  PRD_angle_dep = (input.PRD_angle_dep  &&  atmos.NPRDactive > 0);
  XRD           = (input.XRD  &&  atmos.NPRDactive > 0);

  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL,"Magneto_optical",NC_UBYTE, 1,
               (unsigned char *) &input.magneto_optical ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "PRD_angle_dep", NC_UBYTE, 1,
               (unsigned char *) &PRD_angle_dep )))    ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "XRD",           NC_UBYTE, 1,
               (unsigned char *) &XRD )))              ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "Background_polarization",
    NC_UBYTE,1,(unsigned char *) &input.backgr_pol ))) ERR(ierror,routineName);
  
  switch (input.startJ) {
  case UNKNOWN:
    strcpy(startJ, "Unknown");
    break;
  case LTE_POPULATIONS:
    strcpy(startJ, "LTE_POPULATIONS");
    break;
  case ZERO_RADIATION:
    strcpy(startJ, "ZERO_RADIATION");
    break;
  case OLD_POPULATIONS:
    strcpy(startJ, "OLD_POPULATIONS");
    break;
  case NEW_J:
    strcpy(startJ, "NEW_J");
    break;
  case OLD_J:
    strcpy(startJ, "OLD_J");
    break;
  }
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL, "Start_J", strlen(startJ),
			      startJ )))  ERR(ierror,routineName);

  switch (input.StokesMode) {
  case NO_STOKES:
    strcpy(StokesMode, "NO_STOKES");
    break;
  case FIELD_FREE:
    strcpy(StokesMode, "FIELD_FREE");
    break;
  case POLARIZATION_FREE:
    strcpy(StokesMode, "POLARIZATION_FREE");
    break;
  case FULL_STOKES:
    strcpy(StokesMode, "FULL_STOKES");
    break;
  }
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Stokes_mode",strlen(StokesMode),
			      StokesMode )))  ERR(ierror,routineName);

  switch (atmos.angleSet.set) {
  case SET_VERTICAL:
    strcpy(angleSet, "SET_VERTICAL");
    break;
  case SET_GL:
    strcpy(angleSet, "SET_GL");
    break;
  case SET_A2:
    strcpy(angleSet, "SET_A2");
    break;
  case SET_A4:
    strcpy(angleSet, "SET_A4");
    break;
  case SET_A6:
    strcpy(angleSet, "SET_A6");
    break;
  case SET_A8:
    strcpy(angleSet, "SET_A8");
    break;
  case SET_B4:
    strcpy(angleSet, "SET_B4");
    break;
  case SET_B6:
    strcpy(angleSet, "SET_B6");
    break;
  case SET_B8:
    strcpy(angleSet, "SET_B8");
    break;
  case NO_SET:
    strcpy(angleSet, "NO_SET");
    break;
  }
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Angle_set",strlen(angleSet),
			      angleSet )))  ERR(ierror,routineName);

  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Atmos_file",
       strlen(input.atmos_input), input.atmos_input))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Abundances_file",
       strlen(input.abund_input), input.abund_input))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Kurucz_PF_data",
       strlen(input.pfData),      input.KuruczData)))  ERR(ierror,routineName);

  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "N_max_iter",    NC_INT, 1,
			      &input.NmaxIter )))     ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_delay",      NC_INT, 1,
			      &input.Ngdelay )))      ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_order",      NC_INT, 1,
			      &input.Ngorder )))      ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_period",     NC_INT, 1,
			      &input.Ngperiod )))     ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_N_max_iter",NC_INT, 1,
			      &input.PRD_NmaxIter ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_delay",  NC_INT, 1,
			      &input.PRD_Ngdelay )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_order",  NC_INT, 1,
			      &input.PRD_Ngorder )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_period", NC_INT, 1,
			      &input.PRD_Ngperiod ))) ERR(ierror,routineName);

  if ((ierror=nc_put_att_double(ncid_input, NC_GLOBAL, "Metallicity", NC_DOUBLE, 1,
				&input.metallicity ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_double(ncid_input, NC_GLOBAL, "Lambda_reference",
                   NC_DOUBLE, 1, &atmos.lambda_ref ))) ERR(ierror,routineName);

  // Need to write names of output files


  /* --- Definitions for the ATMOS group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid_atmos, "nx",        mpi.nx,          &nx_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "ny",        mpi.ny,          &ny_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nz",        atmos.Nspace,    &nspace_id))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nhydr",     atmos.H->Nlevel, &nhydr_id )))
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nelements", atmos.Nelem,     &nelem_id ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nrays",     geometry.Nrays,  &nrays_id ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "strlen",   ATOM_ID_WIDTH+1, &atstr_id ))) 
    ERR(ierror,routineName);

  /* variables*/
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspace_id;
  if ((ierror = nc_def_var(ncid_atmos, "temperature",        NC_FLOAT, 3, dimids,
			   &temp_var)))  ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_atmos, "electron_density",   NC_DOUBLE, 3, dimids,
			   &ne_var)))    ERR(ierror,routineName); 
  if ((ierror = nc_def_var(ncid_atmos, "velocity_z",         NC_FLOAT, 3, dimids,
			   &vz_var)))    ERR(ierror,routineName); 
  if ((ierror = nc_def_var(ncid_atmos, "velocity_turbulent", NC_FLOAT, 3, dimids,
			   &vturb_var))) ERR(ierror,routineName); 
  if (atmos.Stokes) {
    if ((ierror=nc_def_var(ncid_atmos, "B",                  NC_FLOAT, 3, dimids,
			   &B_var)))     ERR(ierror,routineName); 
    if ((ierror=nc_def_var(ncid_atmos, "gamma_B",            NC_FLOAT, 3, dimids,
			   &gB_var)))    ERR(ierror,routineName); 
    if ((ierror=nc_def_var(ncid_atmos, "chi_B",              NC_FLOAT, 3, dimids,
			   &chiB_var)))  ERR(ierror,routineName); 
  }
  dimids[0] = nhydr_id;
  dimids[1] = nx_id;
  dimids[2] = ny_id;
  dimids[3] = nspace_id; 
  if ((ierror = nc_def_var(ncid_atmos, "hydrogen_populations",NC_FLOAT, 4, dimids,
			   &nh_var)))    ERR(ierror,routineName); 
  dimids[0] = nelem_id;
  if ((ierror = nc_def_var(ncid_atmos, "element_weight",      NC_DOUBLE,1, dimids,
			   &ew_var)))    ERR(ierror,routineName);  
  if ((ierror = nc_def_var(ncid_atmos, "element_abundance",   NC_DOUBLE,1, dimids,
			   &ab_var)))    ERR(ierror,routineName); 
  dimids[1] = atstr_id;
  if ((ierror = nc_def_var(ncid_atmos, "element_id",          NC_CHAR,  2, dimids,
                           &eid_var)))   ERR(ierror,routineName); 
  /* next 3 came from geometry */
  dimids[0] = nrays_id;
  if ((ierror = nc_def_var(ncid_atmos, "muz",                 NC_DOUBLE,1, dimids,
			   &mu_var)))    ERR(ierror,routineName);  
  if ((ierror = nc_def_var(ncid_atmos, "wmu",                 NC_DOUBLE,1, dimids,
			   &wmu_var)))   ERR(ierror,routineName);  
  dimids[0] = nspace_id;
  if ((ierror = nc_def_var(ncid_atmos, "height",              NC_FLOAT, 1, dimids,
			   &height_var))) ERR(ierror,routineName);
  dimids[0] = nx_id;
  if ((ierror = nc_def_var(ncid_atmos, "x",                   NC_FLOAT, 1, dimids,
			   &x_var)))      ERR(ierror,routineName);
  dimids[0] = ny_id;
  if ((ierror = nc_def_var(ncid_atmos, "y",                   NC_FLOAT, 1, dimids,
			   &y_var)))      ERR(ierror,routineName);

  /* attributes */
  if ((ierror = nc_put_att_ubyte(ncid_atmos, NC_GLOBAL, "moving", NC_UBYTE, 1,
		     (unsigned char *) &atmos.moving ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_ubyte(ncid_atmos, NC_GLOBAL, "stokes", NC_UBYTE, 1,
	             (unsigned char *) &atmos.Stokes ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, temp_var,  "units",  1,
				 "K" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, ne_var,    "units",  4,
				 "m^-3" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, vz_var,    "units",  6,
				 "m s^-1" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, vturb_var, "units",  6,
				 "m s^-1" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, nh_var,    "units",  4,
				 "m^-3" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, ew_var,    "units", 17,
				 "atomic_mass_units" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, height_var,"units",  1,
				 "m" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, x_var,     "units",  1,
				 "m" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, y_var,     "units",  1,
				 "m" ))) ERR(ierror,routineName);
  
  /* --- Definitions for the MPI group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid_mpi, "nx",          mpi.nx,   &nx_id   ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_mpi, "ny",          mpi.ny,   &ny_id   ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_mpi, "nprocesses",  mpi.size, &nproc_id))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_def_dim(ncid_mpi, "nactive_atom",atmos.Nactiveatom,&naa_id)))
    ERR(ierror,routineName);  
  if ((ierror = nc_def_dim(ncid_mpi, "strlen",   ARR_STRLEN, &atstr_id ))) 
    ERR(ierror,routineName);

  /* variables*/
  dimids[0] = nx_id;
  if ((ierror = nc_def_var(ncid_mpi, "xnum",        NC_USHORT, 1, dimids,
			   &xnum_var))) ERR(ierror,routineName);
  dimids[0] = ny_id;
  if ((ierror = nc_def_var(ncid_mpi, "ynum",        NC_USHORT, 1, dimids,
			   &ynum_var))) ERR(ierror,routineName);
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  /* this one is a matrix with floats: XX.YY, where XX = mpi.rank, YY=task number */
  if ((ierror = nc_def_var(ncid_mpi, "taskmap",     NC_FLOAT,  2, dimids,
			   &tm_var)))   ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, "iterations",  NC_USHORT, 2, dimids,
			   &it_var)))   ERR(ierror,routineName);

  dimids[0] = naa_id;
  dimids[1] = nx_id;
  dimids[2] = ny_id;
  if ((ierror = nc_def_var(ncid_mpi, "convergence", NC_UBYTE,  3, dimids,
			   &conv_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, "delta_max",   NC_FLOAT,  3, dimids,
			   &dm_var)))   ERR(ierror,routineName);
  dimids[0] = nproc_id;
  dimids[1] = atstr_id;
  if ((ierror = nc_def_var(ncid_mpi, "ntasks",      NC_LONG,   1, dimids,
			   &ntsk_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, "hostname",    NC_CHAR,   2, dimids,
			   &host_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, "finish_time", NC_CHAR,   2, dimids,
			   &ft_var))) ERR(ierror,routineName);

  // mpi.tasks(size), mpi.taskmap (but in xnum, ynum terms!), mpi.taskmap for each process
  // number of iterations for each mpi.ix, mpi.iy, and also if they converged, convergence
  // limit.

  /* attributes */
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_start", NC_INT, 1,
			      &input.p15d_x0 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_end",   NC_INT, 1,
			      &input.p15d_x1 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_step",  NC_INT, 1,
			      &input.p15d_xst ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_start", NC_INT, 1,
			      &input.p15d_y0 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_end",   NC_INT, 1,
			      &input.p15d_y1 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_step",  NC_INT, 1,
			      &input.p15d_yst ))) ERR(ierror,routineName);

  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);

  /* Copy stuff to the IO data struct */
  io.in_ncid       = ncid;
  io.in_input_ncid = ncid_input;
  io.in_atmos_ncid = ncid_atmos;
  io.in_mpi_ncid   = ncid_mpi;
  
  io.in_atmos_T    = temp_var;
  io.in_atmos_ne   = ne_var;
  io.in_atmos_vz   = vz_var;
  io.in_atmos_vt   = vturb_var;
  io.in_atmos_B    = B_var;
  io.in_atmos_gB   = gB_var;
  io.in_atmos_chiB = chiB_var;
  io.in_atmos_nh   = nh_var;
  io.in_atmos_ew   = ew_var;
  io.in_atmos_ab   = ab_var;
  io.in_atmos_eid  = eid_var;
  io.in_atmos_mu   = mu_var;
  io.in_atmos_wmu  = wmu_var;
  io.in_atmos_z    = height_var;
  io.in_atmos_x    = x_var;
  io.in_atmos_y    = y_var;
  
  io.in_mpi_xnum   = xnum_var;
  io.in_mpi_ynum   = ynum_var;
  io.in_mpi_tm     = tm_var;
  io.in_mpi_it     = it_var;
  io.in_mpi_conv   = conv_var;
  io.in_mpi_dm     = dm_var;
  io.in_mpi_ntsk   = ntsk_var;
  io.in_mpi_host   = host_var;


  return;
}
/* ------- end   --------------------------   init_ncdf_indata.c  --- */

/* ------- begin --------------------------   close_ncdf_indata.c --- */
void close_ncdf_indata(void)
/* Closes the spec netCDF file */ 
{
  const char routineName[] = "close_ncdf_indata";
  int        ierror;

  if ((ierror = nc_close(io.in_ncid))) ERR(ierror,routineName);
  return; 
}
/* ------- end   --------------------------   close_ncdf_indata.c --- */

/* ------- begin --------------------------   writeAtmos_p.c --- */
void writeAtmos_p(void)
{
  const char routineName[] = "writeAtmos_p";
  int     ierror, ncid, i;
  size_t  start[] = {0, 0, 0, 0};
  size_t  count[] = {1, 1, 1, 1};

  ncid = io.in_atmos_ncid;

  /* put atmosphere variables */
  start[0] = mpi.ix; 
  start[1] = mpi.iy;
  count[2] = atmos.Nspace;


  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_T, start, count,
				   atmos.T ))) ERR(ierror,routineName);
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_ne, start, count,
				   atmos.ne ))) ERR(ierror,routineName);
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_vz, start, count,
				   geometry.vel ))) ERR(ierror,routineName);
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_vt, start, count,
				   atmos.vturb ))) ERR(ierror,routineName);  
  if (atmos.Stokes) {
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_B, start, count,
				     atmos.B ))) ERR(ierror,routineName); 
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_gB, start, count,
				     atmos.gamma_B ))) ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_chiB, start, count,
				     atmos.chi_B ))) ERR(ierror,routineName);
  }

  /* put hydrogen populations */
  start[1] = mpi.ix;
  start[2] = mpi.iy;
  count[2] = 1;
  count[3] = atmos.Nspace;

  for (i=0; i < atmos.H->Nlevel; i++) {
    start[0] = (size_t) i;
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_nh, start, count,
				     atmos.H->n[i] ))) ERR(ierror,routineName);  
  }

  /* arrays of number of elements */
  start[1] = 0;
  count[1] = ATOM_ID_WIDTH+1;
  for (i=0; i < atmos.Nelem; i++) {
    start[0] = (size_t) i;
    if ((ierror = nc_put_var1_double(ncid, io.in_atmos_ew,  &start[0], 
		     &atmos.elements[i].weight ))) ERR(ierror,routineName);  
    if ((ierror = nc_put_var1_double(ncid, io.in_atmos_ab,  &start[0], 
		     &atmos.elements[i].abund )))  ERR(ierror,routineName); 
    if ((ierror = nc_put_vara_text(ncid, io.in_atmos_eid, start, count,
       (const char *)&atmos.elements[i].ID )))     ERR(ierror,routineName); 
  }

  /* arrays from geometry */
  if ((ierror = nc_put_var_double(ncid, io.in_atmos_mu,  geometry.muz )))
    ERR(ierror,routineName);   
  if ((ierror = nc_put_var_double(ncid, io.in_atmos_wmu, geometry.wmu )))
    ERR(ierror,routineName);   
  if ((ierror = nc_put_var_double(ncid, io.in_atmos_z,   geometry.height )))
    ERR(ierror,routineName);   

  /* write x and y, convert from MPI indices to actual values */
  for (i=0; i < mpi.nx; i++) {
    start[0] = (size_t) i;
    if ((ierror = nc_put_var1_double(ncid, io.in_atmos_x,  &start[0], 
	     &infile.x[mpi.xnum[i]]))) ERR(ierror,routineName);  
  }

  for (i=0; i < mpi.ny; i++) {
    start[0] = (size_t) i;
    if ((ierror = nc_put_var1_double(ncid, io.in_atmos_y,  &start[0], 
	     &infile.y[mpi.ynum[i]]))) ERR(ierror,routineName);  
  }


  return;
}
/* ------- end   --------------------------   writeAtmos_p.c --- */

/* ------- begin --------------------------   writeMPI_p.c --- */
/* ------- end   --------------------------   writeMPI_p.c --- */

// Now needed: 
// * routine to write indata.atmos stuff;
// * routine to write indata.mpi   stuff;
