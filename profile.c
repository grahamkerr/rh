/* ------- file: -------------------------- profile.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jan 12 10:23:58 2010 --

       --------------------------                      ----------RH-- */

/* --- Absorption profile for bound-bound transitions.

 Note: Since the transfer code in general allows overlapping
       transitions we cannot use the symmetry property of the Gauss and
       Voigt functions in the case of Zeeman shifts. This is
       because the frequency points are not in general symmetrically
       arranged around zero when lines overlap with other lines or
       continua. Instead we have to store the line profile for both the
       up- and down directions along a ray mu.

 Note: If the option input.limit_memory is set profiles are written
       to file for each line seperately. This prevents huge memory
       allocations in the case of multi-dimensional geometry, at the
       cost of more I/O overhead.

       Naming convention for the profile functions (see for instance
       J. Stenflo 1994, in "Solar Magnetic Fields", p. 108 & 115):

         - phi & psi_pi refer to the unshifted \pi component.
 
        - {phi,psi}_sp refer to the \sigma^+ redshifted component
           that is right-handed circularly polarized (i.e., the
           rotation of the electric vector is clockwise in a fixed
           plane, as seen by the observer).

         - {phi,psi}_sm refer to the sigma^- blueshifted component
           that is left-handed circularly polarized.

       --                                              -------------- */

 
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */

void freeZeeman(ZeemanMultiplet *zm);
double ConvStarkVoigt(AtomicLine *line, int k, double a, double v, double *F, enum VoigtAlgorithm algorithm);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Profile.c --------------- */

void Profile(AtomicLine *line)
{
  const char routineName[] = "Profile";
  register int la, k, mu, n, to_obs, nz;

  char    filename[MAX_LINE_SIZE];
  int     lamu, Nlamu, NrecStokes;
  double *adamp = NULL, **v = NULL, **v_los = NULL, *vB = NULL, *sv = NULL, 
         *vbroad = NULL, Larmor, H, F, wlamu, vk, phi_pi, phi_sm, phi_sp, 
         phi_delta, phi_sigma, psi_pi, psi_sm, psi_sp, psi_delta, psi_sigma, 
         sign, sin2_gamma, *phi = NULL, *phi_Q = NULL, *phi_U = NULL, 
         *phi_V = NULL, *psi_Q = NULL, *psi_U = NULL, *psi_V = NULL;

  Atom *atom = line->atom;
  ZeemanMultiplet *zm = NULL;

  if (line->doVCS_Stark==TRUE){
      InitLemkeStark(line);
    }

  if (!line->Voigt) {
    sprintf(messageStr,
	    "Magnetic lines cannot have GAUSSian profiles. Line %d -> %d",
	    line->j, line->i);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  getCPU(3, TIME_START, NULL);

  /* --- Initialize the ratio of PRD to CRD profiles to 1.0, and
         allocate memory for the elastic collisions in Qelast. This
         array is filled in Damping.

         If Profile() is called from adjustStokes() then the PRD
         profile ratio rho should be kept, and should not be 
         reinitialized. --   
                                   -------------- */
  // printf("\n ||||| inside Profile");
  if (line->PRD && line->rho_prd == NULL) {
    if (input.PRD_angle_dep == PRD_ANGLE_DEP) 
      Nlamu = 2*atmos.Nrays * line->Nlambda;
    else
      Nlamu = line->Nlambda;

    line->rho_prd = matrix_double(Nlamu, atmos.Nspace);
    for (la = 0;  la < Nlamu;  la++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	line->rho_prd[la][k] = 1.0;
    }
    line->Qelast = (double *) malloc(atmos.Nspace * sizeof(double));
  }

  vbroad = atom->vbroad;
  // adamp  = (double *) calloc(atmos.Nspace * sizeof(double));
  adamp  = (double *) calloc(atmos.Nspace, sizeof(double));
  if (line->Voigt) {
    // printf("\n ....... calling Damping for %2s %d -> %d", atom->ID, line->j, line->i);
    Damping(line, adamp);
  }
  line->wphi = (double *) calloc(atmos.Nspace, sizeof(double));

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
    Larmor = (Q_ELECTRON / (4.0*PI*M_ELECTRON)) * (line->lambda0*NM_TO_M);

    zm = Zeeman(line);
    sprintf(messageStr,
	    " -- Atom %2s, line %3d -> %3d has %2d Zeeman components\n",
	    atom->ID, line->j, line->i, zm->Ncomponent);
    Error(MESSAGE, routineName, messageStr);
  }

  /* --- Initialize permanent storage for line profiles -- ---------- */

  if (input.limit_memory) {
    sprintf(filename, (atom->ID[1] == ' ') ?
	    "profile.%.1s_%d-%d.dat" : "profile.%.2s_%d-%d.dat", atom->ID,
	    line->j, line->i);
    if ((line->fd_profile =
	 open(filename, O_RDWR | O_CREAT, PERMISSIONS)) == -1) {
      sprintf(messageStr, "Unable to open profile file %s", filename);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
      NrecStokes = (input.magneto_optical) ? 7 : 4;
      phi = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));

      /* --- Assign pointers to subarrays of phi --    -------------- */

      phi_Q = phi + atmos.Nspace;
      phi_U = phi + 2*atmos.Nspace;
      phi_V = phi + 3*atmos.Nspace;
	      
      if (input.magneto_optical) {
	psi_Q = phi + 4*atmos.Nspace;
	psi_U = phi + 5*atmos.Nspace;
	psi_V = phi + 6*atmos.Nspace;
      }
    } else {
      NrecStokes = 1;
      phi = (double *) malloc(atmos.Nspace * sizeof(double));
    }
  } else {
    if (atmos.moving || 
	(line->polarizable && (input.StokesMode > FIELD_FREE))) {
      Nlamu = 2*atmos.Nrays*line->Nlambda;
      line->phi = matrix_double(Nlamu, atmos.Nspace);

      if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
	line->phi_Q = matrix_double(Nlamu, atmos.Nspace);
	line->phi_U = matrix_double(Nlamu, atmos.Nspace);
	line->phi_V = matrix_double(Nlamu, atmos.Nspace);
	     		      
	if (input.magneto_optical) {
	  line->psi_Q = matrix_double(Nlamu, atmos.Nspace);
	  line->psi_U = matrix_double(Nlamu, atmos.Nspace);
	  line->psi_V = matrix_double(Nlamu, atmos.Nspace);
	}
      }
    } else
      line->phi = matrix_double(line->Nlambda, atmos.Nspace);
  }

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {

    /* --- Temporary storage for inner loop variables, vB is the
           Zeeman splitting due to the local magnetic field -- ------ */  

    vB = (double *) malloc(atmos.Nspace * sizeof(double));
    sv = (double *) malloc(atmos.Nspace * sizeof(double));

    for (k = 0;  k < atmos.Nspace;  k++) {
      vB[k] = Larmor * atmos.B[k] / vbroad[k];
      sv[k] = 1.0 / (SQRTPI * vbroad[k]);
    }
  }

  /* --- Calculate the absorption profile and store for each line -- */

  if (atmos.moving ||
      (line->polarizable && (input.StokesMode > FIELD_FREE))) {

    v_los = matrix_double(atmos.Nrays, atmos.Nspace);
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	v_los[mu][k] = vproject(k, mu) / vbroad[k];
      }
    }
    v = matrix_double(atmos.Nspace, line->Ncomponent);

    for (la = 0;  la < line->Nlambda;  la++) {
      for (n = 0;  n < line->Ncomponent;  n++) {
	for (k = 0;  k < atmos.Nspace;  k++) {
	  v[k][n] = (line->lambda[la] - line->lambda0 - line->c_shift[n]) *
	    CLIGHT / (vbroad[k] * line->lambda0);
	}
      }

      for (mu = 0;  mu < atmos.Nrays;  mu++) {
	wlamu = getwlambda_line(line, la) * 0.5*atmos.wmu[mu];

	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  sign = (to_obs) ? 1.0 : -1.0;
	  lamu = 2*(atmos.Nrays*la + mu) + to_obs;

	  /* --- Assign pointers to the proper phi and psi arrays and
	     zero the profiles in case of conservative memory
	     option. In the normal case the call matrix_double
	     initializes the whole array to zero -- ------------- */

	  if (input.limit_memory) {
	    for (k = 0;  k< NrecStokes*atmos.Nspace;  k++) phi[k] = 0.0;
	  } else {
	    phi = line->phi[lamu];
	    if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
	      phi_Q = line->phi_Q[lamu];
	      phi_U = line->phi_U[lamu];
	      phi_V = line->phi_V[lamu];

	      if (input.magneto_optical) {
		psi_Q = line->psi_Q[lamu];
		psi_U = line->psi_U[lamu];
		psi_V = line->psi_V[lamu];
	      }
	    }
	  }
	    
	  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      sin2_gamma = 1.0 - SQ(atmos.cos_gamma[mu][k]);

	      /* --- For the sign conventions to the phi and psi
		 contributions depending on the direction along the ray

		 See:
		 -- A. van Ballegooijen: "Radiation in Strong Magnetic
		    Fields", in Numerical Radiative Transfer, W. Kalkofen
		    1987, p. 285 --                    -------------- */
            
              /* --- Sum over isotopes --              -------------- */

	      for (n = 0;  n < line->Ncomponent;  n++) {
		vk = v[k][n] + sign * v_los[mu][k];

		phi_sm = phi_pi = phi_sp = 0.0;
		psi_sm = psi_pi = psi_sp = 0.0;

                /* --- Sum over Zeeman sub-levels --   -------------- */

		for (nz = 0;  nz < zm->Ncomponent;  nz++) {
		  H = Voigt(adamp[k], vk - zm->shift[nz]*vB[k],
			    &F, HUMLICEK);

		  switch (zm->q[nz]) {
		  case -1:
		    phi_sm += zm->strength[nz] * H;
		    psi_sm += zm->strength[nz] * F;
		    break;
		  case  0:
		    phi_pi += zm->strength[nz] * H;
		    psi_pi += zm->strength[nz] * F;
		    break;
		  case  1:
		    phi_sp += zm->strength[nz] * H;
		    psi_sp += zm->strength[nz] * F;
		  }
		}
		phi_sigma = (phi_sp + phi_sm) * line->c_fraction[n];
		phi_delta = 0.5*phi_pi * line->c_fraction[n] - 0.25*phi_sigma;

		phi[k]   += (phi_delta*sin2_gamma + 0.5*phi_sigma) * sv[k];
		phi_Q[k] += sign *
		  phi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv[k];
		phi_U[k] +=
		  phi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv[k];
		phi_V[k] += sign *
		  0.5*(phi_sp - phi_sm) * atmos.cos_gamma[mu][k] * sv[k];

		if (input.magneto_optical) {
		  psi_sigma = (psi_sp + psi_sm) * line->c_fraction[n];
		  psi_delta = 0.5*psi_pi * line->c_fraction[n] -
		    0.25*psi_sigma;

		  psi_Q[k] += sign *
		    psi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv[k];
		  psi_U[k] +=
		    psi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv[k];
		  psi_V[k] += sign *
		    0.5 * (psi_sp - psi_sm) * atmos.cos_gamma[mu][k] * sv[k];
		}
	      }
              /* --- Ensure proper normalization of the profile -- -- */

	      line->wphi[k] += wlamu * phi[k];
	    }
	  } else {

	    /* --- Field-free case --                  -------------- */

            for (k = 0;  k < atmos.Nspace;  k++) {
	      for (n = 0;  n < line->Ncomponent;  n++) {
		vk = v[k][n] + sign * v_los[mu][k];
	              if (line->doVCS_Stark) {
                 phi[k] = ConvStarkVoigt(line, k, adamp[k], vk, NULL, ARMSTRONG) *
                   line->c_fraction[n] / (SQRTPI * atom->vbroad[k]);
                } else {
        phi[k] += Voigt(adamp[k], vk, NULL, ARMSTRONG) *
        line->c_fraction[n] / (SQRTPI * atom->vbroad[k]);
        }
      }
        line->wphi[k] += phi[k] * wlamu;
      }
    }
	  if (input.limit_memory) writeProfile(line, lamu, phi);
	}
      }
    }
  } else {

    /* --- Angle-independent case --                   -------------- */   

    for (la = 0;  la < line->Nlambda;  la++) {
      wlamu = getwlambda_line(line, la);
      
      if (input.limit_memory)
	for (k = 0;  k < atmos.Nspace;  k++) phi[k] = 0.0;
      else
	phi = line->phi[la];
      
      for (k = 0;  k < atmos.Nspace;  k++) {
  for (n = 0;  n < line->Ncomponent;  n++) {
    vk = (line->lambda[la] - line->lambda0 - line->c_shift[n]) *
      CLIGHT / (line->lambda0 * atom->vbroad[k]);
    if (line->doVCS_Stark) {
      phi[k] = ConvStarkVoigt(line, k, adamp[k], vk, NULL, ARMSTRONG) * 
        line->c_fraction[n] / (SQRTPI * atom->vbroad[k]);
    } else {
      phi[k] += Voigt(adamp[k], vk, NULL, ARMSTRONG) *
        line->c_fraction[n] / (SQRTPI * atom->vbroad[k]);
    }
  }
  line->wphi[k] += phi[k] * wlamu;
      }
      if (input.limit_memory) writeProfile(line, la, phi);
    }
  }
  /* --- Store the inverse of the profile normalization -- ---------- */

  for (k = 0;  k < atmos.Nspace;  k++) line->wphi[k] = 1.0 / line->wphi[k];

  /* --- Clean up --                                     ------------ */

  free(adamp);
  if (input.limit_memory) free(phi);

  if (atmos.moving || (line->polarizable && (input.StokesMode > FIELD_FREE))) {
    if (zm)
      freeZeeman(zm);
    free(zm);
    free(vB);
    free(sv);

    freeMatrix((void **) v);
    freeMatrix((void **) v_los);

    sprintf(messageStr, "Stokes prof %7.1f", line->lambda0);
  } else
    sprintf(messageStr, "Profile %7.1f", line->lambda0);

  getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- Profile.c --------------- */

/* ------- begin -------------------------- MolecularProfile.c ------ */

void MolecularProfile(MolecularLine *mrt)
{
  register int la, k, mu;

  double *adamp, *v, vk, **v_los, wlamu, *sv, *phi_down, *phi_up,
          wlambda;
  Molecule *molecule = mrt->molecule;

  getCPU(3, TIME_START, NULL);

  mrt->wphi = (double *) calloc(atmos.Nspace, sizeof(double));
  adamp = (double *) calloc(atmos.Nspace, sizeof(double));
  if (mrt->Voigt) MolecularDamping(mrt, adamp);

  if (atmos.moving) {
    v = (double *) malloc(atmos.Nspace * sizeof(double));

    /* --- Store line-of-sight velocity in Doppler units to avoid
           having to recompute it for every wavelength -- ----------- */

    v_los = matrix_double(atmos.Nrays, atmos.Nspace);
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	v_los[mu][k] = vproject(k, mu) / molecule->vbroad[k];
    }
    /* --- Initialize permanent storage for line profiles -- -------- */

    mrt->phi = matrix_double(2*atmos.Nrays * mrt->Nlambda, atmos.Nspace);
  } else
    mrt->phi = matrix_double(mrt->Nlambda, atmos.Nspace);

  /* --- Temporary storage for Doppler width scaling -- ------------- */  

  sv = (double *) malloc(atmos.Nspace * sizeof(double));
  for (k = 0;  k < atmos.Nspace;  k++)
    sv[k] = 1.0 / (SQRTPI * molecule->vbroad[k]);

  /* --- Calculate the absorption profile and store for each line - - */

  if (atmos.moving) {
    for (la = 0;  la < mrt->Nlambda;  la++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	v[k] = (mrt->lambda[la] - mrt->lambda0) * CLIGHT /
	  (molecule->vbroad[k] * mrt->lambda0);

      for (mu = 0;  mu < atmos.Nrays;  mu++) {
	wlamu = getwlambda_mrt(mrt, la) * 0.5*atmos.wmu[mu];

	/* --- First the downward, then the upward direction -- ----- */

	phi_down = mrt->phi[2*(atmos.Nrays*la + mu)];
	phi_up   = mrt->phi[2*(atmos.Nrays*la + mu) + 1];

	for (k = 0;  k < atmos.Nspace;  k++) {
	  phi_down[k] = Voigt(adamp[k], v[k] - v_los[mu][k],
			      NULL, ARMSTRONG) * sv[k];
	  phi_up[k]   = Voigt(adamp[k], v[k] + v_los[mu][k],
			      NULL, ARMSTRONG) * sv[k];
	  mrt->wphi[k] += (phi_down[k] + phi_up[k]) * wlamu;
	}
      }
    }
  } else {
 
    /* --- Stationary case --                          -------------- */   

    for (la = 0;  la < mrt->Nlambda;  la++) {
      wlambda = getwlambda_mrt(mrt, la);

      for (k = 0;  k < atmos.Nspace;  k++) {
	vk = (mrt->lambda[la] - mrt->lambda0) * CLIGHT /
	  (mrt->lambda0 * molecule->vbroad[k]);
	mrt->phi[la][k] = Voigt(adamp[k], vk, NULL, ARMSTRONG) * sv[k];
	mrt->wphi[k] += mrt->phi[la][k] * wlambda;
      }
    }
  }
  for (k = 0;  k < atmos.Nspace;  k++) mrt->wphi[k] = 1.0 / mrt->wphi[k];

  /* --- Clean up --                                     ------------ */

  if (atmos.moving) {
    free(v);
    freeMatrix((void **) v_los);
  }
  free(adamp);  free(sv);

  sprintf(messageStr, "Profile %7.1f", mrt->lambda0);
  getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- MolecularProfile.c ------ */

/* ------- begin -------------------------- getProfiles.c ----------- */

void getProfiles(void)
{
  register int nact, kr;

  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  MolecularLine *mrt;

  /* --- Calculate profiles for Non-LTE after all necessary ingredients
         like Hydrogen poulation fro broadening are available -- ---- */

  getCPU(2, TIME_START, NULL);

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      Profile(line);
    }
  }

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    for (kr = 0;  kr < molecule->Nrt;  kr++) {
      mrt = &molecule->mrt[kr];
      MolecularProfile(mrt);
    }
  }

  getCPU(2, TIME_POLL, "Profiles");
}
/* ------- end ---------------------------- getProfiles.c ----------- */

/* ------- begin -------------------------- InitLemkeStark.c -------- */

#define NALPHA 66
#define NT 8
#define NNE 17
#define NLINES 21
#define NSERIES 4
#define NTOTAL NSERIES*NLINES*NALPHA*NNE*NT
#define NLNS NSERIES*NLINES
#define NWL 1501
void InitLemkeStark(AtomicLine *line)
{
  int i, j, k, series, lne;
  int mp[NSERIES][NLINES], mne[NSERIES][NLINES], nline[NSERIES];
  double svcs[NSERIES][NLINES][NALPHA][NNE][NT], alpha[NALPHA], lgT[NT], lgNe[NNE];
  double log_alpha0[NSERIES][NLINES], log_ne0[NSERIES][NLINES];
  double log_alpha_inc[NSERIES][NLINES], log_ne_inc[NSERIES][NLINES];
  double sum,***arr3d, *tempWL, *tempWLmp, wlmp;
  FILE *fplemkestark;
  const char routineName[] = "InitLemkeStark";

  if (strstr(line->atom->ID, "H ") == NULL) {
     sprintf(messageStr, "Model is not a hydrogen atom: %s", line->atom->ID);
     Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* GSK: Changed this to have the .dat file as an input via keyword input, similar to emistab */
  // if ((fplemkestark = fopen("LemkeStark.dat","rb")) == NULL)
  // {
  //    sprintf(messageStr, "Error opening LemkeStark.dat");
  //    Error(ERROR_LEVEL_2, routineName, messageStr);
  // }
  if ((fplemkestark = fopen(input.lemke_file,"rb")) == NULL)
  {
     sprintf(messageStr, "Error opening LemkeStark file");
     Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  printf("\n>>>> Reading LemkeStark file ('%s') for calculating Hydrogen broadening\n",input.lemke_file);

  if (fread(svcs[0][0][0][0],sizeof(double),NTOTAL,fplemkestark) != NTOTAL 
     || fread(lgT, sizeof(double), NT, fplemkestark) != NT
     || fread(log_alpha0[0], sizeof(double), NLNS, fplemkestark) != NLNS 
     || fread(log_ne0[0], sizeof(double), NLNS, fplemkestark) != NLNS
     || fread(log_alpha_inc[0], sizeof(double), NLNS, fplemkestark) != NLNS
     || fread(log_ne_inc[0], sizeof(double), NLNS, fplemkestark) != NLNS
     || fread(mp[0], sizeof(int), NLNS, fplemkestark) != NLNS
     || fread(mne[0], sizeof(int), NLNS, fplemkestark) != NLNS
     || fread(nline, sizeof(int), NSERIES, fplemkestark) != NSERIES)
  { 
     sprintf(messageStr, "Error reading LemkeStark data");
     Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  fclose(fplemkestark);
  
  
  series = line->i; 
  lne = line->j-series-1;
  if (series >= NSERIES || lne >= NLINES || lne >= nline[series])
  {
    
     sprintf(messageStr, "Error: VCS_STARK requested for unknown line");
     Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  for (i = 0; i< mp[series][lne]; i++){
    // alpha[i] = exp10(log_alpha0[series][lne] + log_alpha_inc[series][lne]*(double)(i));
    alpha[i] = pow((double)10,log_alpha0[series][lne] + log_alpha_inc[series][lne]*(double)(i));
  }
  for (i = 0; i< mne[series][lne]; i++){
    lgNe[i] = log_ne0[series][lne] + log_ne_inc[series][lne]*i;
  }
  
  line->VCS_stark = malloc( sizeof(VCS_Stark));
  line->VCS_stark->DopplerWL = matrix_double(atmos.Nspace, NWL);
  line->VCS_stark->S = matrix_double(atmos.Nspace, NWL-1);
  line->VCS_stark->N = NWL;

  arr3d = malloc(sizeof(double**) * (mp[series][lne]-1));
  for (i = 0; i<mp[series][lne]-1; i++) arr3d[i] = malloc(sizeof(double*) * mne[series][lne]);
  tempWL = malloc(sizeof(double) * mp[series][lne]);
  tempWLmp = malloc(sizeof(double) * (mp[series][lne]-1));
  for (k = 0;  k < atmos.Nspace;  k++)
  {
     tempWL[0] = (1.25e-9 * pow(atmos.ne[k]*1e-6,2./3.) * alpha[0])/10 * CLIGHT / (line->lambda0 * line->atom->vbroad[k]);
     for (i = 1; i < mp[series][lne]; i++) // start at i=1 to throw out the initial svcs which is always 0.0
     {
       tempWL[i]= (1.25e-9 * pow(atmos.ne[k]*1e-6,2./3.) * alpha[i])/10 * CLIGHT / (line->lambda0 * line->atom->vbroad[k]); /* WL in Angstroms then converted to voigt WL units) */
       tempWLmp[i-1] = (tempWL[i] + tempWL[i-1])*.5;
       for (j=0; j<mne[series][lne];j++) arr3d[i-1][j] = svcs[series][lne][i][j];
     }
     line->VCS_stark->DopplerWL[k][0] = tempWL[0];
     for (i = 0; i <NWL-1; i++)
     {
        // line->VCS_stark->DopplerWL[k][i+1] = exp10((double)(i)/(NWL-1) * (log10(tempWL[mp[series][lne]-1]) - log10(tempWL[0])) + log10(tempWL[0]));
        line->VCS_stark->DopplerWL[k][i+1] = pow((double)10,(double)(i)/(NWL-1) * (log10(tempWL[mp[series][lne]-1]) - log10(tempWL[0])) + log10(tempWL[0]));
        wlmp = (line->VCS_stark->DopplerWL[k][i+1] + line->VCS_stark->DopplerWL[k][i])*.5;
        if (wlmp >= tempWLmp[mp[series][lne]-2]) wlmp = tempWLmp[mp[series][lne]-2]*.999;
        // line->VCS_stark->S[k][i] = exp10(TriLinear(mp[series][lne]-1, tempWLmp, wlmp, mne[series][lne], lgNe, log10(atmos.ne[k]*1e-6),NT, lgT, log10(atmos.T[k]), arr3d, TRUE)); /* 1e-6 converts m^-3 to cm^-3 */
        line->VCS_stark->S[k][i] = pow((double)10,TriLinear(mp[series][lne]-1, tempWLmp, wlmp, mne[series][lne], lgNe, log10(atmos.ne[k]*1e-6),NT, lgT, log10(atmos.T[k]), arr3d, TRUE)); /* 1e-6 converts m^-3 to cm^-3 */
     }
     
     /* Normalize to 0.5 */
     sum = 0;
     for (i = 0; i< NWL-1; i++)
     {
       sum += (line->VCS_stark->DopplerWL[k][i+1] - line->VCS_stark->DopplerWL[k][i]) * (line->VCS_stark->S[k][i]);
     }
     for (i=0; i < NWL-1; i++) line->VCS_stark->S[k][i] /= (2*sum);
  }
     k = 160;
     // printf(">>> H %d -> %d: S[%d][750] = %e ; n_elec[%d] = %e ; temp[%d] = %e ; lambda0 = %e ; DopplerWL[750] = %e ; vbroad[%d] = %e \n\n", line->j, line->i, k, line->VCS_stark->S[k][750], k, atmos.ne[k], k, atmos.T[k], line->lambda0, line->VCS_stark->DopplerWL[k][750], k, line->atom->vbroad[k]);

  for (i = 0; i<mp[series][lne]-1; i++) free(arr3d[i]);
  free(arr3d);
  free(tempWL);
  free(tempWLmp);
}

/* ------- end ---------------------------- InitLemkeStark.c -------- */
