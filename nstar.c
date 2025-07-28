#include "pluto.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#ifndef ORBIT_INTEGRATION
  #define ORBIT_INTEGRATION YES
#endif

#ifndef UNIT_TIME
  #define  UNIT_TIME   (UNIT_LENGTH/UNIT_VELOCITY)
#endif

#define ORBIT_PARALLEL YES
#define STARSEARCH  YES

static double rk_timestep=1.e-4;
static double eps_abs=1e-6, eps_rel=1e-6;
static double h_init = 1e-8;
Nstar g_nstar;
static int nprocs;

#if NSTAR == SCHWAR
  static int opinit = 0;
  static const int glorder=100;
  static double toffset=0.0;
  // static const gsl_odeiv2_step_type * integrate_stepper = gsl_odeiv2_step_rk8pd;  
  OrbitParams orbitparam;
#elif NSTAR == FOKPLA
  static int opinit = 0;
  static double toffset=0.0;
  FPParams orbitparam; /* orbitparam is the same */
#endif

static double const_g = 6.6743e-08, const_msun = 1.988409870698051e+33;
static double UNIT_GM = 1.3271244e+26/(UNIT_VELOCITY*UNIT_VELOCITY*UNIT_LENGTH);

void SetupNstar (Nstar *ns, Grid *grid)
/* . Read orbit.ini and initialize Nstar stucture. 
   . Called in main.c                       
   . 2022 Apr 27 Zhao Su  */
/* ********************** */
{
    FILE *fp;
    char orbit_file[256], line[256], *token;
    int nline, ch, i, j, nstar;
    double lifetime;
    double *orb;

    /* 1. Read orbital initial file */
    #if ORBIT_INTEGRATION
      printLog("======= USE COORD&VELOCITY PHASE AS INITIAL CONDITION =======\n");
    #endif
    #if NSTAR == SIMPLE
      sprintf(orbit_file, "%s/orbit.ini", RuntimeGet()->output_dir);
      fp = fopen(orbit_file, "r");
      if (fp == NULL){
          printLog ("!nstar.c: reading orbit.ini fails!\n");
          QUIT_PLUTO(1);
      }
      nline=0;
      while (fgets(line, 256, fp) != NULL){
        nline++;
      }
      ns->nstar = nline;
      ns->lifetime = RuntimeGet()->tstop; //just set lifetime to simulation duration
      fclose(fp);
      nstar = ns->nstar;
    #elif NSTAR == SCHWAR || NSTAR == FOKPLA
      printLog ("!nstar.c: schwar.ini->");
      sprintf(orbit_file, "%s/schwar.ini", RuntimeGet()->output_dir);
      fp = fopen(orbit_file, "r");
      if (fp == NULL){
          printLog ("!nstar.c: reading schwar.ini fails!\n");
          QUIT_PLUTO(1);
      }
      fscanf(fp, "%d %d %lf\n", &nline, &nstar, &lifetime); //first line
      ns->nstar = nstar;
      ns->lifetime = lifetime;
      printLog ("%d %d %lf\n", nline, nstar, lifetime);
      fclose(fp);
    #endif

    /* 1.1 initialize for Nstar */
    #if ORBIT_INTEGRATION
      ns->time = 0.0;
    #endif
    ns->evotime = 0.0;
    /* . stellar wind and accretion radius initialization
       . ns->rstar is half length of one cell for AMR 
       . ns->racc is half of diagonal length for level 0 grid cell */
    #ifdef CHOMBO
      ns->rstar = 1./2.*grid->dl_min[IDIR]*pow(2., (double) grid->level);
      ns->racc = 1.2*sqrt((double) DIMENSIONS)*ns->rstar; //1.2 is just a approximation
    #else
      ns->rstar = g_inputParam[WIND_RADIUS]*grid->dl_min[IDIR];
      #ifdef ACC_RADIUS
        ns->racc = g_inputParam[ACC_RADIUS]*grid->dl_min[IDIR]; //setup accretion radius
      #else
        ns->racc = 0.6*sqrt((double) DIMENSIONS)*grid->dl_min[IDIR]; //use inner 2^dim cells if accretion radius not set up. 1.2 is just a approximation
      #endif
    #endif
    printLog ("!nstar.c: accretion radius %12.6e. minimum stellar radius %12.6e.\n", ns->racc, ns->rstar);

    ns->coord = ARRAY_2D(nstar, 3, double);    
    ns->v = ARRAY_2D(nstar, 3, double);    
    ns->phase = ARRAY_2D(nstar, 6, double);
    ns->r = ARRAY_1D(nstar, double);  
    #if NSTAR == SIMPLE
      ns->orbit = ARRAY_2D(nstar, 6, double);
      ns->P = ARRAY_2D(nstar, 3, double);
      ns->Q = ARRAY_2D(nstar, 3, double);
      ns->omega = ARRAY_1D(nstar, double);
      ns->M = ARRAY_1D(nstar, double);
      ns->E = ARRAY_1D(nstar, double);

      /* ************************ */
      printLog ("> NSTAR: start reading orbit initial file...\n");
      printLog ("-------------------------------\n");
      #if ORBIT_INTEGRATION
        printLog ("X    Y    Z    VX    VY    VZ\n");
      #else
        printLog ("SMA   ECC   INC   ASC   AP   MA\n");
      #endif

      fp = fopen(orbit_file, "r");
      nline = 0;
      while (fgets(line, 256, fp) != NULL){
        j=0;
        token = strtok(line, " ");
        while (token != NULL){
          ns->orbit[nline][j] = atof(token);
          printLog ("%.3f  ", atof(token));
          token = strtok(NULL, " ");
          j++;
        }
        nline++;
        printLog ("\n");
      }
      printLog ("--------star number: %i---------\n", ns->nstar);
      fclose(fp);
      /* ************************ */

      /* 2. Calculate coordinates from orbital parameters */
      for (j=0;j<ns->nstar;j++){
        #if ORBIT_INTEGRATION == NO
          #if DIMENSIONS == 2
          ns->orbit[j][INC] = 0.;
          ns->orbit[j][ASC] = 0.;
          #endif
          orb = ns->orbit[j];
          ns->P[j][IDIR] = cos(orb[ASC])*cos(orb[AP]) - sin(orb[ASC])*sin(orb[AP])*cos(orb[INC]);
          ns->P[j][JDIR] = sin(orb[ASC])*cos(orb[AP]) + cos(orb[ASC])*sin(orb[AP])*cos(orb[INC]);
          ns->P[j][KDIR] = sin(orb[AP])*sin(orb[INC]);
          ns->Q[j][IDIR] = -cos(orb[ASC])*sin(orb[AP]) - sin(orb[ASC])*cos(orb[AP])*cos(orb[INC]);
          ns->Q[j][JDIR] = -sin(orb[ASC])*sin(orb[AP]) + cos(orb[ASC])*cos(orb[AP])*cos(orb[INC]);
          ns->Q[j][KDIR] = cos(orb[AP])*sin(orb[INC]);
          ns->omega[j] = sqrt(g_inputParam[G_CONS]*g_inputParam[M_BH]/pow(orb[SMA], 3.));
          ns->M[j] = orb[MA] + ns->omega[j]*g_time;
          ns->E[j] = KeplerSolver(ns->M[j], orb[ECC]);
          // printLog ("P: %12.6e %12.6e %12.6e \n Q: %12.6e %12.6e %12.6e \n %12.6e %12.6e %12.6e \n",
          //             ns->P[j][IDIR], ns->P[j][KDIR], ns->P[j][KDIR], 
          //             ns->Q[j][IDIR], ns->Q[j][KDIR], ns->Q[j][KDIR],
          //             ns->omega[j], ns->M[j], ns->E[j]
          //             );
          for (i=0;i<3;i++){
            ns->coord[j][i] = orb[SMA]*(cos(ns->E[j]) - orb[ECC])*ns->P[j][i]
                             + orb[SMA]*sqrt(1-orb[ECC]*orb[ECC])*sin(ns->E[j])*ns->Q[j][i];
          }
          ns->r[j] = sqrt(DIM_EXPAND(ns->coord[j][IDIR]*ns->coord[j][IDIR], 
                                      +ns->coord[j][JDIR]*ns->coord[j][JDIR],
                                      +ns->coord[j][KDIR]*ns->coord[j][KDIR]));
          for (i=0;i<3;i++){
            ns->v[j][i] = orb[SMA]*orb[SMA]*ns->omega[j]/ns->r[j]*(
              -sin(ns->E[j])*ns->P[j][i] + sqrt(1-orb[ECC]*orb[ECC])*cos(ns->E[j])*ns->Q[j][i]
            );
          }
        #else 
        /* For ORBIT_INTEGRATION case, just setup coordinate and velocity. 
         * The orbit initialization data file contains xyz (column 123) and Vxyz (column 456). */
          orb = ns->orbit[j];
          for (i=0;i<3;i++){
            ns->coord[j][i] = orb[i];
            ns->v[j][i] = orb[i+3];
            ns->phase[j][i] = orb[i];
            ns->phase[j][i+3] = orb[i+3];
            }
      }
      UpdateNstar (ns, grid);
      ns->time = g_time;
      #endif
      
    #elif NSTAR == SCHWAR || NSTAR == FOKPLA
      ns->Orbtype = ARRAY_1D(nstar, int);
      ns->Orbindex = ARRAY_1D(nstar, int);
      #ifdef OFFSET_TIME
        toffset = g_inputParam[OFFSET_TIME];
        printLog ("!nstar.c: orbit offset time %12.6e.\n", toffset);
      #endif
      ns->start = (int) ((int) ((g_time+toffset)/ns->lifetime) * ns->nstar);
      LoadSchwarzschildStars(ns->start+1, nstar, ns, 0);
      ns->evotime = fmod((g_time+toffset), ns->lifetime);
      ns->time = (double) ((int) ((g_time+toffset)/ns->lifetime) * ns->lifetime);
      OrbitParamInit();
      if (opinit == 1) printLog ("!nstar.c: orbitparam initialization succeed.\n");
    #endif
      
    
      /* Nstar Structure initialized */
      ns->init = 1;
      #ifdef PARALLEL
        printLog ("!nstar.c: MPI_Comm_size: ");
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        printLog ("%d.\n", nprocs);
        printLog ("!nstar.c: STARSEARCH %d. ORBIT_PARALLEL %d.\n", STARSEARCH, ORBIT_PARALLEL);
      #endif
      /* ************************* */
}

void UpdateNstar (Nstar *ns, Grid *grid)
/* . Update Nstar structure */
{   
    // printLog ("!nstar.c: UpdateNstar called.\n");
    // SIMPLE condition for Keperian orbit or simple orbit integration
    #if NSTAR == SIMPLE
      int i, n;
      double omega, *orb;
      #if ORBIT_INTEGRATION == YES
        for (n=0;n<ns->nstar;n++){
            double tstop, tstep, time;
            double *y, yout[6], yerr[6];

            time = ns->time;
            y = ns->phase[n];
            //printLog ("g_time:%.6e time:%.6e", g_time, time);
            /* Update to g_time */
            while (time < g_time){
                if (g_time >= time+rk_timestep) tstep = rk_timestep;
                else tstep = g_time-time;
                //printLog ("integration step: %.6e \n", tstep);
                rkck(y, time, tstep, y, yerr, HamiltonianDerivs);
                time += tstep;
            }
            
            for (i=0;i<3;i++){
                ns->coord[n][i] = ns->phase[n][i];
                ns->v[n][i] = ns->phase[n][i+3];
            }
        }
        ns->time = g_time;
      #else //Kepler orbit case
        //printLog ("!nstar.c: UpdateNstar called.\n");
        for (n=0;n<ns->nstar;n++){
            orb = ns->orbit[n];
            omega = ns->omega[n];
            //Update for mean anomaly, eccentric anomaly, coordinates, and velocity.
            ns->M[n] = orb[MA] + omega*g_time;
            ns->E[n] = KeplerSolver(ns->M[n], orb[ECC]);
            for (i=0;i<3;i++){
              ns->coord[n][i] = orb[SMA]*(cos(ns->E[n]) - orb[ECC])*ns->P[n][i] 
                              + orb[SMA]*sqrt(1-orb[ECC]*orb[ECC])*sin(ns->E[n])*ns->Q[n][i];
            }
            /* AMR case: move star to cell center */
            #ifdef CHOMBO
                int ilo, imid, ihi, idim;
                DIM_LOOP (idim) {
                ilo = 0;
                ihi = grid->np_int_glob[idim];
                while (ilo != ihi-1){
                  imid = (ilo + ihi)/2;
                  xmid = g_domBeg[idim] + imid*grid->dl_min[idim]; //uniform grid and Cartesian geometry
                  if (xmid <= ns->coord[n][idim]) {
                    ilo = imid;
                  }
                  else if (xmid > x_star) {
                    ihi = imid;
                  }
                }
                ns->coord[n][idim] = g_domBeg[idim] + (0.5 + ilo)*grid->dl_min[idim];
                }
          #endif
          /* **************************** */
          ns->r[n] = sqrt(DIM_EXPAND(ns->coord[n][IDIR]*ns->coord[n][IDIR], +ns->coord[n][JDIR]*ns->coord[n][JDIR], +ns->coord[n][KDIR]*ns->coord[n][KDIR]));
          for (i=0;i<3;i++){
            ns->v[n][i] = orb[SMA]*orb[SMA]*omega/ns->r[n]*(
              -sin(ns->E[n])*ns->P[n][i] + sqrt(1-orb[ECC]*orb[ECC])*cos(ns->E[n])*ns->Q[n][i]
            );
          }
          //printLog ("!nstar.c: star %i: x=%12.6e, y=%12.6e, vx=%12.6e, vy=%12.6e\n", n, ns->coord[n][IDIR],
                  //ns->coord[n][JDIR], ns->v[n][IDIR], ns->v[n][JDIR]);
        }
        ns->time = g_time;
      #endif //end of interation and Keplerian for SIMPLE

    //case 1: SCHWAR - Schwarzschild orbits
    //case 2: FOKPLA - Fokker-Planck isotropic model
    #elif NSTAR == SCHWAR || NSTAR == FOKPLA
      int i, n, j;
      static int first_call = 1, ndriver, start_idx;
      static gsl_odeiv2_system sys = {orbit_ode,
                                      NULL,
                                      6,
                                      &orbitparam
                                     };
      static gsl_odeiv2_driver **driver;
      #if ORBIT_PARALLEL == YES
        static int *recvcounts, *displs;
        static double **allphase;
      #endif
      
      
      if (driver == NULL){
        #if ORBIT_PARALLEL == NO
          ndriver = ns->nstar;
          start_idx = 0;
        #else
          int base = ns->nstar / nprocs;
          int rem  = ns->nstar % nprocs;
          
          if (prank < rem){ndriver = base + 1; start_idx = prank*ndriver;}
          else{ndriver = base; start_idx = rem*(base+1)+(prank-rem)*base;}

          recvcounts = malloc(nprocs*sizeof(int));
          displs = malloc(nprocs*sizeof(int));
          for (i=0;i<nprocs;i++){
            int ln_r, st_r;
            if (i < rem) {
                ln_r = base + 1;
                st_r = i * ln_r;
            } else {
                ln_r = base;
                st_r = rem*(base+1) + (i - rem)*base;
            }
            recvcounts[i] = ln_r * 6;    /* 每行 6 个 double */
            displs[i]     = st_r * 6;
          }
          
          allphase = ARRAY_2D(ndriver, 6, double);
          printLog ("Driver: %d. Start index: %d. \n", ndriver, start_idx);
        #endif

        driver = malloc(ndriver * sizeof(gsl_odeiv2_driver *));
        for (i=0;i<ndriver;i++) {
          driver[i] = gsl_odeiv2_driver_alloc_y_new(&sys,
                                    gsl_odeiv2_step_rkck,
                                    h_init,
                                    eps_abs, 
                                    eps_rel);
        }
        printLog ("!nstar.c: initialize ode driver\n");
      }

      if (ns->time < g_time+toffset){
        // load new stars when evotime exceeds lifetime, also alloc new driver
        double dt =  g_time+toffset - ns->time;
        if (ns->evotime + dt> ns->lifetime){
            ns->start += ns->nstar;
            LoadSchwarzschildStars(ns->start+1, ns->nstar, ns, 0);
            for (i=0;i<ndriver;i++) {
              gsl_odeiv2_driver_reset(driver[i]);
            }
            ns->time += ns->lifetime - ns->evotime;
            ns->evotime = 0.0;
            dt = g_time + toffset - ns->time;
        }
        //step forward
        for (i=0;i<ndriver;i++) {
          double tc = ns->time;
          int status = gsl_odeiv2_driver_apply(driver[i], &tc, g_time+toffset, ns->phase[i+start_idx]);
          if(status != GSL_SUCCESS)
          {
            printLog ("error at t=%.3f: %s\n", tc, gsl_strerror(status));
            QUIT_PLUTO(1);
          }
          else{
            for (j=0;j<3;j++){
              ns->coord[i+start_idx][j] = ns->phase[i+start_idx][j];
              ns->v[i+start_idx][j] = ns->phase[i+start_idx][j+3];
              #if ORBIT_PARALLEL
                allphase[i][j] = ns->phase[i+start_idx][j];
                allphase[i][j+3] = ns->phase[i+start_idx][j+3];
              #endif
            }
          }
        }
        ns->time = g_time+toffset;
        ns->evotime += dt;

        //Gather at rank0 and send to each processor if integrate in parallel
        #if ORBIT_PARALLEL == YES
          MPI_Gatherv(
                      &allphase[0][0],  /* 本 rank 发送缓冲 */
                      ndriver * 6,               /* 本 rank 发送元素数 */
                      MPI_DOUBLE,
                      &ns->phase[0][0],          /* root 接收缓冲 */
                      recvcounts, 
                      displs,
                      MPI_DOUBLE,
                      0, MPI_COMM_WORLD
                  );
          MPI_Bcast(
                    &ns->phase[0][0],  /* 从 phase[0][0] 开始 */
                    ns->nstar * 6,         /* 总共 nstar*6 个 double */
                    MPI_DOUBLE,
                    0, MPI_COMM_WORLD
                );

          for (i=0;i<ns->nstar;i++){
            for (j=0;j<3;j++){
              ns->coord[i][j] = ns->phase[i][j];
              ns->v[i][j] = ns->phase[i][j+3];
              //if (i==0) {
              //  int k;
              //  printLog("prank %d: %f ", prank, g_time);
              //  for (k=0;k<6;k++) printLog("%f ", ns->phase[i][k]);
              //  printLog("\n");
              //}
            }
          }
        #endif
      }
    #endif
}

/* *************************** */
void UpdateAGBwind (const Data *d, double dt, Grid *grid)
{
  int   i, j, k, nv, nstar, ns, idim;
  double *x1, *x2,*x3, *dx1, *dx2, *dx3;
  double  r, r0, cs_amb, r1, r_bh, r_acc, T, mu, temp;
  double  Vwind, rho, rho_amb, distance, distance1, dM, dmdt, omega, omega1;
  double  *cs, *vs;
  double  unit_mass, rho_plus, rho_mswind, rho_mswind_add, new_prim[NVAR], old_prim[NVAR];
  double  kin, kinw, kinnew;

  // printLog ("!nstar.c: UpdateAGBwind called.\n");
  /* initialize and update g_nstar */
  if (g_nstar.coord == NULL) SetupNstar(&g_nstar, grid);
  #if NSTAR == SIMPLE
    #ifndef CHOMBO
      if (g_nstar.time < g_time) 
    #else
      if ((g_nstar.time < g_time) && (grid->level == 0))
    #endif
  #elif NSTAR == SCHWAR || NSTAR == FOKPLA
    if (g_nstar.time < g_time+toffset) 
  #endif
  {
    UpdateNstar(&g_nstar, grid);
  }

  unit_mass = UNIT_DENSITY*pow(UNIT_LENGTH, 3.);
  x1 = grid->xgc[IDIR];
  x2 = grid->xgc[JDIR];
  x3 = grid->xgc[KDIR];

  dx1 = grid->dx[IDIR];
  dx2 = grid->dx[JDIR];
  dx3 = grid->dx[KDIR];

  nstar = g_nstar.nstar;
  r0 = g_nstar.rstar;
  T = g_inputParam[T_WIND]; //wind temperature
  Vwind = g_inputParam[V_WIND]; //wind speed (at simulation surface)
  dM = g_inputParam[dM_RATE]; //mass loss rate, in solar mass per year
  dmdt = dM*(CONST_Msun/unit_mass)*(UNIT_TIME/CONST_yr);
  #if DIMENSIONS == 2
    rho = dmdt/(2.*CONST_PI*r0*(g_domEnd[KDIR]-g_domBeg[KDIR]))/Vwind;
    rho_plus = dmdt*dt/(CONST_PI*r0*r0*(g_domEnd[KDIR] - g_domBeg[KDIR]));
  #elif DIMENSIONS == 3
    rho_plus = dmdt*dt/(4./3.*CONST_PI*pow(r0, 3.));
  #endif
  mu = 0.615752;
  /* ****** N-star coordinates update ******* */
  //need to distinguish static and AMR grid?
  //check like CH_MPI or other definitions
  //static grid still has grid->level 0
  //only update for g_intStage 0 and 1
  /* **************************************** */
    #if STARSEARCH == NO
      DOM_LOOP(k,j,i){   
       double radius_scale;
       radius_scale = MAX(MAX(dx1[i], dx2[j]), dx3[k])/grid->dl_min[IDIR];
       rho_plus = dmdt*dt/(4./3.*CONST_PI*pow(r0, 3.))/pow(radius_scale, 3.0);
       
       for (ns=0;ns<g_nstar.nstar;ns++) {        
         cs = g_nstar.coord[ns]; //coordinates of star
         r  = sqrt(DIM_EXPAND((x1[i]-cs[IDIR])*(x1[i]-cs[IDIR]), 
                            + (x2[j]-cs[JDIR])*(x2[j]-cs[JDIR]),
                            + (x3[k]-cs[KDIR])*(x3[k]-cs[KDIR])
                              )); //distance to star
        //printLog ("x=%12.6e, y=%12.6e, r=%12.6e\n", x1[i], x2[j], r);
        if (r <= g_nstar.rstar*radius_scale){
          vs = g_nstar.v[ns]; //velocity of star
          if (Vwind*UNIT_VELOCITY > 1.e7){ 
            //fast wind with v > 100 km/s, set rho proportional to r-2.
            //in this case, rho_plus is M/(4*pi*rstar^3)*(rstar/r)^2
            double scale;
            scale = pow(g_nstar.rstar*radius_scale/MAX(r, MAX(MAX(dx1[i], dx2[j]), dx3[k])/2.), 2.);
            rho_plus = dmdt*dt/(4.*CONST_PI*pow(g_nstar.rstar*radius_scale, 3.))*scale;
          }

          for (nv=0;nv<NVAR;nv++) {
              new_prim[nv] = old_prim[nv] = d->Vc[nv][k][j][i];
          }
          new_prim[RHO] += rho_plus;
          DIM_EXPAND(
            new_prim[VX1] = (old_prim[RHO]*old_prim[VX1] + rho_plus*(Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR]))/new_prim[RHO];, 
            new_prim[VX2] = (old_prim[RHO]*old_prim[VX2] + rho_plus*(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR]))/new_prim[RHO];,
            new_prim[VX3] = (old_prim[RHO]*old_prim[VX3] + rho_plus*(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR]))/new_prim[RHO];)
          kin = 0.5*old_prim[RHO]*(DIM_EXPAND(old_prim[VX1]*old_prim[VX1], +old_prim[VX2]*old_prim[VX2], +old_prim[VX3]*old_prim[VX3]));
          kinw = 0.5*rho_plus*(DIM_EXPAND((Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR])*(Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR]), 
                  +(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR])*(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR]),
                 +(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR])*(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR])));
          kinnew = 0.5*new_prim[RHO]*(DIM_EXPAND(new_prim[VX1]*new_prim[VX1], +new_prim[VX2]*new_prim[VX2], +new_prim[VX3]*new_prim[VX3]));
          new_prim[PRS] = (kin+old_prim[PRS]/(g_gamma - 1.)+kinw+rho_plus*T/(KELVIN*mu*(g_gamma-1.))-kinnew)*(g_gamma - 1.);
            new_prim[TRC] = (old_prim[TRC]*old_prim[RHO] + rho_plus)/new_prim[RHO];
          for (nv=0;nv<NVAR;nv++) {
              d->Vc[nv][k][j][i] = new_prim[nv];
          }
        }
      }
    }
      #elif STARSEARCH == YES
    for (ns=0;ns<g_nstar.nstar;ns++) { 
      // determine whether the star locates within the box
      int ilo, ihi, imid, idim, is, js, ks;
      double xs, ys, zs, tmp, rstar, sl[3], sr[3], il[3], ir[3];

      cs = g_nstar.coord[ns];
      xs = g_nstar.coord[ns][0];
      ys = g_nstar.coord[ns][1];
      zs = g_nstar.coord[ns][2];

      // bisearch along x-axis 
      ilo=grid->gbeg[0]; ihi=grid->gend[0]+1;
      if ((xs<grid->xl_glob[0][ilo])||(xs>grid->xl_glob[0][ihi])) continue;
      else{
        while (ilo != ihi-1)
        {
          imid = (ilo + ihi)/2;
          tmp = grid->xl_glob[0][imid];
          if (xs < tmp) ihi=imid;
          else ilo=imid;
        }
        is = ilo;
      }

      // bisearch along y-axis
      ilo=grid->gbeg[1]; ihi=grid->gend[1]+1;
      if ((ys<grid->xl_glob[1][ilo])||(ys>grid->xl_glob[1][ihi])) continue;
      else{
        while (ilo != ihi-1)
        {
          imid = (ilo + ihi)/2;
          tmp = grid->xl_glob[1][imid];
          if (ys < tmp) ihi=imid;
          else ilo=imid;
        }
        js = ilo;
      }

      // bisearch along z-axis
      ilo=grid->gbeg[2]; ihi=grid->gend[2]+1;
      if ((zs<grid->xl_glob[2][ilo])||(zs>grid->xl_glob[2][ihi])) continue;
      else{
        while (ilo != ihi-1)
        {
          imid = (ilo + ihi)/2;
          tmp = grid->xl_glob[2][imid];
          if (zs < tmp) ihi=imid;
          else ilo=imid;
        }
        ks = ilo;
      }

      //stellar radius and density augment adaptively determined by star position
      rstar = g_inputParam[WIND_RADIUS]*MAX(MAX(grid->dx_glob[IDIR][is], grid->dx_glob[JDIR][js]), grid->dx_glob[KDIR][ks]);
      rho_plus = dmdt*dt/(4./3.*CONST_PI*pow(rstar, 3.0));

      for (idim=0;idim<3;idim++){
        sl[idim] = cs[idim] - rstar;
        sr[idim] = cs[idim] + rstar;
      }

      //ignore stars outside the box
      if (((xs-rstar) > x1[grid->lend[IDIR]]) || ((xs+rstar) < x1[grid->lbeg[IDIR]])) continue;
      if (((ys-rstar) > x2[grid->lend[JDIR]]) || ((ys+rstar) < x2[grid->lbeg[JDIR]])) continue;
      if (((zs-rstar) > x3[grid->lend[KDIR]]) || ((zs+rstar) < x3[grid->lbeg[KDIR]])) continue;

      //bisearch for loop boundary
      for (idim=0;idim<3;idim++){
        ilo=grid->beg[idim]; ihi=grid->end[idim]+1;
        if ((sl[idim]<grid->xl_glob[idim][ilo])) il[idim] = grid->lbeg[idim];
        else{
          while (ilo != ihi-1)
          {
            imid = (ilo + ihi)/2;
            tmp = grid->xl_glob[idim][imid];
            if (sl[idim] < tmp) ihi=imid;
            else ilo=imid;
          }
          il[idim] = ilo - (grid->beg[idim] - grid->lbeg[idim]); //transfer to local index
        }
      }

      for (idim=0;idim<3;idim++){
        ilo=grid->beg[idim]; ihi=grid->end[idim]+1;
        if ((sr[idim]>grid->xl_glob[idim][ihi])) ir[idim] = grid->lend[idim];
        else{
          while (ilo != ihi-1)
          {
            imid = (ilo + ihi)/2;
            tmp = grid->xl_glob[idim][imid];
            if (sr[idim] < tmp) ihi=imid;
            else ilo=imid;
          }
          ir[idim] = ilo - (grid->beg[idim] - grid->lbeg[idim]); //transfer to local index
        }
      }

      //loop in the local grid. il and ir are for local index.
      for (i=il[IDIR];i<ir[IDIR]+1;i++){
        for (j=il[JDIR];j<ir[JDIR]+1;j++){
          for (k=il[KDIR];k<ir[KDIR]+1;k++){
            r  = sqrt(DIM_EXPAND((x1[i]-cs[IDIR])*(x1[i]-cs[IDIR]), 
                                + (x2[j]-cs[JDIR])*(x2[j]-cs[JDIR]),
                                + (x3[k]-cs[KDIR])*(x3[k]-cs[KDIR])
                                  )); //distance to star
            //printLog ("x=%12.6e, y=%12.6e, r=%12.6e\n", x1[i], x2[j], r);
            if (r <= rstar){
              vs = g_nstar.v[ns]; //velocity of star
              if (Vwind*UNIT_VELOCITY > 1.e7){ 
                //fast wind with v > 100 km/s, set rho proportional to r-2.
                //in this case, rho_plus is M/(4*pi*rstar^3)*(rstar/r)^2
                double scale;
                scale = pow(rstar/MAX(r, MAX(MAX(dx1[i], dx2[j]), dx3[k])/2.), 2.);
                rho_plus = dmdt*dt/(4.*CONST_PI*pow(rstar, 3.))*scale;
              }

              for (nv=0;nv<NVAR;nv++) {
                  new_prim[nv] = old_prim[nv] = d->Vc[nv][k][j][i];
              }
              new_prim[RHO] += rho_plus;
              DIM_EXPAND(
                new_prim[VX1] = (old_prim[RHO]*old_prim[VX1] + rho_plus*(Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR]))/new_prim[RHO];, 
                new_prim[VX2] = (old_prim[RHO]*old_prim[VX2] + rho_plus*(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR]))/new_prim[RHO];,
                new_prim[VX3] = (old_prim[RHO]*old_prim[VX3] + rho_plus*(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR]))/new_prim[RHO];)
              kin = 0.5*old_prim[RHO]*(DIM_EXPAND(old_prim[VX1]*old_prim[VX1], +old_prim[VX2]*old_prim[VX2], +old_prim[VX3]*old_prim[VX3]));
              kinw = 0.5*rho_plus*(DIM_EXPAND((Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR])*(Vwind*(x1[i]-cs[IDIR])/MAX(1.e-12, r) + vs[IDIR]), 
                      +(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR])*(Vwind*(x2[j]-cs[JDIR])/MAX(1.e-12, r) + vs[JDIR]),
                    +(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR])*(Vwind*(x3[k]-cs[KDIR])/MAX(1.e-12, r) + vs[KDIR])));
              kinnew = 0.5*new_prim[RHO]*(DIM_EXPAND(new_prim[VX1]*new_prim[VX1], +new_prim[VX2]*new_prim[VX2], +new_prim[VX3]*new_prim[VX3]));
              new_prim[PRS] = (kin+old_prim[PRS]/(g_gamma - 1.)+kinw+rho_plus*T/(KELVIN*mu*(g_gamma-1.))-kinnew)*(g_gamma - 1.);
              #if NTRACER > 0
                new_prim[TRC] = (old_prim[TRC]*old_prim[RHO] + rho_plus)/new_prim[RHO];
              #endif
              for (nv=0;nv<NVAR;nv++) {
                  d->Vc[nv][k][j][i] = new_prim[nv];
              }
            }
          }
        }
      }
    }
  #endif
}

void LoadSchwarzschildStars(int nskip, int nstar, Nstar *ns, int nsbegin)
/* Skip first "start" lines and read nstart lines
 * Usually this is nstar*round+1
 * ***************************************** */
{ 
    FILE *fp;
    char line[256], orbit_file[256];
    int i, j, count = 0;
    
    sprintf(orbit_file, "%s/schwar.ini", RuntimeGet()->output_dir);
    fp = fopen(orbit_file, "r");
    // 跳过前 start 行
    for (i = 0; i < nskip; i++) {
        if (!fgets(line, sizeof(line), fp)) {
            printLog("Error: Failed to skip line %d\n", i);
            QUIT_PLUTO(1);
        }
    }

    // 从第 start 行开始，读取 nstar 行
    for (i = nsbegin; i < nsbegin+nstar; i++) {
        if (!fgets(line, sizeof(line), fp)) {
            printLog("Error: Failed to read line %d\n", nskip + i);
            break;
        }

        int orbindex, orbtype;
        double x, y, z, vx, vy, vz;
        int parsed = sscanf(line, "%d %lf %lf %lf %lf %lf %lf %d",
                            &orbindex, &x, &y, &z, &vx, &vy, &vz, &orbtype);

        if (parsed != 8) {
            printLog("Warning: Malformed line at %d, parsed %d values\n", nskip + i - nsbegin, parsed);
            continue;
        }

        ns->Orbindex[i] = orbindex;
        ns->coord[i][0] = x;
        ns->coord[i][1] = y;
        ns->coord[i][2] = z;
        ns->v[i][0] = vx;
        ns->v[i][1] = vy;
        ns->v[i][2] = vz;
        ns->Orbtype[i] = orbtype;
        for (j=0; j<3; j++){
          ns->phase[i][j] = ns->coord[i][j];
          ns->phase[i][j+3] = ns->v[i][j];
        }
    }
    printLog ("Load %d stars from line %d\n", nstar, nskip);
    fclose(fp);
}

// double F(double, double, double, double, double *)
double potfunc(double x, double y, double z, double tau, double *args){
    double m, sigma, b, c, delta, epsilon;
    double tau2;

    tau2 = tau*tau;
    m = args[0];
    sigma = args[1];
    b = args[2];
    c = args[3];

    
    delta = 1-b*b;
    epsilon = 1-c*c;
    return exp(-tau2/(2*sigma*sigma)*(x*x+y*y/(1-delta*tau2)+z*z/(1-epsilon*tau2)))/sqrt((1-delta*tau2)*(1-epsilon*tau2));
}

void forcefunc(double x, double y, double z, double tau, double *args, double *accel){
    double sigma2=args[1]*args[1], tau2=tau*tau;
    double b, c, delta, epsilon, pot;

    b = args[2];
    c = args[3];
    delta = 1-b*b;
    epsilon = 1-c*c;
    pot = potfunc(x, y, z, tau, args);

    accel[0] = -pot*(-tau2/sigma2*x);
    accel[1] = -pot*(-tau2/sigma2*y/(1-delta*tau2));
    accel[2] = -pot*(-tau2/sigma2*z/(1-epsilon*tau2));
}

double calcTriPot(double x, double y, double z, double *args, int glorder, double *glx, double *glw){
   double pot;
   int i;

   pot = 0.0;
   for (i=0; i<glorder; i++){
       pot += potfunc(x, y, z, glx[i], args)*glw[i]; 
   }
   pot *= -args[0]/args[1]*sqrt(2./CONST_PI);
   return pot;
}

void calcTriAccel(double x, double y, double z, double *args, int glorder, double *glx, double *glw, double *accel){
   double accel_tmp[3];
   int i, j;
   
   for (j=0; j<3; j++) accel[j]=0.0;
   for (i=0; i<glorder; i++){ 
       forcefunc(x, y, z, glx[i], args, accel_tmp);
       for (j=0; j<3; j++){
            accel[j] += accel_tmp[j]*glw[i];
       }
   }
   for (j=0; j<3; j++){
        accel[j] *= -args[0]/args[1]*sqrt(2./CONST_PI);
   }
}

void calcBHAccel(double x, double y, double z, double gmbh, double *ax, double *ay, double *az){
    double r;
    r = sqrt(x*x + y*y + z*z);

    *ax += -gmbh/pow(r, 3)*x;
    *ay += -gmbh/pow(r, 3)*y;
    *az += -gmbh/pow(r, 3)*z;
}

double calcTotalPot(double x, double y, double z, OrbitParams *p){
  int i;
  double pot=0.0;

  for(i=0; i<p->Nmge; i++){
    pot += calcTriPot(x, y, z, p->args[i], p->glorder, p->glx, p->glw);
  }
  pot += -p->gmbh/sqrt(x*x + y*y + z*z);
  // print ("calcTotalPot: %d %f %d %f %f %f\n", 
  //       p->Nmge, p->args[0][0], p->glorder, p->glx[0], p->glw[0], p->gmbh);
  // printLog ("calcTotalPot: %f %f %f %f", x, y, z, pot);
  return pot;
}

void calcTotalAccel(double x, double y, double z, OrbitParams *p, double *accel){
  double ax=0, ay=0, az=0;
  int i;
  double a_i[3];

  for(i=0; i<p->Nmge; i++){
      calcTriAccel(x, y, z,
                  p->args[i],
                  p->glorder, p->glx, p->glw,
                  a_i);
      ax += a_i[0];
      ay += a_i[1];
      az += a_i[2];
  }
  calcBHAccel(x, y, z, p->gmbh, &ax, &ay, &az);
  accel[0] = ax;
  accel[1] = ay;
  accel[2] = az;
  // printLog ("calcTotalAccel: %f %f %f %f %f %f \n", x, y, z, ax, ay, az);
}

void calcFPAccel(double x, double y, double z, FPParams *p, double *accel)
/* calculate acceleration of a FP potential */
{
  double logr = log10(sqrt(x*x + y*y + z*z));
  double r3 = pow(x*x + y*y + z*z, 3.0/2.0);
  int klo, kmid, khi, rmid;
  double gmencl, gmtotal;
  /* ----------------------------------------------
        Binary search for radius
   ---------------------------------------------- */

  klo = 0;
  khi = p->Ngrid - 1;
  //printLog ("binary search start...");
  if (logr < p->logr[klo]){
    gmencl = 0.0;
  }
  else{
    while (klo != (khi - 1)){
      kmid = (klo + khi)/2;
      rmid = p->logr[kmid];
      //printLog ("Tmid = %.3f \n", Tmid);
      if (logr <= rmid){
        khi = kmid;
      }else if (logr > rmid){
        klo = kmid;
      }
    }
    gmencl = pow(10, ((logr - p->logr[klo])*p->loggmstar[khi] + (p->logr[khi]-logr)*p->loggmstar[klo])/(p->logr[khi] - p->logr[klo]));
  }
  gmtotal = gmencl + p->gmbh;
  accel[0] = -gmtotal/r3*x;
  accel[1] = -gmtotal/r3*y;
  accel[2] = -gmtotal/r3*z;
}

/* ODE 系统： state = {x,y,z,vx,vy,vz} */
#if NSTAR == SCHWAR
  int orbit_ode(double t, const double y[], double dydt[], void *params){
      OrbitParams *p = (OrbitParams*) params;

      /* 位置导数 = 速度 */
      dydt[0] = y[3];
      dydt[1] = y[4];
      dydt[2] = y[5];

      calcTotalAccel(y[0], y[1], y[2], p, dydt+3);
      return GSL_SUCCESS;
  }
#elif NSTAR == FOKPLA
  int orbit_ode(double t, const double y[], double dydt[], void *params){
      FPParams *p = (FPParams*) params;

      /* 位置导数 = 速度 */
      dydt[0] = y[3];
      dydt[1] = y[4];
      dydt[2] = y[5];

      calcFPAccel(y[0], y[1], y[2], p, dydt+3);
      return GSL_SUCCESS;
  }
#endif

#if NSTAR == SCHWAR
  void OrbitParamInit()
  /* Initialize orbitparam for potential and acceleration calculation*/
  {
    if (opinit == 0){
        // ============ read mge data ==========
      FILE *fp;
      char orbit_file[256], line[256], *token;
      int Nmge=0, i, j;
      gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(glorder);
      double *glx, *glw;

      fp = fopen("mge.dat", "r");
      while (fgets(line, 256, fp) != NULL){
          Nmge++;
      }
      rewind(fp);
      printLog ("MGE components: %i\n", Nmge);

      double **args = malloc(Nmge * sizeof(double*));
      for (i=0;i<Nmge;i++){
        *(args+i) = malloc(4*sizeof(double));
      }
      i = 0;
      while (fgets(line, 256, fp) != NULL){
        j=0;
        token = strtok(line, " ");
        while (token != NULL){
          if (j == 0) args[i][j] = atof(token)*UNIT_GM; //times unitGM
          else args[i][j] = atof(token);
          printLog ("%.6f  ", args[i][j]);
          token = strtok(NULL, " ");
          j++;
        }
        printLog ("\n");
        i++;
      }
      fclose(fp);

      
      orbitparam.args = args;
      orbitparam.Nmge = Nmge;
      orbitparam.glorder = glorder;
      orbitparam.glx = malloc(glorder * sizeof(double));
      orbitparam.glw = malloc(glorder * sizeof(double));
      for(i=0; i<glorder; i++) gsl_integration_glfixed_point(0.0, 1.0, i, (orbitparam.glx)+i, (orbitparam.glw)+i, t);
      orbitparam.gmbh = g_inputParam[M_BH]*UNIT_GM;
      
      // output log
      printLog("orbitparam glorder: %d\n", glorder);
      printLog("orbitparam glx[0]: %f\n", orbitparam.glx[0]);
      printLog("orbitparam glw[0]: %f\n", orbitparam.glw[0]);
      printLog("orbitparam GMBH: %f\n", orbitparam.gmbh);

      opinit = 1;
    }
  }
#elif NSTAR == FOKPLA
  void OrbitParamInit()
  /* Initialize Fokker-Planck parameters for potential and acceleration calculation*/
  {
    if (opinit == 0){
      // ============ read enclosed mass ==========
      FILE *fp;
      char line[256], *token;
      int Nr=0, i, j;

      fp = fopen("FPprofile.dat", "r");
      while (fgets(line, 256, fp) != NULL){
          Nr++;
      }
      rewind(fp);
      printLog ("Fokker-Planck grid: %i\n", Nr);

      double *logr = malloc(Nr * sizeof(double)), *loggmstar = malloc(Nr * sizeof(double));
      i = 0;
      while (fgets(line, 256, fp) != NULL){
        j=0;
        token = strtok(line, " ");
        while (token != NULL){
          if (j == 0) logr[i] = log10(atof(token));
          if (j == 1) loggmstar[i] = log10(atof(token)*UNIT_GM);   
          token = strtok(NULL, " ");
          j++;
        }
        if (i < 10) printLog ("%.6e  %.6e \n", pow(10, logr[i]), pow(10, loggmstar[i])/UNIT_GM);
        i++;
      }
      fclose(fp);

      orbitparam.logr = logr;
      orbitparam.loggmstar = loggmstar;
      orbitparam.gmbh = g_inputParam[M_BH]*UNIT_GM;
      orbitparam.Ngrid = Nr;
      
      // output log
      printLog("orbitparam GMBH: %f\n", orbitparam.gmbh);
      opinit = 1;
    }
  }
#endif

void rkck(double *y, double x, double h, double *yout, double *yerr, void (*derivs)(double, double *, double *))
{
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2, b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
                  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
                  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0, b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0, c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
                  dc5 = -277.00/14336.0;

    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc6=c6-0.25; 
    double *dydx, *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;
    int i, n=6;
    
    //printLog ("rkck called \n");
    dydx = (double *) malloc(n*sizeof(double));
    ak2 = (double *) malloc(n*sizeof(double));
    ak3 = (double *) malloc(n*sizeof(double));
    ak4 = (double *) malloc(n*sizeof(double));
    ak5 = (double *) malloc(n*sizeof(double));
    ak6 = (double *) malloc(n*sizeof(double));
    ytemp = (double *) malloc(n*sizeof(double));
    
    (*derivs)(x,y,dydx);
    for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i]; 
    (*derivs)(x+a2*h,ytemp,ak2);

    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]); 
    (*derivs)(x+a3*h,ytemp,ak3);

    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]); 
    (*derivs)(x+a4*h,ytemp,ak4);

    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]); 
    (*derivs)(x+a5*h,ytemp,ak5);

    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]); 
    (*derivs)(x+a6*h,ytemp,ak6);

    for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); 

    for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    
    free(dydx);
    free(ak2);
    free(ak3);
    free(ak4);
    free(ak5);
    free(ak6);
    free(ytemp);
}

void HamiltonianDerivs(double t, double *y, double *dydx)
/* Calculate right hand side of motion equation for a Plummer NSC + SMBH
 * This function can be substituded with the orbit_ode function
 * ************************************************************** */
{
    static double G, Mbh, r0, Mnsc;
    static int first_call=1;
    double r;

    if (first_call){
        // G=g_inputParam[G_CONS]; 
        // Mbh=g_inputParam[M_BH];
        // r0=g_inputParam[R0];
        // Mnsc=g_inputParam[M_NSC];
        first_call = 0;
    }

    r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    dydx[0] = y[3];
    dydx[1] = y[4];
    dydx[2] = y[5];
    dydx[3] = -G*Mnsc*y[0]/pow(r0*r0 + r*r, 3./2.) - G*Mbh*y[0]/(r*r*r);
    dydx[4] = -G*Mnsc*y[1]/pow(r0*r0 + r*r, 3./2.) - G*Mbh*y[1]/(r*r*r);
    dydx[5] = -G*Mnsc*y[2]/pow(r0*r0 + r*r, 3./2.) - G*Mbh*y[2]/(r*r*r);
}
