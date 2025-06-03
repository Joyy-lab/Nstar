typedef struct OrbitParams_ {
        double **args;      /* 传给 calcAccel 的 args 数组 */
        int Nmge;
        int     glorder;   /* Gauss–Legendre 节点数 */
        double *glx, *glw; /* 节点和权重数组 */
        double gmbh;
    } OrbitParams;

typedef struct Nstar_{
    //constants
    int init; //whether stucture is initialzied
    int nstar; //number of stars implemented in the simulation
    int nsample; //number of star sample
    double rstar; //stellar wind injection radius
    double racc; //accretion radius
    double lifetime; //lifetime of the stars
    
    
    #if NSTAR == SIMPLE
      //constants for Keplerian orbits
      double **orbit; //six orbital parameters
      double **P; //first vector for coordinates and velocity calculation
      double **Q; //second vector for coordinates and velocity calculation
      double *omega;
      //variables for Keplerian orbits
      double *M; //mean anomaly
      double *E; //eccentric anomaly

    #elif NSTAR == SCHWAR
      //constants for Schwarzschild orbits
      int *Orbtype; //orbit type, 1-3 is x, y, z tube orbit, 4 is box orbit, 5 is other
      int *Orbindex; //orbit index for each star
      int start;
    #endif 

    //shared variables
    double time;
    double evotime; //evolutionary time of the stars
    double **coord; //coordinates
    double **v; //velocity
    double **phase; //coordinates&velocity
    double *r; //distance 
  } Nstar;

void UpdateNstar (Nstar *ns, Grid *grid);
void SetupNstar (Nstar *ns, Grid *grid);
void UpdateAGBwind (const Data *d, double dt, Grid *grid);


extern Nstar g_nstar;
extern OrbitParams orbitparam;
    
double F(double x, double y, double z, double tau, double *args);
void force(double x, double y, double z, double tau, double *args, double *accel);
double calcTriPot(double x, double y, double z, double *args, int glorder, double *glx, double *glw);
void calcTriAccel(double x, double y, double z, double *args, int glorder, double *glx, double *glw, double *accel);
void calcBHAccel(double x, double y, double z, double gmbh, double *ax, double *ay, double *az);
int orbit_ode(double t, const double y[], double dydt[], void *params);
void OrbitParamInit();

double   KeplerSolver (double, double); //init.c

#if NSTAR == SIMPLE
    void rkck(double *y, double x, double h, double *yout, double *yerr, void (*derivs)(double, double *, double *));
    void HamiltonianDerivs(double, double, double *);
#elif NSTAR ==  SCHWAR
    void LoadSchwarzschildStars(int nskip, int nstar, Nstar *ns, int nsbegin);
    double potfunc(double x, double y, double z, double tau, double *args);
    void forcefunc(double x, double y, double z, double tau, double *args, double *accel);
    double calcTriPot(double x, double y, double z, double *args, int glorder, double *glx, double *glw);
    void calcTriAccel(double x, double y, double z, double *args, int glorder, double *glx, double *glw, double *accel);
    void calcBHAccel(double x, double y, double z, double gmbh, double *ax, double *ay, double *az);
    double calcTotalPot(double x, double y, double z, OrbitParams *p);
    void calcTotalAccel(double x, double y, double z, OrbitParams *p, double *accel);
    int orbit_ode(double t, const double y[], double dydt[], void *params);
    void OrbitParamInit();
#endif

