/* Constants */
extern double he3_gyro_;
extern double he3_amass_;
extern double he3_mmass_;
extern double const_na_;
extern double const_kb_;
extern double const_r_;
extern double const_h_;
extern double const_hbar_;
extern double const_pi_;
extern double he3_pcr_;
extern double he3_tcr_;
extern double he3_pabn_;
extern double he3_tabn_;
extern double he3_psmin_;
extern double he3_tsmin_;
extern double he3_pa_;
extern double he3_ta_;
extern double he3_pb_;
extern double he3_tb_;
extern double he3_pneel_;
extern double he3_tneel_;

/* Phase diagram functions */
double he3_pvap_(double *T);   /* Vapor pressure [bar] vs T [K] */
double he3_pmelt_(double *T);  /* Melting pressure [bars] vs T [K] */
double he3_pmelt_gr_(double *T);  /* Greywall-86 melting pressure [bars] vs T [K] */
double he3_tc_(double *P);     /* T_c [mK] vs P [bar] */
double he3_tab_(double *P);    /* T_ab [mK] vs P [bar] */

/* Fermi-liquid values, measured */
double he3_vm_(double *P);     /* Molar Volume [cm**3/mole] vs P [bar] */
double he3_gammaf_(double *P); /* R-Gas constant GAMMA=C/RT [1/(K*mol)] vs P [bar] */
double he3_c1_(double *P);     /* First sound velosity [cm/s] vs P [bar] */
double he3_tmag_(double *P);   /* Magnetic temperature [K] vs P [bar] */

/* Fermi-liquid values, derived */
double he3_rho_(double *P);    /* Density [g/cm^3] vs P [bar] */
double he3_2n0_(double *P);    /* 2N0 vs P [bar] */
double he3_mm_(double *P);     /* Effective mass / atom mass vs P [bar] */
double he3_meff_(double *P);   /* Effective mass [g] vs P [bar] */
double he3_pf_(double *P);     /* Fermi momentum [sgs] vs P [bar] */
double he3_vf_(double *P);     /* Fermi velocity [cm/s] vs P [bar] */
double he3_chi_n_(double *P);  /* Normal liquid susceptibility vs P [bar]  */
double he3_f0a_(double *P);    /* F0a vs P [bar] (Z0/4)*/
double he3_f0s_(double *P);    /* F0s vs P [bar] */
double he3_f1a_(double *P);    /* F1a vs P [bar] */
double he3_f1s_(double *P);    /* F1s vs P [bar] */
double he3_a_(double *P);      /* average atomic spacing, angstr. */
double he3_gdk_(double *P);    /* average dipolar coupling energy, K */
double he3_tfeff_(double *P);  /* effective fermi temperature, K */

/* Gap */
double he3_bcsgap_(double *ttc);             /* BCS energy gap */
double he3_bcsgap_fast_(double *ttc);        /* BCS energy gap approximation */
double he3_trivgap_(double *ttc, double *p); /* Trivial strong-coupling correction to the BCS gap*/
double he3_yosida_(double *ttc, double *gap, double *n); /* Yosida functions */

double he3_chi_b_(double *ttc, double *P);  /* B-phase susceptibility vs P [bar]  */
double he3_nu_b_(double *ttc, double *P);   /* B-phase Leggett frequency, Hz */

double he3_z3_(double *ttc, double *gap);
double he3_z5_(double *ttc, double *gap);
double he3_z7_(double *ttc, double *gap);
double he3_yosida0_(double *ttc, double *gap); /* Y0 function -- old */
double he3_yosida0_fast_(double *ttc, double *gap);      /* Yosida0 approximation*/

/* Transport */
double he3_scatt_l1a_(double *P); /* Scattering parameters */
double he3_scatt_g0_(double *P);
double he3_scatt_d0_(double *P);
double he3_coll_int_(double *xi, double *ttc, double *gap,
            double *g0, double *d0); /* Collision integral in Einzel approximation */
double he3_coll_int_lt_(double *xi, double *ttc, double *gap,
            double *g0, double *d0); /* Collision integral for low temp (good for < 0.7Tc) */
double he3_coll_int_ht_(double *xi, double *ttc, double *gap,
            double *g0, double *d0); /* Collision integral for high temp */

double he3_tau_n0_(double *ttc, double *P);    /* Normal state quasiparticle lifetime at the Fermi level, s */
double he3_tau_n_av_(double *ttc, double *P);  /* Thermal average quasiparticle lifetime, s */
double he3_tau0_(double *ttc, double *p);      /* Bogoliubov quasiparticle lifetime at the Fermi level, s */
double he3_tau_av_(double *ttc, double *p);    /* Thermal average quasiparticle lifetime, s */
double he3_fpath_(double *ttc, double *p);     /* Mean free path of Bogoliubov quasiparticles */
double he3_tau_dperp_(double *ttc, double *p); /* Spin diffusion perp transport time, s */
double he3_tau_dpar_(double *ttc, double *p);  /* Spin diffusion parallel transport time, s */
double he3_sdiff_(double *ttc, double *p, double *nu0);  /* Spin diffusion D_perp */

/* Other */
double he3_exp_nu_b_(double *ttc, double *P);  /* B-phase Leggett frequency, Hz */
double he3_flegg_(double *P, double *ttc); /* Legget freq^2, [Hz^2] vs P [bar], T/Tc */

double he3_swvel_(double *P, double *ttc);     /* Osheroff's spin wave vel. [cm/s] vs P [bar], T [mK] */
double he3_swvel_par_(double *P, double *ttc); /* Perp Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */
double he3_swvel_per_(double *P, double *ttc); /* Parallel Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */

double he3_ds_exp_(double *P, double *ttc); /* spin diffusion coeff. in superfluid He3 (measured) */
double he3_dn_exp_(double *P, double *ttc); /* spin diffusion in normal He3 */
double he3_d_exp_(double *P, double *ttc);  /* combined normal + superfluid spin diffusion */

double he3_susept_(double *P, double *ttc); /* Suseptibility [sgs] vs P [bar], T [mK] */

double he3_tau_r_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */
double he3_tau_f_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */

