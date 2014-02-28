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
extern double const_2pi_;
extern double he3_pcr_;
extern double he3_tcr_;
extern double he3_pabn_;
extern double he3_tabn_;
extern double he3_pm_;
extern double he3_tm_;
extern double he3_pa_;
extern double he3_ta_;
extern double he3_pb_;
extern double he3_tb_;
extern double he3_ps_;
extern double he3_ts_;

extern double he3_pa_plts_;
extern double he3_ta_plts_;
extern double he3_pb_plts_;
extern double he3_tb_plts_;
extern double he3_ps_plts_;
extern double he3_ts_plts_;

/* Phase diagram functions */
double he3_pvap_(double *T);   /* Vapor pressure [bar] vs T [K] */
double he3_pmelt_(double *T);  /* Melting pressure [bars] vs T [K] (Greywall86) */
double he3_tc_(double *P);     /* T_c [mK] vs P [bar] (Greywall86) */
double he3_tab_(double *P);    /* T_ab [mK] vs P [bar] (Greywall86) */

double he3_pmelt_plts_(double *T);  /* PLTS2000 melting pressure [bars] vs T [K] */
double he3_gr2plts_(double *T);     /* Greywall86 -> PLTS2000 temperature conversion */
double he3_plts2gr_(double *T);     /* PLTS2000 -> Greywall86 temperature conversion */

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
double he3_yosida_par_(double *ttc, double *gap);
double he3_yosida_perp_(double *ttc, double *gap);
double he3_z3_(double *ttc, double *gap);
double he3_z5_(double *ttc, double *gap);
double he3_z7_(double *ttc, double *gap);
double he3_lambda_(double *ttc, double *gap);
double he3_rho_nb_(double *ttc, double *P); /* B-phase Normal component density \rho_nb/\rho_0 */
double he3_chi_b_(double *ttc, double *P);  /* B-phase susceptibility chi_b/chi_0  */

/* Dipole */
double he3_gd_(double *P); /* Experimental value of Dipolar coefficient g_d, [1/(erg cm3)] */
double he3_ld_(double *ttc, double *P);     /* lambda_D = Delta^2 g_d, [erg/cm3] */
double he3_nu_b_(double *ttc, double *P);   /* B-phase Leggett frequency, Hz */
double he3_nu_b1_(double *ttc, double *P);  /* B-phase Leggett frequency, Hz */

/* Gradient */
double he3_text_lg1_(double *ttc, double *p);   /* Textural parameter lambda_{G1}, erg/cm */
double he3_text_lg2_(double *ttc, double *p);   /* Textural parameter lambda_{G2}, erg/cm */
double he3_text_delta_(double *ttc, double *p); /* Textural parameter delta */
double he3_text_cperp_(double *ttc, double *p); /* perpendicular spin wave velocity, cm/s */
double he3_text_c_(double *ttc, double *p); /* bending stiffness coefficient c, erg/cm */

/* Texture */
double he3_text_a_(double *ttc, double *p);     /* Textural parameter a [erg/cm^3 1/G^2] */
double he3_text_ldv_(double *ttc, double *p);   /* Textural parameter lambda_{DV} */
double he3_text_lhv_(double *ttc, double *p);   /* Textural parameter lambda_{HV} */

double he3_text_d_(double *ttc, double *p); /* Textural parameter d, erg/(cm^2 G^2) */
double he3_text_llh_(double *ttc, double *p, double *omega); /* Textural parameter lambda_LH */
double he3_text_lo_(double *ttc, double *p, double *omega); /* lambda/omega */

double he3_text_vd_(double *ttc, double *p);    /* Dipole velocity v_d, cm/s */
double he3_text_xid_(double *ttc, double *p);  /* Dipole length xi_d, cm */
double he3_text_xih_(double *ttc, double *p, double *h);  /* Magnetic length xi_h, cm */

/* Transport in the normal phase */
double he3_crsect_w_(double *P); /* Crossections */
double he3_crsect_wi_(double *P);
double he3_crsect_wd_(double *P);
double he3_crsect_wl_(double *P);
double he3_scatt_l1a_(double *P); /* Scattering parameters */
double he3_scatt_g0_(double *P);
double he3_scatt_d0_(double *P);
double he3_scatt_w0_(double *P); /* w0 = 1 - 2/3 g0 + d0 */
double he3_tau_n0_(double *ttc, double *P);    /* Normal state quasiparticle lifetime at the Fermi level, s */
double he3_tau_n_av_(double *ttc, double *P);  /* Thermal average quasiparticle lifetime, s */
double he3_tau_nd_(double *ttc, double *P);    /* Spin diffusion transport time for a normal Fermi-liquid, s */
double he3_diffn_hydr_(double *ttc, double *p);   /* Hydrodynamic spin diffusion in normal liquid, cm2/s */
double he3_diffn_perp_(double *ttc, double *p, double *nu0);  /* Spin diffusion D_perp in normal liquid, cm2/s */

/* Transport in the B phase */
double he3_coll_int_(double *xi, double *ttc, double *gap, double *g0, double *d0); /* Collision integral approximation ttc=0..1*/
double he3_coll_int_lt_(double *xi, double *ttc, double *gap, double *g0, double *d0); /* Collision integral for low temp (good for < 0.7Tc) */
double he3_coll_int_ht_(double *xi, double *ttc, double *gap, double *g0, double *d0); /* Collision integral for high temp */
double he3_tau0_(double *ttc, double *p);      /* Bogoliubov quasiparticle lifetime at the Fermi level, s */
double he3_tau0lt_(double *ttc, double *p);    /* he3_tau0_ at T->0 */
double he3_tau_av_(double *ttc, double *p);    /* Thermal average quasiparticle lifetime, s */
double he3_fpath_(double *ttc, double *p);     /* Mean free path of Bogoliubov quasiparticles */
double he3_tau_dperp_(double *ttc, double *p); /* Spin diffusion perp transport time, s */
double he3_tau_dpar_(double *ttc, double *p);  /* Spin diffusion parallel transport time, s */
double he3_diff_hperp_zz_(double *ttc, double *p);    /* Hydrodynamic spin diffusion D_perp_zz, cm2/s*/
double he3_diff_hpar_zz_(double *ttc, double *p);     /* Hydrodynamic spin diffusion D_par_zz, cm2/s*/
double he3_diff_perp_xx_(double *ttc, double *p, double *nu0);     /* Spin diffusion Re D_perp_xx, cm2/s */
double he3_diff_perp_xx_im_(double *ttc, double *p, double *nu0);  /* Spin diffusion Im D_perp_xx, cm2/s */
double he3_diff_perp_zz_(double *ttc, double *p, double *nu0);     /* Spin diffusion Re D_perp_zz, cm2/s */
double he3_diff_perp_zz_im_(double *ttc, double *p, double *nu0);  /* Spin diffusion Im D_perp_zz, cm2/s */
double he3_diff_par_xx_(double *ttc, double *p, double *nu0);      /* Spin diffusion D_par_xx, cm2/s */
double he3_diff_par_zz_(double *ttc, double *p, double *nu0);      /* Spin diffusion D_par_zz, cm2/s */

/* Other */
double he3_xigl_(double *ttc, double *p); /* Extrapolated GL coherence length, cm*/
double he3_vneq_(double *ttc, double *p, double *omega, double *r); /* Equilibrium vortex number */

double he3_swvel_(double *P, double *ttc);     /* Osheroff's spin wave vel. [cm/s] vs P [bar], T [mK] */
double he3_swvel_par_(double *P, double *ttc); /* Perp Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */
double he3_swvel_per_(double *P, double *ttc); /* Parallel Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */

double he3_ds_exp_(double *P, double *ttc); /* spin diffusion coeff. in superfluid He3 (measured) */
double he3_dn_exp_(double *P, double *ttc); /* spin diffusion in normal He3 */
double he3_d_exp_(double *P, double *ttc);  /* combined normal + superfluid spin diffusion */

double he3_tau_r_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */
double he3_tau_f_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */

/* Normal */
double he3_cv_n_(double *t, double *v);    /* C_v */

/* ROTA */

double rota_c_ns_(double *t, double *i); /* Nuclear stage heat capacity [J/K] vs T[K] and I[A] */
double rota_fork_cal_(double *w, double *p, double *n); /* Calibration of fork N, T/Tc, vs width (Hz) and P (bar) */

