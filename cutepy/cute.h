void runfull_xcorr(double *output, int ngal1,
                   double *ra1, double *dec1, double *z1, double *w1,
                   int ngal2,
                   double *ra2, double *dec2, double *z2, double *w2,
                   int nrands1,
                   double *randr1, double *randec1, double *randz1, double *randw1,
                   int nrands2,
                   double *randr2, double *randec2, double *randz2, double *randw2,
                   double r_p_max, int pi_max,
                   int r_p_nbins, int ndecades,
                   int logtrue, double O_M, double O_L);

void run_auto_3d_ps(double *output, int ngals,
                     double *ra, double *dec, double *z, double *w,
                     double r_p_max, int pi_max,
                     int r_p_nbins, int ndecades,
                     int logtrue, double O_M, double O_L);

void run_cross_3d_ps(double *output, int ngal1,
                     double *ra1, double *dec1, double *z1, double *w1,
                     int ngal2,
                     double *ra2, double *dec2, double *z2, double *w2,
                     double r_p_max, int pi_max,
                     int r_p_nbins, int ndecades,
                     int logtrue, double O_M, double O_L);
