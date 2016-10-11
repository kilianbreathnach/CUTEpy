#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "src/define.h"
#include "src/common.h"
#include "cute.h"


typedef struct {
  double dim1_max;
  int dim1_nbin;
  double dim2_max;
  int dim2_nbin;
  double dim3_min;
  double dim3_max;
  int dim3_nbin;
  int logbin;
  int n_logint;
} Binner;


void process_binner(Binner binner)
{
  //////
  // Check that binning options make sense
  if(binner.logbin<0) {
    fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
    exit(1);
  }
  if(binner.logbin) {
    if(binner.n_logint<0) {
      fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
      exit(1);
    }
  }
  if(binner.dim1_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim1 %d\n",binner.dim1_nbin);
    exit(1);
  }
  if(binner.dim2_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim2 %d\n",binner.dim2_nbin);
    exit(1);
  }
  if(binner.dim3_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim3 %d\n",binner.dim3_nbin);
    exit(1);
  }
  if(binner.dim1_max<=0) {
    fprintf(stderr,"CUTE: wrong dim1_max %lf\n",binner.dim1_max);
    exit(1);
  }
  if(binner.dim2_max<=0) {
    fprintf(stderr,"CUTE: wrong dim2_max %lf\n",binner.dim2_max);
    exit(1);
  }
  if((binner.dim3_max<=0)||(binner.dim3_min<0)||
     (binner.dim3_max<=binner.dim3_min)) {
    fprintf(stderr,"CUTE: wrong boundaries for dim3 (%lf , %lf)\n",
	    binner.dim3_min,binner.dim3_max);
    exit(1);
  }

  logbin=binner.logbin;
  n_logint=binner.n_logint;
  if(corr_type==0) {
    nb_dz=binner.dim1_nbin;
    i_dz_max=1./binner.dim1_max;
  }
  else  if(corr_type==1) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
  }
  else if(corr_type==2) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./binner.dim1_max;
    log_r_max=log10(binner.dim1_max);
  }
  else if((corr_type==3)||(corr_type==8)||(corr_type==9)) {
    nb_rl=binner.dim2_nbin;
    i_rl_max=1./binner.dim2_max;
    nb_rt=binner.dim1_nbin;
    i_rt_max=1./binner.dim1_max;
    if(logbin)
      log_rt_max=log10(binner.dim1_max);
  }
  else if(corr_type==4) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./binner.dim1_max;
    log_r_max=log10(binner.dim1_max);
    nb_mu=binner.dim2_nbin;
  }
  else if(corr_type==5) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
    nb_dz=binner.dim2_nbin;
    i_dz_max=1./binner.dim2_max;
    nb_red=binner.dim3_nbin;
    i_red_interval=1./(binner.dim3_max-binner.dim3_min);
    red_0=binner.dim3_min;
  }
  else if(corr_type==6) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
    nb_red=binner.dim3_nbin;
    i_red_interval=1./(binner.dim3_max-binner.dim3_min);
    red_0=binner.dim3_min;
  }
  else {
    fprintf(stderr,"WTF!?\n");
    exit(1);
  }
}


Catalog make_cat(double *ras, double *decs, double *zs, double *ws,
                 int ngals, np_t *sum_w, np_t *sum_w2) {
    /*
     * Makes a CUTE catalogue with the data
     */

    int ii;
    Catalog cat;

    cat.np=ngals;
    cat.red=(double *)my_malloc(cat.np*sizeof(double));
    cat.cth=(double *)my_malloc(cat.np*sizeof(double));
    cat.phi=(double *)my_malloc(cat.np*sizeof(double));
#ifdef _WITH_WEIGHTS
    cat.weight=(double *)my_malloc(cat.np*sizeof(double));
#endif //_WITH_WEIGHTS


    //Read galaxies in mask
    *sum_w=0;
    *sum_w2=0;
    for(ii=0;ii<ngals;ii++) {
      double zz,cth,phi,weight;

      zz = zs[ii];
      if(zz<0) {
        fprintf(stderr,"Wrong redshift = %lf %d\n",zz,ii+1);
        exit(1);
      }

      cth = cos(DTORAD*(90-ras[ii]));
      if((cth>1)||(cth<-1)) {
        fprintf(stderr,"Wrong cos(theta) = %lf %d\n",cth,ii+1);
        exit(1);
      }

      phi=DTORAD*decs[ii];
      phi=wrap_phi(phi);

      weight = ws[ii];

      cat.red[ii]=zz;
      cat.cth[ii]=cth;
      cat.phi[ii]=phi;
#ifdef _WITH_WEIGHTS
      cat.weight[ii]=weight;
      (*sum_w)+=weight;
      (*sum_w2)+=weight*weight;
#else //_WITH_WEIGHTS
      (*sum_w)++;
      (*sum_w2)++;
#endif //_WITH_WEIGHTS
    }

    printf("  Effective n. of particles: %lf\n",(*sum_w));
    printf("\n");
    return cat;
}


void make_CF_cross(histo_t D1D2,histo_t D1R2,
                   histo_t R1D2,histo_t R1R2,
                   np_t sum_wd,np_t sum_wd2,
                   np_t sum_wr,np_t sum_wr2,
                   np_t sum_wd_2,np_t sum_wd2_2,
                   np_t sum_wr_2,np_t sum_wr2_2,
		   double *corr,double *ercorr)
{
    //////
    // Creates correlation function and poisson errors
    // from pair counts DD, DR and RR
    double ed1d2,ed1r2,er1d2,er1r2;
    double dd1d2,dd1r2,dr1d2,dr1r2;
    double norm_d1d2=((double)sum_wd)*sum_wd_2;
    double norm_d1r2=((double)sum_wd)*sum_wr_2;
    double norm_r1d2=((double)sum_wr)*sum_wd_2;
    double norm_r1r2=((double)sum_wr)*sum_wr_2;

    ed1d2=1./sqrt((double)D1D2);
    ed1r2=1./sqrt((double)D1R2);
    er1d2=1./sqrt((double)R1D2);
    er1r2=1./sqrt((double)R1R2);
    dd1d2=(double)(D1D2/norm_d1d2);
    dd1r2=(double)(D1R2/norm_d1r2);
    dr1d2=(double)(R1D2/norm_r1d2);
    dr1r2=(double)(R1R2/norm_r1r2);

    double c,ec;
    if((D1D2==0)||(D1R2==0)||(R1D2==0)||(R1R2==0)) {
      c=0;
      ec=0;
    }
    else {
      c=(dd1d2-dd1r2-dr1d2+dr1r2)/dr1r2;
      ec=(1+c)*(ed1d2+ed1r2+er1d2+er1r2);
    }

    *corr=c;
    *ercorr=ec;
}

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
                   double O_M, double O_L) {
    /*
     * This computes xi(r_p, pi) for two data sets with two
     * sets of random points. It does so using a modified version
     * of CUTE. The need for a param file is done away with and
     * most parameters can be set with the python call. The code
     * assumes that bins are to be linear in pi and logspace in
     * r_p. The code returns a set of histograms for the pair
     * counts with correlation and error estimates.
     */

    // set the cosmology from given values
    omega_M = O_M;
    omega_L = O_L;

    // CUTE variables that we need for the calculation
    np_t sum_wd,sum_wd2,sum_wr,sum_wr2,sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
    Catalog cat_dat_1, cat_dat_2, cat_ran_1, cat_ran_2;

    Box3D *boxes_dat_1, *boxes_ran_1, *boxes_dat_2, *boxes_ran_2;
    int *indices_dat_1, *indices_ran_1, *indices_dat_2, *indices_ran_2;
    int nfull_dat_1, nfull_ran_1, nfull_dat_2, nfull_ran_2;
    int redo_boxing;

    // assume that pi_max is a whole number and we take 1 Mpc bins

    // Set up a binner object
    Binner binner;
    binner.logbin = 1;
    binner.n_logint = r_p_nbins / ndecades;
    binner.dim1_nbin = r_p_nbins;
    binner.dim2_nbin = pi_max;
    binner.dim3_nbin = 1;
    binner.dim1_max = r_p_max;
    binner.dim2_max = pi_max;
    binner.dim3_min = 0.08;
    binner.dim3_max = 0.45;  // this is a dummy bin

    // consistency check and setting more variables
    process_binner(binner);

    // Set up the histograms for the measurements.
    histo_t *D1D2=(histo_t *)my_calloc(nb_rt*nb_rl,sizeof(histo_t));
    histo_t *D1R2=(histo_t *)my_calloc(nb_rt*nb_rl,sizeof(histo_t));
    histo_t *R1D2=(histo_t *)my_calloc(nb_rt*nb_rl,sizeof(histo_t));
    histo_t *R1R2=(histo_t *)my_calloc(nb_rt*nb_rl,sizeof(histo_t));

    timer(4);

    set_r_z();

    // Now to create CUTE catalogues with the arrays
    cat_dat_1 = make_cat(ra1, dec1, z1, w1, ngal1,
                         &sum_wd, &sum_wd_2);
    cat_dat_2 = make_cat(ra2, dec2, z2, w2, ngal2,
                         &sum_wd2, &sum_wd2_2);
    cat_ran_1 = make_cat(randr1, randec1, randz1, randw1, nrands1,
                         &sum_wr, &sum_wr_2);
    cat_ran_2 = make_cat(randr2, randec2, randz2, randw2, nrands2,
                         &sum_wr2, &sum_wr2_2);

    // set up the boxes for pair-counting
    init_3D_params_cross(cat_dat_1,cat_ran_1,cat_dat_2,cat_ran_2,3);

    boxes_ran_1=mk_Boxes3D_from_Catalog(cat_ran_1,&indices_ran_1,&nfull_ran_1);
    boxes_ran_2=mk_Boxes3D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    boxes_dat_2=mk_Boxes3D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_ran_1);

    // Where the meat of the crunching happens.
    // This is the same code as my original edit of CUTE.
    printf("\n");

    printf(" - Cross-correlating randoms \n");
    timer(0);
    cross_3d_ps_bf(nfull_ran_1,indices_ran_1,
  		 boxes_ran_1,boxes_ran_2,R1R2);
    timer(2);
    printf(" - Cross-correlating R1D2 \n");
    cross_3d_ps_bf(nfull_dat_2,indices_dat_2,
                   boxes_dat_2,boxes_ran_1,R1D2);
    timer(2);

    redo_boxing = cross_followup(cat_dat_1, 3);

    boxes_dat_1=mk_Boxes3D_from_Catalog(cat_dat_1,&indices_dat_1,&nfull_dat_1);
    free_Catalog(cat_dat_1);

    if (redo_boxing)
    {
      printf(" - ** Have to rebox :( \n");
      boxes_ran_2=mk_Boxes3D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
      free_Catalog(cat_ran_2);
      boxes_dat_2=mk_Boxes3D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
      free_Catalog(cat_dat_2);
    }
    else
    {
      free_Catalog(cat_ran_2);
      free_Catalog(cat_dat_2);
    }

    timer(2);
    printf(" - Cross-correlating data \n");
    cross_3d_ps_bf(nfull_dat_1,indices_dat_1,
                   boxes_dat_1,boxes_dat_2,D1D2);
    timer(2);
    printf(" - Cross-correlating D1R2 \n");
    cross_3d_ps_bf(nfull_dat_1,indices_dat_1,
                   boxes_dat_1,boxes_ran_2,D1R2);
    timer(1);

    printf("\n");

    // Fill the output array that the python wrapper gets.
    // (write_CF used to be called here)
    int ii;
    for(ii=0;ii<nb_rt;ii++) {
      int jj;
      double rt=pow(10,((ii+0.5)-nb_rt)/n_logint+log_rt_max);
      for(jj=0;jj<nb_rl;jj++) {
        double corr,ercorr;
        double rl=(jj+0.5)/(nb_rl*i_rl_max);
        int ind=jj+nb_rl*ii;
        make_CF_cross(D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind],
                      sum_wd,sum_wd2,sum_wr,sum_wr2,
                      sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2,
                      &corr,&ercorr);
        output[4 * (ii * pi_max + jj)] = rl;
        output[4 * (ii * pi_max + jj) + 1] = rt;
        output[4 * (ii * pi_max + jj) + 2] = corr;
        output[4 * (ii * pi_max + jj) + 3] = ercorr;
      }
    }

    printf("*** Cleaning up\n");
    free_Boxes3D(n_boxes3D,boxes_dat_1);
    free_Boxes3D(n_boxes3D,boxes_dat_2);
    free_Boxes3D(n_boxes3D,boxes_ran_1);
    free_Boxes3D(n_boxes3D,boxes_ran_2);
    free(indices_dat_1);
    free(indices_ran_1);
    free(indices_dat_2);
    free(indices_ran_2);
    end_r_z();
    free(D1D2);
    free(D1R2);
    free(R1D2);
    free(R1R2);

    // ...and return the arrays
}
