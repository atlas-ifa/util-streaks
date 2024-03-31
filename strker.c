/* Cross-correlate an image with a tilted streak kernel, find streaks */
/* E.g.
  strker /atlas/diff/02a/58071/02a58071o0503c.diff.fz \
         /atlas/red/02a/58071/02a58071o0503c.fits.fz \
         /tmp/foo.s -len 10 -nlam 1 -bin 2 | tee /tmp/det.all

  fm 1,2,10,11,12 /tmp/det.all -tol 2,g,g,g,g -grp -grpfile /tmp/foo.grp -grpfmt "%8.1f" | sort -k 3g > /tmp/foo
  awk '$2>800 && $2<930 && $3>4200 && $3<4500{print $0}' /tmp/foo


*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "/atlas/src/trunk/util/sortlib/tsort.h"

#include "fitsio.h"

#include "psf2d.h"

#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))
// #define ABS(x)  ((x)<0?(-(x)):(x))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define LN10 (2.302585092994)	/* ln(10) */

/****Declared parameters and functions needed by Ari's code****/
#define PI 3.14159265358979
#define MAXMEM 40000000000.0 /*Refuse to allocate a matrix or vector that will
                             use more than 40GB of memory.*/
#define SIGFWHM (double)2.354820045 /*Factor by which the FWHM of a 
                                      Gaussian exceeds sigma.*/
#define TRAILCENGFAC 3.0 /*Factor by which centering half-widths should
                   exceed Gaussian sigma*/
#define TRAILCENRLEV 0.2
#define TRAILCENPAPIXSAMP 0.5 /*Sampling interval in pixels around circumference
                              for position angle finding*/
#define TRAILCENPAPIXRANGE 4.0 /*Range in pixels around circumference
                              for position angle finding*/
#define TRAILCENITNUM 3 /*Number of iterations for trail centroiding*/

#define TRAILWIDTH_FWHMSCALE 5.0 /*Factor by which the trail extraction
                                   half-width, and the maximum FWHM probed,
                                   exceeds the input FWHM.*/
#define TRAILWIDTH_FWHMSAMP 0.1 /*Trail width sampling interval in units of the FWHM*/
#define MAXTRAILWIDTH 1.5 /*Maximum trail width in units of the FWHM*/
#define TRAILMAXFADE 0.1 /*We give up tracing a trail if its brightness
                           drops to less than this fraction of the nominal
                           central flux.*/

#define MAXTRAILLEN 1000.0 /*Maximum trail length probed by endfind*/
#define MAXNTPH 46000 /*Maximum number of total triggers, to avoid
                        exceeding the limits of calloc() when allocating
                        the pair arrays with size ntph*ntph*/
#define MAXNWHISK 30000 /*Maximum number of detections fed into whiskers*/

#define PAFIND_DEBUG 0
#define ENDFIND_DEBUG 0
#define PREMATCH_DEBUG 0
#define DEBUGWHISK 0

#define ARRAYERR(x) {fprintf(stderr,"ERROR allocating array called \'%s\': ABORTING\n",x); fflush(stderr); return(4);}
#define DSQUARE(x) ((double)(x)*(double)(x))
/****End parameters and functions needed by Ari's code****/

static int VERBOSE=0;		/* verbosity level */

/* Timing diagnostics */
#define BLAB(msg) { gettimeofday(&tv1, NULL); telapse = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec); if(VERBOSE>0) fprintf(stderr, "%7.3f - %s\n", telapse, msg); }


/* Macros to manage FITS error handling */
#define FRE { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); } ; status = 0; }

#define FFATAL { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); exit(1); }  }


/* Structure to manage detections */
typedef struct {
   int hlen;		/* correlate with (2*hlen+1) */
   double phicorr;	/* [deg] correlation angle */
   int ndet;		/* number of detections */
   PSF2D_PARAM *wpar;	/* detection fit parameters */
   PSF2D_FLUX *wflux;	/* detection fit flux */
} XCDET;

/* Structure for stars in the field */
typedef struct {
   double ra;		/* [rad] RA */
   double dec;		/* [rad] Dec */
   double x;		/* [pix] x */
   double y;		/* [pix] y */
   double m;		/* [AB] magnitude */
} STAR;

/* Probability vector entries: 0-999 */
#define Ptr 0	// Ptr = "transient" real object, possibly moving
#define Pvr 1	// Pvr = "variable" coincident with a bright star
#define Pmv 2	// Pmv = "moving" assesses reality of trail
#define Psc 3	// Psc = "scar" clutter from bad star subtration or diffraction
#define Pbn 4	// Pbn = "burn" is a persistence trail away from serial register
#define Pcr 5	// Pcr = "cosmic ray"


/* Structure to characterize detections */
typedef struct {
   int ndet;		/* how many (duplicate) detections? */
   int lencorr;		/* correlate with (2*hlen+1) */
   double phicorr;	/* [deg] correlation angle */
   PSF2D_PARAM *wpar;	/* detection fit parameters */
   PSF2D_FLUX *wflux;	/* detection fit flux */
   double x0;		/* fitted center */
   double y0;		/* fitted center */
   double len;		/* fitted length = wpar->major */
   double phi;		/* fitted angle = wpar->phi */
   double dphi;		/* fitted angle uncertainty */
   double th1;		/* intercept #1 around circle */
   double th2;		/* intercept #2 around circle */
   double dtheta;	/* intercept around circle uncertainty */
   double lam1;		/* streak-intercept #1 distance */
   double lam2;		/* streak-intercept #2 distance */
   double dlambda;	/* streak-intercept distance uncertainty */
   double flux;	        /* Median flux measurement over a group*/
} WHISKER;

/* Structure to group detections */
typedef struct {
   int n;		/* number of whiskers */
   int ndet;		/* number of detections */
   int P[6];		/* Classification probabilities */
   WHISKER **w;		/* array of grouped whisker pointers */
   double x0;		/* ave center */
   double y0;		/* ave center */
   double len;		/* ave length */
   double phi;		/* ave angle */
   double flux;		/* ave flux */
} STREAK;

/* Test pixel i,j to see if it passes the trigger */
void one_trigger(int nx, int NX, int ny, float *data,
		 int *npix, int ymax, double *skymed, double *skyrms);

/* Do the work on just one object */
int one_object(int nx, int NX, int ny, float *data, float *wgt,
	       int i, int j, int invert, int trail, int force, int okfit,
	       PSF2D_PARAM *wpar, PSF2D_FLUX *wflux);
/* Report just one object */
void one_report(int nx, int NX, int ny, float *data,
		int i, int j, int trail, PSF2D_PARAM *wpar, PSF2D_FLUX *wflux,
		FILE *fp, FILE *fpmom);

/* Set trigger parameters */
void tpctrl_search(int dx0, int dx1, int dy0, int dy1, float bad);

/* Set trigger parameters */
void tpctrl_trigger(double min, double max, double rad, double sig);

/* Set acceptance parameters */
void tpctrl_accept(double move, double wmax, double wmin, double hgt,
		   double chin, double sig);

/* Set fit parameters */
void tpctrl_fit(int npar, int ap, int sky);

/* Look for peaks in sigma image */
int sigma_srch(int nx, int ny, float *data, float *diff, int bin,
	       double sigma, int hlen, double fwhm, double phi/*rad*/, 
               double maxtraillen, int fixpsf, int checkwidth, int checklength, 
               PSF2D_PARAM **wpar, PSF2D_FLUX **wflux);

/* Ari's new version of sigma_srch, which runs some
   diagnostics on the unbinned difference image as well*/
int sigma_srch_ari(int nx, int ny, float *data, float *diff, int bin,
		   double sigma, int hlen, double fwhm,double phi/*rad*/,
                   double maxtraillen, int fixpsf, int checkwidth, int checklength, 
                   PSF2D_PARAM **wpar, PSF2D_FLUX **wflux);

/* Consolidate all the whisker entries according to nearly dup tph entries */
int whisk_cull(int bin, int nphi, XCDET *xcdet, WHISKER **whisk);


/* Group all overlapping detections */
int whiskers(int n, WHISKER *w, int nx, int ny,
	     double dmax, double tmax, STREAK **grp, double dcirc);

/* Fit each streak using initial conditions */
int streak_fit(int nx, int ny, float *data, float *wgt, float badata,
	       int bin, double len/*unbinned*/, double phi/*rad*/, 
               double fwhm/*unbinned*/, int fixpsf,
	       int ngrp, STREAK *strk, int *sortind, double starttime, double timeout,
	       PSF2D_PARAM **wpar, PSF2D_FLUX **wflux, FILE *fp);

/* Classify each streak */
void streak_class(int nx, int ny, int x0, int y0,
		  double zp, double fw, double sat,
		  int nstar, STAR *star,
		  int nstrk, STREAK *streak, PSF2D_PARAM *wpar);


/* Convolve with 1/cos/sin with nlam wavelength at phi and length +/-hlen */
void xdrift_tcs(double phi, int hlen, double nlam, 
	    int nx, int ny, float *H, float *C, float *S);
void ydrift_tcs(double phi, int hlen, double nlam,
	    int nx, int ny, float *H, float *C, float *S);

/* Convolve image with top hat at phi and length +/-hlen */
void xdrift_t(double phi, int hlen, int nx, int ny, float *D);
void ydrift_t(double phi, int hlen, int nx, int ny, float *D);



/* Convolve a 1D array by a Gaussian with sigma: 14nsec/pix */
void gconv(int n, float *a, double sig);

/* Fast median. Does not alter input array. */
float quickmedian(int n, float *key);

/* Median of an array of doubles */
double median(int n, double *buf, int *idx);

/* Count of factors that divide n, from 1 to maxdivisor: nfact[0:maxdivisor] */
void factors(int n, int maxdivisor, int *nfact);

/* Linearly interpolate src(nx/xbin,ny/ybin) -> dest(nx,ny) in ncx,ncy chunks */
/* ncx,ncy must divide the src and dest dimensions evenly */
void linterp_chunk(int nx, int ny, int xbin, int ybin, int ncx, int ncy,
		   float *src, float *dest);
/* Linearly interpolate an array into a bigger one */
void linterp(int mx, int my, int MX, float *src,
	     int nx, int ny, int NX, float *dest);

/* Pick binning factors for background subtraction */
void pick_bin(int nx, int ny, int chx, int chy, int *chunkok, int *bx, int *by);

/* Median bin an image, ignore pixels at badata */
void medianbin(int nx, int ny, float *D, float badata, 
	       int bx, int by, float **B);

/* Clip an image using Laplacian and saturation */
void mowgrass(int nx, int ny, float *R, float badata, int ds,
              double saturate, float **G, double *sky);

/* Sort an array, carry a second */
int qsort2(int n, double *x, int *idx);

/* 2D smooth an image by sig */
void smooth(double sig, int nx, int ny, float *D);

/* Clip and image and bin */
void clipbin(double clipmin, double clipmax, int bin, float badata, int nx, int ny, float *Din, float *Dout);

/* Read a file of stars: ra, dec, x, y, m */
int read_star(char *fname, STAR **star);

/* Print out cfitsio error messages and exit program */
void printerror(int status);

/*Ari's array allocation routines*/
float *vec(int np);
float **mat(int nx,int ny);
double *dvec(int np);
double **dmat(int nx,int ny);
int free_mat(float **mt,int nx,int ny);
int free_dmat(double **mt,int nx,int ny);

/*Ari's trim mean*/
int creepvec01ds_S(double *invec,int pnum,double rejfac,double *output,double *error);

/*Ari's trail reality check*/
int trailsquish_S(float *DIFF,int bin,int nxbin,int nybin,double tpx,double tpy,int hlen,double phi,double fwhm,double *trailflux1,double *trailflux2,double *skylev,double *skyrms);

/*Ari's position angle refinement code*/
int phifind(float *image,int nx,int ny,double tpx,double tpy,double traillen,double *phi,double phirange,double fwhm,double gfac,float pixsamp,float RLEV);

/*Ari's cross-trail centroiding code*/
int crosstrailcen01_S(float *image,int nx,int ny,double *tpx,double *tpy,double traillen,double phi,double fwhm,double gfac,int itnum,float RLEV);

/*Ari's code for estimating the width of a trail,
in order to reject broad, fuzzy trails that cannot have
been made by asteroids.*/
int trailwidth(float *DIFF,int nx,int ny,double tpx,double tpy,int halflen,double phi,double fwhm,double *trailwidth,double *flux);

/*Ari's code for trying to find the ends of the trails*/
int trailendfind(float *DIFF,int nx,int ny,double tpx,double tpy,int halflen,int halfsearchlen,int maxhalflen,double phi,double fwhm,double *endx,double *endy);

int writeemptyfile(char strkfile[]);

/* Append a FITS file as an extension at the end of another one */
int main(int argc, char **argv)
{
   int i, j, k, bitpix, status=0, naxis, anynull, nstrk, nstrk0, nstar, *idx;
   int hlen, dcs, nx, ny, bx, by, chx, chy, bin, bkgsub, nphi, fixpsf;
   int nxbin, nybin, nxsub, nysub, x0sub, y0sub, chunkok, nsquash;
   char *diffile, *redfile, *newfile, *newbang, *starfile;
   char line[1024];
   fitsfile *fnew, *ff;   /* pointers to the FITS files */
   long int naxes[3], fpixel[3];
//   size_t nbyte;
   LONGLONG nelements;
   double len, phi, xcphi, nlam, csamp, pi=4*atan(1.0), dr=atan(1.0)/45;
   float *DIFF, *IMG, *R, *D, *T, *C, *S, *B, *G, badata_dif, badata_img;
   double difclipmin, difclipmax, sig, redclipmin, redclipmax, badata_replacement;
   double eadu, rtn, scsum, rms, rdnoise;
   double majmin, majmax, minmin, minmax, chinmax, *buf;
   double srchsig, imgsky, xtra,maxtraillen;
   double *sortflux,starttime,timeout=1e30;
   int *sortind,maxdets=1000000;
   char *tphfile, *whiskfile, *grpfile, *strkfile;
   FILE *fp;

   XCDET *xcdet;
   int ntph, nwhisk, dowhisk,do_oldsigma_srch=0;
   int checkwidth=1;
   int checklength=1;
   int freepsf=0;
   double dmax, tmax, unsquashable, unsquashed;
   WHISKER *whisk;
   STREAK *streak;
   STAR *star;
   
   PSF2D_PARAM *wpar;	/* detection fit parameters */
   PSF2D_FLUX *wflux;	/* detection fit flux */

   struct timeval tv0, tv1;
   double telapse;

   gettimeofday(&tv0, NULL);
   starttime =  tv0.tv_sec + 1e-6*tv0.tv_usec;

   if(argc < 3) {
      printf("strker: insufficient arguments\n");
      printf("syntax: strker diffile redfile [options]\n");
      exit(1);
   }
   
   diffile = argv[1];
   redfile = argv[2];
   newfile = NULL;
   starfile = NULL;		/* input file of stars */
   nstar = 0;
   tphfile = NULL;		/* tphot results from sigma images */
   whiskfile = NULL;		/* whiskers of all detections */
   grpfile = NULL;		/* grouped whiskers initial conditions */
   strkfile = "-";		/* final streak fit */
   badata_img = 0.0;		/* Bad data value in image */
   badata_dif = -31415.0;	/* Bad data value in difference */
//   badata_replacement = 1e10;	/* Bad data replacement value for difference */
   badata_replacement = 0;	/* Bad data replacement value for difference */

/* Parse the options */
   bin = 2;		/* Bin input image by bin */
   difclipmin = -1000;	/* Clip diff output, post-bin */
   difclipmax =  80000;	/* Clip diff output, post-bin */
   redclipmin = -1000;	/* Clip red output, post-bin */
   redclipmax = 80000;	/* Clip red output, post-bin */
/* These chunks are here to compensate for hotpants bad sky sub of sta1600 */
/* Alter if hotpants works in different chunk sizes */
   chx = 8;		/* Respect chx x chy chunks for linterp */
   chy = 8;		/* Respect chx x chy chunks for linterp */
   sig = 1.0;		/* Smooth by sig, post-bin */
   srchsig = 5.0;	/* Search sigma */

   nxsub = nysub = 0;	/* Subarray size */
   x0sub = y0sub = 0;	/* Subarray origin (0 based) */

   dowhisk = 1;		/* Group detections into whiskers? */

   nsquash = 4;         /* Number of extra factors of image divided into sigma */
   unsquashable = 50.0; /* Sigma limit above which sources are immune from squashing */

   phi = -1000;		/* Streak direction wrt x [rad] */
   csamp = 1.0;		/* Sample circumference every csamp */
   len = -1;		/* Full correlation kernel length */
   hlen = -1;		/* Correlation length is bin*(2*hlen+1) */
   fixpsf = -1;		/* Fix the tphot PSF fit to a given len,phi? */
   nlam = 1;		/* Number of cos/sin wavelengths in (2*hlen+1) */
   dcs = 0;		/* dcs=0,1,2,3,4,5 for sig,imgtop,difftop,cos,img,dif */
   bkgsub = 1;		/* Subtract the background? */
   eadu = 2.5;		/* (effective) electrons per ADU */
   rdnoise = 12;	/* [e-] read noise */

   minmin = 0.5;	/* [pix] minimum fwmin to report streak */
   minmax = 10.0;	/* [pix] maximum fwmin to report streak */
   majmin = 0.0;	/* [pix] minimum fwmax to report streak */
   majmax = 2000;	/* [pix] maximum fwmax to report streak */
   chinmax = 1e5;	/* [pix] maximum chin to report streak */

   dmax = 1.0;          /*along-streak duplicate-matching distance,
                          in units of quadrature-sum streak length,
                          approximately (it's complicated).*/
   maxtraillen = MAXTRAILLEN;

/* Parse arguments */
   for(i=3; i<argc; i++) {
      if(strcmp(argv[i], "-phi") == 0) {   // [deg] Streak direction wrt x
	 sscanf(argv[++i], "%lf", &phi);
	 if(phi < -180) phi = phi + 360;
	 if(phi >  180) phi = phi - 360;
	 if(phi < -45) phi = phi + 180;
	 if(phi > 135) phi = phi - 180;
	 phi *= pi / 180;	// -pi/4 < phi < 3*pi/4

      } else if(strcmp(argv[i], "-cstep") == 0) { // Circumference sampling
	 sscanf(argv[++i], "%lf", &csamp);

      } else if(strcmp(argv[i], "-len") == 0) {   // Full corr kernel length
	 sscanf(argv[++i], "%lf", &len);

      } else if(strcmp(argv[i], "-hlen") == 0) {  // Corr ker half len after bin
	 sscanf(argv[++i], "%d", &hlen);

      } else if(strcmp(argv[i], "-nlam") == 0) {  // number of sin/cos periods
	 sscanf(argv[++i], "%lf", &nlam);

      } else if(strcmp(argv[i], "-bin") == 0) {   // bin factor
	 sscanf(argv[++i], "%d", &bin);

      } else if(strcmp(argv[i], "-smooth") == 0) {   // post-bin smooth
	 sscanf(argv[++i], "%lf", &sig);

      } else if(strcmp(argv[i], "-nsquash") == 0) {   // image division factors
	 sscanf(argv[++i], "%d", &nsquash);

      } else if(strcmp(argv[i], "-srchsig") == 0) {   // search depth
	 sscanf(argv[++i], "%lf", &srchsig);

      } else if(strcmp(argv[i], "-sig") == 0) {       // search depth
	 sscanf(argv[++i], "%lf", &srchsig);

      } else if(strcmp(argv[i], "-badata") == 0) {   // diff badata replacement value
	 sscanf(argv[++i], "%lf", &badata_replacement);

      } else if(strcmp(argv[i], "-fixpsf") == 0) {   // fixed PSF?
	 fixpsf = 1;
	 freepsf = 0;
      } else if(strcmp(argv[i], "-freepsf") == 0) {   // fixed PSF?
	 fixpsf = 0;
	 freepsf = 1;
      } else if(strcmp(argv[i], "-subarray") == 0) {   // x0,nx,y0,ny
	 if(sscanf(argv[++i], "%d,%d,%d,%d", &x0sub, &nxsub, &y0sub, &nysub) != 4) {
	    printf("Error parsing x0,nx,y0,ny from %s\n", argv[i]);
	    exit(1);
	 }

      } else if(strcmp(argv[i], "-clip") == 0) {  // clip levels for image
	 sscanf(argv[++i], "%lf,%lf", &redclipmin, &redclipmax);

      } else if(strcmp(argv[i], "-CLIP") == 0) {  // clip levels for diff
	 sscanf(argv[++i], "%lf,%lf", &difclipmin, &difclipmax);

      } else if(strcmp(argv[i], "-dmax") == 0) {  // along-streak duplicate-matching distance
	 sscanf(argv[++i], "%lf", &dmax);

      } else if(strcmp(argv[i], "-maxtrail") == 0) {  // maximum trail length probed by endfind
	 sscanf(argv[++i], "%lf", &maxtraillen);

      } else if(strcmp(argv[i], "-timeout") == 0) {  // maximum run time in seconds
	 sscanf(argv[++i], "%lf", &timeout);

      } else if(strcmp(argv[i], "-maxdets") == 0) {  // maximum number of streak detections
	 sscanf(argv[++i], "%d", &maxdets);

      } else if(strcmp(argv[i], "-dcs") == 0) {   // dcs=0,1,2,3,4,5 for
	 sscanf(argv[++i], "%d", &dcs);           // sig,imgtop,diftop,cos,img,dif

      } else if(strcmp(argv[i], "-dfile") == 0) { // diagnostic file
	 newfile = argv[++i];

      } else if(strcmp(argv[i], "-tph") == 0) {   // tphot analysis of sigmas
	 tphfile = argv[++i];

      } else if(strcmp(argv[i], "-nowhisk") == 0) {   // all whiskers
	 dowhisk = 0;

      } else if(strcmp(argv[i], "-whisk") == 0) {   // all whiskers
	 whiskfile = argv[++i];

      } else if(strcmp(argv[i], "-grp") == 0) {   // grouped whiskers
	 grpfile = argv[++i];

      } else if(strcmp(argv[i], "-strk") == 0) {   // final streaks
	 strkfile = argv[++i];

      } else if(strcmp(argv[i], "-star") == 0) {   // star file
	 starfile = argv[++i];

      } else if(strcmp(argv[i], "-bkg") == 0) {	  // Don't subtract background?
	 bkgsub = 0;

      } else if(strcmp(argv[i], "-oldsigsrch") == 0) {	  // Use old version sigma_srch
	 do_oldsigma_srch = 1;

      } else if(strcmp(argv[i], "-checkwidthoff") == 0) {	  // Don't check for and reject broad, fuzzy trails
	 checkwidth = 0;

      } else if(strcmp(argv[i], "-checklengthoff") == 0) {	  // Don't attempt to refine trail lengths.
	 checklength = 0;

      } else if(strcmp(argv[i], "-verb") == 0) {
	 VERBOSE = 1;

      } else if(strcmp(argv[i], "-VERB") == 0) {
	 VERBOSE = 2;

      } else if(strcmp(argv[i], "-verbose") == 0) {
	 sscanf(argv[++i], "%d", &VERBOSE);

      } else {
	 printf("Unrecognized argument `%s'\n", argv[i]);
      }
   }

/* Fixed psf?  Implied by specifying both len and phi on the command line */
   if(len > 0 && phi > -100 && fixpsf < 0) fixpsf = 1;
   if(fixpsf < 0) fixpsf = 0;
/* Enable user to override fixpsf using the -freepsf flag*/
   if(freepsf==1) fixpsf=0; /*This will apply even if len and phi are fixed*/

/* Fill in len and/or hlen if not specified */
   if(len < 0 && hlen < 0) {
      hlen = 8;					/* default hlen=8 */
      len = bin * (2*hlen+1);			/* hlen is integer */
   } else {
      if(hlen > 0) len = bin * (2*hlen+1);	/* hlen is integer */
      if(len > 0) hlen = NINT((len/bin-1)/2);	/* len is real */
   }

   if(VERBOSE > 1) {
      printf("Starting analysis of %s at hlen %d\n", diffile, hlen);
   }
   BLAB("Args parsed");

   if(telapse>timeout) {
     /*Out of time already and nothing done. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) already before doing any processing\n",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }

   if(starfile != NULL && (nstar = read_star(starfile, &star)) < 0) exit(1);
   if(VERBOSE > 1) {
      printf("Read %d stars for classification\n", nstar);
   }


/////////////////////////////
// READ WRITE BOILER PLATE //
/////////////////////////////
/* Read the reduced file */
   status = 0;
   
   if(nxsub<=0 || nysub<=0) {
      strcpy(line, redfile);
   } else {
      sprintf(line, "%s[%d:%d,%d:%d]", redfile, x0sub+1,x0sub+nxsub,y0sub+1,y0sub+nysub);
   }
   if ( fits_open_image(&ff, line, READONLY, &status) ) FFATAL;

/* What about it? */
   if( fits_get_img_param(ff, 2, &bitpix, &naxis, naxes, &status) ) FFATAL;

   if(VERBOSE > 1) {
      printf("Input reduced (variance) file has %d axes %ld x %ld and bitpix %d\n", 
	     naxis, naxes[0], naxes[1], bitpix);
   }

   nx = naxes[0];
   ny = naxes[1];

   fpixel[0] = fpixel[1] = fpixel[2] = 1;    /* first pixel to read */
   nelements = naxes[0] * naxes[1];          /* number of pixels to read */

/* Allocate space to read data */
   IMG = (float *)calloc(naxes[0]*naxes[1], sizeof(float));

   if( fits_read_pix(ff, TFLOAT, fpixel, nelements, NULL,
		     IMG, &anynull, &status) ) {
      if(status != NUM_OVERFLOW) FFATAL;
      if(VERBOSE > 1) FRE;
      status = 0;
   }
   BLAB("Input reduced file read");

   if(telapse>timeout) {
     /*Out of time already and nothing done. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) after reading reduced image\n",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }

/* Read the difference file: use open_image vs open_file in case compressed */
   status = 0;
   if(nxsub<=0 || nysub<=0) {
      strcpy(line, diffile);
   } else {
      sprintf(line, "%s[%d:%d,%d:%d]", diffile, x0sub+1,x0sub+nxsub,y0sub+1,y0sub+nysub);
   }
   if ( fits_open_image(&ff, line, READONLY, &status) ) FFATAL;

/* What about it? */
   if( fits_get_img_param(ff, 2, &bitpix, &naxis, naxes, &status) ) FFATAL;

   if(VERBOSE > 1) {
      printf("Input diff file has %d axes %ld x %ld and bitpix %d\n", 
	     naxis, naxes[0], naxes[1], bitpix);
   }

   nx = naxes[0];
   ny = naxes[1];

   fpixel[0] = fpixel[1] = 1;                /* first pixel to read */
   nelements = naxes[0] * naxes[1];          /* number of pixels to read */


/* Allocate space to read data */
   DIFF = (float *)calloc(naxes[0]*naxes[1], sizeof(float));

   if( fits_read_pix(ff, TFLOAT, fpixel, nelements, NULL,
		     DIFF, &anynull, &status) ) {
      if(status != NUM_OVERFLOW) FFATAL;
      if(VERBOSE > 1) FRE;
      status = 0;
   }
   BLAB("Input difference file read");

   if(telapse>timeout) {
     /*Out of time already and nothing done. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) after reading difference image\n",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }

/* Get the bad data value */
   if( fits_read_key(ff, TFLOAT, "MASKVAL", &badata_dif, NULL, &status) ) FRE;
   
   double magzpt, etime, saturate, fwhm;
   /*Reasonable defaults:*/
   magzpt = 20.0;
   etime = 30.0;
   saturate = 65000.0;
   fwhm = 2.0;

/* Get the exposure time */
   if( fits_read_key(ff, TDOUBLE, "EXPTIME", &etime, NULL, &status) ) FRE;
   
/* Get the zeropoint */
   if( fits_read_key(ff, TDOUBLE, "MAGZPT", &magzpt, NULL, &status) ) FRE;
   
/* Adjust zeropoint to 1 ADU per exposure instead of 1 ADU per second */
   magzpt += 2.5/log(10) * log(MAX(1,etime));

/* Get the saturation value */
   if( fits_read_key(ff, TDOUBLE, "SATURATE", &saturate, NULL, &status) ) FRE;
   
/* Get the PSF */
   if( fits_read_key(ff, TDOUBLE, "FWHM", &fwhm, NULL, &status) ) FRE;
   
/* Get the sky background */
//   if( fits_read_key(ff, TDOUBLE, "BCKGND", &bckgnd, NULL, &status) ) FRE;
   


/* Close the input image */
   if( fits_close_file(ff, &status) ) FRE;
   BLAB("Input closed");
////////////////////////////////////////////////////////////////



/* Arrays
 * IMG = difference image
 * DIFF = difference image
 * R = reduced image, binned by bin, smoothed by sig
 * D = difference image, binned by bin
 * G = 'grass level', varying sky match to reduced image
 (temporary)
 * B = difference background, median binned by bx,by
 * C = difference background, linterp back to D
 (then)
 * D = difference image, binned by bin, C subtracted, smoothed by sig
 * B = copy of R, tophat correlated reduced image
 * T = copy of D, tophat correlated difference image
 * C = cosine correlated difference image
 * S = sine correlated difference image (temporary)
 (then)
 * S = sigma image: S = [T-sqrt(C*C+S*S)]/[B*eadu/(bin*bin*len)]
 */


/* Mash all bad difference pixels to badata_replacement? */
   if(badata_replacement < 1e9) {
      for(i=0; i<nx*ny; i++) if(DIFF[i] == badata_dif) DIFF[i] = badata_replacement;
   }

/* Clip and bin as requested */
   R = (float *)calloc(nx/bin*ny/bin, sizeof(float));
   clipbin(redclipmin, redclipmax, bin, badata_img, nx, ny, IMG, R);
   BLAB("Clip and bin reduced");

   D = (float *)calloc(nx/bin*ny/bin, sizeof(float));
   clipbin(difclipmin, difclipmax, bin, badata_dif, nx, ny, DIFF, D);
   BLAB("Clip and bin difference");

   if(telapse>timeout) {
     /*Out of time already and nothing done. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) after clipping and binning\n",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }


/* Update nx,ny to mean the binned resulting image dimensions */
   nxbin = nx / bin;
   nybin = ny / bin;

/* Select binning factors for sky subtraction */
   pick_bin(nxbin, nybin, chx, chy, &chunkok, &bx, &by);


/* Allocate space for cosine, and sine arrays (tophat is in D) */
   C = (float *)calloc(nxbin*nybin, sizeof(float));
   S = (float *)calloc(nxbin*nybin, sizeof(float));

/* Make a median-bin version for difference background */
   if(bkgsub) {
      medianbin(nxbin, nybin, D, badata_dif, bx, by, &B);
      BLAB("Median bin diff for background");

#if 0
      int j;
      for(j=0; j<nybin/by; j++) {
	 for(i=0; i<10; i++) printf(" %7.1f", B[i+j*(nxbin/bx)]);
	 printf("\n");
      }
#endif

/* Linearly interpolate back to full size in C array, respecting 8x8 chunks */
      if(chunkok) {
	 linterp_chunk(nxbin, nybin, bx, by, chx, chy, B, C);
/* Linearly interpolate back to full size in C array, no chunks */
      } else {
	 linterp((nxbin+bx-1)/bx, (nybin+by-1)/by, (nxbin+bx-1)/bx, B, nxbin, nybin, nxbin, C);
      }
      BLAB("Linterp diff back to full size");

/* Subtract the background */
      for(i=0; i<nxbin*nybin; i++) D[i] -= C[i];
      BLAB("Background subtracted difference");

      free(B);
   }
 

   mowgrass(nxbin, nybin, R, badata_img, 2, saturate, &G, &imgsky);
   medianbin(nxbin, nybin, G, badata_img, bx, by, &B);
   BLAB("Median red diff for background");
   linterp((nxbin+bx-1)/bx, (nybin+by-1)/by, (nxbin+bx-1)/bx, B, nxbin, nybin, nxbin, G);
   BLAB("Linterp red back to full size");
   free(B);

/* Smooth entire image by sig */
   smooth(sig, nxbin, nybin, R);
   BLAB("Reduced smoothed");

   smooth(sig, nxbin, nybin, D);
   BLAB("Difference smoothed");

   if(telapse>timeout) {
     /*Out of time already and nothing done. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) after smoothing\n",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }

/* Make a copy of reduced in case of muliple hlen and phi, B = 'blur' */
   B = (float *)calloc(nxbin*nybin, sizeof(float));

/* Make a copy of diff for muliple hlen and phi, T = 'tophat' */
   T = (float *)calloc(nxbin*nybin, sizeof(float));

///////////////////////////////////////////
// RETURN HERE FOR MULTIPLE HLEN AND PHI //
///////////////////////////////////////////
   
   nphi = 1;
   if(phi < -pi) nphi = hlen*pi / csamp + 1;
   xcdet = (XCDET *)calloc(nphi, sizeof(XCDET));

/* Total number of detections */
   ntph = 0;

   for(k=0; k<nphi; k++) {
      if(nphi>1) phi = pi * k/(hlen*pi/csamp);

/* Ensure that -45 <= xcphi <= 135 */
      xcphi = phi;
      if(xcphi < -pi/4)   xcphi = xcphi + pi;
      if(xcphi > 0.75*pi) xcphi = xcphi - pi;

/* Make a copy of red and diff */
      memcpy(B, R, nxbin*nybin*sizeof(float));
      memcpy(T, D, nxbin*nybin*sizeof(float));

/* Convolve the red with tophat then diff image with tophat,cos,sin */
      if(ABS(xcphi) < pi/4) {
	 xdrift_t(xcphi, hlen, nxbin, nybin, B);
	 xdrift_tcs(xcphi, hlen, nlam, nxbin, nybin, T, C, S);
      } else {
	 ydrift_t(xcphi, hlen, nxbin, nybin, B);
	 ydrift_tcs(xcphi, hlen, nlam, nxbin, nybin, T, C, S);
      }
      sprintf(line, "Drifted at phi %.1f (%.1f)", phi/dr, xcphi/dr);
      BLAB(line);
      if(telapse>timeout) {
	/*Out of time already. Print warning
	  message and then write empty streak file*/
	fprintf(stderr,"strker timeout (%lf>%lf) with phi %.1f\n",telapse,timeout, phi/dr);
	fprintf(stderr,"writing output file with no data\n");
	writeemptyfile(strkfile);
	return(1);
      }

/* Make the desired combination: [T-sqrt(C*C+S*S)]/RMS -> S = 'sigma' */
      rtn = 1 / (eadu*bin*bin*pi*sig);
      for(i=0; i<nxbin*nybin; i++) {
	 scsum = S[i]*S[i] + C[i]*C[i];
	 if(scsum > 0) scsum = sqrt(scsum);
	 rms = rtn * B[i];
	 if(rms > 0) rms = sqrt(rms);
	 else        rms = 1;
	 S[i] = (T[i] - scsum) / rms;
	 
/* Extra squash of things above background in B[i] */
         xtra = G[i] * (2*hlen+1) / B[i];
         if(xtra > 1) xtra = 1;
	 unsquashed = S[i]; /*preserves pre-squash value*/
         for(j=0; j<nsquash; j++) S[i] *= xtra; /*perform squashing*/
	 /*reverse squashing in the case of very bright sources*/
	 if(unsquashed > unsquashable) S[i]=unsquashed;
      }
      BLAB("Sigma calculated");

      if(telapse>timeout) {
	/*Out of time already. Print warning
	  message and then write empty streak file*/
	fprintf(stderr,"strker timeout (%lf>%lf) after sigma calculation\n",telapse,timeout);
	fprintf(stderr,"writing output file with no data\n");
	writeemptyfile(strkfile);
	return(1);
      }

/* Detect significant peaks */
      xcdet[k].hlen = hlen;
      xcdet[k].phicorr = phi;
      if(do_oldsigma_srch) {
	xcdet[k].ndet = sigma_srch(nxbin, nybin, S, DIFF, bin, srchsig, hlen, fwhm,
				   phi, maxtraillen, fixpsf, checkwidth, checklength, &xcdet[k].wpar, &xcdet[k].wflux);
      } else {
	xcdet[k].ndet = sigma_srch_ari(nxbin, nybin, S, DIFF, bin, srchsig, hlen, fwhm,
				       phi, maxtraillen, fixpsf, checkwidth, checklength, &xcdet[k].wpar, &xcdet[k].wflux);
      }

      ntph += xcdet[k].ndet;

      printf("Found %d detections at phi %.1f total %d\n", xcdet[k].ndet, phi/dr, ntph);
      fflush(stdout);

   }

   if(ntph>MAXNTPH)
     {
       printf("Warning: maximum number %d of detections exceeded\n",MAXNTPH);
       ntph = MAXNTPH;
       printf("Only the first %d detections will be considered\n",ntph);
     }

/* Consolidate all the tphot entries into whiskers, consolidating dups */
   if(DEBUGWHISK>0) printf("about to run whisk_cull\n");
   fflush(stdout);
   nwhisk = whisk_cull(bin, nphi, xcdet, &whisk);

   printf("%d sigma detections reduced to %d whiskers\n", ntph, nwhisk);
   fflush(stdout);
   if(nwhisk>=MAXNWHISK)
     {
       nwhisk = MAXNWHISK;
       printf("Since this is more than the maximum allowed\n");
       printf("by whiskers(), only %d will be analyzed\n",nwhisk);
     }

/* Group the whiskers into streaks to be subjected to a trail tphot fit */
   if(dowhisk) {
     // dmax = 1.0; dmax is now set by the user.
      tmax = 0.1;
   } else {
      dmax = tmax = 0.0;
   }
   nstrk0 = whiskers(nwhisk, whisk, bin*nxbin, bin*nybin,
		    dmax, tmax, &streak, csamp*bin);
      
//   if(VERBOSE > 1) printf("%d groups found\n", nstrk0);

   
   sprintf(line, "%d whiskers reduced to %d grouped streaks", nwhisk, nstrk0);
   BLAB(line);

   if(telapse>timeout) {
     /*Out of time already. Print warning
       message and then write empty streak file*/
     fprintf(stderr,"strker timeout (%lf>%lf) after whiskers run",telapse,timeout);
     fprintf(stderr,"writing output file with no data\n");
     writeemptyfile(strkfile);
     return(1);
   }

/* Create a weight file in IMG */
// FIXME: good idea to wipe out environs of streak with wgt=0?
   for(i=0; i<nx*ny; i++) {
      IMG[i] = (IMG[i] > 0) ? 1/(IMG[i]+rdnoise*rdnoise/eadu/eadu) : 0;
   }

   /*Sort streak array output by whiskers according to decreasing
     flux, so we can analyze the brightest ones first.*/
   /*allocate vectors*/
   sortflux = (double *)calloc(nstrk0, sizeof(double));
   sortind = (int *)calloc(nstrk0, sizeof(int));
   /*Load vectors*/
   for(k=0; k<nstrk0; k++) {
     sortflux[k] = -streak[k].flux;
     sortind[k] = k;
   }
   /*Sort*/
   tsort_d(nstrk0,sortflux,sortind);

/* Fit each streak using initial conditions */
   if(nstrk0>maxdets) nstrk0=maxdets; /*Impose limit on the number of detections to be analyzed.*/
   nstrk = streak_fit(nx, ny, DIFF, IMG, badata_dif,
		      bin, len, phi, fwhm, fixpsf, nstrk0, streak,
		      sortind, starttime, timeout, &wpar, &wflux, NULL);

   sprintf(line, "%d streaks produce %d distinct trail fits", nstrk0, nstrk);
   BLAB(line);


/* Classify each streak */
   streak_class(nx, ny, x0sub, y0sub, magzpt, fwhm, saturate,
		nstar, star, nstrk, streak, wpar);

   BLAB("streaks classified");

   

//////////////////////////////////////////////////
// Write the tphot detections... ordered by phi //
//////////////////////////////////////////////////
   if(tphfile != NULL) {
      if(strcmp(tphfile, "-") == 0) fp = stdout;
      else if( (fp=fopen(tphfile, "w")) == NULL) {
	 fprintf(stderr, "Cannot open tphfile `%s'\n", tphfile);
	 exit(1);
      }
      for(k=0; k<nphi; k++) {
	 for(i=0; i<xcdet[k].ndet; i++) {

/* Correct coords for (possible) subarray offset */
	    xcdet[k].wpar[i].x0 += x0sub;
	    xcdet[k].wpar[i].y0 += y0sub;
/* Tell us about it */
	    one_report(nxbin, nxbin, nybin, S,
		       (int)xcdet[k].wpar[i].x0+x0sub,
		       (int)xcdet[k].wpar[i].y0+y0sub,
		       0, xcdet[k].wpar+i, xcdet[k].wflux+i, fp, NULL);
	 }
      }
      if(strcmp(tphfile, "-") != 0) fclose(fp);
   }
//////////////////////////////////////////////////




/////////////////////////////////////////////////
// Write the group properties: maj,min,len,phi //
/////////////////////////////////////////////////
   if(grpfile != NULL) {
      if(strcmp(grpfile, "-") == 0) fp = stdout;
      else if( (fp=fopen(grpfile, "w")) == NULL) {
	 fprintf(stderr, "Cannot open grpfile `%s'\n", grpfile);
	 exit(1);
      }
      fprintf(fp, "# Grp Ndet x0     y0    len   phi\n");
      for(k=0; k<nstrk; k++) {
	 fprintf(fp, "%3d %3d %6.0f %6.0f %5.0f %5.0f\n",
		k, streak[k].n,
		streak[k].x0+x0sub, streak[k].y0+y0sub,
		 streak[k].len, streak[k].phi/dr);
      }
      if(strcmp(grpfile, "-") != 0) fclose(fp);
   }
/////////////////////////////////////////////////
   


/////////////////////////////////
// Write a whisker detail file //
/////////////////////////////////
   if(whiskfile != NULL) {
      if(strcmp(whiskfile, "-") == 0) fp = stdout;
      else if( (fp=fopen(whiskfile, "w")) == NULL) {
	 fprintf(stderr, "Cannot open whisker file `%s'\n", whiskfile);
	 exit(1);
      }
      fprintf(fp, "# Grp Ndet phicor    x0      y0      len      phi     th1      lam1     th2     lam2\n");

      for(k=0; k<nstrk; k++) {
	 for(i=0; i<streak[k].n; i++) {
	    fprintf(fp, "%3d %3d %7.1f  %7.1f %7.1f  %7.1f %7.1f  %7.1f %7.0f  %7.1f %7.0f\n",
		   k, i, streak[k].w[i]->phicorr/dr,
		   streak[k].w[i]->x0+x0sub, streak[k].w[i]->y0+y0sub,
		   streak[k].w[i]->len, streak[k].w[i]->phi/dr,
		   streak[k].w[i]->th1/dr, streak[k].w[i]->lam1,
		   streak[k].w[i]->th2/dr, streak[k].w[i]->lam2);

	 }
      }
      if(strcmp(whiskfile, "-") != 0) fclose(fp);
   }
/////////////////////////////////


   

/////////////////////////////////
// Write the streak fit output //
/////////////////////////////////
   if(strkfile != NULL) {
      if(strcmp(strkfile, "-") == 0) fp = stdout;
      else if( (fp=fopen(strkfile, "w")) == NULL) {
	 fprintf(stderr, "Cannot open streak file `%s'\n", strkfile);
	 exit(1);
      }
   
/* Put out these in order of increasing y */
      idx = (int *)calloc(nstrk, sizeof(int));
      buf = (double *)calloc(nstrk, sizeof(double));
      for(k=i=0; k<nstrk; k++) {
/* Rejection criteria */
#if 0
	 if(wpar[k].x0<0 || wpar[k].x0>nx ||
	    wpar[k].y0<0 || wpar[k].y0>ny || 
	    wpar[k].major < majmin || wpar[k].major > majmax ||
	    wpar[k].minor < minmin || wpar[k].minor > minmax) continue;
#else
	 if(wpar[k].x0<0 || wpar[k].x0>nx ||
	    wpar[k].y0<0 || wpar[k].y0>ny || 
	    wpar[k].major < majmin || wpar[k].major > majmax ||
	    wpar[k].minor < minmin || wpar[k].minor > minmax ||
	    wpar[k].chin > chinmax) continue;
#endif
	 buf[i] = wpar[k].y0;
	 idx[i] = k;
	 i++;
      }
      j = qsort2(i, buf, idx);
      
      fprintf(fp, "#  x0       y0       peak   dpeak  sky    flux    dflux   length  minor    phi  err nit ndet chin Ptr Pvr Pcr Pbn Psc Pmv\n");

      for(j=0; j<i; j++) {
	 k = idx[j];
// Suppress output if there is a fit error?  Probably not a good idea.
//	 if(wpar[k].bayer == 0 || VERBOSE > 0) {
	 /*Convert phi convention from -90..90 to 0..180*/
	 if(wpar[k].phi<0.0) wpar[k].phi += pi;
	 /*Write output file*/
	    fprintf(fp, "%8.2f %8.2f  %7.1f %6.1f %5.1f  %8.0f %7.0f  %6.2f %6.2f %6.2f %3x %2d %2d %8.1f %3d %3d %3d %3d %3d %3d\n",
		    wpar[k].x0+x0sub, wpar[k].y0+y0sub,
		    wpar[k].peak, wpar[k].dpeak, wpar[k].sky,
		    wflux[k].flux, wflux[k].dflux, 
		    wpar[k].major, wpar[k].minor, wpar[k].phi*180/pi,
		    wpar[k].bayer, wpar[k].niter, streak[k].ndet, wpar[k].chin,
		    streak[k].P[Ptr], streak[k].P[Pvr], streak[k].P[Pcr],
		    streak[k].P[Pbn], streak[k].P[Psc], streak[k].P[Pmv]);
//	 }
      }
      if(strcmp(strkfile, "-") != 0) fclose(fp);
      free(idx);
      free(buf);
   }
//////////////////////////////////




////////////////////////////////////
// Write a diagnostic image file? //
////////////////////////////////////
   if(newfile != NULL) {

/* We are *creating* the new file and overwriting if it exists */
/* In cfitsio parlance this means that it must start with an ! */
      newbang = malloc(strlen(newfile) + 2);
      if(newfile[0] != '!') {
	 strcpy(newbang+1, newfile);
	 newbang[0] = '!';
      } else {
	 strcpy(newbang, newfile);
      }

/* Initialize the new, output file */
      status = 0;
      if( fits_create_file(&fnew, newbang, &status) ) FRE;

      BLAB("Output file created");

/* Initialize the new, output image */
      bitpix = FLOAT_IMG;		/* signed shorts */
      naxis = 2;
      naxes[0] = nxbin;
      naxes[1] = nybin;
      if( fits_create_img(fnew, bitpix, naxis, naxes, &status) ) FFATAL;
      BLAB("Output image created");

      fpixel[0] = fpixel[1] = 1;                /* first pixel to write      */
      nelements = naxes[0] * naxes[1];          /* number of pixels to write */
      
/* Diagnostic: write an array to the FITS file */
      if(dcs == 0) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, S, &status) ) FRE;
      } else if(dcs == 1) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, B, &status) ) FRE;
      } else if(dcs == 2) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, T, &status) ) FRE;
      } else if(dcs == 3) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, C, &status) ) FRE;
      } else if(dcs == 4) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, R, &status) ) FRE;
      } else if(dcs == 5) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, D, &status) ) FRE;
      } else if(dcs == 6) {
	 if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, G, &status) ) FRE;
      }

      BLAB("Pixels written");

/* Finally close the output image */
      if( fits_close_file(fnew, &status) ) FRE;
   }
////////////////////////////////////////////////////////////////

   BLAB("Done");

   free(D);
   free(C);
   free(B);
   free(T);
   free(R);
   free(S);
   free(IMG);
   free(DIFF);

   return(0);
}


/* Local maximum */
void sigma_tune_localmax(int nx, int ny, float *data, int i0, int j0, int dr,
			 int *imax, int *jmax, double *max, double *min)
{
   int i, j;
   *imax = i0;
   *jmax = j0;
   *max = data[i0+j0*nx];
   *min = data[i0+j0*nx];
   for(j=-dr; j<=dr; j++) {
      for(i=-dr; i<=dr; i++) {
	 if(data[i0+i+(j0+j)*nx] > *max) {
	    *imax = i + i0;
	    *jmax = j + j0;
	    *max = data[i0+i+(j0+j)*nx];
	 }		  
	 if(data[i0+i+(j0+j)*nx] < *min) {
	    *min = data[i0+i+(j0+j)*nx];
	 }		  
      }
   }
   return;
}

/* Tune up the center, length, and phi of a sigma_srch() detection */
void sigma_tune(int nx, int ny, float *data, double sigma,
		PSF2D_PARAM *wpar, double *sigmax)
{
   int r, dr, i, j, i0, i1, j0, j1, ibig, jbig, done0, done1;
   double phi, big0, big1, min, prev0, prev1;

   *sigmax = 0;
   i0 = i1 = wpar->x0;
   j0 = j1 = wpar->y0;
   prev0 = prev1 = 0;
   dr = 2;

   phi = wpar->phi;

//   printf("Tune:  %8.1f %8.1f  %8.1f %8.1f %8.1f\n",
//	  wpar->x0, wpar->y0,
//	  wpar->major, wpar->minor, wpar->phi*45/atan(1));

   done0 = done1 = 0;

/* Move in each direction from the center point looking for a loss of signal */
   for(r=dr; r<MAX(nx,ny); r+=dr) {

/* Move counter to phi */
      if(!done0) {
	 i = wpar->x0 - r * cos(phi);
	 j = wpar->y0 - r * sin(phi);
	 if(i>=dr && i<nx-dr && j>=dr && j<ny-dr) {
	    sigma_tune_localmax(nx, ny, data, i, j, dr,
				&ibig, &jbig, &big0, &min);
	 } else {
	    done0 = 1;
	 }
      }

/* Advance i0,j0 if signal detected */
      if(!done0 && big0 >= sigma) {
	 i0 = ibig;
	 j0 = jbig;
	 prev0 = big0;
	 *sigmax = MAX(*sigmax, big0);
//	 printf("%d %d %7.1f\n", i0, j0, big0);
      }

/* Stop i0,j0 advance if low signal detected (but not masked area) */
      if(big0 < sigma && min > -sigma) done0 = 1;

/* Move along phi */
      if(!done1) {
	 i = wpar->x0 + r * cos(phi);
	 j = wpar->y0 + r * sin(phi);
	 if(i>=dr && i<nx-dr && j>=dr && j<ny-dr) {
	    sigma_tune_localmax(nx, ny, data, i, j, dr,
				&ibig, &jbig, &big1, &min);
	 } else {
	    done1 = 1;
	 }
      }
/* Advance i1,j1 if signal detected */
      if(!done1 && big1 >= sigma) {
	 i1 = ibig;
	 j1 = jbig;
	 prev1 = big1;
	 *sigmax = MAX(*sigmax, big1);
//	 printf("%d %d %7.1f %7.1f\n", i1, j1, big1, phi*45/atan(1));
      }

/* Stop i1,j1 advance if low signal detected (but not masked area) */
      if(big1 < sigma && min > -sigma) done1 = 1;

//      printf("%4d %4d %4d %8.1f %8.1f %d  %4d %4d %8.1f %8.1f %d\n", r, i0, j0, big0, min, done0, i1, j1, big1, min, done1);

/* Ends reached? */
      if(done0 && done1) break;

/* Adjust phi */
      phi = atan2(j1-j0, i1-i0);
   }

/* Update center, length, phi */
   wpar->x0 = 0.5 * (i0 + i1);
   wpar->y0 = 0.5 * (j0 + j1);
   wpar->major = sqrt((i1-i0)*(i1-i0)+(j1-j0)*(j1-j0));
   wpar->phi = phi;

   if(VERBOSE > 2) {
      printf("Tuned: %7.1f %7.1f  %6.1f %4.1f %6.1f  %4d %4d %5.1f  %4d %4d %5.1f  %5.1f\n",
	     wpar->x0, wpar->y0,
	     wpar->major, wpar->minor, wpar->phi*45/atan(1),
	     i0,j0, prev0, i1,j1, prev1, *sigmax);
   }

   return;
}



#define NPSF_ALLOC 8192		/* Number to allocate in a hunk */

extern int PSF2D_DEBUG;
extern int TPFN_TEST;
extern int PSF2D_TRAILITER;
extern int nwpar;

/* Look for peaks in sigma image */
int sigma_srch(int nx, int ny, float *data, float *diff, int bin,
		   double sigma, int hlen, double fwhm,double phi/*rad*/,
                   double maxtraillen, int fixpsf, int checkwidth, int checklength, 
                   PSF2D_PARAM **wpar, PSF2D_FLUX **wflux)
{
   int i, j, NX=nx, npix, nobj, nfit, err, nalloc, edge;
   int invert, trail, force, okfit;
   double sigmax;
   double len = 2*hlen+1;       // full length on binned image
   double srchmin = 5.0;	// Minimum value above 0 to trigger
   double srchmax = 1e5;	// Maximum peak to trigger
   double srchrad = 5.0;	// Local max required over radius r to trigger
   double srchsig = 5.0;	// Minimum value = sky+sig*rms to trigger

   double maxmove = 10;		// Fit ctr must coincide with peak this closely
   double fwmax = 200.0;	// Maximum FW to accept
   double fwmin = 1.0;		// Minimum FW to accept
//   double hgtmin = sigma;	// Minimum peak-sky to keep
// but a trigger at sigma may fit to less than sigma...
   double hgtmin = sigma/2;	// Minimum peak-sky to keep
   double chinmax = 1000.0;	// Maximum chi/N to keep
   double sig = -1.0;		// Minimum SNR in the flux determination
//   double minormax = 6.0;	// Maximum minor axis
// truly the minor should not exceed the streak length, but it might
   double minormax = 2*len;	// Maximum minor axis   	

// JTJT
   float *wgt=(float *)calloc(nx*ny, sizeof(float));
   for(i=0; i<nx*ny; i++) {
      wgt[i] = (data[i] > -sigma) ? 1 : 0;
   }

/* Set search parameters */
   invert = 0;		// Invert the sign of the image?
   trail = 0;		// NOT a trailed fit.
//   okfit = 1;		// Does the fit have to be OK to keep?
   okfit = 0;		// Does the fit have to be OK to keep?
   force = 0;		// Disregard all errors that inhibit output

/* Tell tphot about trigger and acceptance parameters */
   srchrad = MAX(srchrad, 0.5*len);	// Local max required over radius r to trigger
   srchsig = sigma;	// Minimum value = sky+sig*rms to trigger
   tpctrl_trigger(srchmin, srchmax, srchrad, srchsig);

/* Set the borders to N times the length */
   edge = 3 * len;
   tpctrl_search(edge, edge, edge, edge, 0.0);

   maxmove = MAX(maxmove, 0.5*len);	// Fit ctr must coincide with peak this closely
   tpctrl_accept(maxmove, fwmax, fwmin, hgtmin, chinmax, sig);

/* Use just a 4 parameter fit */
   if(fixpsf) tpctrl_fit(4/*npar*/, 15/*ap*/, 40/*sky*/);

/* Allocate some space for the results */
   *wpar = (PSF2D_PARAM *)calloc(NPSF_ALLOC, sizeof(PSF2D_PARAM));
   *wflux = (PSF2D_FLUX *)calloc(NPSF_ALLOC, sizeof(PSF2D_FLUX));
   nalloc = NPSF_ALLOC;

/* Find the objects: npix is the starting address in the image */
   nfit = 0;		// How many triggers were there?
   nobj = 0;		// How many fits passed acceptance?

//   printf("nx,y= %d %d sigma= %.2f len= %.2f phi=%.5f fixpsf= %d  data= %.3f\n", nx, ny, sigma, len, phi, fixpsf, data[3460+3624*nx]);

   npix = 0;		// Initial pixel to search
//   int iii=0;	// JT
   while(1) {
//      iii++;	// JT
//      if( (iii%1) == 0) {
//	 printf( "%6d %6d %6d\n", iii, i, j);
//	 fflush(stdout);
//      }

/* Find a pixel at address npix->i,j that passes the trigger criteria */
      one_trigger(nx, NX, ny, data, &npix, ny, NULL, NULL);
      if(npix >= nx*ny) break;
      i = npix % nx;
      j = npix / nx;
      nfit++;

///////
#if 0
      int inside=ABS(i-3460)<5 && ABS(j-3624)<5;
      if(inside) printf("Trigger %d %d %.3f\n", i, j, data[i+j*nx]);
      PSF2D_DEBUG = TPFN_TEST = 2 * inside;
#endif
///////

/* Ensure we still have enough space */
      if(nobj >= nalloc) {
	 *wpar = (PSF2D_PARAM *)realloc(*wpar,
				(nalloc+NPSF_ALLOC)*sizeof(PSF2D_PARAM));
	 *wflux = (PSF2D_FLUX *)realloc(*wflux,
				(nalloc+NPSF_ALLOC)*sizeof(PSF2D_FLUX));
	 nalloc += NPSF_ALLOC;
      }
#ifdef TPHOT_DEBUG
      printf("New det at x= %5d y= %5d  peak= %8.1f\n", i, j, data[i+j*NX]);
#endif

//      fprintf(stderr, "nalloc %d nobj %d at %6d %6d  0x%08x %d %d\n", nalloc, nobj, i, j, *wpar, ((long int *)(*wpar))[-1], ((long int *)(*wpar))[-1]/sizeof(PSF2D_PARAM));


/* OK, got one... analyze it */
      (*wpar)[nobj].x0 = i;
      (*wpar)[nobj].y0 = j;
      if(fixpsf) {
//	 (*wpar)[nobj].major = len/3;	/* Use length/3 for sigma image */
	 (*wpar)[nobj].major = len;	/* Use full length for sigma image */
	 (*wpar)[nobj].minor = 2;	/* Don't cheap out on minor axis! */
	 (*wpar)[nobj].phi = phi;
	 (*wpar)[nobj].ninit = 5;
      } else {
	 (*wpar)[nobj].major = -1;
	 (*wpar)[nobj].minor = -1;
	 (*wpar)[nobj].phi = -100;
	 (*wpar)[nobj].ninit = 0;
      }

//      PSF2D_DEBUG = TPFN_TEST = i==744 && j==4560;
//      PSF2D_DEBUG = TPFN_TEST = 2 * (i==4105 && j==1684);
//      PSF2D_DEBUG = TPFN_TEST = 2 * (i==149 && j==150);

// JTJT
      err = one_object(nx, NX, ny, data, NULL, i, j,
		       invert, trail, force, okfit,
		       (*wpar)+nobj, (*wflux)+nobj);

///////
#if 0
/* Iterate again if we've hit the iteration maximum */
      if((*wpar)[nobj].niter == PSF2D_TRAILITER) {
         printf("Iteration= %d try again\n", (*wpar)[nobj].niter);
	 (*wpar)[nobj].ninit = 5;
         err = one_object(nx, NX, ny, data, NULL, i, j,
		       invert, trail, force, okfit,
		       (*wpar)+nobj, (*wflux)+nobj);
	 (*wpar)[nobj].niter += PSF2D_TRAILITER;
      }
#endif
///////


      if(PSF2D_DEBUG) {
	 printf("  err %d  i,j %d %d minor %6.2f  chin %6.1f\n",
		err, i, j, (*wpar)[nobj].minor, (*wpar)[nobj].chin);

	 one_report(nx, NX, ny, data, i, j, trail,
		    (*wpar)+nobj, (*wflux)+nobj, stdout, NULL);
      }
      
/* Keep it? */
      if(err == 0 && (*wpar)[nobj].minor <= minormax) {
#ifdef TPHOT_DEBUG
	 one_report(nx, NX, ny, data, i, j, trail,
		    (*wpar)+nobj, (*wflux)+nobj, stdout, NULL);
#endif

//	 one_report(nx, NX, ny, data, i, j, trail, (*wpar)+nobj, (*wflux)+nobj, stdout, NULL);

/* Tune up the streak parameters if they are not fixed */
	 if(!fixpsf) {
	    sigma_tune(nx, ny, data, MAX(3.0,0.5*sigma), *wpar+nobj, &sigmax);
//	    sigma_tune(nx, ny, data, sigma, *wpar+nobj, &sigmax);
/* Trail in sigma image will be a trapezoid to ~0.5 sigma *BUT* cos/sin also */
//	    if((*wpar)[nobj].major > len) (*wpar)[nobj].major -= len;
//	    if((*wpar)[nobj].major > 2*len)
//	       (*wpar)[nobj].major -= 2*len*(1-0.5*sigma/sigmax);
	 }

	 nobj++;
      }
   }
// JTJT
//   printf("ntrig = %d\n", nfit);
//   tpctrl_pareport();
   
/* Correct the parameters according to bin factor */	 
   for(i=0; i<nobj; i++) {
     (*wpar)[i].x0 *= bin;
     (*wpar)[i].y0 *= bin;
     (*wpar)[i].major *= bin;
     (*wpar)[i].minor *= bin;
     (*wpar)[i].sky *= bin*bin;
     (*wpar)[i].peak *= bin*bin;
     (*wpar)[i].dpeak *= bin*bin;
     (*wflux)[i].sky *= bin*bin;
     (*wflux)[i].flux *= bin*bin;
     (*wflux)[i].dflux *= bin*bin;
   }

   free(wgt);
   
   return(nobj);
}

#define LONGBRIGHTMAXNUM 1000 /*Max number of trails in internal Long Bright catalog*/
#define LONGBRIGHTLEN 100.0 /*Minimum length for inclusion in Long Bright catalog*/
#define LONGBRIGHTSNR 20.0  /*Minimum SNR for inclusion in Long Bright catalog*/

/* Look for peaks in sigma image. This is
   Ari's new version of sigma_srch, which runs some
   diagnostics on the unbinned difference image as well*/
int sigma_srch_ari(int nx, int ny, float *data, float *DIFF, int bin,
		   double sigma, int hlen, double fwhm,double phi/*rad*/,
                   double maxtraillen, int fixpsf, int checkwidth, int checklength, 
                   PSF2D_PARAM **wpar, PSF2D_FLUX **wflux)
{
  int i, j, NX=nx, npix, nobj, nfit, nalloc, edge,firstct,maxct;

   double srchmin = 5.0;	// Minimum value above 0 to trigger
   double srchmax = 1e5;	// Maximum peak to trigger
   double srchrad = 1.5;	// Local max required over radius r to trigger
   double srchsig = 5.0;	// Minimum value = sky+sig*rms to trigger

   double maxmove = 10;		// Fit ctr must coincide with peak this closely
   double fwmax = 200.0;	// Maximum FW to accept
   double fwmin = 1.0;		// Minimum FW to accept
//   double hgtmin = sigma;	// Minimum peak-sky to keep
// but a trigger at sigma may fit to less than sigma...
   double hgtmin = sigma/2;	// Minimum peak-sky to keep
   double chinmax = 1000.0;	// Maximum chi/N to keep
   double sig = -1.0;		// Minimum SNR in the flux determination
//   double minormax = 6.0;	// Maximum minor axis

   double fluxsig,pixarea,len,tf1,tf2,skylev,skynoise,phirange,traillen;
   double tpx,tpy,phitemp,flux,tw;
   int trailendcase,trailstartcase,endfound,startfound,trailwidthOK=1;
   double endx,endy,startx,starty,te_tpx,te_tpy,te_phitemp,te_traillen;
   int temphlen,halfsearchlen,longbrightnum=0,longbrightct,longbrightmatch=0;
   double **longbrightmat,alongdist,acrossdist,dx,dy,cv,sv;
   int maxhalflen;

   maxhalflen = ceil(maxtraillen/2.0);

   trailendcase=trailstartcase=endfound=startfound=0;

   len = 2*hlen+1;
   phirange = 0.3;
   traillen = 2*hlen*bin+1;

/* Tell tphot about trigger and acceptance parameters */
   // srchrad = MAX(srchrad, 0.5*len);	// Local max required over radius r to trigger
   srchsig = sigma;	// Minimum value = sky+sig*rms to trigger
   tpctrl_trigger(srchmin, srchmax, srchrad, srchsig);

/* Set the borders to N times the length */
   edge = 3 * len;
   tpctrl_search(edge, edge, edge, edge, 0.0);

   maxmove = MAX(maxmove, 0.5*len);	// Fit ctr must coincide with peak this closely
   tpctrl_accept(maxmove, fwmax, fwmin, hgtmin, chinmax, sig);

/* Use just a 4 parameter fit */
   if(fixpsf) tpctrl_fit(4/*npar*/, 15/*ap*/, 40/*sky*/);

/* Allocate some space for the results */
   *wpar = (PSF2D_PARAM *)calloc(NPSF_ALLOC, sizeof(PSF2D_PARAM));
   *wflux = (PSF2D_FLUX *)calloc(NPSF_ALLOC, sizeof(PSF2D_FLUX));
   nalloc = NPSF_ALLOC;
   longbrightmat = dmat(LONGBRIGHTMAXNUM,5); /*txp, tpy, phi(rad), len, SNR*/

/* Find the objects: npix is the starting address in the image */
   nfit = 0;		// How many triggers were there?
   nobj = 0;		// How many fits passed acceptance?
   firstct = maxct = 0;
   longbrightnum = 0;

//   printf("nx,y= %d %d sigma= %.2f len= %.2f phi=%.5f fixpsf= %d  data= %.3f\n", nx, ny, sigma, len, phi, fixpsf, data[3460+3624*nx]);

   npix = 0;		// Initial pixel to search
//   int iii=0;	// JT
   while(1) {
//      iii++;	// JT
//      if( (iii%1) == 0) {
//	 printf( "%6d %6d %6d\n", iii, i, j);
//	 fflush(stdout);
//      }

/* Find a pixel at address npix->i,j that passes the trigger criteria */
      one_trigger(nx, NX, ny, data, &npix, ny, NULL, NULL);
      if(npix >= nx*ny) break;
      i = npix % nx;
      j = npix / nx;
      nfit++;

///////
#if 0
      int inside=ABS(i-3460)<5 && ABS(j-3624)<5;
      if(inside) printf("Trigger %d %d %.3f\n", i, j, data[i+j*nx]);
      PSF2D_DEBUG = TPFN_TEST = 2 * inside;
#endif
///////

/* Ensure we still have enough space */
      if(nobj >= nalloc) {
	 *wpar = (PSF2D_PARAM *)realloc(*wpar,
				(nalloc+NPSF_ALLOC)*sizeof(PSF2D_PARAM));
	 *wflux = (PSF2D_FLUX *)realloc(*wflux,
				(nalloc+NPSF_ALLOC)*sizeof(PSF2D_FLUX));
	 nalloc += NPSF_ALLOC;
      }
#ifdef TPHOT_DEBUG
      printf("New det at x= %5d y= %5d  peak= %8.1f\n", i, j, data[i+j*NX]);
#endif

//      fprintf(stderr, "nalloc %d nobj %d at %6d %6d  0x%08x %d %d\n", nalloc, nobj, i, j, *wpar, ((long int *)(*wpar))[-1], ((long int *)(*wpar))[-1]/sizeof(PSF2D_PARAM));


/* OK, got one... analyze it */
      /*Get location of local max in unbinned coordinates*/
      tpx = (double)(i*bin) + (double)bin/2.0;
      tpy = (double)(j*bin) + (double)bin/2.0;
      traillen = 2*hlen*bin+1;
      temphlen = hlen*bin;

      longbrightmatch=0;
      /*See if it matches something in the long bright trail matrix*/
      longbrightct=1;
      while(longbrightct<=longbrightnum && longbrightmatch==0)
	{
	  /*txp, tpy, phi(rad), len, SNR*/
	  dx = tpx-longbrightmat[longbrightct][1];
	  dy = tpy-longbrightmat[longbrightct][2];
	  if(dx*dx+dy*dy<DSQUARE(longbrightmat[longbrightct][4]/2.0))
	    {
	      sv = sin(longbrightmat[longbrightct][3]);
	      cv = cos(longbrightmat[longbrightct][3]);
	      /*dot product with along-trail unit vector*/
	      alongdist = cv*dx + sv*dy;
	      /*Subtract along-trail component from dx,dy to get cross trail distances*/
	      dx -= alongdist*cv;
	      dy -= alongdist*sv;
	      acrossdist = sqrt(dx*dx + dy*dy);
	      /*Account for possibility that alongdist is negative*/
	      alongdist = fabs(alongdist);
	      if(alongdist<=longbrightmat[longbrightct][4]/2.0 && acrossdist<=ceil(fwhm+1.0)) 
		{
		  longbrightmatch=1;
		  if(PREMATCH_DEBUG) printf("Matched source at %.2f,%.2f with bright long trail %d: %.2f,%.2f %.3f %.2f %.3f\n",tpx,tpy,longbrightct,longbrightmat[longbrightct][1],longbrightmat[longbrightct][2],longbrightmat[longbrightct][3]*180.0/PI,longbrightmat[longbrightct][4],longbrightmat[longbrightct][5]);
		}
	    }
	  longbrightct++;
	}

      if(PAFIND_DEBUG) printf("START START: %f, %f\n",tpx,tpy);

      if(ENDFIND_DEBUG) printf("nobj = %d, pos = %f,%f, flux=%f\n",nobj,tpx,tpy,data[i+j*NX]);

      maxct+=1;
      trailsquish_S(DIFF,bin,nx,ny,tpx,tpy,hlen,phi,fwhm,&tf1,&tf2,&skylev,&skynoise);
      /*tf1 and tf2 are the total trail fluxes based on trim-means
        over a collapsed, 1-D trail vector with 10% and 50% rejection, respectively.
        skynoise is trailsquish estimate of the true single-pixel
        sky rms near this trail on the original difference image.*/
      pixarea = (2*hlen*bin+1)*(2*ceil(fwhm/2.0)+1.0); /*pixel area over which trailsquish
                                                         calculated the flux*/
      fluxsig = skynoise*sqrt(pixarea); /*uncertainty on trail flux due to sky noise*/

      /* Keep it if the 10% trim  mean flux indicates a 
         greater than srchsig sigma detection above sky noise, 
         and the 50% trim mean flux is at least half as bright as the
         10% trim mean flux.*/
      if(ENDFIND_DEBUG) printf("tf1=%f, tf2=%f, fluxsig=%f: ratios %f %f %f\n",tf1,tf2,fluxsig,tf1/fluxsig,tf2/fluxsig,tf2/tf1);
      if(tf1/fluxsig>srchsig && tf2>=0.5*tf1 && !longbrightmatch)
	{
	  firstct++;
	  /*Refine phi and cross-trail centroid*/
	  /*obtain celestial position angle in degrees*/
	  phitemp=phi;
	  phirange = TRAILCENPAPIXRANGE/((double)(hlen*bin));
	  /*Two iterations refining PA and cross-trail centroid*/
	  if(ENDFIND_DEBUG) printf("PASS1: %f, %f, tlen=%f ,phi=%f, phirng=%f\n",tpx,tpy,traillen,phitemp,phirange);
	  phifind(DIFF,nx*bin,ny*bin,tpx,tpy,traillen,&phitemp,phirange,fwhm,TRAILCENGFAC,TRAILCENPAPIXSAMP,TRAILCENRLEV);
	  crosstrailcen01_S(DIFF,nx*bin,ny*bin,&tpx,&tpy,traillen,phitemp,fwhm,TRAILCENGFAC,TRAILCENITNUM,TRAILCENRLEV);
	  phifind(DIFF,nx*bin,ny*bin,tpx,tpy,traillen,&phitemp,phirange,fwhm,TRAILCENGFAC,TRAILCENPAPIXSAMP,TRAILCENRLEV);
	  crosstrailcen01_S(DIFF,nx*bin,ny*bin,&tpx,&tpy,traillen,phitemp,fwhm,TRAILCENGFAC,TRAILCENITNUM,TRAILCENRLEV);
	  while(phitemp>=PI) phitemp -= PI;
	  if(phitemp<0) phitemp += PI;

	  if(ENDFIND_DEBUG) printf("REFINE1: %f, %f, phitemp = %f\n",tpx,tpy,phitemp*180.0/PI);

	  /*See if we still pass the trailsquish cut*/
	  trailsquish_S(DIFF,bin,nx,ny,tpx,tpy,hlen,phitemp,fwhm,&tf1,&tf2,&skylev,&skynoise);
	  pixarea = (2*hlen*bin+1)*(2*ceil(fwhm/2.0)+1.0);
	  fluxsig = skynoise*sqrt(pixarea);
	  if(ENDFIND_DEBUG) printf("tf1=%f, tf2=%f, fluxsig=%f: ratios %f %f %f\n",tf1,tf2,fluxsig,tf1/fluxsig,tf2/fluxsig,tf2/tf1);
	  if(tf1/fluxsig>srchsig && tf2>=0.5*tf1)
	    {
	      if(ENDFIND_DEBUG) printf("PASS2: %f, %f\n",tpx,tpy);
	      /*Seems to be bright enough, and with a reasonably consistent
                flux. Check if the trail is too fuzzy.*/
	      trailwidthOK=1;
	      if(checkwidth) {
		tw=fwhm;
		trailwidth(DIFF,nx*bin,ny*bin,tpx,tpy,hlen*bin,phitemp,fwhm,&tw,&flux);
		// printf("trailwidth result: %.2f,%2.f phi %.2f width is %f sig %f\n",tpx,tpy,phitemp*180.0/PI,tw/fwhm,flux/fluxsig);
		trailwidthOK=1;
		/*adding 1 pixel in quadrature accounts for interpolation blurring*/
		if(ENDFIND_DEBUG) printf("tw=%f, comp=%f, flux=%f, fluxsig=%f, ratio=%f\n",tw,sqrt(DSQUARE(MAXTRAILWIDTH*fwhm)+1.0),flux,fluxsig,flux/fluxsig);
		if(tw>sqrt(DSQUARE(MAXTRAILWIDTH*fwhm)+1.0))
		  {
		    if(flux/fluxsig>3.0) trailwidthOK=0;
		    /*Trail gets classified as too wide only if trailwidth
		      reported a sufficiently significant detection that
		      we can trust the width measurement*/
		  }
	      }
	      if(trailwidthOK)
		{
		  temphlen = hlen*bin;
		  /*Refind trail length if requested*/
		  if(checklength) {
		    /*Try to refine the length of the trail*/
		    endfound=startfound=trailendcase=trailstartcase=0;
		    halfsearchlen = temphlen*3;
		    te_tpx = tpx;
		    te_tpy = tpy;
		    te_phitemp = phitemp;
		    te_traillen = traillen;
		    while((!endfound || !startfound) && temphlen<sqrt(nx*nx*bin*bin+ny*ny*bin*bin))
		      {
			temphlen = hlen*bin;
			/*Trail in forward direction*/
			if(!endfound)
			  {
			    trailendcase = trailendfind(DIFF,nx*bin,ny*bin,te_tpx,te_tpy,temphlen,halfsearchlen,maxhalflen,te_phitemp,fwhm,&endx,&endy);
			    if(ENDFIND_DEBUG) printf("endx=%f, endy=%f\n",endx,endy);
			    if(trailendcase==-1) break; /*Too faint to measure, break out of loop*/
			    else if(trailendcase==2) endfound=1; /*Trail faded away; no point searching further*/
			    else if(trailendcase==1) endfound=0; /*Trail seems to continue further in this direction*/
			    else if(trailendcase==0) endfound=1; /*Definitive end-of-trail found*/
			  }
			/*Trail in backward direction*/
			if(!startfound)
			  {
			    trailstartcase = trailendfind(DIFF,nx*bin,ny*bin,te_tpx,te_tpy,temphlen,halfsearchlen,maxhalflen,te_phitemp+PI,fwhm,&startx,&starty);
			    if(ENDFIND_DEBUG) printf("startx=%f, starty=%f\n",startx,starty);
			    if(trailstartcase==-1) break; /*Too faint to measure, break out of loop*/
			    else if(trailstartcase==2) startfound=1; /*Trail faded away; no point searching further*/
			    else if(trailstartcase==1) startfound=0; /*Trail seems to continue further in this direction*/
			    else if(trailstartcase==0) startfound=1; /*Definitive end-of-trail found*/
			  }
			if(trailendcase>=0 && trailstartcase>=0)
			  {
			    if((0.5*endx + 0.5*startx)>0.0 && (0.5*endx + 0.5*startx)<nx*bin && (0.5*endy + 0.5*starty)>0 && (0.5*endy + 0.5*starty)<ny*bin)
			      {
				te_tpx = 0.5*endx + 0.5*startx;
				te_tpy = 0.5*endy + 0.5*starty;
				te_traillen = sqrt(DSQUARE(endx-startx)+DSQUARE(endy-starty));
				halfsearchlen = te_traillen*1.5;
				if(halfsearchlen<temphlen*3) halfsearchlen=temphlen*3;
				/*Refine phi and cross-trail centroid*/
				phirange = TRAILCENPAPIXRANGE/(te_traillen*0.5);
				/*Two iterations refining phi and cross-trail centroid*/
				if(ENDFIND_DEBUG) printf("phi = %f, phirange=%f, traillen=%f, tp=%f,%f\n",te_phitemp*180.0/PI,phirange*180.0/PI,te_traillen,te_tpx,te_tpy);
				phifind(DIFF,nx*bin,ny*bin,te_tpx,te_tpy,te_traillen,&te_phitemp,phirange,fwhm,TRAILCENGFAC,TRAILCENPAPIXSAMP,TRAILCENRLEV);
				if(ENDFIND_DEBUG) printf("phi = %f, phirange=%f, traillen=%f, tp=%f,%f\n",te_phitemp*180.0/PI,phirange*180.0/PI,te_traillen,te_tpx,te_tpy);
				crosstrailcen01_S(DIFF,nx*bin,ny*bin,&te_tpx,&te_tpy,te_traillen,te_phitemp,fwhm,TRAILCENGFAC,TRAILCENITNUM,TRAILCENRLEV);
				if(ENDFIND_DEBUG) printf("phi = %f, phirange=%f, traillen=%f, tp=%f,%f\n",te_phitemp*180.0/PI,phirange*180.0/PI,te_traillen,te_tpx,te_tpy);
				phifind(DIFF,nx*bin,ny*bin,te_tpx,te_tpy,te_traillen,&te_phitemp,phirange,fwhm,TRAILCENGFAC,TRAILCENPAPIXSAMP,TRAILCENRLEV);
				if(ENDFIND_DEBUG) printf("phi = %f, phirange=%f, traillen=%f, tp=%f,%f\n",te_phitemp*180.0/PI,phirange*180.0/PI,te_traillen,te_tpx,te_tpy);
				crosstrailcen01_S(DIFF,nx*bin,ny*bin,&te_tpx,&te_tpy,te_traillen,te_phitemp,fwhm,TRAILCENGFAC,TRAILCENITNUM,TRAILCENRLEV);
				if(ENDFIND_DEBUG) printf("phi = %f, phirange=%f, traillen=%f, tp=%f,%f\n",te_phitemp*180.0/PI,phirange*180.0/PI,te_traillen,te_tpx,te_tpy);
				if(ENDFIND_DEBUG) printf("TEF %f,%f len=%f, phi=%f\n",te_tpx,te_tpy,te_traillen,te_phitemp*180.0/PI);
			      }
			    else 
			      {
				/*Calculated trail center is off the image. Print warning
                                  and then exit while loop, leaving trail parameters
                                  at their last good values.*/
				fprintf(stderr, "WARNING: endfind produce trail center off image\nstart %f,%f to end %f,%f. Center NOT updated\n",startx,starty,endx,endy);
				endfound = startfound = 1; /*This gets us out of the while loop*/
			      }
			  }
			if(ENDFIND_DEBUG) printf("Bottom of while loop with trailecase %d,%d, foundcase %d,%d, TEF %f,%f len=%f, phi=%f\n",trailstartcase,trailendcase,startfound,endfound,te_tpx,te_tpy,te_traillen,te_phitemp*180.0/PI);
		      }
		    endx = te_tpx + te_traillen/2.0*cos(te_phitemp);
		    startx = te_tpx - te_traillen/2.0*cos(te_phitemp);
		    endy = te_tpy + te_traillen/2.0*sin(te_phitemp);
		    starty = te_tpy - te_traillen/2.0*sin(te_phitemp);
		    if(ENDFIND_DEBUG) printf("Finished trail length analysis %f %f: start %f,%f, end %f,%f\n",te_tpx,te_tpy,startx,starty,endx,endy);
		    /*Did it work?*/
		    if(trailendcase>=0 && trailstartcase>=0 && te_traillen>traillen)
		      {
			/*trailendfind reported success in both directions, and
                          found a trail length longer than the previous nominal value.
                          Hence, we trust the outputs and will copy them to the
                          operational values that will be written to the output struct*/
			tpx = te_tpx;
			tpy = te_tpy;
			while(te_phitemp>=PI) te_phitemp -= PI;
			if(te_phitemp<0) te_phitemp += PI;
			phitemp = te_phitemp;
			traillen = te_traillen;
			/*Also, recalculate fluxes using trailsquish*/
			temphlen = ceil(traillen/2.0)/bin;
			trailsquish_S(DIFF,bin,nx,ny,tpx,tpy,temphlen,phitemp,fwhm,&tf1,&tf2,&skylev,&skynoise);
			pixarea = (2*temphlen*bin+1)*(2*ceil(fwhm/2.0)+1.0);
			fluxsig = skynoise*sqrt(pixarea);
		      }
		  }
		  if(ENDFIND_DEBUG) printf("phi=%f,phitemp=%f\n",phi*180.0/PI,phitemp*180.0/PI);
		  nobj++;
		  (*wpar)[nobj-1].phi = phitemp;
		  (*wpar)[nobj-1].x0 = tpx;
		  (*wpar)[nobj-1].y0 = tpy;
		  (*wpar)[nobj-1].major = traillen;
		  (*wpar)[nobj-1].minor = fwhm;
		  (*wpar)[nobj-1].sky = skylev;
		  (*wpar)[nobj-1].peak = tf1/(2.0*temphlen*bin+1.0);
		  (*wpar)[nobj-1].dpeak = fluxsig/(2.0*temphlen*bin+1.0);
		  (*wflux)[nobj-1].sky = skylev;
		  (*wflux)[nobj-1].flux = tf1;
		  (*wflux)[nobj-1].dflux = fluxsig;
		  /*Load into longbrightmat if appropriate*/
		  if(tf1/fluxsig > LONGBRIGHTSNR && traillen > LONGBRIGHTLEN && longbrightnum < LONGBRIGHTMAXNUM)
		    {
		      longbrightnum+=1;
		      /*txp, tpy, phi(rad), len, SNR*/
		      longbrightmat[longbrightnum][1] = tpx;
		      longbrightmat[longbrightnum][2] = tpy;
		      longbrightmat[longbrightnum][3] = phitemp;
		      longbrightmat[longbrightnum][4] = traillen;
		      longbrightmat[longbrightnum][5] = tf1/fluxsig;
		      if(PREMATCH_DEBUG) printf("Loaded long bright trail %d: %.2f,%.2f, phi=%.3f, len=%.2f, SNR=%.3f\n",longbrightnum,longbrightmat[longbrightnum][1],longbrightmat[longbrightnum][2],longbrightmat[longbrightnum][3]*180.0/PI,longbrightmat[longbrightnum][4],longbrightmat[longbrightnum][5]);
		    }
		}
	    }
	}
   }
   free_dmat(longbrightmat,LONGBRIGHTMAXNUM,5);
   return(nobj);
}

#undef LONGBRIGHTMAXNUM
#undef LONGBRIGHTLEN
#undef LONGBRIGHTSNR


/* Consolidate all the tphot entries into whiskers, consolidating dups */
int whisk_cull(int bin, int nphi, XCDET *xcdet, WHISKER **whisk)
{
   int i, j, k, l, n, ntph, maxdet, ngrp, nlink, *link;
   double dmax, dx, dy, *buf;
   PSF2D_PARAM **wpar;	/* pointers to all detection fits */
   int *idx;
   char *pairs, *used;

   if(DEBUGWHISK>0) printf("inside whisk_cull\n");
   fflush(stdout);

/* Friends of friends linkage distance between x,y */
   dmax = 3;

/* Total number of detections over all phi */
   for(k=ntph=maxdet=0; k<nphi; k++) {
      ntph += xcdet[k].ndet;
      maxdet = MAX(maxdet, xcdet[k].ndet);
      if(DEBUGWHISK>0) printf("k=%d, xcdet[k].ndet=%d, maxdet=%d, ntph=%d\n",k,xcdet[k].ndet,maxdet,ntph);
      fflush(stdout);
   }

   if(ntph>MAXNTPH)
     {
       printf("Warning: actual detection count %d found by whisk_cull\n",ntph);
       printf("exceeds calloc-imposed maximum of %d for the pair array\n",MAXNTPH);
       ntph = MAXNTPH;
       printf("Only the first %d detections will be considered\n",ntph);
     }

/* Filled array of tph pointers and an index to get back to XCDET array */
   idx = (int *)calloc(ntph, sizeof(int));
   if(DEBUGWHISK>0) printf("Allocated idx OK\n");
   fflush(stdout);
   wpar = (PSF2D_PARAM **)calloc(ntph, sizeof(PSF2D_PARAM *));
   if(DEBUGWHISK>0) printf("Allocated wpar OK\n");
   fflush(stdout);
   for(k=n=0; k<nphi; k++) {
      for(i=0; i<xcdet[k].ndet; i++) {
	if(n<ntph)
	  {
	    idx[n] = i + k*maxdet;
	    wpar[n] = xcdet[k].wpar + i;
	    if(DEBUGWHISK>0) printf("k=%d, i=%d, idx=%d, n=%d\n",k,i,idx[n],n);
	    fflush(stdout);
	    n++;
	  }
      }
   }

/* Allocate whisker array that is definitely big enough */
   *whisk = (WHISKER *)calloc(ntph, sizeof(WHISKER));

// ORIGINAL CODE: ACCEPT ALL STREAKS
#if 0
/* Assemble all the detections: N(phi) for nphi */
   j = 0;
   for(k=0; k<nphi; k++) {
      for(i=0; i<xcdet[k].ndet; i++) {
/* Create whisker information for all detections */
	 (*whisk)[j].wpar = xcdet[k].wpar+i;
	 (*whisk)[j].wflux = xcdet[k].wflux+i;
	 (*whisk)[j].lencorr = bin*(2*xcdet[k].hlen+1);
	 (*whisk)[j].phicorr = xcdet[k].phicorr;
	 (*whisk)[j].x0 = xcdet[k].wpar[i].x0;
	 (*whisk)[j].y0 = xcdet[k].wpar[i].y0;
	 (*whisk)[j].len = xcdet[k].wpar[i].major;
	 (*whisk)[j].phi = xcdet[k].wpar[i].phi;
	 j++;
      }
   }
   return(ntph);
#endif
   if(DEBUGWHISK>0) printf("About to allocate pair array");
   fflush(stdout);

   pairs = (char *)calloc(ntph*ntph, sizeof(char));
   if(DEBUGWHISK>0) printf("Allocated pair array");
   fflush(stdout);
   used = (char *)calloc(ntph, sizeof(char));
   if(DEBUGWHISK>0) printf("Allocated used array");
   fflush(stdout);

   for(k=0; k<ntph; k++) {
     if(DEBUGWHISK>0) printf("k=%d\n",k);
     fflush(stdout);
      pairs[k+k*ntph] = 1;		// a detection always pairs with itself
      used[k] = 0;
      for(j=k+1; j<ntph; j++) {
	 dx = wpar[j]->x0 - wpar[k]->x0;
	 dy = wpar[j]->y0 - wpar[k]->y0;
	 pairs[j+k*n] = pairs[k+j*n] = ABS(dx) <= dmax && ABS(dy) <= dmax;
      }
   }

   if(DEBUGWHISK>0) printf("Pair stuff done\n");
   fflush(stdout);

/* Create groups from the pairs */
   ngrp = 0;

/* Allocate a large enough list to contain any single group */
   link = (int *)calloc(ntph, sizeof(int));

/* Allocate space for 6 medians */
   buf = (double *)calloc(6*ntph, sizeof(double));

   if(DEBUGWHISK>0) printf("Pairing arrays allocated\n");
   fflush(stdout);

   for(k=0; k<ntph; k++) {
      if(used[k]) continue;

/* Start a new group with unused whisker k */
      nlink = 1;
      link[0] = k;
      used[k] = 1;
/* Chase down all the friends of friends to the end of the list */
      for(j=0; j<nlink; j++) {
	 for(i=0; i<ntph; i++) {
	    if(!used[i] && pairs[i+link[j]*ntph]) {
	       link[nlink++] = i;
	       used[i] = 1;
	       if(DEBUGWHISK>0) printf("j=%d k=%d link[j]=%d nlink=%d\n", j, i, link[j], nlink);
	       fflush(stdout);
	    }
	 }
      }

#if 0
      printf("k= %d nlink= %d\n", k, nlink);
      for(i=0; i<nlink; i++) printf(" %d", link[i]);
      printf("\n");
#endif

   
/* Save this new group */

/* Representative detection */
      i = idx[link[0]] % maxdet;
      j = idx[link[0]] / maxdet;
      (*whisk)[ngrp].wpar = xcdet[j].wpar + i;
      (*whisk)[ngrp].wflux = xcdet[j].wflux + i;
      (*whisk)[ngrp].lencorr = bin*(2*xcdet[j].hlen+1);
      (*whisk)[ngrp].phicorr = xcdet[j].phicorr;

/* Median properties of the group */
      for(l=0; l<nlink; l++) {
	 i = idx[link[l]] % maxdet;
	 j = idx[link[l]] / maxdet;
	 buf[l+0*nlink] = xcdet[j].wpar[i].x0;
	 buf[l+1*nlink] = xcdet[j].wpar[i].y0;
	 buf[l+2*nlink] = xcdet[j].wpar[i].major;
	 buf[l+3*nlink] = cos(2*xcdet[j].wpar[i].phi);
	 buf[l+4*nlink] = sin(2*xcdet[j].wpar[i].phi);
	 buf[l+5*nlink] = xcdet[j].wflux[i].flux/xcdet[j].wpar[i].major; /*Flux per unit trail length*/
      }

      (*whisk)[ngrp].ndet = nlink;
      (*whisk)[ngrp].x0 = median(nlink, buf+0*nlink, NULL);
      (*whisk)[ngrp].y0 = median(nlink, buf+1*nlink, NULL);
      (*whisk)[ngrp].len = median(nlink, buf+2*nlink, NULL);
      (*whisk)[ngrp].phi = 0.5 * atan2(median(nlink, buf+4*nlink, NULL),
				       median(nlink, buf+3*nlink, NULL));
      (*whisk)[ngrp].flux = median(nlink, buf+5*nlink, NULL);
#if 0
      printf("Grp: %4d %4d %8.1f %8.1f %8.1f %8.1f\n",
	     ngrp, nlink, (*whisk)[ngrp].x0, (*whisk)[ngrp].y0,
	     (*whisk)[ngrp].len, (*whisk)[ngrp].phi*57.296);
      for(l=0; l<nlink; l++) {
	 i = idx[link[l]] % maxdet;
	 j = idx[link[l]] / maxdet;
	 printf("     %4d %4d %8.1f %8.1f %8.1f %8.1f\n",
		i, j, xcdet[j].wpar[i].x0, xcdet[j].wpar[i].y0,
		xcdet[j].wpar[i].major, xcdet[j].wpar[i].phi*57.296);
      }
#endif

      ngrp++;
   }


//   printf("whisk_cull: ntph %d ngrp %d\n", ntph, ngrp);


   free(wpar);
   free(idx);
   free(pairs);
   free(used);
   free(link);
   free(buf);
   return(ngrp);
}


/* Group all overlapping detections */
int whiskers(int n, WHISKER *w, int nx, int ny,
	     double dmax, double tmax, STREAK **grp, double dcirc)
{
   int i, j, k, swap, ngrp, *link, nlink;
   char *pairs, *used;
   double pi=4*atan(1.0);
   double avelen,aveflux;
   double theta, lam, R, c, s, dx, dy, det, xc, yc, dth, dth0, dthx, dsig;

#ifdef DUMP_THE_DIST
   FILE *fp = fopen("/tmp/foo.dist", "w");
#endif

   if(n > 30000) {
      fprintf(stderr, "Calling wiskers() with n= %d?  Fagedabodit\n", n);
      exit(1);
   }

// FIXME: tmax depends on length, csamp, and SNR, should be ~csamp/(len/bin/2)
//   dmax = 4.0;		// how many 'sigma' for circle distance pair?
//   tmax = 0.10;		// [rad] angular tolerance to pair detections

   pairs = (char *)calloc(n*n, sizeof(char));
   used = (char *)calloc(n, sizeof(char));

/* Circumscribing circle */
   R = sqrt(0.25*(nx*ny+ny*ny));
   
#ifdef SHOW_WHISKERS
   printf(" det phicor     x       y        len     phi      th1     lam1    th2     lam2\n");
#endif

/* Chord intersections with circumscribed circle */
   for(i=0; i<n; i++) {
      c = cos(w[i].phi);
      s = sin(w[i].phi);
      dx = w[i].x0 - nx/2;
      dy = w[i].y0 - ny/2;
      det = (c*dx+s*dy)*(c*dx+s*dy) - (dx*dx+dy*dy-R*R);
      if(det<0) {
	 printf("%7.1f %7.1f %7.3f %7.3f %7.1f\n", dx, dy, c, s, R);
      }
      lam = -(c*dx+s*dy) + sqrt(det);
      xc = dx + lam*c;
      yc = dy + lam*s;
      theta = atan2(yc, xc);
      w[i].lam1 = lam;
      w[i].th1 = theta;
      lam = -(c*dx+s*dy) - sqrt(det);
      w[i].lam2 = lam;
      w[i].th2 = atan2(dy+lam*s, dx+lam*c);
#ifdef SHOW_WHISKERS
      printf("%3d %7.1f  %7.1f %7.1f  %7.1f %7.1f  %7.1f %7.0f  %7.1f %7.0f\n",
	     i, w[i].phicorr*180/pi,
	     w[i].x0, w[i].y0, w[i].len, w[i].phi*180/pi,
	     w[i].th1*180/pi, w[i].lam1, w[i].th2*180/pi, w[i].lam2);
#endif

   }

/* Distance between all pairs */
/*   dsig is the difference between det to circle, divided by streak length */
/*   dth is the difference between chord intercepts on circle */

   for(k=0; k<n; k++) {
      pairs[k+k*n] = 1;		// a whisker always pairs with itself
      for(j=k+1; j<n; j++) {
	 dth0 = fmod(w[j].th1-w[k].th1 + 2*pi, 2*pi);
	 if(dth0 > pi) dth0 -= 2*pi;
	 dthx = fmod(w[j].th1-w[k].th2 + 2*pi, 2*pi);
	 if(dthx > pi) dthx -= 2*pi;

/* Are the two chord ends swapped? */
	 swap = ABS(dthx) < ABS(dth0);
	    
/* Angular and distance agreement */
	 if(! swap) {
	    dth = fmod(w[j].th2-w[k].th2 + 2*pi, 2*pi);
	    if(dth >  pi) dth -= 2*pi;
	    dth = sqrt(dth*dth + dth0*dth0);
	    dsig = sqrt( (w[j].lam1-w[k].lam1)*(w[j].lam1-w[k].lam1) + 
			 (w[j].lam2-w[k].lam2)*(w[j].lam2-w[k].lam2) );
	    dsig /= sqrt(w[j].len*w[j].len + w[k].len*w[k].len);
	 } else {
	    dth = fmod(w[j].th2-w[k].th1 + 2*pi, 2*pi);
	    if(dth >  pi) dth -= 2*pi;
	    dth = sqrt(dth*dth + dthx*dthx);
	    dsig = sqrt( (w[j].lam1+w[k].lam2)*(w[j].lam1+w[k].lam2) + 
			 (w[j].lam2+w[k].lam1)*(w[j].lam2+w[k].lam1) );
	    dsig /= sqrt(w[j].len*w[j].len + w[k].len*w[k].len);
	 }

/* Are these two paired? */
	 pairs[j+k*n] = pairs[k+j*n] = dsig <= dmax && dth <= tmax;
	 
#ifdef DUMP_THE_DIST
   fprintf(fp, "%3d %3d %d %d %7.1f %7.1f %7.2f %7.2f %7.0f %7.0f %7.0f %7.0f\n",
		 k, j, swap, pairs[j+k*n], dth0*180/pi, dthx*180/pi, dth*180/pi, dsig,
		 w[k].lam1, w[k].lam2, w[j].lam1, w[j].lam2);
#endif
      }
   }

#ifdef SHOW_THE_MATRIX
   printf("        ");
   for(k=0; k<n; k+=10) printf("%2d        ", k/10);
   printf("\n");
   printf("        ");
   for(k=0; k<n; k++) printf("%d", k%10);
   printf("\n");
   for(k=0; k<n; k++) {
      printf("k= %3d  ", k);
      for(i=0; i<n; i++) printf("%c", pairs[i+k*n]?'x':'.');
      printf("\n");
   }
#endif

/* Create groups from the pairs */
   ngrp = 0;
/* Allocate a large enough array for all groups */
   *grp = (STREAK *)calloc(n, sizeof(STREAK));

/* Allocate a large enough list a single group */
   link = (int *)calloc(n, sizeof(int));

   for(k=0; k<n; k++) {
      if(used[k]) continue;
      
/* Start a new group with unused whisker k */
      nlink = 1;
      link[0] = k;
      used[k] = 1;
/* Chase down all the friends of friends to the end of the list */
      for(j=0; j<nlink; j++) {
	 for(i=0; i<n; i++) {
	    if(!used[i] && pairs[i+link[j]*n]) {
	       link[nlink++] = i;
	       used[i] = 1;
//	       printf("%3d %3d %3d %3d\n", j, i, link[j], nlink);
	    }
	 }
      }

#if 0
      printf("k= %d nlink= %d\n", k, nlink);
      for(i=0; i<nlink; i++) printf(" %d", link[i]);
      printf("\n");
#endif


/* Save this new group */
      (*grp)[ngrp].n = nlink;
      (*grp)[ngrp].w = (WHISKER **)calloc(nlink, sizeof(WHISKER *));
      xc = yc = dx = dy = avelen = aveflux = 0;
      (*grp)[ngrp].ndet = 0;
      for(i=0; i<nlink; i++) {
	 (*grp)[ngrp].w[i] = w + link[i];
	 (*grp)[ngrp].ndet += w[link[i]].ndet;
/* Average properties */
	 xc += w[link[i]].x0;
	 yc += w[link[i]].y0;
	 dx += cos(2*w[link[i]].phi);
	 dy += sin(2*w[link[i]].phi);
	 avelen += w[link[i]].len;
	 aveflux += w[link[i]].flux*w[link[i]].len;
      }
      (*grp)[ngrp].x0 = xc / nlink;
      (*grp)[ngrp].y0 = yc / nlink;
      (*grp)[ngrp].phi = 0.5 * atan2(dy, dx);
      (*grp)[ngrp].len = avelen / nlink;
      (*grp)[ngrp].flux = aveflux / nlink;
      ngrp++;

   }
   free(pairs);
   free(used);
   free(link);
   return(ngrp);
}


/* Fit each streak using initial conditions */
int streak_fit(int nx, int ny, float *data, float *wgt, float badata,
	       int bin, double len/*unbinned*/, double phi/*rad*/, double fwhm/*unbinned*/, int fixpsf,
	       int nstrk, STREAK *strk, int *sortind, double starttime, double timeout,
	       PSF2D_PARAM **wpar, PSF2D_FLUX **wflux, FILE *fp)
{
   int i, j, k, NX=nx, err;
   int invert, trail, force, okfit;
//   double dr=atan(1.0)/45;

   double maxmove = 50;		// Fit ctr must coincide with peak this closely
   double fwmax = 1500.0;	// Maximum FW to accept
   double fwmin = 1.0;		// Minimum FW to accept
   double hgtmin = 8.0;		// Minimum peak-sky to keep
   double chinmax = 30000.0;	// Maximum chi/N to keep
   double sig = -1.0;		// Minimum SNR in the flux determination

   double minormax = 6.0;	// Maximum minor axis

   struct timeval tv1;
   double telapse;
   int ancount=0;
   gettimeofday(&tv1, NULL);
   telapse =  tv1.tv_sec + 1e-6*tv1.tv_usec - starttime;


/* Set search parameters */
   invert = 0;		// Invert the sign of the image?
   trail = 1;		// Fit a trailed PSF?
   okfit = 1;		// Does the fit have to be OK to keep?
   force = 1;		// Disregard all errors that inhibit output
   force = 0;		// Disregard all errors that inhibit output

/* Use just a 4 paramter fit */
   if(fixpsf) tpctrl_fit(4/*npar*/, 15/*ap*/, 40/*sky*/);

   tpctrl_search((int)len, (int)len, (int)len, (int)len, badata);

/* Allocate some space for the results */
   *wpar = (PSF2D_PARAM *)calloc(nstrk, sizeof(PSF2D_PARAM));
   *wflux = (PSF2D_FLUX *)calloc(nstrk, sizeof(PSF2D_FLUX));

/* Allow 40 iterations */
   PSF2D_TRAILITER = 40;


/* Fit the objects */
   ancount=0;
   for(k=0; k<nstrk; k++) {

     gettimeofday(&tv1, NULL);
     telapse =  tv1.tv_sec + 1e-6*tv1.tv_usec - starttime;
     if(telapse>timeout) continue; /*Don't analyze any more streaks*/

      (*wpar)[k].x0 = strk[sortind[k]].x0;
      (*wpar)[k].y0 = strk[sortind[k]].y0;

      if(fixpsf) {
	 (*wpar)[k].major = len;
	 (*wpar)[k].minor = fwhm;
	 (*wpar)[k].phi = phi;
      } else {
// FIXME!  Want initialized value for major at least half the length of result
	 (*wpar)[k].major = strk[sortind[k]].len;
	 (*wpar)[k].minor = 1.5;
	 (*wpar)[k].phi = strk[sortind[k]].phi;
      }
      (*wpar)[k].ninit = 5;

//      printf("%6d %d %8.1f %8.1f  %8.1f %8.1f %8.1f\n",		// JT
//	     k, fixpsf, strk[k].x0, strk[k].y0,
//	     (*wpar)[k].major, (*wpar)[k].minor, (*wpar)[k].phi*45/atan(1));


      i = NINT(strk[sortind[k]].x0);
      j = NINT(strk[sortind[k]].y0);

/* Set acceptance criteria */
      maxmove = strk[sortind[k]].len;
      tpctrl_accept(maxmove, fwmax, fwmin, hgtmin, chinmax, sig);
      
/* Enable very detailed dump for group k */
//      PSF2D_DEBUG = TPFN_TEST = 2*(k==21);
//      PSF2D_DEBUG = k==6;
//      TPFN_TEST = 2*PSF2D_DEBUG;

/* Fit one object */
      err = one_object(nx, NX, ny, data, wgt, i, j,
		       invert, trail, force, okfit,
		       (*wpar)+k, (*wflux)+k);

/* Has this gone seriously astray?  If so don't iterate */
      if((*wpar)[k].minor <= minormax && (*wpar)[k].major <= fwmax &&
	 (*wpar)[k].niter == PSF2D_TRAILITER) {
	 (*wpar)[k].minor = 2.5;	// In case it wandered off...
	 err = one_object(nx, NX, ny, data, wgt, i, j,
			  invert, trail, force, okfit,
			  (*wpar)+k, (*wflux)+k);
	 (*wpar)[k].niter += PSF2D_TRAILITER;
      }


/* Coopt bayer structure entry to save the err */
      (*wpar)[k].bayer = err;


#if 0
      printf("%2d %2d err %d init %6.0f -> %7.1f %7.1f %7.1f %6.2f %6.1f  chin %6.1f\n",
	     k, strk[sortind[k]].n, err, 2*strk[sortind[k]].len,
	     (*wpar)[k].x0, (*wpar)[k].y0,
	     (*wpar)[k].major, (*wpar)[k].minor,
	     (*wpar)[k].phi*180/atan(1), (*wpar)[k].chin);
#endif
      
/* Is this a decent looking streak? */
      if(err == 0 && (*wpar)[k].minor <= minormax
	          && (*wpar)[k].major <= fwmax
	          && (*wpar)[k].chin <= chinmax) {
/* Blabber about it? */
	 if(fp != NULL) {
	    one_report(nx, NX, ny, data, i, j, trail,
		       (*wpar)+k, (*wflux)+k, fp, NULL);
	 }
      }

/* Correct for tphot trail fit length bias of about -0.2 PSF */
      if((*wpar)[k].minor > 1.0 && (*wpar)[k].minor < minormax) {
	 (*wpar)[k].major += 0.16 * (*wpar)[k].minor;
      }
      
/* Try to help a non-converged fit with a stupid sky level */
      if((*wpar)[k].sky == badata) {
	 (*wpar)[k].sky = 0;
	 (*wpar)[k].peak += badata;
      }

/* Correct a stupid major axis to something less stupid (but still deletable) */
      (*wpar)[k].major = MIN((*wpar)[k].major, sqrt(nx*nx+ny*ny));       
      /*keep count of the number of streaks actually analyzed,
        in case the timeout prevents us from analyzing all of them*/
      ancount+=1;
   }

   if(ancount<nstrk)
     {
       fprintf(stderr,"Warning, due to timeout, streak_fit could only analyze\n");
       fprintf(stderr,"%d out of the total %d input streaks\n",ancount,nstrk);
     }
   nstrk=ancount;

/* Squeeze out dups signaled by nearly identical x,y */

   char *pairs, *used;
   int *link, ngrp, l, lmin, nlink;
   double dmax=3, chinmin, dx, dy;

   pairs = (char *)calloc(nstrk*nstrk, sizeof(char));
   used = (char *)calloc(nstrk, sizeof(char));

   for(k=0; k<nstrk; k++) {
      pairs[k+k*nstrk] = 1;		// a detection always pairs with itself
      used[k] = 0;
      for(j=k+1; j<nstrk; j++) {
	 dx = (*wpar)[j].x0 - (*wpar)[k].x0;
	 dy = (*wpar)[j].y0 - (*wpar)[k].y0;
	 pairs[j+k*nstrk] = pairs[k+j*nstrk] = ABS(dx) <= dmax && ABS(dy) <= dmax;
      }
   }


/* Create groups from the pairs */
   ngrp = 0;

/* Allocate a large enough list to contain any single group */
   link = (int *)calloc(nstrk, sizeof(int));

   for(k=0; k<nstrk; k++) {
      if(used[k]) continue;

/* Start a new group with unused streak fit k */
      nlink = 1;
      link[0] = k;
      used[k] = 1;
/* Chase down all the friends of friends to the end of the list */
      for(j=0; j<nlink; j++) {
	 for(i=0; i<nstrk; i++) {
	    if(!used[i] && pairs[i+link[j]*nstrk]) {
	       link[nlink++] = i;
	       used[i] = 1;
	    }
	 }
      }

/* If a singleton and nothing yet saved, move on */
      if(nlink == 1 && ngrp == k) {
	 ngrp++;
	 continue;
      }

/* Among this group find the one with the lowest chin */      
      chinmin = (*wpar)[link[0]].chin;
      lmin = 0;
      for(l=1; l<nlink; l++) {
	 if((*wpar)[link[l]].chin < chinmin) {
	    chinmin = (*wpar)[link[l]].chin;
	    lmin = l;
	 }
      }

/* Overwrite current grp with this one */
      (*wpar)[ngrp] = (*wpar)[link[lmin]];
      (*wflux)[ngrp] = (*wflux)[link[lmin]];
      ngrp++;
   }

   free(pairs);
   free(link);
   free(used);

   return(ngrp);
}



// #define WRITE_CLASS_VARS

/* Classify each streak */
void streak_class(int nx, int ny, int x0, int y0,
		  double zp, double fw, double sat,
		  int nstar, STAR *star,
		  int nstrk, STREAK *streak, PSF2D_PARAM *wpar)
{
   int i, k, istar;
   double pi=4*atan(1);
   double flux, fstar, dx, dy, d2, prox;
   double phi, daz, proxmin, proxsum, dmin, Ptot;
   double MAJMIN, MAJMAX, MINMIN, MINMAX, CHINMAX;
   double PMAX, DSPIKE, BURNTHRESH, BURNPHI;

   PMAX = 999;		/* Maximum probability */

   DSPIKE = 3.0;	/* Augmentation in diffraction spike directions */
   BURNTHRESH = 1;	/* Burn threshold: m < zp-2.5log(sat*fw*fw*BURNTHRESH) */
   BURNPHI = 0.2;	/* Streak angle alignment with up-down for burns */

   MINMIN = 1.0;	/* [pix] minimum fwmin to report streak */
   MINMAX = 6.0;	/* [pix] maximum fwmin to report streak */
   MAJMIN = 0.0;	/* [pix] minimum fwmax to report streak */
   MAJMAX = 2000.0;	/* [pix] maximum fwmax to report streak */
   CHINMAX = 1E5;	/* [pix] maximum chin to report streak */ 

#ifdef WRITE_CLASS_VARS
   FILE *fp;
   fp = fopen("/tmp/strker.vars", "w");
#endif

   for(k=0; k<nstrk; k++) {
      streak[k].P[Ptr] = streak[k].P[Pvr] = streak[k].P[Pmv] =
	 streak[k].P[Psc] = streak[k].P[Pbn] = streak[k].P[Pcr] = 0;

/* Off image?  All probability falls to "other" */
      if(wpar[k].x0<0 || wpar[k].x0>nx || wpar[k].y0<0 || wpar[k].y0>ny) continue;

/* Evaluate "closest star", and evaluate the burn probability for each */
      proxmin = proxsum = 1e-10;
      dmin = 1e20;
      istar = -1;
      for(i=0; i<nstar; i++) {
	 fstar = exp(-0.4*log(10)*(star[i].m-zp));
	 dx = wpar[k].x0+x0 - star[i].x;
	 dy = wpar[k].y0+y0 - star[i].y;
	 d2 = dx*dx + dy*dy;

	 prox = fstar/sat / pow(d2+fw*fw, 1.5);

/* Augment at +/-45 deg for diffraction spikes (rather STA1600-centric) */
	 phi = atan2(dy, dx);
	 phi = ABS(phi);		// range reduce det loc angle wrt star
	 if(phi > pi/2) phi = pi - phi;	// range reduce

/* Sigma offset in azimuthal direction from +/-45 deg */
	 daz = exp(-(phi-pi/4)*(phi-pi/4)*d2/(2*fw*fw));
/* Augment also in the side-to-side direction from diffraction */
	 daz += exp(-dy*dy/(2*fw*fw));
/* Augment also in the up-down direction in case of bleed */
	 daz += exp(-dx*dx/(2*fw*fw));
	 prox *= (1+ DSPIKE * daz);

	 proxsum += prox;
	 
	 if(prox > proxmin) {
	    proxmin = prox;
	    istar = i;
	 }

	 if(d2 < dmin) dmin = d2;

#if 0
	 printf("%3d %7.1f %7.1f  %3d %7.1f %7.1f  %9.1f %7.1f %7.3f  %10.1ef %10.1e\n",
		k, wpar[k].x0+x0, wpar[k].y0+y0, i, star[i].x, star[i].y,
		fstar, phi*57.296, daz, prox, proxmin);
#endif


/* Burns only occur towards the chip center from saturated stars (STA-centric) */
	 phi = ABS(wpar[k].phi);	// range reduce streak direction
	 if(phi > pi/2) phi = pi - phi;	// range reduce

#if 0
//	 if(ABS(dx)<1 && NINT(wpar[k].x0) == 8820 && ABS(dy)<60) {
//	 if(NINT(wpar[k].x0) == 1618 && NINT(wpar[k].y0) == 8892 &&
	 if(NINT(wpar[k].x0) == 839 && NINT(wpar[k].y0) == 8203 &&
	    fstar > BURNTHRESH*sat*fw*fw && ABS(dx) < 3*fw) {
	    printf("%5.2f %6.0f %6.0f  %6.1f %6.1f %10.0f %10.0f %7.3f\n",
		   star[i].m, star[i].x, star[i].y, dx, dy, fstar, BURNTHRESH*sat*fw*fw, pi/2-phi);
	 }
#endif
	 
	 if(fstar > BURNTHRESH*sat*fw*fw &&		// bright enough to burn
	    ABS(dx) < 3*fw &&				// underneath
	    (star[i].y-ny/2)*(wpar[k].y0+y0-ny/2) > 0 && // same side of ccd
	    dy*(star[i].y-ny/2) < 0 &&			// burn towards center
	    (pi/2-phi) < BURNPHI) {			// trail upwards
	    prox = fstar/sat / sqrt(d2+fw*fw);		// burns go like 1/r (?)
	    streak[k].P[Pbn] += 100 * PMAX * prox * exp(-dx*dx/(2*fw*fw));

#if 0
	    printf("%7.0f %7.0f %5.2f %6.1f %8.1f %10.0f %10.0f %7.3f %10.6f %8.0f %8d\n",
		   wpar[k].x0, wpar[k].y0, star[i].m, dx, dy, fstar, BURNTHRESH*sat*fw*fw, pi/2-phi,
		   prox, 100 * PMAX * prox * exp(-dx*dx/(2*fw*fw)), streak[k].P[Pbn]);
#endif

	 }
      }
      dmin = sqrt(dmin);

#ifdef WRITE_CLASS_VARS
      if(istar >= 0) {
	 dx = wpar[k].x0+x0 - star[istar].x;
	 dy = wpar[k].y0+y0 - star[istar].y;
	 fprintf(fp, "%7.1f %7.1f  %5.2f %6.1f %6.1f %5.2f  %10.7f\n",
		 wpar[k].x0+x0, wpar[k].y0+y0, star[istar].m, dx, dy, dmin,
		 proxmin);
      }
#endif


/* Estimate of flux in this streak */
      flux = wpar[k].peak * wpar[k].major * wpar[k].minor;
/* Flux from this star */
      fstar = 1;
      if(istar >= 0) fstar = exp(-0.4*log(10)*(star[istar].m-zp));

/* variable probability rises fast if we're really close to a star */
      streak[k].P[Pvr] = PMAX * MIN(1,fw/MAX(0.1,100*dmin));
      Ptot = streak[k].P[Pvr];

/* scar probability rises fast with prox, ~0.005 is likely a scar */
      streak[k].P[Psc] = PMAX * MIN(1, 100*proxsum);
      streak[k].P[Psc] /= (1+exp(MIN(50, 10*(flux/fstar-1))));
      streak[k].P[Psc] = MIN(PMAX-Ptot, streak[k].P[Psc]);
      Ptot += streak[k].P[Psc];

/* If it's a scar, don't call it a burn */
      streak[k].P[Pbn] = MIN(PMAX, streak[k].P[Pbn]);
      streak[k].P[Pbn] = MIN(PMAX-Ptot, streak[k].P[Pbn]);
      Ptot += streak[k].P[Pbn];

/* CR probability rises fast as minor/fw approaches 0.6 */
      streak[k].P[Pcr] = PMAX / (1+exp(MIN(50, 10*(wpar[k].minor/fw-0.6))));
/* But it shouldn't be too crazy long */
      streak[k].P[Pcr] /= (1+exp(MIN(50, 10*(wpar[k].major/fw-30))));
      streak[k].P[Pcr] = MIN(PMAX-Ptot, streak[k].P[Pcr]);
      Ptot += streak[k].P[Pcr];

/* transient probability uses a number of criteria */
      streak[k].P[Ptr] =  PMAX *
	 (wpar[k].major > MAJMIN && wpar[k].major < MAJMAX &&
	  wpar[k].minor > MINMIN && wpar[k].minor < MINMAX &&
	  wpar[k].chin < CHINMAX && wpar[k].bayer == 0);
      streak[k].P[Ptr] = MIN(PMAX-Ptot, streak[k].P[Ptr]);
      Ptot += streak[k].P[Ptr];

/* But if nobody claims this thing, disregard fit error and give it up to 80% */
      if(Ptot < PMAX) {
	 Ptot -= streak[k].P[Ptr];
	 streak[k].P[Ptr] =  (0.8*PMAX) *
	    (wpar[k].major > MAJMIN && wpar[k].major < MAJMAX &&
	     wpar[k].minor > MINMIN && wpar[k].minor < MINMAX &&
	     wpar[k].chin < CHINMAX);
	 streak[k].P[Ptr] = MIN(PMAX-Ptot, streak[k].P[Ptr]);
	 Ptot += streak[k].P[Ptr];
      }

/* mover probability rises fast as length/fw rises above 1 */
      streak[k].P[Pmv] = PMAX / (1+exp(MIN(50, 10*(-wpar[k].major/fw+1.0))));
      if(streak[k].P[Pmv] > MAJMAX) streak[k].P[Pmv] = 0;
      streak[k].P[Pmv] = MIN(PMAX - streak[k].P[Pcr], streak[k].P[Pmv]);
   }

#ifdef WRITE_CLASS_VARS
   fclose(fp);
#endif

   return;
}






/* 2D smooth an image by sig */
void smooth(double sig, int nx, int ny, float *D)
{
   int i, j;
   
   float *buf;
   buf = (float *)calloc(ny, sizeof(float));

/* x smooth */
   for(j=0; j<ny; j++) gconv(nx, D+j*nx, sig);

/* y smooth */
   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) buf[j] = D[i+j*nx];
      gconv(ny, buf, sig);
      for(j=0; j<ny; j++) D[i+j*nx] = buf[j];
   }

   free(buf);
   return;
}


/* Clip an image and bin */
void old_clipbin(double clipmin, double clipmax, int bin, float badata,
	     int nx, int ny, float *Din, float *Dout)
{
   int i, j, k, l, n;
   double sum, d;
// Why did I impose this?
//   if(bin > 1) {
      for(j=0; j<ny/bin; j++) {
	 for(i=0; i<nx/bin; i++) {
	    sum = 0.0;
	    n = 0;
	    for(l=j*bin; l<(j+1)*bin; l++) {
	       for(k=i*bin; k<(i+1)*bin; k++) {
		  d = Din[k+l*nx];
		  if(Din[k+l*nx] < clipmin) d = clipmin;
		  if(Din[k+l*nx] > clipmax) d = clipmax;
		  sum += d;
		  n++;
	       }
	    }
	    Dout[i+j*(nx/bin)] = sum / MAX(1,n);
	 }
      }
//   }
   return;
}

/* Clip an image and bin */
void clipbin(double clipmin, double clipmax, int bin, float badata,
	     int nx, int ny, float *Din, float *Dout)
{
   int i, j, k, l, n, allbad;
   double sum, d;
   for(j=0; j<ny/bin; j++) {
      for(i=0; i<nx/bin; i++) {
	 sum = 0.0;
	 n = 0;
	 allbad = 1;
	 for(l=j*bin; l<(j+1)*bin; l++) {
	    for(k=i*bin; k<(i+1)*bin; k++) {
	       if(Din[k+l*nx] != badata) {
		  d = Din[k+l*nx];
		  if(Din[k+l*nx] < clipmin) d = clipmin;
		  if(Din[k+l*nx] > clipmax) d = clipmax;
		  sum += d;
		  allbad = 0;
		  n++;
	       }
	    }
	 }
	 if(allbad) {
	    Dout[i+j*(nx/bin)] = badata;
	 } else {
	    Dout[i+j*(nx/bin)] = sum / MAX(1,n);
	 }
      }
   }
   return;
}


/* Linearly interpolate an array into a bigger one */
void linterp(int mx, int my, int MX, float *src,
	     int nx, int ny, int NX, float *dest)
{
   int i, j, k, l;
   double xrat, yrat, xfrac, yfrac;

   xrat = ((double)mx) / nx;
   yrat = ((double)my) / ny;
   for(j=0; j<ny; j++) {
      l = (j+0.5)*yrat - 0.49999999;
      yfrac = (j+0.5)*yrat - 0.5 - l;
      if(l < 0) { l++; yfrac -= 1.0; }
      if(l > my-2) { l--; yfrac += 1.0; }
      for(i=0; i<nx; i++) {
	 k = (i+0.5)*xrat - 0.49999999;
	 xfrac = (i+0.5)*xrat - 0.5 - k;
	 if(k < 0) { k++; xfrac -= 1.0; }
	 if(k > mx-2) { k--; xfrac += 1.0; }

	 dest[i+j*NX] = ( (1-xfrac)*(1-yfrac)*src[k+l*MX] +
			     xfrac *(1-yfrac)*src[(k+1)+l*MX] +
			  (1-xfrac)*   yfrac *src[k+(l+1)*MX] +
			     xfrac *   yfrac *src[(k+1)+(l+1)*MX]);
      }
   }
   return;
}

/* Linearly interpolate src(nx/xbin,ny/ybin) -> dest(nx,ny) in ncx,ncy chunks */
/* ncx,ncy must divide the src and dest dimensions evenly */
void linterp_chunk(int nx, int ny, int xbin, int ybin, int ncx, int ncy,
		   float *src, float *dest)
{
   int i, j, sx=(nx+xbin-1)/xbin, sy=(ny+ybin-1)/ybin;

   for(j=0; j<ncy; j++) {
      for(i=0; i<ncx; i++) {
	 linterp(sx/ncx, sy/ncy, sx, src+i*(sx/ncx)+j*(sy/ncy)*sx,
		 nx/ncx, ny/ncy, nx, dest+i*(nx/ncx)+j*(ny/ncy)*nx);
      }
   }
   return;
}


#define LAP_LIM 0.01
#define SAT_LIM 0.1
#define SKY_LIM 0.2

/* Clip an image using Laplacian and saturation */
void mowgrass(int nx, int ny, float *R, float badata, int ds,
              double saturate, float **G, double *sky)
{
   int i, j;
   long long int r;
   
/* Allocate an output image */
   *G = (float *)calloc(nx*ny, sizeof(float));

   r = 5121953;
   for(i=0; i<1000; i++) {
      r = (1664525L * r + 1013904223L) % 4294967296L;
      (*G)[i] = R[r%(nx*ny)];
   }
   *sky = quickmedian(1000, *G);

/* Fill G with the Laplacian of R and clip */
   for(j=ds; j<ny-ds; j++) {
      for(i=ds; i<nx-ds; i++) {
         (*G)[i+j*nx] = R[i+j*nx] -
            0.25*(R[i-ds+j*nx]+R[i+ds+j*nx]+R[i+(j-ds)*nx]+R[i+(j+ds)*nx]);
         if((*G)[i+j*nx] >  LAP_LIM*saturate ||
            (*G)[i+j*nx] < -LAP_LIM*saturate ||
            R[i+j*nx] > (*sky)+SAT_LIM*saturate ||
            R[i+j*nx] < SKY_LIM*(*sky)) {
            (*G)[i+j*nx] = badata;
         } else {
            (*G)[i+j*nx] = R[i+j*nx];
         }
      }
   }

   for(j=0; j<ny; j++) {
      for(i=0; i<ds; i++) (*G)[i+j*nx] = (*G)[nx-1-i+j*nx] = badata;
   }

   for(j=0; j<ds; j++) {
      for(i=0; i<nx; i++) (*G)[i+j*nx] = (*G)[i+(ny-1-j)*nx] = badata;
   }

   return;
}


/* Median bin D(nx,ny) by bx,by, return B, ignore pixels at badata */
void medianbin(int nx, int ny, float *D, float badata, int bx, int by, float **B)
{
   int i, j, k, l, n, nbx=(nx+bx-1)/bx, nby=(ny+by-1)/by;
   float *buf;

   buf = (float *)calloc(bx*by, sizeof(float));
   *B = (float *)calloc(nbx*nby, sizeof(float));

   for(j=0; j<nby; j++) {
      for(i=0; i<nbx; i++) {
	 n = 0;
	 for(l=j*by; l<MIN(ny-1,(j+1)*by); l++) {
	    for(k=i*bx; k<MIN(nx-1,(i+1)*bx); k++) {
	       if(D[k+l*nx] != badata) buf[n++] = D[k+l*nx];
	    }
	 }
	 (*B)[i+j*nbx] = quickmedian(n, buf);
      }
   }

   free(buf);

   return;
}


/* y value for at 0-based column i from start row j, slope dy */
#define SLOPE(i,j,slope,n,k) { k = (int)((i+0.5)*slope+j) ; \
				       if(k<0) k = 0; if(k>n-1) k = n-1; }

/* Convolve with 1/cos/sin with nlam wavelength at phi and length +/-hlen */
/* xdrift requires -pi/4 <= phi <= pi/4 */
void xdrift_tcs(double phi, int hlen, double nlam, 
	    int nx, int ny, float *D, float *C, float *S)
{
   int i, j, k, kp, ix0, ix1, im, ip;
   float *d;
   double Dp, Cp, Sp, tp, lam;
   double y0, y1, x0, x1, dy, c1, s1, cL, sL, cL1, sL1, pi=4*atan(1.0);

   d = (float *)calloc(nx, sizeof(float));

   dy = tan(phi);
   if(ABS(dy*nx) < 0.5) dy = 0;

   lam = (2*hlen+1) / nlam;
   c1 = cos(2*pi/lam);
   s1 = sin(2*pi/lam);
   cL = cos(2*pi*hlen/lam);
   sL = sin(2*pi*hlen/lam);
   cL1 = cos(2*pi*(hlen+1)/lam);
   sL1 = sin(2*pi*(hlen+1)/lam);

/* Min and max y launch point: index is int(TSK) */
   y0 = 0.5;
   y1 = ny - 0.5;
   if(dy > 0) y0 = y0 - (nx - 0.5) * dy;
   if(dy < 0) y1 = y1 - (nx - 0.5) * dy;

/* Loop over launch points for each y row */
   for(j=NINT(floor(y0)); j<(int)y1; j++) {

/* The first and last x pixel */
      x0 = 0.5;
      x1 = nx - 0.5;
      if(j < 0) x0 = -j / dy;
      if(j > ny-1) x0 = -(j-(ny-1)) / dy;
      if(j+nx*dy < 0) x1 = -j / dy;
      if(j+nx*dy > ny-1) x1 = ((ny-1)-j) / dy;

      ix0 = MIN( MAX( (int)x0, 0), nx-1);
      ix1 = MIN( MAX( (int)x1, 0), nx-1);

/* Load initial values up to forward reach of first pixel */
      for(i=MAX(ix0-hlen-1,0); i<ix0; i++) d[i] = 0.0;
      Dp = Cp = Sp = 0.0;
      for(i=ix0; i<MIN(ix0+hlen,nx-1); i++) {
	 SLOPE(i, j, dy, ny, k);
	 d[i] = D[i+nx*k];
	 Dp += d[i];
	 tp = Cp;
	 Cp = Cp*c1 + Sp*s1 + d[i]*cL;
	 Sp = Sp*c1 - tp*s1 + d[i]*sL;
      }

/* Loop over all x for this trajectory */
      for(i=ix0; i<=ix1; i++) {
	 SLOPE(i, j, dy, ny, k);
	 im = MAX(i-hlen-1, 0);
	 ip = MIN(i+hlen, nx-1);
	 SLOPE(ip, j, dy, ny, kp);
/* New original values */
	 d[ip] = D[ip+nx*kp];
	 if(i-ix0 > hlen) {
	    D[i+nx*k] = Dp + d[ip] - d[im];	 
	    C[i+nx*k] = Cp*c1 + Sp*s1 + d[ip]*cL - d[im]*cL1;
	    S[i+nx*k] = Sp*c1 - Cp*s1 + d[ip]*sL + d[im]*sL1;
	 } else {
	    D[i+nx*k] = Dp + d[ip];
	    C[i+nx*k] = Cp*c1 + Sp*s1 + d[ip]*cL;
	    S[i+nx*k] = Sp*c1 - Cp*s1 + d[ip]*sL;
	 }
	 Dp = D[i+nx*k];
	 Cp = C[i+nx*k];
	 Sp = S[i+nx*k];
      }
   }

   free(d);
   return;
}

/* Convolve with 1/cos/sin with nlam wavelength at phi and length +/-hlen */
/* ydrift requires pi/4 <= phi <= 3*pi/4 */
void ydrift_tcs(double phi, int hlen, double nlam, 
	    int nx, int ny, float *D, float *C, float *S)
{
   int i, j, k, kp, iy0, iy1, im, ip;
   float *d;
   double Dp, Cp, Sp, tp, lam;
   double y0, y1, x0, x1, dx, c1, s1, cL, sL, cL1, sL1, pi=4*atan(1.0);

   d = (float *)calloc(ny, sizeof(float));

   if(ABS(cos(phi))*ny < 0.5) {
      dx = 0;
   } else {
      dx = 1/tan(phi);
   }

   lam = (2*hlen+1) / nlam;
   c1 = cos(2*pi/lam);
   s1 = sin(2*pi/lam);
   cL = cos(2*pi*hlen/lam);
   sL = sin(2*pi*hlen/lam);
   cL1 = cos(2*pi*(hlen+1)/lam);
   sL1 = sin(2*pi*(hlen+1)/lam);

/* Min and max y launch point: index is int(TSK) */
   x0 = 0.5;
   x1 = nx - 0.5;
   if(dx > 0) x0 = x0 - (ny - 0.5) * dx;
   if(dx < 0) x1 = x1 - (ny - 0.5) * dx;

/* Loop over launch points for each x column */
   for(j=NINT(floor(x0)); j<(int)x1; j++) {

/* The first and last y pixel */
      y0 = 0.5;
      y1 = ny - 0.5;
      if(j < 0) y0 = -j / dx;
      if(j > nx-1) y0 = -(j-(nx-1)) / dx;
      if(j+ny*dx < 0) y1 = -j / dx;
      if(j+ny*dx > nx-1) y1 = ((nx-1)-j) / dx;

      iy0 = MIN( MAX( (int)y0, 0), ny-1);
      iy1 = MIN( MAX( (int)y1, 0), ny-1);

/* Load initial values up to forward reach of first pixel */
      for(i=MAX(iy0-hlen-1,0); i<iy0; i++) d[i] = 0.0;

      Dp = Cp = Sp = 0.0;
      for(i=iy0; i<MIN(iy0+hlen,ny-1); i++) {
	 SLOPE(i, j, dx, nx, k);
	 d[i] = D[k+nx*i];
	 Dp += d[i];
	 tp = Cp;
	 Cp = Cp*c1 + Sp*s1 + d[i]*cL;
	 Sp = Sp*c1 - tp*s1 + d[i]*sL;
      }

/* Loop over all y for this trajectory */
      for(i=iy0; i<=iy1; i++) {
	 SLOPE(i, j, dx, nx, k);
	 im = MAX(i-hlen-1, 0);
	 ip = MIN(i+hlen, ny-1);
	 SLOPE(ip, j, dx, nx, kp);
/* New original values */
	 d[ip] = D[kp+nx*ip];
	 if(i-iy0 > hlen) {
	    D[k+nx*i] = Dp + d[ip] - d[im];	 
	    C[k+nx*i] = Cp*c1 + Sp*s1 + d[ip]*cL - d[im]*cL1;
	    S[k+nx*i] = Sp*c1 - Cp*s1 + d[ip]*sL + d[im]*sL1;
	 } else {
	    D[k+nx*i] = Dp + d[ip];
	    C[k+nx*i] = Cp*c1 + Sp*s1 + d[ip]*cL;
	    S[k+nx*i] = Sp*c1 - Cp*s1 + d[ip]*sL;
	 }
	 Dp = D[k+nx*i];
	 Cp = C[k+nx*i];
	 Sp = S[k+nx*i];
      }
   }

   free(d);
   return;
}

/* Convolve image with top hat at phi and length +/-hlen */
/* xdrift requires -pi/4 <= phi <= pi/4 */
void xdrift_t(double phi, int hlen, int nx, int ny, float *D)
{
   int i, j, k, kp, ix0, ix1, im, ip;
   float *d;
   double Dp, y0, y1, x0, x1, dy;

   d = (float *)calloc(nx, sizeof(float));

   dy = tan(phi);
   if(ABS(dy*nx) < 0.5) dy = 0;

/* Min and max y launch point: index is int(TSK) */
   y0 = 0.5;
   y1 = ny - 0.5;
   if(dy > 0) y0 = y0 - (nx - 0.5) * dy;
   if(dy < 0) y1 = y1 - (nx - 0.5) * dy;

/* Loop over launch points for each y row */
   for(j=NINT(floor(y0)); j<(int)y1; j++) {

/* The first and last x pixel */
      x0 = 0.5;
      x1 = nx - 0.5;
      if(j < 0) x0 = -j / dy;
      if(j > ny-1) x0 = -(j-(ny-1)) / dy;
      if(j+nx*dy < 0) x1 = -j / dy;
      if(j+nx*dy > ny-1) x1 = ((ny-1)-j) / dy;

      ix0 = MIN( MAX( (int)x0, 0), nx-1);
      ix1 = MIN( MAX( (int)x1, 0), nx-1);

/* Load initial values up to forward reach of first pixel */
      for(i=MAX(ix0-hlen-1,0); i<ix0; i++) d[i] = 0.0;
      Dp = 0.0;
      for(i=ix0; i<MIN(ix0+hlen,nx-1); i++) {
	 SLOPE(i, j, dy, ny, k);
	 d[i] = D[i+nx*k];
	 Dp += d[i];
      }

/* Loop over all x for this trajectory */
      for(i=ix0; i<=ix1; i++) {
	 SLOPE(i, j, dy, ny, k);
	 im = MAX(i-hlen-1, 0);
	 ip = MIN(i+hlen, nx-1);
	 SLOPE(ip, j, dy, ny, kp);
/* New original values */
	 d[ip] = D[ip+nx*kp];
	 if(i-ix0 > hlen) {
	    D[i+nx*k] = Dp + d[ip] - d[im];	 
	 } else {
	    D[i+nx*k] = Dp + d[ip];
	 }
	 Dp = D[i+nx*k];
      }
   }

   free(d);
   return;
}

/* Convolve image with top hat at phi and length +/-hlen */
/* ydrift requires pi/4 <= phi <= 3*pi/4 */
void ydrift_t(double phi, int hlen, int nx, int ny, float *D)
{
   int i, j, k, kp, iy0, iy1, im, ip;
   float *d;
   double Dp, y0, y1, x0, x1, dx;

   d = (float *)calloc(ny, sizeof(float));

   if(ABS(cos(phi))*ny < 0.5) {
      dx = 0;
   } else {
      dx = 1/tan(phi);
   }

/* Min and max y launch point: index is int(TSK) */
   x0 = 0.5;
   x1 = nx - 0.5;
   if(dx > 0) x0 = x0 - (ny - 0.5) * dx;
   if(dx < 0) x1 = x1 - (ny - 0.5) * dx;

/* Loop over launch points for each x column */
   for(j=NINT(floor(x0)); j<(int)x1; j++) {

/* The first and last y pixel */
      y0 = 0.5;
      y1 = ny - 0.5;
      if(j < 0) y0 = -j / dx;
      if(j > nx-1) y0 = -(j-(nx-1)) / dx;
      if(j+ny*dx < 0) y1 = -j / dx;
      if(j+ny*dx > nx-1) y1 = ((nx-1)-j) / dx;

      iy0 = MIN( MAX( (int)y0, 0), ny-1);
      iy1 = MIN( MAX( (int)y1, 0), ny-1);

/* Load initial values up to forward reach of first pixel */
      for(i=MAX(iy0-hlen-1,0); i<iy0; i++) d[i] = 0.0;

      Dp = 0.0;
      for(i=iy0; i<MIN(iy0+hlen,ny-1); i++) {
	 SLOPE(i, j, dx, nx, k);
	 d[i] = D[k+nx*i];
	 Dp += d[i];
      }

/* Loop over all y for this trajectory */
      for(i=iy0; i<=iy1; i++) {
	 SLOPE(i, j, dx, nx, k);
	 im = MAX(i-hlen-1, 0);
	 ip = MIN(i+hlen, ny-1);
	 SLOPE(ip, j, dx, nx, kp);
/* New original values */
	 d[ip] = D[kp+nx*ip];
	 if(i-iy0 > hlen) {
	    D[k+nx*i] = Dp + d[ip] - d[im];	 
	 } else {
	    D[k+nx*i] = Dp + d[ip];
	 }
	 Dp = D[k+nx*i];
      }
   }

   free(d);
   return;
}

#define EDGETAPER 3.0	/* taper on smoothing closer than this to edge */

/* Convert Gaussian sigma to Y&vV coefficients */
void sig2q(double sig, double *a1, double *a2, double *a3, double *B)
{
   double q, q1, q2;
   double b0, b1, b2, b3;

   if(sig <= 0.0) {
      *a1 = *a2 = *a3 = 0.0;
      *B = 1.0;
      return;
   }

/* Ensure continuity and q>0, despite Y&vV */
   if(sig >= 3) {
      q = 0.98711*sig - 0.96330;
   } else if(sig >= 2) {
      q1 = 0.98711*sig - 0.96330;
      q2 = 3.97156 - 4.14554*sqrt(1-0.26891*sig);
      q = (sig-2) * q1 + (3-sig) * q2;
   } else if(sig >= 0.5) {
      q = 3.97156 - 4.14554*sqrt(1-0.26891*sig);
   } else {
      q2 = 3.97156 - 4.14554*sqrt(1-0.26891*0.5);
      q = 2*sig * q2;
   }

/* Recursion coefficients from Y&vV */
   b0 = 1.57825 + 2.44413*q + 1.42810*q*q + 0.422205*q*q*q;
   b1 =           2.44413*q + 2.85619*q*q + 1.266610*q*q*q;
   b2 =                      -1.42810*q*q - 1.266610*q*q*q;
   b3 =                                     0.422205*q*q*q;

   *a1 = b1 / b0;
   *a2 = b2 / b0;
   *a3 = b3 / b0;
   *B = 1 - (*a1+*a2+*a3);
}

/* Recursion from Young & van Vliet Signal Processing 44,139, 1995 */
/* Convolve a 1D array by a Gaussian with sigma: 14nsec/pix */
void gconv(int n, float *a, double sig)
{
   int i, edge=1;
   double B, a1, a2, a3;

   if(sig <= 0.0) return;
   if(n < 4) return;

/* Forward filter pass */
   for(i=3; i<n; i++) {
/* Edge pixels are indeterminate and require special care */
      if(i < EDGETAPER*sig) {
	 sig2q(i/EDGETAPER, &a1, &a2, &a3, &B);
      } else if(edge) {
	 sig2q(sig, &a1, &a2, &a3, &B);
	 edge = 0;
      }
      a[i] = B*a[i] + a1*a[i-1] + a2*a[i-2] + a3*a[i-3];
   }

/* Backward filter pass */
   for(i=n-4; i>=0; i--) a[i] = B*a[i] + a1*a[i+1] + a2*a[i+2] + a3*a[i+3];
   return;
}

/* Median using Heapsort algorithm (based on Num.Rec, ex sextractor). */
float quickmedian(int n, float *ra)
{
   int		l, j, ir, i;
   float	rra;
   if (n<2) return *ra;
   ra--;	/* Bleah, fake up fortran indexing */
   for (l = ((ir=n)>>1)+1;;) {
      if (l>1) {
	 rra = ra[--l];
      } else {
	 rra = ra[ir];
	 ra[ir] = ra[1];
	 if (--ir == 1) {
	    ra[1] = rra;
	    return n&1? ra[n/2+1] : 0.5*(ra[n/2]+ra[n/2+1]);
	 }
      }
      for (j = (i=l)<<1; j <= ir;) {
	 if (j < ir && ra[j] < ra[j+1]) ++j;
	 if (rra < ra[j]) {
	    ra[i] = ra[j];
	    j += (i=j);
	 } else {
	    j = ir + 1;
	 }
      }
      ra[i] = rra;
   }
}

#define NFACT 12	/* max divisor factor for sky sub */

/* Pick binning factors for background subtraction */
void pick_bin(int nx, int ny, int chx, int chy, int *chunkok, int *bx, int *by)
{
   int i;
   int nfact[NFACT+1], nxfact[NFACT+1], nyfact[NFACT+1];

/* Find an appropriate divisor for nx,ny for sky subtraction */
   factors(nx, NFACT, nxfact);
   factors(ny, NFACT, nyfact);
   if(VERBOSE > 1) {
      printf("x factors 2-%d from %d:", NFACT, nx);
      for(i=2; i<=NFACT; i++) printf(" %d", nxfact[i]);
      printf("\n");
      printf("y factors 2-%d from %d:", NFACT, ny);
      for(i=2; i<=NFACT; i++) printf(" %d", nyfact[i]);
      printf("\n");
   }

/* Can we respect a chx x chy (8x8) chunkification? */
   if( chx<NFACT && nxfact[chx] >= 1 && chy<NFACT && nyfact[chy] >= 1) {
      *chunkok = 1;

/* We want ~5-10 median samples over a chunk */
      factors(chx, NFACT, nfact);
      for(i=2; i<=NFACT; i++) nxfact[i] -= nfact[i];
/* Test possibilities 2*3, 5, 2*2*2, 2*2, 2*10 then give up */
      if( nxfact[6]>0 ) {
	 *bx = nx / chx / 6;
      } else if( nxfact[5]>0 ) {
	 *bx = nx / chx / 5;
      } else if( nxfact[4]>0 ) {
	 *bx = nx / chx / 4;
      } else if( nxfact[8]>0 ) {
	 *bx = nx / chx / 8;
      } else if( nxfact[10]>0 ) {
	 *bx = nx / chx / 10;
      } else {
	 *bx = nx / 6;	/* Median bin by bx for background, post-bin */
	 *chunkok = 0;
      }

/* We want ~5-10 median samples over a chunk */
      factors(chy, NFACT, nfact);
      for(i=2; i<=NFACT; i++) nyfact[i] -= nfact[i];
/* Test possibilities 2*3, 5, 2*2*2, 2*2, 2*10 then give up */
      if( nyfact[6]>0 ) {
	 *by = ny / chy / 6;
      } else if( nyfact[5]>0 ) {
	 *by = ny / chy / 5;
      } else if( nyfact[4]>0 ) {
	 *by = ny / chy / 4;
      } else if( nyfact[8]>0 ) {
	 *by = ny / chy / 8;
      } else if( nyfact[10]>0 ) {
	 *by = ny / chy / 10;
      } else {
	 *by = ny / 6;	/* Median bin by by for background, post-bin */
	 *chunkok = 0;
      }

      if(VERBOSE > 1) {
	 printf("chunk OK %d, bkg xbin %d x %d ybin %d x %d\n", *chunkok, *bx, chx, *by, chy);
      }

   } else {
      *chunkok = 0;
      *bx = nx / 6;	/* Median bin by bx for background, post-bin */
      *by = ny / 6;	/* Median bin by by for background, post-bin */
      if(VERBOSE > 1) {
	 printf("chunk not OK, bkg xbin %d ybin %d\n", *bx, *by);
      }
   }
   return;
}

/* Least divisor of num */
int least(int num)
{
   int i;
   double rt;
   rt = sqrt(num+1.0);
   for(i=2; i<=(long long int)rt; i++) if((num/i)*i == num) return(i);
   return(num);
}

/* Count of factors that divide n, from 1 to maxdivisor: nfact[0:maxdivisor] */
void factors(int n, int maxdivisor, int *nfact)
{
   int i, j, ntest;
   nfact[0] = nfact[1] = 0;	/* sort of... */
   for(i=2; i<=maxdivisor; i++) {
      ntest = n;
      nfact[i] = 0;
      for(j=0; j<n; j++) {
	 if(i*(ntest/i) == ntest) {
	    nfact[i] += 1;
	    ntest /= i;
	 } else {
	    break;
	 }
      }
   }
   return;
}

/* Open a file and return the count of lines */
FILE *fopencnt(char *path, char *mode, int *nline)
{
   int i, prev=0;
   FILE *fp;
   if( (fp=fopen(path, mode)) == NULL) return(fp);
   *nline = 0;
   while( (i=fgetc(fp)) != EOF) {
      if(i == '\n') *nline += 1;
      prev = i;
   }
   if(prev != '\n') *nline += 1;	/* Didn't end with newline */
   rewind(fp);
   return(fp);
}

/* Read a file of stars: ra, dec, x, y, m */
int read_star(char *fname, STAR **star)
{
   int n, nline;
   double dr=atan(1)/45;
   char line[1024];
   FILE *fp;
   if( (fp=fopencnt(fname, "r", &nline)) == NULL) {
      fprintf(stderr, "Cannot open star file %s\n", fname);
      return(-1);
   }

   *star = (STAR *)calloc(nline, sizeof(STAR));

   n = 0;
   while(fgets(line, 1024, fp) != NULL) {
      if(line[0] == '#') continue;
      if(sscanf(line, "%lf %lf %lf %lf %lf", &(*star)[n].ra, &(*star)[n].dec, 
	       &(*star)[n].x, &(*star)[n].y, &(*star)[n].m) != 5) {
	 fprintf(stderr, "Cannot parse 'ra dec x y m' from %s\n", line);
	 fclose(fp);
	 return(-1);
      }
      (*star)[n].ra *= dr;
      (*star)[n].dec *= dr;
      n++;
   }
   
   fclose(fp);
   return(n);
}


/* Print out cfitsio error messages and exit program */
void printerror(int status)
{
   if (status) {
      fits_report_error(stderr, status); /* print error report */
      exit( status );    /* terminate the program, returning error status */
   }
   return;
   
}


#define SKYX1 2.0 /*Distance of inner edge of sky subtraction
                    band, measured perpendicular to the trail,
                    in units of the FWHM.*/
#define SKYX2 5.0 /*Distance of outer edge of sky subtraction
                    band, measured perpendicular to the trail,
                    in units of the FWHM.*/

#define REJFRACS1 0.1 /*Trim mean rejection fraction number 1.*/
#define REJFRACS2 0.5 /*Trim mean rejection fraction number 2.
                        For a real trail, both rejection fractions
                        should give about the same mean. In particular,
                        the higher rejection fraction should not
                        lead to a much fainter mean.*/

#define ROBUSTREJ 0.1 /*Standard trim mean rejection factor*/
#define ROBUSTSIGSCALE 1.267 /*Corrects trim RMS to Gaussian sigma
                               if underlying distribution is Gaussian*/

/*trailsquish_S: January 28, 2021:
Given an input image, the central coordinates, celestial
position angle (RADIANS), half-width, and half-length for a trail,
and inner and other widths for a sky subtraction band,
collapse the trail to 1-D sum, and calculate the local
sky noise.

Pixel convention stuff: assumes input coordinates tpx,tpy
are in the tphot convention, and input image DIFF has lower
left pixel as 0,0. Uses FITS convention internally, but
input and output uses tphot convention.

Bin convention: nxbin, nybin are the pixel dimensions of
a BINNED image, but DIFF is assumed unbinned. hlen is
BINNED half-length, but fwhm is assumed to be unbinned.*/

int trailsquish_S(float *DIFF,int bin,int nxbin,int nybin,double tpx,double tpy,int hlen,double phi,double fwhm,double *trailflux1,double *trailflux2,double *skylev,double *skyrms)
{

  double parad,dx,dy,cv,sv,tx,ty,xlo,xhi;
  int trailnx,trailny,i,j,ix,iy;
  float **trailim;
  double *jvec,*skyvec,error,*skyall,sky,x,y,*squishvec,tf1,tf2;
  int skynum,skynumall,nx,ny,halflen,halfwidth,skyx1,skyx2;

  /*convert phi to PA, both in radians*/
  parad = phi - PI/2.0;
  if(parad<0.0) parad += 2.0*PI;

  /*Convert binned quantities to unbinned*/
  nx = nxbin*bin;
  ny = nybin*bin;
  halflen = hlen*bin;

  /*Convert location of local max from tph-convention to FITS convention*/
  x = tpx + 0.5;
  y = tpy + 0.5;

  /*Note well that fwhm is from the original image,
    prior to binning, and needs no correction.*/
  /*Derive photometric widths and sky-subtraction
    regions from fwhm*/
  halfwidth = ceil(fwhm/2.0);
  skyx1 = ceil(fwhm*SKYX1);
  skyx2 = ceil(fwhm*SKYX2);

  trailnx = 2*skyx2+1;
  trailny = 2*halflen+1;

  /*Trail will go vertically up through the center of the trail image*/
  trailim = mat(trailnx,trailny);
  squishvec = dvec(trailny);
  skyvec = dvec(trailny);
  skyall = dvec(trailny*trailnx); /*This vector is for estimating the background noise*/
  jvec = dvec(trailnx);

  for(i=1;i<=trailnx;i++)
    {
      for(j=1;j<=trailny;j++) trailim[i][j]=0.0;
    }

  sv = sin(parad);
  cv = cos(parad);

  /*Load trail image*/
  for(i=1;i<=trailnx;i++)
    {
      /*offset from center of trailim*/
      dx = (double)i - (double)(skyx2+1);
      for(j=1;j<=trailny;j++)
	{
	  /*offset from center of trailim*/
	  dy = (double)j - (double)(halflen+1);
	  /*Rotation to match trail on image*/
	  tx = x + dx*cv - dy*sv;
	  ty = y + dy*cv + dx*sv;
	  /*pixel at lower left on original image*/
	  ix = tx;
	  iy = ty;
	  if(ix>=1 && ix<nx && iy>=1 && iy<ny)
	    {
	      /*Interpolate pixels out of tph-convention 1-D vector DIFF into
		FITS-convention 2-D matrix trailim*/
	      xlo = (1.0-tx+(double)ix)*DIFF[nx*(iy-1) + ix-1] + (tx-(double)ix)*DIFF[nx*(iy-1) + ix];
	      xhi = (1.0-tx+(double)ix)*DIFF[nx*iy + ix-1] + (tx-(double)ix)*DIFF[nx*iy + ix];
	      trailim[i][j] = (1.0-ty+(double)iy)*xlo + (ty-(double)iy)*xhi;
	    }
	}
    }

  /*Create sky subtraction vector, and also estimate
    the background noise using skyall*/
  skynumall=0;
  for(j=1;j<=trailny;j++)
    {
      skynum=0;
      for(i=1;i<=skyx2+1-skyx1;i++)
	{
	  if(isnormal(trailim[i][j]))
	    {
	      skynum+=1;
	      jvec[skynum] = trailim[i][j];
	      skynumall+=1;
	      skyall[skynumall] = trailim[i][j];
	    }
	}
      for(i=skyx2+1+skyx1;i<=trailnx;i++)
	{
	  if(isnormal(trailim[i][j]))
	    {
	      skynum+=1;
	      jvec[skynum] = trailim[i][j];
	      skynumall+=1;
	      skyall[skynumall] = trailim[i][j];
	    }
	}
      if(skynum>0) creepvec01ds_S(jvec,skynum,ROBUSTREJ,skyvec+j,&error);
      else skyvec[j]=0.0;
    }

  /*Calculate overall background noise*/
  if(skynumall>0) creepvec01ds_S(skyall,skynumall,ROBUSTREJ,&sky,&error);
  else error=0.0;
  *skyrms = error*ROBUSTSIGSCALE;
  *skylev = sky;

  /*Collapse trail vector*/
  for(j=1;j<=trailny;j++)
    {
      squishvec[j] = 0.0;
      for(i=skyx2+1-halfwidth;i<=skyx2+1+halfwidth;i++) squishvec[j] += trailim[i][j] - skyvec[j];
      // printf("TSQvec: %d %f\n",j-halflen-1,squishvec[j]);
    }

  /*Find trim means over trail vectors for two different rejection fractions*/
  creepvec01ds_S(squishvec,trailny,REJFRACS1,&tf1,&error);
  creepvec01ds_S(squishvec,trailny,REJFRACS2,&tf2,&error);
  *trailflux1 = tf1*(double)trailny;
  *trailflux2 = tf2*(double)trailny;

  free(skyvec);
  free(skyall);
  free(jvec);
  free(squishvec);
  free_mat(trailim,trailnx,trailny);
  return(0);
}


#define BADPIXVAL -30000.0

/*trailendfind: February 05, 2021:
Given an input image, the central coordinates, celestial
position angle, half-width, and half-length for a trail,
and inner and outer widths for a sky subtraction band,
search for an end of the trail at any location from
the nominal center up to twice the half-length from this
center. Return 0 if the end is found, 1 if the trail
seems to continue beyond the search length specified
by the calling center; 2 if the trail continued at
least out to the nominal half-length but then faded
away; and -1 if the trail end cannot be satisfactorily 
measured -- e.g. the central flux was too faint OR
the trail faded to nothing at a distance less than
the nominal half-length.*/

int trailendfind(float *DIFF,int nx,int ny,double tpx,double tpy,int halflen,int halfsearchlen,int maxhalflen,double phi,double fwhm,double *endx,double *endy)
{
  
  double xcorn,ycorn,parad,dx,dy,cv,sv,tx,ty;
  int trailnx,trailny,i,j,l,xl,xh,yl,yh,ix,iy;
  float **trailim;
  double *jvec,*skyvec,error,*skyall,x,y,*squishvec;
  int skynum,skynumall,halfwidth,skyx1,skyx2;
  double skyrms,skylev,*posvec,*negvec,posav,negav,baseflux,*posflux,*negflux,sigflux;
  int peakcond,dropcond,tdone,treturn,lp1=0;
  double endpos,peakflux,*endmetric,*medvec;
  double rlev,pixarea,testx,testy,testl;
  int badoffim=0;
  
  if(tpx<=0.0 || tpx>=nx || tpy<=0.0 || tpy>=ny)
    {
      fprintf(stderr,"ERROR: trailendfind called with input coordinates\noff the image: %f %f\n",tpx,tpy);
      return(1);
    }

  if(halfsearchlen<=halflen*2)
    {
      printf("ERROR: trailendfind called with no room to search!\n");
      return(-1);
    }

  /*convert phi to PA, both in radians*/
  parad = phi - PI/2.0;
  if(parad<0.0) parad += 2.0*PI;

  /*Convert location of local max from tph-convention to FITS convention*/
  x = tpx + 0.5;
  y = tpy + 0.5;

  /*Set creeping mean rejection fraction so we could
    reject a bad region of size one PSF*/
  // rlev = ceil(fwhm)/(double)halflen;
  rlev = 0.0;

  /*Derive photometric widths and sky-subtraction
    regions from fwhm*/
  halfwidth = ceil(fwhm/2.0)+1.0;
  skyx1 = ceil(fwhm*SKYX1);
  skyx2 = ceil(fwhm*SKYX2);

  trailnx = 2*skyx2+1;
  trailny = 2*halfsearchlen+1;

  /*Trail will go vertically up through the center of the trail image*/
  trailim = mat(trailnx,trailny);
  squishvec = dvec(trailny);
  skyvec = dvec(trailny);
  skyall = dvec(trailny*trailnx); /*This vector is for estimating the background noise*/
  jvec = dvec(trailnx);
  posvec=dvec(halflen);
  negvec=dvec(halflen);
  medvec=dvec(2*halflen+1);

  for(i=1;i<=trailnx;i++)
    {
      for(j=1;j<=trailny;j++) trailim[i][j]=0.0;
    }

  sv = sin(parad);
  cv = cos(parad);

  /*Find the corners of a box on the original image that
    is guaranteed to contain the whole trail and associated
    sky subtraction regions*/
  /*Upper right corner*/
  xcorn = x - (double)(halfsearchlen+1)*sv + (double)(skyx2+1)*cv;
  ycorn = y + (double)(halfsearchlen+1)*cv + (double)(skyx2+1)*sv;
  xl=xh=xcorn;
  yl=yh=ycorn;
  /*Upper left corner*/
  xcorn = x - (double)(halfsearchlen+1)*sv - (double)(skyx2+1)*cv;
  ycorn = y + (double)(halfsearchlen+1)*cv - (double)(skyx2+1)*sv;
  if(xl>xcorn) xl=xcorn;
  if(xh<xcorn) xh=xcorn;
  if(yl>ycorn) yl=ycorn;
  if(yh<ycorn) yh=ycorn;
  /*Lower right corner*/
  xcorn = x + (double)(halfsearchlen+1)*sv + (double)(skyx2+1)*cv;
  ycorn = y - (double)(halfsearchlen+1)*cv + (double)(skyx2+1)*sv;
  if(xl>xcorn) xl=xcorn;
  if(xh<xcorn) xh=xcorn;
  if(yl>ycorn) yl=ycorn;
  if(yh<ycorn) yh=ycorn;
  /*Lower left corner*/
  xcorn = x + (double)(halfsearchlen+1)*sv - (double)(skyx2+1)*cv;
  ycorn = y - (double)(halfsearchlen+1)*cv - (double)(skyx2+1)*sv;
  if(xl>xcorn) xl=xcorn;
  if(xh<xcorn) xh=xcorn;
  if(yl>ycorn) yl=ycorn;
  if(yh<ycorn) yh=ycorn;

  /*Fudge by 1 pixel just to be sure*/
  xl--;
  yl--;
  xh++;
  yh++;
  /*Respect image boundaries*/
  if(xl<1) xl=1;
  if(xh>nx) xh=nx;
  if(yl<1) yl=1;
  if(yh>ny) yh=ny;

  /*Load trail image*/
  if(ENDFIND_DEBUG) printf("TEFPAR: xy=%f,%f,sv,cv=%f,%f,PA=%f,nx,ny=%d,%d,halflen=%d,rlev=%f\n",x,y,sv,cv,parad*180.0/PI,nx,ny,halflen,rlev);
  for(i=xl;i<=xh;i++)
    {
      for(j=yl;j<=yh;j++)
	{
	  dx = (double)i-x;
	  dy = (double)j-y;
	  tx = 1;
	  /*dot product with along-trail unit vector*/
	  ty = dy*cv - dx*sv;
	  /*dot product with across-trail unit vector*/
	  tx = dx*cv + dy*sv;
	  /*Add 0.5 so integer conversion rounds correctly*/
	  tx += 0.5;
	  ty += 0.5;
	  ix = skyx2+1 + tx;
	  iy = halfsearchlen+1 + ty;
	  /*Read pixel out of tph-convention 1-D vector DIFF into
            FITS-convention 2-D matrix trailim*/
	  if(ix>=1 && ix<=trailnx && iy>=1 && iy<=trailny)
	    {
	      if(isnormal(DIFF[nx*(j-1)+i-1]) && DIFF[nx*(j-1)+i-1]>BADPIXVAL)
		{
		  trailim[ix][iy] += DIFF[nx*(j-1)+i-1];
		}
	    }
	}
    }

  /*Create sky subtraction vector, and also estimate
    the background noise using skyall*/
  skynumall=0;
  for(j=1;j<=trailny;j++)
    {
      skynum=0;
      for(i=1;i<=skyx2+1-skyx1;i++)
	{
	  if(isnormal(trailim[i][j]))
	    {
	      skynum+=1;
	      jvec[skynum] = trailim[i][j];
	      skynumall+=1;
	      skyall[skynumall] = trailim[i][j];
	    }
	}
      for(i=skyx2+1+skyx1;i<=trailnx;i++)
	{
	  if(isnormal(trailim[i][j]))
	    {
	      skynum+=1;
	      jvec[skynum] = trailim[i][j];
	      skynumall+=1;
	      skyall[skynumall] = trailim[i][j];
	    }
	}
      if(skynum>0) creepvec01ds_S(jvec,skynum,ROBUSTREJ,skyvec+j,&error);
      else skyvec[j]=0.0;
    }

  /*Calculate overall background noise*/
  if(skynumall>0) creepvec01ds_S(skyall,skynumall,ROBUSTREJ,&skylev,&skyrms);
  else skyrms=0.0;
  skyrms *= ROBUSTSIGSCALE;

  /*Collapse trail vector*/
  for(j=1;j<=trailny;j++)
    {
      squishvec[j] = 0.0;
      for(i=skyx2+1-halfwidth;i<=skyx2+1+halfwidth;i++) squishvec[j] += trailim[i][j] - skyvec[j];
      if(ENDFIND_DEBUG) printf("TEFvec: %d %f\n",j-halfsearchlen-1,squishvec[j]);
    }

  /*Find baseline flux value*/
  for(j=halfsearchlen+1-halflen;j<=halfsearchlen+1+halflen;j++)
    {
      medvec[j-(halfsearchlen-halflen)]=squishvec[j];
    }
  pixarea = (2*halflen+1)*(2*halfwidth+1);
  sigflux = skyrms * sqrt(pixarea);
  creepvec01ds_S(medvec,halflen*2+1,TRAILCENRLEV,&baseflux,&error);
  baseflux *= (halflen*2+1);
  if(ENDFIND_DEBUG) printf("skyrms = %f, pixarea = %f, sigflux = %f, baseflux = %f, SNR = %f\n",skyrms,pixarea,sigflux,baseflux,baseflux/sigflux);

  if(baseflux/sigflux<3.0) 
    {
      /*Trail is too faint for accurate measurement*/
      if(ENDFIND_DEBUG) printf("Trail is too faint for accurate measurement, SNR = %f\n",baseflux/sigflux);
      free(skyvec);
      free(skyall);
      free(jvec);
      free(squishvec);
      free(posvec);
      free(negvec);
      free_mat(trailim,trailnx,trailny);
      free(medvec);
      return(-1);
    }

  posflux = dvec(halfsearchlen-halflen+2);
  negflux = dvec(halfsearchlen-halflen+2);
  endmetric = dvec(halfsearchlen-halflen+2);
  for(l=0;l<=halfsearchlen-halflen+1;l++)
    {
      /*Load negative vector, extending forward along the trail starting l+1 pixels from nominal center*/
      for(j=halfsearchlen+1+l;j<=halfsearchlen+l+halflen;j++) negvec[j-(halfsearchlen+l)] = squishvec[j];
      /*Load positive vector, extending backward along the trail starting l pixels from nominal center*/
      for(j=halfsearchlen+l-halflen+1;j<=halfsearchlen+l;j++) posvec[j-(halfsearchlen+l-halflen)] = squishvec[j];
      creepvec01ds_S(negvec,halflen,rlev,&negav,&error);
      creepvec01ds_S(posvec,halflen,rlev,&posav,&error);
      negflux[l+1] = negav*(halflen*2+1);
      posflux[l+1] = posav*(halflen*2+1);
      endmetric[l+1] = posflux[l+1]-negflux[l+1];
    }

  peakflux=0.0;
  peakcond=dropcond=peakflux=endpos=tdone=0;
  treturn=-1;
  l=0;
  while(l<=halfsearchlen-halflen+1 && !tdone)
    {
      if(ENDFIND_DEBUG) printf("l=%d, baseflux=%f, posflux=%f, negflux=%f, endmetric=%f, peakflux=%f, ratios: %f %f %f %f %f\n",l,baseflux,posflux[l+1],negflux[l+1],endmetric[l+1],peakflux,posflux[l+1]/sigflux,negflux[l+1]/sigflux,posflux[l+1]/baseflux,peakflux/sigflux,endmetric[l+1]/peakflux);

      /*Track peak value of endmetric reached so far*/
      if(endmetric[l+1]>peakflux && negflux[l+1]<sigflux*3.0)
	{
	  peakflux = endmetric[l+1];
	  lp1 = l;
	}
      if(posflux[l+1]<sigflux*3.0 || posflux[l+1]<baseflux*TRAILMAXFADE)
        {
          /*Trail has become too faint to trace further*/
          if(l>halflen)
	    {
	      endpos = l+0.5;
	      /*It was traced out to an interesting length before fading*/
	      *endx = tpx + cos(phi)*endpos;
	      *endy = tpy + sin(phi)*endpos;
	      tdone=1;
	      treturn=2;
	      if(ENDFIND_DEBUG) printf("Trail was traced out %f pixels, then faded out\n",endpos);
	    }
	  else
	    {
	      /*The trail faded before we even reached the end
                of the nominal input trail. This is a failed measurement*/
	      tdone=1;
	      treturn=-1;
	      if(ENDFIND_DEBUG) printf("Trail faded out after only %d pixels, less than halflen %d\n",l,halflen);
	    }
	}
      else if(!tdone)
	{
	  /*Trail is still bright enough to trace*/
	  /*Make sure we haven't run off the edge of the image*/
	  badoffim=0;
	  testl=0.0;
	  testx = tpx + cos(phi)*((double)l+0.5);
	  testy = tpy + sin(phi)*((double)l+0.5);
	  if(testx<0.0 && badoffim==0)
	    {
	      /*Trail has run off left edge of image*/
	      if(cos(phi)<0.0) testl = -tpx/cos(phi);
	      else badoffim=1; /*Should be logically impossible for good input*/
	      if(testy<0.0 && badoffim==0)
		{
		  /*Pathological case: ran into lower left corner*/
		  if(sin(phi)<0.0 && -tpy/sin(phi) < testl)
		    {
		      /*It actually hits the bottom of the image before
                        the left edge, so the bottom will define the length*/
		      testl = -tpy/sin(phi);
		    }
		  else if(sin(phi)>=0.0) badoffim=3; /*Should be logically impossible for good input*/
		}
	      else if(testy>ny && badoffim==0)
		{
		  /*Pathological case: ran into uppper left corner*/
		  if(sin(phi)>0.0 && ((double)ny-tpy)/sin(phi) < testl)
		    {
		      /*It actually hits the top of the image before
                        the left edge, so the top will define the length*/
		      testl = ((double)ny-tpy)/sin(phi);
		    }
		  else if(sin(phi)<=0.0) badoffim=4; /*Should be logically impossible for good input*/
		}
	    }
	  else if(testx>nx && badoffim==0)
	    {
	      /*Trail has run off right edge of image*/
	      if(cos(phi)>0.0) testl = ((double)nx - tpx)/cos(phi);
	      else badoffim=2; /*Should be logically impossible for good input*/
	      if(testy<0.0 && badoffim==0)
		{
		  /*Pathological case: ran into lower right corner*/
		  if(sin(phi)<0.0 && -tpy/sin(phi) < testl)
		    {
		      /*It actually hits the bottom of the image before
                        the right edge, so the bottom will define the length*/
		      testl = -tpy/sin(phi);
		    }
		  else if(sin(phi)>=0.0) badoffim=3; /*Should be logically impossible for good input*/
		}
	      else if(testy>ny && badoffim==0)
		{
		  /*Pathological case: ran into uppper right corner*/
		  if(sin(phi)>0.0 && ((double)ny-tpy)/sin(phi) < testl)
		    {
		      /*It actually hits the top of the image before
                        the right edge, so the top will define the length*/
		      testl = ((double)ny-tpy)/sin(phi);
		    }
		  else if(sin(phi)<=0.0) badoffim=4; /*Should be logically impossible for good input*/
		}
	    }
	  else if(testy<0 && badoffim==0)
	    {
	      if(sin(phi)<0.0) testl = -tpy/sin(phi);
	      else badoffim=3;
	      /*We've already covered the pathological corner cases*/
	    }
	  else if(testy>ny && badoffim==0)
	    {
	      if(sin(phi)>0.0) testl = ((double)ny-tpy)/sin(phi);
	      else badoffim=4;
	      /*We've already covered the pathological corner cases*/
	    }
	  if(testl>0.0 && badoffim==0)
	    {
	      /*Trail ran off the image*/
	      if(peakcond==1 && lp1+0.5 < testl)
		{
		  /*Go back to local max of endmetric, assign that as the trail end*/
		  endpos = lp1+0.5;
		  if(ENDFIND_DEBUG) printf("Trail ran off image at after %d pixels: End-of-trail at last local max, %f\n",l,endpos);
		}
	      else
		{
		  /*No prior local max: assign trail end where it leaves the image*/
		  endpos = testl;
		  if(ENDFIND_DEBUG) printf("Trail ran off image after %f pixels: End-of-trail at image edge.\n",endpos);
		}
	      *endx = tpx + cos(phi)*endpos;
	      *endy = tpy + sin(phi)*endpos;
	      tdone=1;
	      treturn=0;
	    }
	  if(badoffim!=0)
	    {
	      tdone=1;
	      treturn=-1;
	      if(ENDFIND_DEBUG) printf("Logically incoherent running-off-image case %d,\nincoo %f,%f, testend %f,%f, phi=%f, testl=%f\n",badoffim,tpx,tpy,testx,testy,phi*180.0/PI,testl);
	    }
	  if(l>=maxhalflen && !tdone) 
	    {
	      /*We've probed the trail out as far as we're allowed to*/
	      if(peakcond==1)
		{
		  /*Go back to local max of endmetric, assign that as the trail end*/
		  endpos = lp1+0.5;
		  if(ENDFIND_DEBUG) printf("Maximum trail length of %d pixels reached: Definitive end-of-trail at last local max, %f\n",l,endpos);
		}
	      else
		{
		  /*No prior local max: assign trail end right here*/
		  endpos = l+0.5;
		  if(ENDFIND_DEBUG) printf("Maximum trail length of %.2f pixels reached: Declaring this to be end-of-trail\n",endpos);
		}
	      *endx = tpx + cos(phi)*endpos;
	      *endy = tpy + sin(phi)*endpos;
	      tdone=1;
	      treturn=0;
	    }
	  if(negflux[l+1]<0.0 && !tdone) 
	    {
	      /*Best possible sign that trail has ended*/
	      if(peakcond==1)
		{
		  /*Go back to local max of endmetric, assign that as the trail end*/
		  endpos = lp1+0.5;
		  if(ENDFIND_DEBUG) printf("Flux dropped to zero after %d pixels: Definitive end-of-trail at last local max, %f\n",l,endpos);
		}
	      else
		{
		  /*No prior local max: assign trail end right here*/
		  endpos = l+0.5;
		  if(ENDFIND_DEBUG) printf("Flux dropped to zero after %f pixels: Definitive end-of-trail found\n",endpos);
		}
	      *endx = tpx + cos(phi)*endpos;
	      *endy = tpy + sin(phi)*endpos;
	      tdone=1;
	      treturn=0;
	    }
	  if(peakflux>sigflux*3.0 && peakflux>baseflux*TRAILMAXFADE) peakcond=1; /*Next local max in peakflux will be trail end*/
	  if(!tdone && peakcond==1 && endmetric[l+1] < peakflux/3.0) 
	    {
	      dropcond=1; /*The flux dropped enough to indicate local max is significant*/
	      /*Since peakflux got set, we know negflux dropped below 3-sigma, so we
                believe the local max marks the end of the trail.*/
	      endpos = lp1+0.5;
	      *endx = tpx + cos(phi)*endpos;
	      *endy = tpy + sin(phi)*endpos;
	      tdone=1;
	      treturn=0;
	      if(ENDFIND_DEBUG) printf("Probable end-of-trail found after %f pixels\n",endpos);
	    }
	  /*What if peakcond never gets met? That means the trail
            stays pretty much uniformly bright, and we haven't found
            the end of it.*/
	  /*What if peakcond gets set but dropcond never does?
            That's a bit weird*/
	  /*What if the trail actually ends but negflux never
            dropped below 3 sigma?*/
	}
      l++;
    }
  if(!tdone)
    {
      /*We never found the end of the trail: it must extend beyond our search range*/
      endpos = l;
      *endx = tpx + cos(phi)*endpos;
      *endy = tpy + sin(phi)*endpos;
      treturn=1;
      if(ENDFIND_DEBUG) printf("Trail continued past end of %d-pixel search region\n",l);
    }

  free(skyvec);
  free(skyall);
  free(jvec);
  free(squishvec);
  free(posvec);
  free(negvec);
  free_mat(trailim,trailnx,trailny);
  free(medvec);
  free(posflux);
  free(negflux);
  free(endmetric);
  return(treturn);
}



#undef SKYX1
#undef SKYX2
#undef REJFRACS1
#undef REJFRACS2

/*vec(): May 01, 2019: Duplicate functionality of Numerical Recipes vector() function
without using copyrighted code.*/
float *vec(int np)
{
  float *vt;
  float memuse;
  
  memuse=(float)np*(float)sizeof(float);
  if(np<1)
    {
      fprintf(stderr,"Error in vec(): called with size %d\n",np);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in vec(): allocated vector would be too big (%d = %.0f bytes\n",np,memuse);
      return(NULL);
    }
  vt = (float*) calloc(1+np,sizeof(float));
  return(vt);
}

/*mat(): May 01, 2019: Duplicate functionality of Numerical Recipes matrix() function
without using copyrighted code.*/
float **mat(int nx,int ny)
{
  float **mt;
  float memuse;
  int i;

  memuse=(float)nx*(float)ny*(float)sizeof(float);
  if(nx<1 || ny<1)
    {
      fprintf(stderr,"Error in mat(): called with dimensions %d %d\n",nx,ny);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in mat(): allocated matrix would be too big (%dx%d = %.0f bytes\n",nx,ny,memuse);
      return(NULL);
    }

  mt = (float**)calloc(nx+1,sizeof(float*));
  for(i=0;i<=nx;i++) 
    {
      mt[i]=(float*)calloc(ny+1,sizeof(float));
      if(mt[i] == NULL) return(NULL);
    }
  return(mt);
}

/*dvec(): May 01, 2019: Duplicate functionality of Numerical Recipes dvector() function
without using copyrighted code.*/
double *dvec(int np)
{
  double *vt;
  float memuse;

  memuse=(float)np*(float)sizeof(double);
  if(np<1)
    {
      fprintf(stderr,"Error in dvec(): called with size %d\n",np);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in dvec(): allocated vector would be too big (%d = %.0f bytes is greater than %.0f\n",np,memuse,MAXMEM);
      return(NULL);
    }
  vt = (double*) calloc(1+np,sizeof(double));
  return(vt);
}

/*dmat(): May 01, 2019: Duplicate functionality of Numerical Recipes dmatrix() function
without using copyrighted code.*/
double **dmat(int nx,int ny)
{
  double **mt;
  float memuse;
  int i;

  memuse=(float)nx*(float)ny*(float)sizeof(double);
  if(nx<1 || ny<1)
    {
      fprintf(stderr,"Error in mat(): called with dimensions %d %d\n",nx,ny);
      return(NULL);
    }
  else if(memuse>MAXMEM)
    {
      fprintf(stderr,"Error in mat(): allocated matrix would be too big (%dx%d = %.0f bytes\n",nx,ny,memuse);
      return(NULL);
    }

  mt = (double**)calloc(nx+1,sizeof(double*));
  for(i=0;i<=nx;i++) 
    {
      mt[i]=(double*)calloc(ny+1,sizeof(double));
      if(mt[i] == NULL) return(NULL);
    }
  return(mt);
}

/*free_mat(): May 01, 2019: free a matrix allocated by the mat() function*/
int free_mat(float **mt,int nx,int ny)
{ 
  int i;

  for(i=0;i<=nx;i++) free(mt[i]);
  free(mt);
  return(0);
}

/*free_dmat(): May 01, 2019: Free a matrix allocated with dmat().*/
int free_dmat(double **mt,int nx,int ny)
{ 
  int i;
  for(i=0;i<=nx;i++) free(mt[i]);
  free(mt);
  return(0);
}


/*creepvec01ds_S: January 28, 2021:
Exactly like creepvec01ds in oholiab.c, but
uses tsort.

Description of ancestor program creepvec01ds
from the oholiab library: Exactly like creepvec01,
but does not print warning messages.

Description of ancestor program creepvec01:
Given a vector, its length, and
a rejection factor, carries out creeping mean
averaging of the vector and outputs the result.
Also outputs the rms of the values that were
retained.  If there was only one value, the rms
is set to -1 to signify this.  Note that the
error returned is indeed the straightforward
RMS of the retained values, not an attempt to
calculate the RMS of the mean.

September 07, 2016: added code to catch the possibility
that the entire input vector is NANs. In this case the code
will return mean zero and rms -1.

April 27, 2017: replaced !isnan check with isnormal check.*/

int creepvec01ds_S(double *invec,int pnum,double rejfac,double *output,double *error)
{
  int rejnum,i,nrej;
  double rms,*invec2;
  double mean1,norm;
  int beginning, end;
  double del1, del2;
  int length;
  

  if(pnum<1)
    {
      /*There was no data. Return with an error.*/
      printf("ERROR: pnum=0: creepvec01 called with no input points!!\n");
     *output=0.0;
      *error=-1.0;
      return(1);
    } 
  else
    {
      /*We have some data. Screen it for NANs*/
      invec2 = dvec(pnum);
      if(invec2==NULL) ARRAYERR("invec2");
      length = 0;
      for(i=1;i<=pnum;i++)
	{
	  if(isnormal(invec[i]) || invec[i]==0.0)
	    {
	      length+=1;
	      invec2[length] = invec[i];
	    }
	}
    }

  /*Determine how many entries will be rejected*/
  rejnum = rejfac*length;
  
  if(length<=0)
    {
      /*We did have some data points, but they were all NANs.*/
      /*printf("ERROR: creepvec01 called on a vector of pure NANs!!\n");*/
      *output=0.0;
      *error=-1.0;
      free(invec2);
      return(1);
    }
  else if(rejnum>=length)
    {
      /*printf("Warning: you have asked creepvec01 to\n");
	printf("reject all of your points!\nAll but one will be rejected.\n");*/
      rejnum = length-1;
    }

  /*If no points are to be rejected, let it work as an ordinary average*/
  if(rejnum<=0)
    {
      norm = mean1 = 0.0;
      for(i=1;i<=length;i++)
	{
	  norm+=1.0;
	  mean1+=invec2[i];
	}
      mean1/=norm;
      rms = 0.0;
      for(i=1;i<=length;i++)
	{
	  rms+=DSQUARE(mean1 - invec2[i]);
	}
      if(norm==1.0)
	{
	  rms = -1.0;
	}
      else {
	rms = sqrt(rms/(norm-1.0));
      }
    }
  else{
    /* Sort the array */
    tsort_d(length,invec2+1,NULL);
   
    /* Magic happens here */
    
    /* Compute mean */
    norm = mean1 = 0.0;
    for(i=1;i<=length;i++)
      {
	norm+=1.0;
	mean1+=invec2[i];
      }
    mean1/=norm;
    
    beginning = 1;
    end = length;
    
    for(nrej=0;nrej<rejnum;nrej++)
      {
	del1 = fabs(invec2[beginning]- mean1);
	del2 = fabs(invec2[end] - mean1);
	
	if(del1 < del2) {
	  mean1 = (mean1*(double)norm - (double)invec2[end])/((double)norm-(double)1);
	  end--;

	}
	else {
	  mean1 = (mean1*(double)norm - (double)invec2[beginning])/((double)norm-(double)1);
	  beginning++;
	}
	norm--;
      }

    /* Recompute Mean */
    norm = mean1 = 0.0;
    for(i=beginning;i<=end;i++)
      {
	norm+=1.0;
	mean1 += invec2[i];
      }
    mean1 /= norm;

    rms = 0.0;
    for(i=beginning;i<=end;i++)
      {
	rms+=DSQUARE(invec2[i]-mean1);
      }
    if(norm==1.0)
      {
	rms = -1.0;
      }
    else{
      rms = sqrt(rms/(norm-1.0));
    }
  }    

  *output = (double)mean1;
    
  if(rms==-1.0)
    {
      *error = -1.0;
    }
  else{
    *error = rms;
  }
  
  free(invec2);
  return(0);
}

/*phifind: Febuary 03, 2020: searches a range of phi
in radians to find the value giving the largest
integrated flux along the trail.*/
int phifind(float *image,int nx,int ny,double tpx,double tpy,double traillen,double *phi,double phirange,double fwhm,double gfac,float pixsamp,float RLEV)
{
  double x,y,cv,sv,phitemp,bestphi,flux,maxflux;
  double dx,dy,fx,fy;
  double phistep,*acrossvec,*alongvec,error,gsigma;
  int halflen,halfwidth,i,j,ix,iy,sgct;
  
  /*Convert to fits convention*/
  x = tpx+0.5;
  y = tpy+0.5;

  halflen = (traillen/2.0)+0.5;
  halfwidth = (gfac*fwhm/SIGFWHM)+0.5;

  acrossvec = dvec(2*halfwidth+1);
  alongvec = dvec(2*halflen+1);

  phistep = pixsamp/(double)halflen;
  gsigma = fwhm/SIGFWHM;

  phitemp=*phi-phirange;
  maxflux=0;
  bestphi=*phi;
  if(PAFIND_DEBUG) printf("phifind called with %d,%d, traillen=%.2f, hlen=%d, %.2f,%.2f, phi=%.2f, range %.2f, step %.2f\n",nx,ny,traillen,halflen,tpx,tpy,180.0*bestphi/PI,180.0*phirange/PI,180.0*phistep/PI);
  while(phitemp<=*phi+phirange)
    {
      sv = sin(phitemp);
      cv = cos(phitemp);
      for(i=-halfwidth;i<=halfwidth;i++)
	{
	  sgct=0;
	  for(j=-halflen;j<=halflen;j++)
	    {
	      /*Calculate pixel location*/
	      dx = (double)j*cv + (double)i*sv;
	      dy = (double)j*sv - (double)i*cv;
	      fx = x+dx;
	      fy = y+dy;
	      /*round to integers*/
	      ix = fx+0.5;
	      iy = fy+0.5;
	      /*convert to lower left 0,0 convention*/
	      ix-=1;
	      iy-=1;
	      if(ix>=0 && ix<nx && iy>=0 && iy<ny) 
		{
		  if((isnormal(image[iy*nx+ix]) || image[iy*nx+ix]==0) && image[iy*nx+ix]>BADPIXVAL)
		    {
		      sgct+=1;
		      alongvec[sgct] = image[iy*nx+ix];
		    }
		}
	    }
	  if(sgct<=0) acrossvec[i+halfwidth+1] = 0.0;
	  else if(sgct==1) acrossvec[i+halfwidth+1] = alongvec[1];
	  else creepvec01ds_S(alongvec,sgct,RLEV,acrossvec+i+halfwidth+1,&error);
	}
      flux=0.0;
      for(i=-halfwidth;i<=halfwidth;i++)
	{
	  flux += exp(-DSQUARE(i)/2.0/DSQUARE(gsigma))*acrossvec[i+halfwidth+1];
	}
      if(flux>maxflux)
	{
	  maxflux=flux;
	  bestphi=phitemp;
	}
      if(PAFIND_DEBUG) printf("phi=%f, flux=%f, bestphi=%f, maxflux=%f\n",phitemp*180.0/PI,flux,bestphi*180.0/PI,maxflux);
      phitemp+=phistep;
    }
  *phi=bestphi;

  free(acrossvec);
  free(alongvec);
  return(0);
}

/*crosstrailcen01_S: January 29, 2021:
Like crosstrailcen01,
but uses a 1-D image and inputs in the tphot
pixel convention, rather than a 2-D image and
inputs in the FITS convention. Also takes
phi in RADIANS, rather than position angle
in degrees.

Description of ancestor program crosstrailcen01:
Given x,y coordinates for the
center of a trailed image; the trail length in pixels; the
image position angle in degrees, and the approximate fwhm
of the PSF, find a more accurate centroid in the cross-trail
direction using one-dimensional Gaussian centroiding.*/
int crosstrailcen01_S(float *image,int nx,int ny,double *tpx,double *tpy,double traillen,double phi,double fwhm,double gfac,int itnum,float RLEV)
{
  int itct,halflen,ilen,halfwidth,iwidth;
  double xnow,ynow,norm,cen,gsigma,x,y;
  double gaussfac,sv,cv;
  double *alongvec,*acrossvec,fx,fy,error;
  int ix,iy,i,j,sgct;

  if(*tpx<0.0 || *tpx>nx || *tpy<0.0 || *tpy>ny)
    {
      fprintf(stderr,"ERROR: crosstrailcen01_S called with input coordinates\noff the image: %f %f\n",*tpx,*tpy);
      return(1);
    }

  halflen = (traillen/2.0)+0.5;
  ilen = halflen*2+1;
  halfwidth = (gfac*fwhm/SIGFWHM)+0.5;
  iwidth = halfwidth*2+1;

  alongvec = dvec(ilen);
  acrossvec = dvec(iwidth);

  x = *tpx+0.5;
  y = *tpy+0.5;

  xnow=x;
  ynow=y;
  gsigma = fwhm/SIGFWHM;
  sv = sin(phi);
  cv = cos(phi);

  for(itct=1;itct<=itnum;itct++)
    {
      cen=norm=0.0;
      for(i=-halfwidth;i<=halfwidth;i++)
	{
	  /*Load pixel vector parallel to the trail*/
	  sgct=0;
	  for(j=-halflen;j<=halflen;j++)
	    {
	      fx = xnow + (double)j*cv + (double)i*sv;
	      fy = ynow + (double)j*sv - (double)i*cv;
	      /*Nearest-neighbor interpolation, plus correction to lower left 0,0 pixel convention.*/
	      ix = fx+0.5-1.0;
	      iy = fy+0.5-1.0;
	      if(ix>=0&&ix<nx&&iy>=0&&iy<ny)
		{
		  if((isnormal(image[iy*nx+ix]) || image[iy*nx+ix]==0.0) && image[iy*nx+iy]>BADPIXVAL)
		    {
		      sgct+=1;
		      alongvec[sgct] = image[iy*nx+ix];
		    }
		}
	    }
	  /*Collapse vector and store the output*/
	  if(sgct>1) creepvec01ds_S(alongvec,sgct,RLEV,acrossvec+i+halfwidth+1,&error);
	  else if (sgct==1) acrossvec[i+halfwidth+1] = alongvec[1];
	  else acrossvec[i+halfwidth+1] = 0.0;
	  /*Calculate Gaussian centroid*/
	  gaussfac = exp(-pow((double)i/gsigma,2.0)/(double)2.0);
	  norm += acrossvec[i+halfwidth+1]*gaussfac;
	  cen += (double)i*acrossvec[i+halfwidth+1]*gaussfac;
	}
      /*Calculate centroid in the cross-trail direction*/
      if(norm>0.0) cen/=norm;
      else 
	{
	  cen=0.0;
	}
      /*Calculate new center of trail*/
      xnow += cen*sv;
      ynow -= cen*cv;
    }
  x=xnow;
  y=ynow;

  /*Convert back to tphot convention*/
  if(x-0.5>0 && x-0.5<nx && y-0.5>0 && y-0.5<ny)
    {
      *tpx = x-0.5;
      *tpy = y-0.5;
    }
  else
    {
      fprintf(stderr,"WARNING: crosstrailcen01_S obtained output coordinates\noff the image: %f %f\n",x-0.5,y-0.5);
      free(alongvec);
      free(acrossvec);
      return(1);
    }

  free(alongvec);
  free(acrossvec);
  return(0);
}

/*trailwidth: February 02, 2021:
Given an input image, the central coordinates, celestial
position angle, half-width, and half-length for a trail,
collapse it in the along trail direction, thus creating
an averaged cross section of the trail, and then calculate
the Gaussian width of the trail.*/
int trailwidth(float *DIFF,int nx,int ny,double tpx,double tpy,int halflen,double phi,double fwhm,double *trailwidth,double *flux)
{
  int trailnx,trailny,halfwidth=ceil(TRAILWIDTH_FWHMSCALE*fwhm);
  float **trailim;
  double sv,cv,dx,dy,tx,ty,error;
  int i,j,ix,iy,sgct=0;
  double parad,x,y,*crossvec,*jvec;
  double test_fwhm,sigma,xal,xsq,yal,xty,gval,delta,a,b;
  double chival,minchi,bfwhm;
  minchi=1.0e30;

  // printf("trailwidth recieved nx,ny = %d,%d, halflen=%d, halfwidth=%d\n",nx,ny,halflen,halfwidth);

  /*convert phi to PA, both in radians*/
  parad = phi - PI/2.0;
  if(parad<0.0) parad += 2.0*PI;
  /*Convert location of local max from tph-convention to FITS convention*/
  x = tpx + 0.5;
  y = tpy + 0.5;
 
  trailnx = 2*halfwidth+1;
  trailny = 2*halflen+1;

  /*Trail will go vertically up through the center of the trail image*/
  trailim = mat(trailnx,trailny);
  crossvec = dvec(trailnx);
  jvec = dvec(trailny);

  for(i=1;i<=trailnx;i++)
    {
      for(j=1;j<=trailny;j++) trailim[i][j]=0.0;
    }

  sv = sin(parad);
  cv = cos(parad);

  /*Load trail image*/
  for(i=1;i<=trailnx;i++)
    {
      for(j=1;j<=trailny;j++)
	{
	  /*Offset from trail center*/
	  dx = (double)(i-halfwidth-1);
	  dy = (double)(j-halflen-1);
	  /*predict floating point position on image*/
	  tx = x - dy*sv + dx*cv;
	  ty = y + dy*cv + dx*sv;
	  /*Add 0.5 so integer conversion rounds correctly*/
	  ix = tx + 0.5;
	  iy = ty + 0.5;
	  /*Read pixel out of tph-convention 1-D vector DIFF into
            FITS-convention 2-D matrix trailim*/
	  if(ix>=1 && ix<=nx && iy>=1 && iy<=ny) trailim[i][j] = DIFF[nx*(iy-1)+ix-1];
	  // printf("i,j = %d,%d, %f %f\n",i,j,DIFF[nx*(iy-1)+ix-1],trailim[i][i]);
	}
    }

  /*Collapse trail vector*/
  for(i=1;i<=trailnx;i++)
    {
      sgct=0;
      for(j=1;j<=trailny;j++)
	{
	  if(isnormal(trailim[i][j]))
	    {
	      sgct+=1;
	      jvec[sgct] = trailim[i][j];
	      // printf("i,j=%d,%d, sgct=%d, %f %f\n",i,j,sgct,trailim[i][j],jvec[sgct]);
	    }
	}
      if(sgct>1) creepvec01ds_S(jvec,sgct,TRAILCENRLEV,crossvec+i,&error);
      else if (sgct==1) crossvec[i] = jvec[1];
      else crossvec[i]=0.0;
      // printf("crossvec = %f\n",crossvec[i]);
    }


  bfwhm = test_fwhm = fwhm;
  while(test_fwhm<=TRAILWIDTH_FWHMSCALE*fwhm)
    {
      sigma = test_fwhm/2.354820045;
      /*Simple least-squares fit: 1-D Gaussian + constant background*/
      xal=xsq=yal=xty=0.0;
      // printf("FWHM = %f\n",test_fwhm);
      for(i=1;i<=trailnx;i++)
	{
	  dx = i-halfwidth-1;
	  gval =  exp(-dx*dx/2.0/sigma/sigma);
	  xal += gval;
	  xsq += gval*gval;
	  xty += gval*crossvec[i];
	  yal += crossvec[i];
	  // printf("%d %f %f\n",i,gval,crossvec[i]);
	}
      delta = (double)trailnx*xsq - xal*xal;
      a = (xsq*yal - xal*xty)/delta;
      b = ((double)trailnx*xty - xal*yal)/delta;
      /*Calculate chi-square value*/
      chival=0;
      for(i=1;i<=trailnx;i++) 
	{
	  dx = i-halfwidth-1;
	  gval =  exp(-dx*dx/2.0/sigma/sigma);
	  chival += DSQUARE(crossvec[i] - a - b*gval);
	}
      if(test_fwhm == fwhm) 
	{
	  minchi=chival;
	  *flux = b*xal;
	}
      else if(chival<minchi)
	{
	  minchi=chival;
	  bfwhm=test_fwhm;
	  *flux = b*xal*(2*halflen+1);
	}
      test_fwhm += fwhm*TRAILWIDTH_FWHMSAMP;
    }
  *trailwidth = bfwhm;

  free_mat(trailim,trailnx,trailny);
  free(crossvec);
  free(jvec);
  return(0);
}

#undef BADPIXVAL

int writeemptyfile(char strkfile[])
{
  FILE *fp;
  if(strkfile != NULL) {
      if(strcmp(strkfile, "-") == 0) fp = stdout;
      else if( (fp=fopen(strkfile, "w")) == NULL) {
	 fprintf(stderr, "Cannot open streak file `%s'\n", strkfile);
	 exit(1);
      }
      fprintf(fp, "#  x0       y0       peak   dpeak  sky    flux    dflux   length  minor    phi  err nit ndet chin Ptr Pvr Pcr Pbn Psc Pmv\n");
      fclose(fp);
  }
  return(0);
}
