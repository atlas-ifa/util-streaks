 /* Read an image, find objects, write table */
/* Usage:
 *
 *        tphot image_file [options]
 *
 */
/* v1.0 120929 John Tonry */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psf2d.h"

/* Failure codes */
#define FAIL_XBORDER  0x0001	/* Failed because outside x border */
#define FAIL_YBORDER  0x0002	/* Failed because outside y border */
#define FAIL_XMOVE    0x0004	/* Failed because x moved too far */
#define FAIL_YMOVE    0x0008	/* Failed because y moved too far */
#define FAIL_MAJBIG   0x0010    /* Failed because major > fwmax */
#define FAIL_MINBIG   0x0020    /* Failed because minor > fwmax */
#define FAIL_MAJSMALL 0x0040    /* Failed because major < fwmin */
#define FAIL_MINSMALL 0x0080    /* Failed because minor < fwmin */
#define FAIL_PEAK     0x0100    /* Failed because peak < netmin */
#define FAIL_CHIN     0x0200    /* Failed because chin > chinmax */
#define FAIL_FLUXNEG  0x0400    /* Failed because flux < 0 */
#define FAIL_FERRNEG  0x0800    /* Failed because dflux < 0 */
#define FAIL_SNR      0x1000    /* Failed because flux/dflux < fitsig */
#define FAIL_OKFIT    0x2000    /* Failed because fit was not OK */

int syntax(char *prog);

/* Print reasons for failure */
int failcode(int failure);

/* Compute median and median RMS for an array */
int imgstat(int x0, int y0, int x1, int y1, int NX, int NY, float badata,
	    float *data, double *med, double *rms, int *nok);

/* Evaluate the 4 point interpolation of sky at image index k */
double skinterp(int i, int j, int nx, int ny,
		int nxsky, int nysky, double *skymed,
		double *skyrms, double *rms);

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
/* Subtract the fit for just one object */
void one_subtract(int nx, int NX, int ny, float *data, float *resid,
		  int i, int j, int trail, int subtract, PSF2D_PARAM *wpar);





/* Quicksort program for doubles: 1980 (John Tonry) */
int qsort8(int n, double *x);

/* Compute moments around a given center and sky level */
int mompsf(int nx, int ny, float *a, float badata,
	   int aprad, double x0, double y0, double sky, PSF2D_MOM *mom);

/* Compute the value of a Waussian fit at (x,y) (leaving out sky) */
double wauss(double x, double y, PSF2D_PARAM *wpar, double *z2);

/* Compute the value of a trailed waussian fit at (x,y) (leaving out sky) */
double wtrail(double x, double y, PSF2D_PARAM *wpar, double *z2);

/* Calculate sx2, sxy, sy2 from sigma's and position angles */
int dosxy(PSF2D_PARAM *par);

int testfits_(char *file, int *mefits, int *bitpix, int filelen);
/* rfitsreal() reads FITS image from disk and returns a header and real array */
int rfitsreal(char **head, int *nx, int *ny, float **data, char *file);
/* wfitsreal() will write a real FITS image onto disk */
int wfitsreal(char *head, void *data, char *file);

int nitfread(char *nitfile, int bitpix,
	     char **fitshead, int *nx, int *ny, float **data, int verbose);
/* Read a FITS file into a floating array */
int cf_readreal(char **head, int *nx, int *ny, float **data, char *file);
/* chfitshead() changes a line in a FITS header */
int chfitshead(int *n, char *header, char *keyword, char *cvalue, int ivalue,
	       double rvalue);
/* ifitshead() reads an integer component of a FITS header */
int ifitshead(char *header, char *keyword, int *ivalue);

int TPFN_TEST=0;
extern int PSF2D_DEBUG;

#define MAXSRT (4093)
#define RESIDRAD (2.5)
#define MAXBIN (10)

/* Global control parameters and defaults */

/* Data file description */
static float badata = 0.0;	// Bad data value
static double sat = 100000;	// Saturated pixel level
static double eadu = 2.0;	// e/ADU for statistics
static double bias = 0.0;	// bias = Image_background - true_sky
static int x0border = 2;	// Image border to avoid
static int x1border = 2;	// Image border to avoid
static int y0border = 2;	// Image border to avoid
static int y1border = 2;	// Image border to avoid
/* Sky specification if non-constant */
static int nxsky = 0;		// Subdivide image in x for sky into M chunks
static int nysky = 0;		// Subdivide image in y for sky into N chunks
static int mosaicsky = 0;	// Treat sky as MxN discontinuous chunks?

/* Trigger specification */
static double srchcut = 100;	// Minimum value above 0 to trigger
static double srchmax = 1e5;	// Maximum peak to trigger
static double srchrad = 5.0;	// Local max required over radius r to trigger
static double srchsig = 5.0;	// Minimum value = sky+sig*rms to trigger

/* Fit specification */
static int nwpar = 7;		// How many Waussian params to fit?
static int aprad = 15;		// Aperture radius for photometry
static int skyrad = 40;		// Sky radius for photometry

/* Acceptance criteria */
static double maxmove = 3.0;	// Fit ctr must coincide with peak this closely
static double fwmax = 100.0;	// Maximum FW to accept
static double fwmin = 0.5;	// Minimum FW to accept
static double netmin = 10.0;	// Minimum peak-sky to keep
static double chinmax = 1000.0;	// Maximum chi/N to keep
static double fitsig = -1.0;	// Minimum SNR in the flux determination

/* Set search parameters */
void tpctrl_search(int dx0, int dx1, int dy0, int dy1, float bad)
{
   x0border = dx0;
   x1border = dx1;
   y0border = dy0;
   y1border = dy1;
   badata = bad;
   return;
}

/* Set trigger parameters */
void tpctrl_trigger(double min, double max, double rad, double sig)
{
   srchcut = min;
   srchmax = max;
   srchrad = rad;
   srchsig = sig;
   return;
}

/* Set fit parameters */
void tpctrl_fit(int npar, int ap, int sky)
{
   nwpar = npar;
   aprad = ap;
   skyrad = sky;
   return;
}

/* Set acceptance parameters */
void tpctrl_accept(double move, double wmax, double wmin, double hgt,
		   double chin, double sig)
{
   maxmove = move;
   fwmax = wmax;
   fwmin = wmin;
   netmin = hgt;
   chinmax = chin;
   fitsig = sig;
   return;
}

/* Report all the internal parameters */
void tpctrl_pareport()
{
   printf("Image border to avoid = %d\n", x0border);
   printf("Image border to avoid = %d\n", x1border);
   printf("Image border to avoid = %d\n", y0border);
   printf("Image border to avoid = %d\n", y1border);
   printf("Subdivide image in x for sky into M chunks = %d\n", nxsky);
   printf("Subdivide image in y for sky into N chunks = %d\n", nysky);
   printf("Treat sky as MxN discontinuous chunks? = %d\n", mosaicsky);
   printf("How many Waussian params to fit? = %d\n", nwpar);
   printf("Aperture radius for photometry = %d\n", aprad);
   printf("Sky radius for photometry = %d\n", skyrad);
   printf("Bad data value = %.1f\n", badata);
   printf("Saturated pixel level = %.1f\n", sat);
   printf("e/ADU for statistics = %.1f\n", eadu);
   printf("bias (Image_background - true_sky) = %.1f\n", bias);
   printf("Minimum value above 0 to trigger = %.1f\n", srchcut);
   printf("Maximum peak to trigger = %.1f\n", srchmax);
   printf("Local max required over radius r to trigger = %.1f\n", srchrad);
   printf("Minimum value = sky+sig*rms to trigger = %.1f\n", srchsig);
   printf("Fit ctr must coincide with peak this closely = %.1f\n", maxmove);
   printf("Maximum FW to accept = %.1f\n", fwmax);
   printf("Minimum FW to accept = %.1f\n", fwmin);
   printf("Minimum peak-sky to keep = %.1f\n", netmin);
   printf("Maximum chi/N to keep = %.1f\n", chinmax);
   printf("Minimum SNR in the flux determination = %.1f\n", fitsig);
   return;
}


int tpfn_main(int argc, char *argv[])
// int main(int argc, char *argv[])
{
   int i, j, npix, nx, NX, ny, sx, sy, wgtnx, wgtny;
   int ii, jj, err, mefits, bitpix, nfind, nfit, trail, force, subtract;
   char *imname, *head, *flatfile, *headflat, *momfile, *wgtfile, *wgthead;
   char *fout, *objinput, *residfile;
   char line[10240];
   int nxflat, nyflat, fieldhdr, rdntf, binfactor, okfit, nwinit;
   float *data, *wgt, *resid, *flat;
   double major, minor, phi;
   double pi=4*atan(1.0), z2;
   int nskyok, neadu, invert, residsky;
   double *skymed, *skyrms, eaduin, sky, rms, flatnorm;
   FILE *fp, *infp=NULL, *fpmom=NULL;
   PSF2D_PARAM wpar;
   PSF2D_FLUX wflux;
   int nthread, maxthread, thread, stripass, stripwidth, ymin, ymax;

/* Default parameters */
   bitpix = -32;	// Want floating point data
   invert = 0;		// Invert the sign of the image?
   fieldhdr = 0;	// Write details of all fields?

   trail = 0;		// Fit a trailed PSF?
   okfit = 1;		// Does the fit have to be OK to keep?
   force = 0;		// Disregard all errors that inhibit output
   subtract = 0;	// Subtract fit from image as we go?

   binfactor = 1;	// Bin the input data?
   okfit = 1;		// Does the fit have to be OK to keep?
   nwinit = 0;		// Initialization params to read if obj input

   fout = "-";		// Output file name
   rdntf = 0;		// Read NITF instead of FITS?

   eaduin = -1.0;	// Electrons per ADU to estimate noise [-1 => auto]

   major = -1.0;	// External specification of major axis?
   minor = -1.0;	// External specification of minor axis?
   phi = -100.0;	// External specification of angle of PSF?

   wgtfile = NULL;	// External weight (inverse variance) file
   wgt = NULL;		// External weight (inverse variance) array

   residsky = 0;	// Subtract sky estimate from residual image?
   objinput = NULL;	// File with x,y values to fit

   residfile = NULL;	// Residual output file
   resid = NULL;	// Residual image

   flatfile = NULL;	// Flatfield file
   flatnorm = 0.0;	// Flatfield normalization factor

   momfile = NULL;	// Moments file

   nthread = 1;		// Number of threads to use [default no MP]
   maxthread = 0;	// Max number of threads available [default no MP]

/* Bomb out if no arguments, now that defaults are set */
   if(argc < 2) {
      syntax(argv[0]);
      exit(0);
   }

   imname = argv[1];	// First argument is the image file

/* Parse the arguments */
   for(i=2; i<argc; i++) {

      if(strcmp(argv[i], "-out") == 0) {
	 fout = argv[++i];

      } else if(strcmp(argv[i], "-obj") == 0) {
	 objinput = argv[++i];
	 nwinit = 2;

      } else if(strcmp(argv[i], "-resid") == 0) {
	 residfile = argv[++i];

      } else if(strcmp(argv[i], "-wgt") == 0) {
	 wgtfile = argv[++i];

      } else if(strcmp(argv[i], "-subtract") == 0) {
	 subtract = 1;

      } else if(strcmp(argv[i], "-residsky") == 0) {
	 residsky = 1;

      } else if(strcmp(argv[i], "-nitf") == 0) {
	 rdntf = 1;

      } else if(strcmp(argv[i], "-flat") == 0) {
	 flatfile = argv[++i];

      } else if(strcmp(argv[i], "-flatnorm") == 0) {
	 sscanf(argv[++i], "%lf", &flatnorm);

      } else if(strcmp(argv[i], "-moment") == 0) {
	 momfile = argv[++i];

      } else if(strcmp(argv[i], "-trail") == 0) {
	 trail = 1;

      } else if(strcmp(argv[i], "-verb") == 0) {
	 TPFN_TEST = 1;

      } else if(strcmp(argv[i], "-VERB") == 0) {
	 TPFN_TEST = 2;

      } else if(strcmp(argv[i], "-debug") == 0) {  /* Undocumented feature */
	 PSF2D_DEBUG = 1;			   /* employ psf2d debug */

      } else if(strcmp(argv[i], "-hdr") == 0) {
	 fieldhdr = 1;

      } else if(strcmp(argv[i], "-bin") == 0) {
	 sscanf(argv[++i], "%d", &binfactor);
	 if(binfactor > MAXBIN || binfactor < 1) {
	    fprintf(stderr, "Error: binfactor %d > MAXBIN %d\n", binfactor, MAXBIN);
	    exit(1);
	 }

      } else if(strcmp(argv[i], "-invert") == 0) {
	 invert = 1;

      } else if(strcmp(argv[i], "-border") == 0) {
	 sscanf(argv[++i], "%d", &x0border);
	 x1border = y0border = y1border = x0border;

      } else if(strcmp(argv[i], "-xlft") == 0) {
	 sscanf(argv[++i], "%d", &x0border);

      } else if(strcmp(argv[i], "-xrgt") == 0) {
	 sscanf(argv[++i], "%d", &x1border);

      } else if(strcmp(argv[i], "-ybot") == 0) {
	 sscanf(argv[++i], "%d", &y0border);

      } else if(strcmp(argv[i], "-ytop") == 0) {
	 sscanf(argv[++i], "%d", &y1border);

      } else if(strcmp(argv[i], "-nxsky") == 0) {
	 sscanf(argv[++i], "%d", &nxsky);

      } else if(strcmp(argv[i], "-nysky") == 0) {
	 sscanf(argv[++i], "%d", &nysky);

      } else if(strcmp(argv[i], "-mosaicsky") == 0) {
	 mosaicsky = 1;

      } else if(strcmp(argv[i], "-rad") == 0) {
	 sscanf(argv[++i], "%lf", &srchrad);

      } else if(strcmp(argv[i], "-badata") == 0) {
	 sscanf(argv[++i], "%f", &badata);

      } else if(strcmp(argv[i], "-sig") == 0) {
	 sscanf(argv[++i], "%lf", &srchsig);

      } else if(strcmp(argv[i], "-min") == 0) {
	 sscanf(argv[++i], "%lf", &srchcut);

      } else if(strcmp(argv[i], "-max") == 0) {
	 sscanf(argv[++i], "%lf", &srchmax);

      } else if(strcmp(argv[i], "-major") == 0) {
	 sscanf(argv[++i], "%lf", &major);

      } else if(strcmp(argv[i], "-minor") == 0) {
	 sscanf(argv[++i], "%lf", &minor);

      } else if(strcmp(argv[i], "-phi") == 0) {
	 sscanf(argv[++i], "%lf", &phi);
	 phi *= pi/180;	// convert to rad

      } else if(strcmp(argv[i], "-npar") == 0) {
	 sscanf(argv[++i], "%d", &nwpar);

      } else if(strcmp(argv[i], "-ninit") == 0) {
	 sscanf(argv[++i], "%d", &nwinit);

      } else if(strcmp(argv[i], "-fwmax") == 0) {
	 sscanf(argv[++i], "%lf", &fwmax);

      } else if(strcmp(argv[i], "-fwmin") == 0) {
	 sscanf(argv[++i], "%lf", &fwmin);

      } else if(strcmp(argv[i], "-net") == 0) {
	 sscanf(argv[++i], "%lf", &netmin);

      } else if(strcmp(argv[i], "-move") == 0) {
	 sscanf(argv[++i], "%lf", &maxmove);

      } else if(strcmp(argv[i], "-okfit") == 0) {
	 sscanf(argv[++i], "%d", &okfit);

      } else if(strcmp(argv[i], "-force") == 0) {
	 force = 1;

      } else if(strcmp(argv[i], "-chin") == 0) {
	 sscanf(argv[++i], "%lf", &chinmax);

      } else if(strcmp(argv[i], "-snr") == 0) {
	 sscanf(argv[++i], "%lf", &fitsig);

      } else if(strcmp(argv[i], "-eadu") == 0) {
	 if(strcmp(argv[++i], "auto") == 0) {
	    eaduin = -1;
	 } else {
	    sscanf(argv[i], "%lf", &eaduin);
	 }

      } else if(strcmp(argv[i], "-bias") == 0) {
	 sscanf(argv[++i], "%lf", &bias);

      } else if(strcmp(argv[i], "-sat") == 0) {
	 sscanf(argv[++i], "%lf", &sat);

      } else if(strcmp(argv[i], "-aprad") == 0) {
	 sscanf(argv[++i], "%d", &aprad);

      } else if(strcmp(argv[i], "-skyrad") == 0) {
	 sscanf(argv[++i], "%d", &skyrad);

      } else {
	 fprintf(stderr, "Unrecognized argument `%s'\n\n", argv[i]);
	 syntax(argv[0]);
	 exit(1);
      }   
   }

/* nwinit signals how many values to read from objinput, default is 2 */
   if(objinput == NULL) nwinit = 0;
   if(nwinit != 0 && nwinit != 2 && nwinit != 5) {
      fprintf(stderr, "Illegal nwinit %d\n", nwinit);
      exit(1);
   }

/* Are we being requested to do a 4 param fit? */
   if(major > 0 && minor > 0 && phi > -pi && 
      (nwinit == 0 || nwinit == 2) && nwpar >= 4) nwpar = 4;

/* Verify that we have shape and external positions for 2 param fit */
   if(nwpar <= 2 && 
      (major < 0 || minor < 0 || phi < -pi || objinput == NULL)) {
      fprintf(stderr, "Findobj: must specify input file and valid major (%.2f), minor (%.2f), phi (%.2f) for npar=2\n", major, minor, phi);
      exit(1);
   }

/* Sky is either constant or subdivided in both dimensions */
/* (Of course the subdivision can be 1 for no subdivision) */
   if(nxsky == 0 || nysky == 0) nxsky = nysky = 0;

/* Open the output file */
   if(strcmp(fout, "-") == 0) {
      fp = stdout;
   } else if( (fp = fopen(fout, "w")) == NULL) {
      fprintf(stderr, "Cannot open output file '%s'\n", fout);
      syntax(argv[0]);
      exit(1);
   }

/* Open an input file? */
   if(objinput != NULL) {
      if(strcmp(objinput, "-") == 0) {
	 infp = stdin;
      } else if( (infp = fopen(objinput, "r")) == NULL) {
	 fprintf(stderr, "Cannot open input file '%s'\n", objinput);
	 syntax(argv[0]);
	 exit(1);
      }
   }

/* Open a moment file? */
   if(momfile != NULL && (fpmom = fopen(momfile, "w")) == NULL) {
      fprintf(stderr, "Cannot open moment file '%s'\n", momfile);
      syntax(argv[0]);
      exit(1);
   }

/* Does the image look like FITS? */
   if(!rdntf) {
      err = testfits_(imname, &mefits, &bitpix, strlen(imname)+1);
      rdntf = mefits != 1;
      bitpix = -32;
   }

/* Read the image */
   if(!rdntf) {		// Read FITS
      data = NULL;
      err = rfitsreal(&head, &nx, &ny, &data, imname);
      if(err) {
#ifdef USECFITSIO
         err = cf_readreal(&head, &nx, &ny, &data, imname);
#else	 
	 fprintf(stderr, "Error from rfitsreal %d reading '%s'\n", err, imname);
	 exit(1);
#endif
      }

   } else {		// Read NITF
      err = nitfread(imname, bitpix, &head, &nx, &ny, &data, TPFN_TEST);
      if(err) {
	 fprintf(stderr, "Error from nitfread %d reading '%s'\n", err, imname);
	 exit(1);
      }
   }

/* Read a weight file if requested */
   if(wgtfile != NULL) {
      if(!rdntf) {		// Read FITS
	 wgt = NULL;
	 err = rfitsreal(&wgthead, &wgtnx, &wgtny, &wgt, wgtfile);
	 if(err) {
#ifdef USECFITSIO
	    err = cf_readreal(&wgthead, &wgtnx, &wgtny, &wgt, wgtfile);
#else	 
	    fprintf(stderr, "Error from rfitsreal %d reading '%s'\n", err, wgtfile);
	    exit(1);
#endif
	 }

      } else {		// Read NITF
	 err = nitfread(wgtfile, bitpix, &wgthead, &wgtnx, &wgtny, &wgt, TPFN_TEST);
	 if(err) {
	    fprintf(stderr, "Error from nitfread %d reading '%s'\n", err, wgtfile);
	    exit(1);
	 }
      }
      if(nx != wgtnx || ny != wgtny) {
	 fprintf(stderr, "Weight file has different dims %d %d than image %d %d\n", wgtnx, wgtny, nx, ny);
	 exit(1);
      }
   }

/* Correct for the bias level */
   for(i=0; i<nx*ny; i++) {
      if(data[i] < 0 || data[i] >= 0) {	/* strip NaNs! */
	 if(data[i] != badata) data[i] -= bias;
      } else {
	 data[i] = badata;
      }
   }

/* Is there a flatfield? */
   if(flatfile != NULL) {
      err = rfitsreal(&headflat, &nxflat, &nyflat, &flat, flatfile);
      if(err) {
	 fprintf(stderr, "Error reading flat file %d\n", err);
	 exit(1);
      }
      if(nxflat != nx || nyflat != ny) {
	 fprintf(stderr, "Error: flat dims %d %d do not match image %d %d\n",
		 nxflat, nyflat, nx, ny);
	 exit(1);
      }
/* Normalize to max if no explicit value given */
      if(flatnorm == 0.0) {
	 for(i=0; i<nx*ny; i++) if(flat[i] > flatnorm) flatnorm = flat[i];
      }
/* Flatten the data */
      for(i=0; i<nx*ny; i++) {
	 if(flat[i] != 0.0) {
	    data[i] /= flat[i] / flatnorm;
	 } else {
	    data[i] = badata;
	 }
      }
   }

/* Bin the image? */
   if(binfactor > 1) {
      for(j=0; j<ny/binfactor; j++) {
	 for(i=0; i<nx/binfactor; i++) {
	    z2 = 0.0;
	    for(jj=0; jj<binfactor; jj++) {
	       for(ii=0; ii<binfactor; ii++) {
		  z2 += data[ii+i*binfactor+(jj+j*binfactor)*nx];
	       }
	    }
	    data[i+j*(nx/binfactor)] = z2;
	 }
      }
      nx /= binfactor;
      ny /= binfactor;
      j = chfitshead(&i, head, "NAXIS1  ", "INTEGER", nx, 0.0);
      j = chfitshead(&i, head, "NAXIS2  ", "INTEGER", ny, 0.0);
   }

   NX = nx;
   if(ifitshead(head, "CNPIX1  ", &sx)) sx = 0;
   if(ifitshead(head, "CNPIX2  ", &sy)) sy = 0;

   if(TPFN_TEST > 0) {
      fprintf(stderr, "Read %s\n  size %d x %d, offset (%d, %d),", 
	     imname, nx, ny, sx, sy);
      fprintf(stderr, "  data %6.1f %6.1f %6.1f\n", data[nx/2-1+ny/2*NX], 
	     data[nx/2+ny/2*NX], data[nx/2+1+ny/2*NX]);
   }

/* Get estimate of sky level and sky rms (over-allocation if mosaicsky) */
   skymed = (double *)calloc((nxsky+1)*(nysky+1), sizeof(double));
   skyrms = (double *)calloc((nxsky+1)*(nysky+1), sizeof(double));
   neadu = 0;
   eadu = 1.0;
   for(j=0; j<=nysky; j++) {
      for(i=0; i<=nxsky; i++) {
	 if(nxsky == 0 && nysky == 0) {		/* Constant sky */
	    imgstat(0, 0, nx-1, ny-1, NX, ny, badata, data,
		    &skymed[0], &skyrms[0], &nskyok);
	 } else {				/* Chunked sky */
	    if(!mosaicsky) {
	       imgstat((int)((i-0.5)*(nx/nxsky)), (int)((j-0.5)*(ny/nysky)), 
		       (int)((i+0.5)*(nx/nxsky)), (int)((j+0.5)*(ny/nysky)),
		       NX, ny, badata, data, 
		       &skymed[i+j*(nxsky+1)], &skyrms[i+j*(nxsky+1)], &nskyok);
	    } else {
	       imgstat((int)((i-0)*(nx/nxsky)), (int)((j-0)*(ny/nysky)), 
		       (int)((i+1)*(nx/nxsky)), (int)((j+1)*(ny/nysky)),
		       NX, ny, badata, data, 
		       &skymed[i+j*(nxsky+1)], &skyrms[i+j*(nxsky+1)], &nskyok);
	    }
	 }

/* Accumulate e/ADU as an average of mean/variance */
	 if(eaduin < 0) {
	    if(skyrms[i+j*(nxsky+1)]>0 && skymed[i+j*(nxsky+1)]>0) {
	       if(neadu == 0) eadu = 0.0;
	       if(skyrms[i+j*(nxsky+1)] > 0) {
		  eadu += skymed[i+j*(nxsky+1)] / 
		     (skyrms[i+j*(nxsky+1)]*skyrms[i+j*(nxsky+1)]);
		  neadu++;
	       }
	    }
	 } else {
	    eadu = eaduin;
	    neadu = 1;
	 }
	 if(TPFN_TEST > 0) {
	    fprintf(stderr, "Sky= %7.1f  RMS= %7.1f  %4d  bias= %7.1f  e/ADU= %7.3f  sat= %7.1f\n", 
		    skymed[i+j*(nxsky+1)], skyrms[i+j*(nxsky+1)], 
		    nskyok, bias, eadu/MAX(1,neadu), sat);
	 }

      }
   }
   if(neadu == 0 && wgtfile == NULL) {
      fprintf(stderr, "Error: tphot sees no RMS in the sky, enter eadu manually\n");
      exit(1);
   } else {
      eadu /= MAX(1,neadu);
   }

/* Invert the sign of the image? */
   if(invert) {
      for(j=0; j<ny; j++) {
	 for(i=0; i<NX; i++) {
	    if(data[i+j*NX] != badata) {
	       sky = skinterp(i, j, nx, ny, nxsky, nysky, skymed, skyrms, &rms);
	       data[i+j*NX] = 2*sky - data[i+j*NX];
	    }
	 }
      }
   }

/* Are we going to write a residual image? */
   if(residfile != NULL) {
      resid = (float *)calloc(NX*ny, sizeof(float));
      for(j=0; j<ny; j++) {
	 for(i=0; i<NX; i++) {
	    resid[i+j*NX] = data[i+j*NX];
	    if(residsky) {
	       sky = skinterp(i, j, nx, ny, nxsky, nysky, skymed, skyrms, &rms);
	       resid[i+j*NX] -= sky;
	    }
	 }
      }
   }

   if(fieldhdr) {
      fprintf(fp, "#  1    xctr    [pix]  %%8.2f  x position of object\n");
      fprintf(fp, "#  2    yctr    [pix]  %%8.2f  y position of object\n");
      fprintf(fp, "#  3    peak    [ADU]  %%7.1f  peak data value on object\n");
      fprintf(fp, "#  4     sky    [ADU]  %%7.1f  sky data value near object\n");
      fprintf(fp, "#  5 peakfit    [ADU]  %%7.1f  peak of fit\n");
      fprintf(fp, "#  6   dpeak    [ADU]  %%7.1f  uncertainty in peak of fit\n");
      fprintf(fp, "#  7  skyfit    [ADU]  %%7.1f  sky background of fit\n");
      fprintf(fp, "#  8    flux    [ADU]  %%9.1f  aperture flux of object\n");
      fprintf(fp, "#  9   dflux    [ADU]  %%8.1f  uncertainty in aperture flux\n");
      if(!trail) {		  
         fprintf(fp, "# 10   major    [pix]  %%6.2f  major axis of object\n");
         fprintf(fp, "# 11   minor    [pix]  %%6.2f  minor axis of object\n");
      } else {
         fprintf(fp, "# 10   trail    [pix]  %%6.2f  trail length of object\n");
         fprintf(fp, "# 11    FWHM    [pix]  %%6.2f  cross-trail width of object\n");
      }
      fprintf(fp, "# 12     phi    [deg]  %%8.3f  CCW angle from x of object\n");
      fprintf(fp, "# 13     err    []     %%2d    fit error code\n");
      fprintf(fp, "# 14    chin    []     %%6.2f  fit chi^2/N\n");

      if(momfile != NULL) {
	 fprintf(fpmom, "#  1    xctr    [pix]  %%8.2f  x position of object\n");
	 fprintf(fpmom, "#  2    yctr    [pix]  %%8.2f  y position of object\n");
	 fprintf(fpmom, "#  3      cx    [pix]  %%8.2f  centroid x position\n");
	 fprintf(fpmom, "#  4      cy    [pix]  %%8.2f  centroid y position\n");
	 fprintf(fpmom, "#  5     xpk    [pix]  %%4d    highest pixel x\n");
	 fprintf(fpmom, "#  6     ypk    [pix]  %%4d    highest pixel y\n");
	 fprintf(fpmom, "#  7   rkron    [pix]  %%7.2f  weighted average radius\n");
	 fprintf(fpmom, "#  8    rvar    [pix]  %%7.2f  weighted RMS radius\n");
	 fprintf(fpmom, "#  9   fkron    [pix]  %%9.1f  flux within 2.5 rkron\n");
	 fprintf(fpmom, "# 10    crad  [pix^2]  %%7.3f  average squared radius wrt highest\n");
	 fprintf(fpmom, "# 11    qrad  [pix^2]  %%7.3f  radial quadrupole\n");
	 fprintf(fpmom, "# 12    trad  [pix^2]  %%7.3f  radial trefoil\n");
	 fprintf(fpmom, "# 13      q2  [pix^2]  %%7.3f  average squared radius wrt cx,cy\n");
	 fprintf(fpmom, "# 14   q2cos  [pix^2]  %%7.3f  + quadrupole\n");
	 fprintf(fpmom, "# 15   q2sin  [pix^2]  %%7.3f  x quadrupole\n");
	 fprintf(fpmom, "# 16   q3cos  [pix^2]  %%7.3f  phi=0 trefoil\n");
	 fprintf(fpmom, "# 17   q3sin  [pix^2]  %%7.3f  phi=30 trefoil\n");
	 fprintf(fpmom, "# 18   q1cos  [pix^2]  %%7.3f  phi=0 monopole\n");
	 fprintf(fpmom, "# 19   q1sin  [pix^2]  %%7.3f  phi=90 monopole\n");
      }
   }

   if(!trail) {
      fprintf(fp, "#   x        y     peakval  skyval  peakfit   dpeak  skyfit     flux     dflux    major  minor    phi  err chi/N\n");
   } else {
      fprintf(fp, "#   x        y     peakval  skyval  peakfit   dpeak  skyfit     flux     dflux    trail   FWHM     phi   err chi/N\n");
   }

   if(momfile != NULL) {
      fprintf(fpmom, "#   x        y        cx       cy     xpk  ypk    rKron    rVar     flux     Crad    Qrad    Trad     Q2     Q2cos   Q2sin   Q3cos   Q3sin   Q1cos   Q1sin\n");
   }

/* Find the objects */
   nfind = nfit = npix = 0;

/* Case 1: external file of triggers */
   if(infp != NULL) {
/* Read an object from the input file */
      while(fgets(line, 10240, infp) != NULL) {
	 if(line[0] == '#') continue;
	 if(nwpar >= 2) {
	    if(nwinit == 2) {
	       if(sscanf(line, "%lf %lf", &wpar.x0, &wpar.y0) != 2) {
		  fprintf(stderr, "Findobj: cannot read x,y from `%s'\n", line);
		  exit(1);
	       }
	       wpar.major = major;
	       wpar.minor = minor;
	       wpar.phi = phi;
	    } else if(nwinit == 5) {
	       if(sscanf(line, "%lf %lf %lf %lf %lf", &wpar.x0, &wpar.y0, 
			 &wpar.major, &wpar.minor, &wpar.phi) != 5) {
		  fprintf(stderr, "Findobj: cannot read x,y,maj,min,phi from `%s'\n", line);
		  exit(1);
	       }
	       wpar.phi *= pi/180;	// convert to rad
	    }
	 } else if(nwpar == 0) {	// No fit: strictly difference image
	    if(sscanf(line, "%lf %lf %lf", &wpar.x0, &wpar.y0, &wpar.peak) != 3) {
	       fprintf(stderr, "Findobj: cannot read x,y,peak from `%s'\n", line);
	       exit(1);
	    }
/* Fill in the rest of wpar since we're not going to fit anything */
	    wpar.beta4 = 1.0;
	    wpar.beta6 = 0.5;
	    if(!trail) {
	       err = dosxy(&wpar);
	    } else {
	       wpar.sxy = phi;
	       wpar.sx2 = pow(2.3548/minor, 2.0);
	       wpar.sy2 = sqrt(major*major+minor*minor) - minor;
	       if(wpar.sy2 > 0) wpar.sy2 = log(0.5*wpar.sy2);
	    }
	 }
	 i = wpar.x0 + 0.5;
	 j = wpar.y0 + 0.5;

	 wpar.ninit = nwinit;

	 if(nwpar > 0) {	// Fit at all or just difference?
	    nfit++;
	    err = one_object(nx, NX, ny, data, wgt, i, j,
			     invert, trail, force, okfit,  &wpar, &wflux);
	    if(err == 0) {
	       nfind++;
	       one_report(nx, NX, ny, data, i, j, trail, &wpar, &wflux, fp, fpmom);
	    } else {
	       continue;	// Don't subtract bad fit
	    }
	 }

/* Subtract from the residual image? */
	 if(residfile != NULL || subtract) {
	    one_subtract(nx, NX, ny, data, resid, i, j, trail, subtract, &wpar);
	 }
      } /* object trigger */
   } /* Input from file */



/* Case 2: search the image for triggers */
   if(infp == NULL) {

/* Two passes with disjoint stripes in y prevents multiple threads from */
/*  interfering if this is multithreaded; otherwise harmless. */
      stripwidth = ny / (2*nthread);
      nfind = nfit = 0;

      if(TPFN_TEST > 0) {
	 printf("Number of threads is %d of maximum %d width %d\n", nthread, maxthread, stripwidth);
      }

      for(stripass=0; stripass<=1; stripass++) {

	 for(thread=0; thread<nthread; thread++) {

	    ymin = (2*thread+stripass) * stripwidth;
	    ymax = ymin + stripwidth;

	    if(TPFN_TEST > 0) {
	       printf("thread %4d : %4d <= y < %d\n", thread, ymin, ymax);
	    }

/* Find the objects: npix is the starting address in the image */
	    npix = nx * ymin;
	    while(1) {

/* Find a pixel at address npix->i,j that passes the trigger criteria */
	       one_trigger(nx, NX, ny, data, &npix, ymax, skymed, skyrms);
	       if(npix >= nx*ymax) break;
	       i = npix % nx;
	       j = npix / nx;

/* OK, got one... analyze it */
	       wpar.x0 = i;
	       wpar.y0 = j;

	       if(TPFN_TEST > 1) {
		  printf("\nNew object at x= %d y= %d  peak= %.1f\n", 
			 i, j, data[i+j*NX]);
	       }

	       wpar.major = major;
	       wpar.minor = minor;
	       wpar.phi = phi;
	       wpar.ninit = nwinit;

	       if(nwpar > 0) {	// Fit at all or just difference?
		  nfit++;
		  err = one_object(nx, NX, ny, data, wgt, i, j,
				   invert, trail, force, okfit,  &wpar, &wflux);
		  if(err == 0) {
		     nfind++;
		     one_report(nx, NX, ny, data, i, j, trail,
				&wpar, &wflux, fp, fpmom);
		  } else {
		     continue;	// Don't subtract bad fit
		  }
	       }

/* Subtract from the residual image? */
	       if(residfile != NULL || subtract) {
		  one_subtract(nx, NX, ny, data, resid, i, j, trail, subtract, &wpar);
	       }
	    } /* object trigger and analysis */

	 } /* thread loop */

      } /* stripass=0,1 */

   } /* trigger search */

   if(TPFN_TEST > 0) printf("Nfit= %d  Nfind= %d\n", nfit, nfind);

/* Write a residual image */
   if(residfile != NULL) wfitsreal(head, resid, residfile);

   exit(0);
}

/* Step forward through image seeking pixels that pass the trigger */
void one_trigger(int nx, int NX, int ny, float *data,
		 int *npix, int ymax, double *skymed, double *skyrms)
{			
   int i, j, ii, jj, nequal, idxeq;
   double sky, rms;

/* Walk forward through the image */
   while((*npix/nx) < ymax) {
      *npix += 1;
      j = *npix / nx;
      i = *npix % nx;

/* Skip pixels too close to the borders */
      if(j < y0border || j > ny-1-y1border) continue;
      if(i < x0border || i > nx-1-x1border) continue;

/* Not bad data, correct? */
      if(data[i+j*NX] == badata) continue;

/* Sky level and RMS */
      if(skymed != NULL && skyrms != NULL) {
	 sky = skinterp(i, j, nx, ny, nxsky, nysky, skymed, skyrms, &rms);
      } else {
/* Or just nominal, if data are already units of sigma */
	 sky = 0.0;
	 rms = 1.0;
      }

/* Pass the cut?  Weirdness with double negatives is in case of NaN */
      if(!(data[i+j*NX] > sky+srchcut && data[i+j*NX] < sky+srchmax)) continue;
      if(!(data[i+j*NX] > sky + srchsig*rms)) continue;

/* Maximum among closest neighbors? */
      if(!(data[i+j*NX] >  data[i+1+j*NX]))   continue;
      if(!(data[i+j*NX] >= data[i-1+j*NX]))   continue;
      if(!(data[i+j*NX] >  data[i+(j+1)*NX])) continue;
      if(!(data[i+j*NX] >= data[i+(j-1)*NX])) continue;

/* Maximum among all 8 neighbors? */
      if(srchrad > 1) {
	 if(!(data[i+j*NX] >  data[i+1+(j+1)*NX])) continue;
	 if(!(data[i+j*NX] >= data[i-1+(j+1)*NX])) continue;
	 if(!(data[i+j*NX] >  data[i+1+(j-1)*NX])) continue;
	 if(!(data[i+j*NX] >= data[i-1+(j-1)*NX])) continue;
      }
	 
/* Bigger search area?  One tie allowed. */
      if(srchrad > sqrt(2.0)) {
	 nequal = 0;
	 idxeq = 0;
	 for(jj=-NINT(srchrad); jj<=NINT(srchrad); jj++) {
	    if(j+jj < 0 || j+jj >= ny) continue;
	    for(ii=-NINT(srchrad); ii<=NINT(srchrad); ii++) {
	       if(i+ii < 0 || i+ii >= nx) continue;
	       if(ii*ii+jj*jj > srchrad*srchrad) continue;
	       if(!(data[i+ii+(j+jj)*NX] != badata)) continue;
	       if(!(data[i+j*NX] >= data[i+ii+(j+jj)*NX])) {
		  nequal = 99;
		  break;
	       }
	       if(!(data[i+j*NX] != data[i+ii+(j+jj)*NX])) {
		  nequal++;
		  idxeq = ii<0 || (ii==0 && jj<0);
	       }
	    }
	    if(nequal > 2) break;
	 }

/* Local maximum?  Tie on the RHS? */
	 if(nequal > 2 || idxeq) continue;
      }

/* OK, got one... analyze it */
      return;
   }
   return;
}


/* Assemble bits that describe various failure modes */
int failbits(int nx, int ny, int i, int j, int trail,
	     PSF2D_PARAM *wpar, PSF2D_FLUX *wflux)
{
   int failure;
   /* Does this pass the remainder of the tests? */
   failure = 0;
   if(wpar->x0<x0border || wpar->x0>nx-1-x1border) failure |= FAIL_XBORDER;
   if(wpar->y0<y0border || wpar->y0>ny-1-y1border) failure |= FAIL_YBORDER;
   if(ABS(wpar->x0-(i+0.5)) > maxmove) failure |= FAIL_XMOVE;
   if(ABS(wpar->y0-(j+0.5)) > maxmove) failure |= FAIL_YMOVE;
   if(wpar->major > fwmax) failure |= FAIL_MAJBIG;
   if(wpar->minor > fwmax) failure |= FAIL_MINBIG;
   if(wpar->major < fwmin && !trail) failure |= FAIL_MAJSMALL;
   if(wpar->minor < fwmin) failure |= FAIL_MINSMALL;
   if(wpar->peak < netmin) failure |= FAIL_PEAK;
   if(wpar->chin > chinmax) failure |= FAIL_CHIN;
   if(wpar->chin == 0.0) failure |= FAIL_CHIN;
   if(fitsig > 0) {
      if(wflux->flux <= 0) failure |= FAIL_FLUXNEG;
      if(wflux->dflux <= 0) failure |= FAIL_FERRNEG;
      if(wflux->dflux > 0 && wflux->flux/wflux->dflux < fitsig) failure |= FAIL_SNR;
   }
   return(failure);
}

/* Do the work on just one object */
int one_object(int nx, int NX, int ny, float *data, float *wgt,
	       int i, int j, int invert, int trail, int force, int okfit,
	       PSF2D_PARAM *wpar, PSF2D_FLUX *wflux)
{
   int failure, err;
   double pi=4*atan(1.0);

   if(!trail) {
      err = wpsf(nx, ny, data, wgt, eadu, sat, badata, aprad, skyrad,
		 nwpar, wpar, wflux);
   } else {
      err = trailpsf(nx, ny, data, wgt, eadu, sat, badata, aprad, skyrad,
		     nwpar, wpar, wflux);
   }

// err = 10 * (number of bad pixels)

/* Does this pass the remainder of the tests? */
   failure = failbits(nx, ny, i, j, trail, wpar, wflux);

#if 0
   if(wpar->x0<x0border || wpar->x0>nx-1-x1border) failure |= FAIL_XBORDER;
   if(wpar->y0<y0border || wpar->y0>ny-1-y1border) failure |= FAIL_YBORDER;
   if(ABS(wpar->x0-(i+0.5)) > maxmove) failure |= FAIL_XMOVE;
   if(ABS(wpar->y0-(j+0.5)) > maxmove) failure |= FAIL_YMOVE;
   if(wpar->major > fwmax) failure |= FAIL_MAJBIG;
   if(wpar->minor > fwmax) failure |= FAIL_MINBIG;
   if(wpar->major < fwmin && !trail) failure |= FAIL_MAJSMALL;
   if(wpar->minor < fwmin) failure |= FAIL_MINSMALL;
   if(wpar->peak < netmin) failure |= FAIL_PEAK;
   if(wpar->chin > chinmax) failure |= FAIL_CHIN;
   if(wpar->chin == 0.0) failure |= FAIL_CHIN;
   if(fitsig > 0) {
      if(wflux->flux <= 0) failure |= FAIL_FLUXNEG;
      if(wflux->dflux <= 0) failure |= FAIL_FERRNEG;
      if(wflux->dflux > 0 && wflux->flux/wflux->dflux < fitsig) failure |= FAIL_SNR;
   }
#endif

/* Debug info for all triggers for objects */
   if(TPFN_TEST > 1) {
      printf("PSF: err= %d x0= %.2f y0= %.2f peak= %.1f sky= %.1f\n",
	     err, wpar->x0, wpar->y0, wpar->peak, wpar->sky);
      printf("  maj= %.2f min= %.2f phi= %.3f chin= %.2f rms= %.2f resid= %.2f %.2f\n",
	     wpar->major, wpar->minor, wpar->phi*180/pi,
	     wpar->chin, wpar->rms, wpar->absresid, wpar->maxresid);
      printf("Uncertainties: dx0= %.2f dy0= %.2f dpeak= %.1f dsky= %.1f\n",
	     wpar->dx0, wpar->dy0, wpar->dpeak, wpar->dsky);
      printf("  dsx2= %.2f dsxy= %.2f dsy2= %.3f\n",
	     wpar->dsx2, wpar->dsxy, wpar->dsy2);
      printf("Apflux: flux= %.1f dflux= %.1f sky=%.1f dsky= %.1f skyrms= %.1f\n",
	     wflux->flux, wflux->dflux, wflux->sky, wflux->dsky, 
	     wflux->skyrms);
      if(failure) {
	 printf("FitFailCode= 0x%04x  ", failure);
	 failcode(failure);
      }
   }

// FIXME: this just doesn't make any sense, and it confuses two searches:
// the initial search in the sigma image which we want to be quite restrictive
// and the final trail fit that we want to be really inclusive.

// err = 10*nignore_pix so this fails objects that happen to lie near bad pixels.
// Bad for initial search in sigma image, good for trail fit
//   if((err%10) && okfit) failure |= FAIL_OKFIT;

// Good for initial search in sigma image, obnoxious bit for trail fit
// I wonder whether a bad pixel nearby should be a failure for initial search...
   if(err && okfit) failure |= FAIL_OKFIT;

// Don't really care about erasing this as long as failure&&!force is out...
/* Allow chin=0 (bad fit) if okfit=0 */
   if(wpar->chin == 0.0 && okfit == 0) {
      failure &= ~FAIL_CHIN;
      err = -1;
   }

/* If invert, flip the sign of the data */
   if(invert) {
      wpar->peak *= -1;
      wflux->flux *= -1;
   }

// I'd like to get the failure code back.  The !force is just a hack to
// get a successful return even though some bits fail.
//   if(failure && !force) return(1);

   return(failure);
}

/* Report just one object */
void one_report(int nx, int NX, int ny, float *data,
		int i, int j, int trail, PSF2D_PARAM *wpar, PSF2D_FLUX *wflux,
		FILE *fp, FILE *fpmom)
{
   int err=0;
   PSF2D_MOM mom;
   double val, pi=4*atan(1.0);

   if(i>=0 && i<=nx && j>=0 && j<=ny) {
      val = data[i+j*NX];
   } else {
      val = 0.0;
   }

/* Tell us about a real object */
   if(!trail) {
      fprintf(fp, "%8.2f %8.2f  %7.1f %7.1f  %7.1f %7.1f %7.1f  %9.1f %8.1f  %6.2f %6.2f %6.1f  %2d %6.2f\n",
	      wpar->x0, wpar->y0, val, wflux->sky, 
	      wpar->peak, wpar->dpeak, wpar->sky, wflux->flux, wflux->dflux, 
	      wpar->major, wpar->minor, wpar->phi*180/pi, err, wpar->chin);
   } else {
      fprintf(fp, "%8.2f %8.2f  %7.1f %7.1f  %7.1f %7.1f %7.1f  %9.1f %8.1f  %6.2f %6.2f %8.3f  %2d %6.2f\n",
	      wpar->x0, wpar->y0, val, wflux->sky, 
	      wpar->peak, wpar->dpeak, wpar->sky, wflux->flux, wflux->dflux, 
	      wpar->major, wpar->minor, wpar->phi*180/pi, err, wpar->chin);
   }


   if(fpmom != NULL) {
/* Compute moments of this object */
      err = mompsf(nx, ny, data, badata, aprad, 
		   wpar->x0, wpar->y0, wflux->sky, &mom);
/* Write the moments of this object */
      fprintf(fpmom, "%8.2f %8.2f  %8.2f %8.2f %4d %4d  %7.2f %7.2f %9.1f  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	      wpar->x0, wpar->y0, mom.cx, mom.cy, mom.xpk, mom.ypk,
	      mom.rkron, mom.rvar, mom.flux, mom.crad, mom.qrad, mom.trad,
	      mom.q2, mom.q2cos, mom.q2sin, mom.q3cos, mom.q3sin,
	      mom.q1cos, mom.q1sin);
   }
   return;
}

/* Subtract the fit for just one object */
void one_subtract(int nx, int NX, int ny, float *data, float *resid,
		  int i, int j, int trail, int subtract, PSF2D_PARAM *wpar)
{
   int ii, jj, rsub;
   double z2;

/* Subtract it from the image or residual image */
   if(resid != NULL || subtract) {
      rsub = sqrt(wpar->major*wpar->major+wpar->minor*wpar->minor) * RESIDRAD + 0.5;
      if(TPFN_TEST > 1) {
	 printf(" Subtract box x= %d %d y= %d %d\n",
		i-rsub,i+rsub, j-rsub, j+rsub);
      }
      for(jj=j-rsub; jj<=j+rsub; jj++) {
	 if(jj < 0 || jj > ny-1) continue;
	 for(ii=i-rsub; ii<=i+rsub; ii++) {
	    if(ii < 0 || ii > nx-1) continue;
	    if(resid != NULL) {
	       if(!trail) {
		  resid[ii+jj*NX] -= wauss(ii+0.5, jj+0.5, wpar, &z2);
	       } else {
		  resid[ii+jj*NX] -= wtrail(ii+0.5, jj+0.5, wpar, &z2);
	       }
	    }
	    if(subtract) {
	       if(!trail) {
		  data[ii+jj*NX] -= wauss(ii+0.5, jj+0.5, wpar, &z2);
	       } else {
		  data[ii+jj*NX] -= wtrail(ii+0.5, jj+0.5, wpar, &z2);
	       }
	    }
	 }
      }
   }
   return;
}

/* Compute median and median RMS for an array */
int imgstat(int x0, int y0, int x1, int y1, int NX, int NY, float badata,
	    float *data, double *med, double *rms, int *nok)
{
   int i, j, k, m, n, nsrt, err;
   double buf[MAXSRT];
   *nok = 0;
   nsrt = MIN((x1-x0+1)*(y1-y0+1), MAXSRT);

//   fprintf(stderr, "%d %d %d %d %d\n", x0, x1, y0, y1, nsrt);
   for(m=n=0; m<nsrt; m++) {
      k = (m/(double)nsrt) * (x1-x0+1)*(y1-y0+1);
      i = (k%(x1-x0+1)) + x0;
      j = (k/(x1-x0+1)) + y0;
      if(i<0 || i>=NX || j<0 || j>=NY) continue;
      if(data[i+j*NX] == badata) continue;
      buf[n++] = data[i+j*NX];
   }
   if(n < 1) {
      *med = *rms = 0.0;
      return(-1);
   }
   if(n < 2) {
      *med = buf[0];
      *rms = 0;
      return(0);
   }
   if( (err = qsort8(n, buf)) ) return(err);
   *med = 0.5 * (buf[(n-1)/2] + buf[n/2]);
   *rms = *med - buf[n/4];
   *nok = n;
   return(0);
}

/* Evaluate a 4 point linear interpolation of sky at image pixel i,j */
double skinterp(int i, int j, int nx, int ny,
	 int nxsky, int nysky, double *skymed, double *skyrms, double *rms)
{
   int is, js;
   double u, v, sky;
   if(nxsky == 0 && nysky == 0) {
      sky = skymed[0];
      *rms = skyrms[0];
   } else {
      is = i / (nx/nxsky);
      js = j / (ny/nysky);
      if(!mosaicsky) {
/* Edge squares are calculated between 0:0.5, not 0:1 => extrapolation */
	 if(is == 0) {
	    u = (((double)i) / (nx/nxsky) - is - 0.25) * (4.0/3.0);
	 } else if(is == nxsky-1) {
	    u = (((double)i) / (nx/nxsky) - is) * (4.0/3.0);
	 } else {
	    u = ((double)i) / (nx/nxsky) - is;
	 }
	 if(js == 0) {
	    v = (((double)j) / (ny/nysky) - js - 0.25) * (4.0/3.0);
	 } else if(js == nysky-1) {
	    v = (((double)j) / (ny/nysky) - js) * (4.0/3.0);
	 } else {
	    v = ((double)j) / (ny/nysky) - js;
	 }
	 sky = (1-u) * (1-v) * skymed[is+js*(nxsky+1)] +
	    u * (1-v) * skymed[is+1+js*(nxsky+1)] +
	    (1-u) * v * skymed[is+(js+1)*(nxsky+1)] +
	    u * v * skymed[is+1+(js+1)*(nxsky+1)];
	 *rms = (1-u) * (1-v) * skyrms[is+js*(nxsky+1)] +
	    u * (1-v) * skyrms[is+1+js*(nxsky+1)] +
	    (1-u) * v * skyrms[is+(js+1)*(nxsky+1)] +
	    u * v * skyrms[is+1+(js+1)*(nxsky+1)];
      } else {
	 sky = skymed[is+js*(nxsky+1)];
	 *rms = skyrms[is+js*(nxsky+1)];
      }
   }
   return(sky);
}

/* Print reasons for failure */
int failcode(int failure)
{
   if(failure & FAIL_XBORDER) printf(" XBORDER");
   if(failure & FAIL_YBORDER) printf(" YBORDER");
   if(failure & FAIL_XMOVE) printf(" XMOVE");
   if(failure & FAIL_YMOVE) printf(" YMOVE");
   if(failure & FAIL_MAJBIG) printf(" MAJBIG");
   if(failure & FAIL_MINBIG) printf(" MINBIG");
   if(failure & FAIL_MAJSMALL) printf(" MAJSMALL");
   if(failure & FAIL_MINSMALL) printf(" MINSMALL");
   if(failure & FAIL_PEAK) printf(" PEAK");
   if(failure & FAIL_CHIN) printf(" CHIN");
   if(failure & FAIL_FLUXNEG) printf(" FLUXNEG");
   if(failure & FAIL_FERRNEG) printf(" FERRNEG");
   if(failure & FAIL_SNR) printf(" SNR");
   if(failure & FAIL_OKFIT) printf(" OKFIT");
   if(failure) printf("\n");
   return(0);
}

int syntax(char *prog)
{
/* Defaults...  Keep consistent with above. */
   int binfactor = 1;		// Bin the input data?
   int okfit = 1;		// Does the fit have to be OK to keep?
   int nwinit = 0;		// Initialization params to read if obj input

   printf("%s image [options]\n\n", prog);
   printf("Finds objects in an image, reports on properties\n");
   printf("[options] include:\n\n");
   printf("  -out fname   Output file name (stdout)\n");
   printf("  -obj fname   Use x,y from fname for triggers (NULL)\n");
   printf("  -resid fname Write a residual FITS image (NULL)\n");
   printf("  -moment fname Write moments for each output object (NULL)\n");
   printf("  -subtract    Subtract fits successively from image (default no)\n");
   printf("  -residsky    Subtract sky model from residual image (default no)\n");
   printf("  -trail       Fit a trailed PSF (default no)\n");
   printf("  -verb        Diagnostic verbosity\n");
   printf("  -VERB        More diagnostic verbosity\n");
   printf("  -hdr         Write details about all fields\n");
   printf("  -bin N       Bin input image by a factor NxN (%d)\n", binfactor);
   printf("  -invert      Invert the sign of the image? (no)\n\n");
   printf("Trigger parameters:\n");
   printf("  -border N    Image border to avoid on all sides (%d)\n", x0border);
   printf("  -xlft N      Left x image border to avoid (%d)\n", x0border);
   printf("  -xrgt N      Right x image border to avoid (%d)\n", x1border);
   printf("  -ybot N      Bottom y image border to avoid (%d)\n", y0border);
   printf("  -ytop N      Top y image border to avoid (%d)\n", y1border);
   printf("  -badata B    Ignore data values of B (%.1f)\n", badata);
   printf("  -rad R       Local max required over radius r to trigger (%.1f)\n", srchrad);
   printf("  -sig S       Minimum value = sky+sig*rms to trigger (%.1f)\n", srchsig);
   printf("  -min Z       Minimum value above sky to trigger (%.1f)\n", srchcut);
   printf("  -max Z       Maximum peak to trigger (%.1f)\n\n", srchmax);
   printf("Fit acceptance parameters:\n");
   printf("  -fwmax W     Maximum FW to keep (%.2f)\n", fwmax);
   printf("  -fwmin W     Minimum FW to keep (%.2f)\n", fwmin);
   printf("  -net Z       Minimum peak-sky to keep (%.1f)\n", netmin);
   printf("  -move M      Max offset between peak and fit (%.1f)\n", maxmove);
   printf("  -chin C      Maximum chi/N to keep (%.2f)\n", chinmax);
   printf("  -snr S       Minimum flux SNR to keep (%.1f)\n", fitsig);
   printf("  -okfit I     Insist (0/1) on successful profile fit (%d) \n", okfit);
   printf("  -force       Output even if some tests failed (default not)\n\n");
   printf("The fit is normally 7 parameters (position, peak, sky, PSF\n");
   printf("but if major, minor, and phi are set a 4 param fit is used:\n");
   printf("  -major M     Set major axis M for 4 param fit [pix]\n");
   printf("  -minor M     Set minor axis M for 4 param fit [pix]\n");
   printf("  -phi P       Set major axis angle to P for 4 param fit [deg]\n");
   printf("  -npar N      How many parameter fit (2/4/7)? (%d)\n", nwpar);
   printf("  -ninit N     Does trigger file have x,y or x,y,maj,min,phi (2/5)? (%d)\n\n", nwinit);
   printf("Noise, signal, and aperture parameters:\n");
   printf("  -eadu E      Electrons per ADU to estimate noise (%.1f)\n", eadu);
   printf("               if E is 'auto' eadu is set to mean/variance\n");
   printf("  -bias B      Image background exceeds sky by B (%.1f)\n", bias);
   printf("  -wgt fname   Weight (inv variance) file name (NULL)\n");
   printf("  -sat S       Saturation level is S (%.1f)\n", sat);
   printf("  -aprad R     Aperture radius for photometry (%d)\n", aprad);
   printf("  -skyrad R    Sky radius for photometry (%d)\n", skyrad);
   printf("  -nxsky M     Subdivide image by M in x for sky (%d)\n", nxsky);
   printf("  -nysky N     Subdivide image by N in y for sky (%d)\n", nysky);
   printf("  -mosiacsky   Use an MxN chunk discontinuous sky model\n");
   printf("  -flat F      Read flat file and divide by it\n");
   printf("  -flatnorm N  Divide flat file by N before application\n");
//   printf("\n  -nitf        Input is NITF, not FITS\n");
   return(0);
}
