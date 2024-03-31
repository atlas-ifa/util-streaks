/* Convolve an image with a tilted streak kernel */
/* E.g.
  strkonv /atlas/diff/02a/58071/02a58071o0503c.diff.fz /tmp/foo.d  -phi 159 -len 16 -nlam 3
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "fitsio.h"

#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))
#define ABS(x)  ((x)<0?(-(x)):(x))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define LN10 (2.302585092994)	/* ln(10) */

#define SHOWTIME		/* define to show time as execution proceeds */

/* Timing diagnostics */
#ifdef SHOWTIME
#define BLAB(msg) { gettimeofday(&tv1, NULL); telapse = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec); fprintf(stderr, "%7.3f - %s\n", telapse, msg); }
#else
#define BLAB(msg)
#endif

static int VERBOSE=0;		/* verbosity level */

// Macro from robospect to do all the boring fits error handling
#define FRE { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); } ; status = 0; }

#define FFATAL { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); exit(1); }  }

/* Convolve image with 1/cos/sin with nlam wavelength at phi and length +/-len */
void xdrift(double phi, int len, double nlam, 
	    int nx, int ny, float *H, float *C, float *S);
void ydrift(double phi, int len, double nlam,
	    int nx, int ny, float *H, float *C, float *S);

/* Convolve a 1D array by a Gaussian with sigma: 14nsec/pix */
void gconv(int n, float *a, double sig);

/* Fast median. Does not alter input array. */
float quickmedian(int n, float *key);

/* Linearly interpolate src(nx/xbin,ny/ybin) -> dest(nx,ny) in ncx,ncy chunks */
/* ncx,ncy must divide the src and dest dimensions evenly */
void linterp_chunk(int nx, int ny, int xbin, int ybin, int ncx, int ncy,
		   float *src, float *dest);
/* Median bin an image, ignore pixels at badata */
void medianbin(int nx, int ny, float *D, float badata, 
	       int bx, int by, float **B);


/* 2D smooth an image by sig */
void smooth(double sig, int nx, int ny, float *D);

/* Clip and image and bin */
void clipbin(double clipmin, double clipmax, int bin, float badata, int nx, int ny, float *D);

/* Print out cfitsio error messages and exit program */
void printerror(int status);

/* Append a FITS file as an extension at the end of another one */
int main(int argc, char **argv)
{
   int i, bitpix, status=0, naxis, anynull;
   int len, dcs, nx, ny, bx, by, chx, chy, bin, bkgsub;
   char *infile, *newfile, *newbang;
   fitsfile *fnew, *ff;   /* pointers to the FITS files */
   long int naxes[2], fpixel[2];
//   size_t nbyte;
   LONGLONG nelements;
   double phi, nlam, pi=4*atan(1.0);
   float *D, *C, *S, *B, badata;
   double clipmin, clipmax, sig;

#ifdef SHOWTIME
   struct timeval tv0, tv1;
   double telapse;
   gettimeofday(&tv0, NULL);
#endif


   if(argc < 3) {
      printf("strkonv: insufficient arguments\n");
      printf("syntax: strkonv infile outfile [options]\n");
      exit(1);
   }

   infile = argv[1];
   newfile = argv[2];
   badata = 0.0;	/* Bad data value */

/* Parse the options */
   bin = 2;		/* Bin input image by bin */
   clipmin = -1000;	/* Clip output, post-bin */
   clipmax = 1000;	/* Clip output, post-bin */
   bx = 110;		/* Median bin by bx for background, post-bin */
   by = 110;		/* Median bin by by for background, post-bin */
   chx = 8;		/* Respect chx x chy chunks for linterp */
   chy = 8;		/* Respect chx x chy chunks for linterp */
   sig = 1.0;		/* Smooth by sig, post-bin */

   phi = 30;		/* Streak direction wrt x [deg] */
   len = 16;		/* Streak length, post-bin */
   nlam = 3;		/* Number of cos/sin wavelengths in (2*len+1) */
   dcs = 0;		/* dcs=0,1,2 for tophat,cos,sin */
   bkgsub = 1;		/* Subtract the background? */

   for(i=3; i<argc; i++) {
      if(strcmp(argv[i], "-phi") == 0) {
	 sscanf(argv[++i], "%lf", &phi);
	 if(phi < -180) phi = phi + 360;
	 if(phi >  180) phi = phi - 360;
	 if(phi < -45) phi = phi + 180;
	 if(phi > 135) phi = phi - 180;
	 phi *= pi / 180;	// -pi/4 < phi < 3*pi/4

      } else if(strcmp(argv[i], "-len") == 0) {
	 sscanf(argv[++i], "%d", &len);

      } else if(strcmp(argv[i], "-nlam") == 0) {
	 sscanf(argv[++i], "%lf", &nlam);

      } else if(strcmp(argv[i], "-bin") == 0) {
	 sscanf(argv[++i], "%d", &bin);

      } else if(strcmp(argv[i], "-sig") == 0) {
	 sscanf(argv[++i], "%lf", &sig);

      } else if(strcmp(argv[i], "-clip") == 0) {
	 sscanf(argv[++i], "%lf,%lf", &clipmin, &clipmax);

      } else if(strcmp(argv[i], "-dcs") == 0) {
	 sscanf(argv[++i], "%d", &dcs);

      } else if(strcmp(argv[i], "-bkg") == 0) {
	 bkgsub = 0;

      } else if(strcmp(argv[i], "-verb") == 0) {
	 VERBOSE = 1;

      } else if(strcmp(argv[i], "-VERB") == 0) {
	 VERBOSE = 2;

      } else {
	 printf("Unrecognized argument `%s'\n", argv[i]);
      }
   }

   BLAB("Args parsed");

/////////////////////////////
// READ WRITE BOILER PLATE //
/////////////////////////////
/* Read the original file: use open_image vs open_file in case compressed */
   status = 0;
   if ( fits_open_image(&ff, infile, READONLY, &status) ) FFATAL;

/* What about it? */
   if( fits_get_img_param(ff, 2, &bitpix, &naxis, naxes, &status) ) FFATAL;

   if(VERBOSE > 1) {
      printf("Input file has %d axes %ld x %ld and bitpix %d\n", 
	     naxis, naxes[0], naxes[1], bitpix);
   }

   nx = naxes[0];
   ny = naxes[1];

   fpixel[0] = fpixel[1] = 1;                /* first pixel to read */
   nelements = naxes[0] * naxes[1];          /* number of pixels to read */

/* Allocate space to read data */
   D = (float *)calloc(naxes[0]*naxes[1], sizeof(float));

   if( fits_read_pix(ff, TFLOAT, fpixel, nelements, NULL,
		     D, &anynull, &status) ) {
      if(status != NUM_OVERFLOW) FFATAL;
      if(VERBOSE > 1) FRE;
      status = 0;
   }
   BLAB("Input file read");


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
   naxes[0] = nx / bin;
   naxes[1] = ny / bin;
   if( fits_create_img(fnew, bitpix, naxis, naxes, &status) ) FFATAL;
   BLAB("Output image created");

/* Close the input image */
   if( fits_close_file(ff, &status) ) FRE;
   BLAB("Header copied");
////////////////////////////////////////////////////////////////





/* Clip and bin as requested */
   clipbin(clipmin, clipmax, bin, badata, nx, ny, D);
   nx /= bin;
   ny /= bin;
   BLAB("Clip and bin");

/* Allocate space for cosine, and sine arrays (tophat is in D) */
   C = (float *)calloc(naxes[0]*naxes[1], sizeof(float));
   S = (float *)calloc(naxes[0]*naxes[1], sizeof(float));

   if(bkgsub) {
/* Make a median-bin version for background */
      medianbin(nx, ny, D, badata, bx, by, &B);
      BLAB("Median bin");

#if 0
      int j;
      for(j=0; j<ny/by; j++) {
	 for(i=0; i<10; i++) printf(" %7.1f", B[i+j*(nx/bx)]);
	 printf("\n");
      }
#endif

/* Linearly interpolate back to full size in C array, respecting 8x8 chunks */
      linterp_chunk(nx, ny, bx, by, chx, chy, B, C);
      BLAB("Linterp");

/* Subtract the background */
      for(i=0; i<nx*ny; i++) D[i] -= C[i];
      BLAB("Background subtracted");
   }

/* Smooth entire image by sig=1 */
   smooth(sig, nx, ny, D);
   BLAB("Smoothed");


/* Convolve the image with tophat 1,c,s */
   if(ABS(phi) < pi/4) {
      xdrift(phi, len, nlam, nx, ny, D, C, S);
   } else {
      ydrift(phi, len, nlam, nx, ny, D, C, S);
   }
   BLAB("Drifted");

/* Make the desired combination: D-sqrt(C*C+S*S) */






/////////////////////////////
// READ WRITE BOILER PLATE //
/////////////////////////////
/* Write the output array */
   fpixel[0] = fpixel[1] = 1;                /* first pixel to write      */
   nelements = naxes[0] * naxes[1];          /* number of pixels to write */

/* write the array of unsigned integers to the FITS file */
   if(dcs == 0) {
      if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, D, &status) ) FRE;
   } else if(dcs == 1) {
      if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, C, &status) ) FRE;
   } else {
      if( fits_write_pix(fnew, TFLOAT, fpixel, nelements, S, &status) ) FRE;
   }


   BLAB("Pixels written");

/* Finally close the output image */
   if( fits_close_file(fnew, &status) ) FRE;
   BLAB("Done");
////////////////////////////////////////////////////////////////


   return(0);
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


/* Clip and image and bin */
void clipbin(double clipmin, double clipmax, int bin, float badata, int nx, int ny, float *D)
{
   int i, j, k, l, n;
   double sum;
   for(i=0; i<nx*ny; i++) {
      if(D[i] < clipmin) D[i] = clipmin;
      if(D[i] > clipmax) D[i] = clipmax;
   }
   if(bin > 1) {
      for(j=0; j<ny/bin; j++) {
	 for(i=0; i<nx/bin; i++) {
	    sum = 0.0;
	    n = 0;
	    for(l=j*bin; l<(j+1)*bin; l++) {
	       for(k=i*bin; k<(i+1)*bin; k++) {
		  if(D[k+l*nx] != badata) {
		     sum += D[k+l*nx];
		     n++;
		  }
	       }
	    }
	    D[i+j*(nx/bin)] = sum / MAX(1,n);
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



/* Median bin D(nx,ny) by bx,by, return B, ignore pixels at badata */
void medianbin(int nx, int ny, float *D, float badata, int bx, int by, float **B)
{
   int i, j, k, l, n, nbx=(nx+bx-1)/bx, nby=(ny+by-1)/by;
   float *buf;

   buf = calloc(bx*by, sizeof(float));
   *B = calloc(nbx*nby, sizeof(float));

   for(j=0; j<nby; j++) {
      for(i=0; i<nbx; i++) {
	 n = 0;
	 for(l=j*by; l<MIN(ny-1,(j+1))*by; l++) {
	    for(k=i*bx; k<MIN(nx-1,(i+1))*bx; k++) {
	       if(D[k+l*nx] != badata) buf[n++] = D[k+l*nx];
	    }
	 }
	 (*B)[i+j*nby] = quickmedian(n, buf);
      }
   }

   free(buf);

   return;
}


/* y value for at 0-based column i from start row j, slope dy */
#define SLOPE(i,j,slope,n,k) { k = (int)((i+0.5)*slope+j) ; \
				       if(k<0) k = 0; if(k>n-1) k = n-1; }

/* Convolve image with 1/cos/sin with nlam wavelength at phi and length +/-len */
void xdrift(double phi, int len, double nlam, 
	    int nx, int ny, float *D, float *C, float *S)
{
   int i, j, k, kp, ix0, ix1, im, ip;
   float *d;
   double Dp, Cp, Sp, tp, lam;
   double y0, y1, x0, x1, dy, c1, s1, cL, sL, cL1, sL1, pi=4*atan(1.0);

   d = (float *)calloc(nx, sizeof(float));

   dy = atan(phi);
   if(ABS(dy*nx) < 0.5) dy = 0;

   lam = (2*len+1) / nlam;
   c1 = cos(2*pi/lam);
   s1 = sin(2*pi/lam);
   cL = cos(2*pi*len/lam);
   sL = sin(2*pi*len/lam);
   cL1 = cos(2*pi*(len+1)/lam);
   sL1 = sin(2*pi*(len+1)/lam);

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
      for(i=MAX(ix0-len-1,0); i<ix0; i++) d[i] = 0.0;
      Dp = Cp = Sp = 0.0;
      for(i=ix0; i<MIN(ix0+len,nx-1); i++) {
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
	 im = MAX(i-len-1, 0);
	 ip = MIN(i+len, nx-1);
	 SLOPE(ip, j, dy, ny, kp);
/* New original values */
	 d[ip] = D[ip+nx*kp];
	 if(i-ix0 > len) {
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

/* Convolve image with 1/cos/sin with nlam wavelength at phi and length +/-len */
void ydrift(double phi, int len, double nlam, 
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

   lam = (2*len+1) / nlam;
   c1 = cos(2*pi/lam);
   s1 = sin(2*pi/lam);
   cL = cos(2*pi*len/lam);
   sL = sin(2*pi*len/lam);
   cL1 = cos(2*pi*(len+1)/lam);
   sL1 = sin(2*pi*(len+1)/lam);

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
      for(i=MAX(iy0-len-1,0); i<iy0; i++) d[i] = 0.0;

      Dp = Cp = Sp = 0.0;
      for(i=iy0; i<MIN(iy0+len,ny-1); i++) {
	 SLOPE(i, j, dx, nx, k);
	 d[i] = D[k+ny*i];
	 Dp += d[i];
	 tp = Cp;
	 Cp = Cp*c1 + Sp*s1 + d[i]*cL;
	 Sp = Sp*c1 - tp*s1 + d[i]*sL;
      }

/* Loop over all y for this trajectory */
      for(i=iy0; i<=iy1; i++) {
	 SLOPE(i, j, dx, nx, k);
	 im = MAX(i-len-1, 0);
	 ip = MIN(i+len, ny-1);
	 SLOPE(ip, j, dx, nx, kp);
/* New original values */
	 d[ip] = D[kp+nx*ip];
	 if(i-iy0 > len) {
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

/* Print out cfitsio error messages and exit program */
void printerror(int status)
{
   if (status) {
      fits_report_error(stderr, status); /* print error report */
      exit( status );    /* terminate the program, returning error status */
   }
   return;
}
