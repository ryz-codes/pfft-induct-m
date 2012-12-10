/* mutual.c
 * Calculates the mutual inductance
 *
 *
 *************************************************************/

#include <math.h>
#include <stdlib.h>
#include "mex.h"

/* OpenMP includes and compilation warnings */
#ifdef _OPENMP
#include <omp.h> /* Needed for VC++ linker */
#pragma message("Parallelism enabled in compilation!")
#else
#pragma message("Warning: Compile with the OpenMP flags to enable parallelism.")
#pragma message("\tSuggested command in GCC\n\t\tmex -v mutual.c CFLAGS=\"-fopenmp -Wall -O3 \\$CFLAGS\" LDFLAGS=\"-fopenmp \\$LDFLAGS\"")
#pragma message("\tSuggested command in VC++\n\t\tmex -v mutual.c COMPFLAGS=\"/openmp /O2 /Wall $COMPFLAGS\"")
#pragma message("\n\n")
#endif

#define PI 3.14159265358979323846
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13
#define MAXsubfils 6          /* maximum subdivisions for subdividing filaments */

/* for deg_mutual, how small a dimension must be to be degenerate */
#define DEG_TOL 1.1e-4

/*-------------------------------------------------------------------------
 * Basic macro functions
 *------------------------------------------------------------------------*/

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

#define nearzero(x) (fabs(x) < EPS)

#define TRUE (1==1)
#define FALSE (1==0)

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef MIN
#define MIN(A,B)  ( (A) < (B) ? (A) : (B) )
#endif

#define compare(x,y,eps) (  (((x)==0 && (y)==0) || (fabs((x) - (y)) < eps*((x) + (y)) )) \
  ? 0 : (  ((x) > (y)) ? 1 : -1 )  )

/* These are missing in some math.h files */
#define atanh(x) 0.5*log((1+(x))/(1-(x)))
#define asinh(x) log((x)+sqrt((x)*(x)+1)) 
#define finite(x) ((x) != HUGE_VAL && (x) != -HUGE_VAL)
  
/*---------------------------------------------------------------------------------------- */
/* Utility functions (some of which were turned into macros) */
/*---------------------------------------------------------------------------------------- */
#define magdiff2(fil1, node1, fil2, node2) \
          ( SQUARE(fil1->x[node1] - fil2->x[node2]) \
	      +SQUARE(fil1->y[node1] - fil2->y[node2]) \
	      +SQUARE(fil1->z[node1] - fil2->z[node2]))
           
#define dotprod(fil1,fil2) \
          (  (fil1->x[1] - fil1->x[0])*(fil2->x[1] - fil2->x[0]) \
          + (fil1->y[1] - fil1->y[0])*(fil2->y[1] - fil2->y[0]) \
          + (fil1->z[1] - fil1->z[0])*(fil2->z[1] - fil2->z[0]) )
   
#define vdotp(v1,v2) \
          (v1[XX]*v2[XX] + v1[YY]*v2[YY] + v1[ZZ]*v2[ZZ])

#define dotp(x1,y1,z1,x2,y2,z2) \
          ((x1)*(x2) + (y1)*(y2) + (z1)*(z2))
     
#define getD(fil,D) \
          D[XX] = fil->x[1] - fil->x[0]; \
          D[YY] = fil->y[1] - fil->y[0]; \
          D[ZZ] = fil->z[1] - fil->z[0]
          
#define getr(x,y,z,s,t,D) \
          *x = s[XX] + t*D[XX]; \
          *y = s[YY] + t*D[YY]; \
          *z = s[ZZ] + t*D[ZZ]

#define mag(x1,y1,z1) \
          sqrt( (x1)*(x1) + (y1)*(y1) + (z1)*(z1) )

#define magsq(x1,y1,z1) \
          ( (x1)*(x1) + (y1)*(y1) + (z1)*(z1) )
          
#define fill_4(vec,E,a,d) \
          vec[0] = (E) - (a); \
          vec[1] = (E) + (d) - (a); \
          vec[2] = (E) + (d); \
          vec[3] = (E)

#define dist_between(x1,y1,z1,x2,y2,z2) \
          sqrt(SQUARE((x1) - (x2)) + SQUARE((y1) - (y2)) + SQUARE((z1) - (z2)))

#define aspectratio(fil) \
          ((fil->height >= fil->width) ? \
                  (fil->height/fil->width) :\
                  (fil->width/fil->height))  
                  
/* COUNTERS */
int num_mutualfil = 0;
int num_perp = 0;
int num_found = 0;
int num_fourfil = 0;
int num_exact_mutual = 0;
int num_self = 0;
int num_quadFil = 0;

typedef struct Filament {
  double x[2], y[2], z[2];  /* endpoints */
  double length, area, width, height;
  double lenvect[3];        /* vector along the length of filament */
  double widthdir[3];
} FILAMENT;

enum degen_type {brick = 0, flat = 1, skinny = 2, too_long = 3, too_short = 4, 
		   short_flat = 5, short_skinny = 6, impossible = 7};


/* Filament utility functions */
void get_wid();
void get_height();
double min_endpt_sep();
double dist_betw_pt_and_fil();

/* Degenerate cases */
double brick_to_brick();
double parallel_fils();
double flat_to_flat_tape();
double flat_to_skinny_tape();

/* Exact terms: close-distance brick to brick */
double exact_mutual();
double parallel_fils(FILAMENT *fil_j, FILAMENT *fil_m,
    int whperp, double *x_j, double *y_j, double dist);

/* Fourfil: mid-distance brick to brick */
void findfourfils();
double fourfil(FILAMENT *fil_j, FILAMENT *fil_m);

/* Selfterm: brick self to brick self */
double self();
double selfterm(FILAMENT *fil);

/* Mutualfil: filament to filament */
double mutualfil(FILAMENT *fil1, FILAMENT *fil2);

/* Main working function mutual() */
double mutual();

/* Matlab interface */
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]);


/*---------------------------------------------------------------------------------------- */
/*  Filament Utility Functions */
/*---------------------------------------------------------------------------------------- */
/* calculates direction of width if not specified */
/* CANDIDATE FOR INLINE EXPANSION */
void get_wid(fil, wid)
FILAMENT *fil;
double *wid;
{

  double wx,wy,wz;
  double mag;

  if (fil->widthdir != NULL) {
    wx = fil->widthdir[XX];
    wy = fil->widthdir[YY];
    wz = fil->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if (fabs(wx/fil->length) < EPS && fabs(wy/fil->length) < EPS) {
      /* if all of x-y is perpendic to length, then choose x direction */
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }
  wid[XX] = wx;
  wid[YY] = wy;
  wid[ZZ] = wz;
}

/* calculates direction of height */
/* CANDIDATE FOR INLINE EXPANSION */
void get_height(fil, wid, height)
FILAMENT *fil;
double *wid, *height;
{
  double wx = wid[XX];
  double wy = wid[YY];
  double wz = wid[ZZ];
  double hx,hy,hz, mag;

  hx = -wy*(fil->z[1] - fil->z[0]) + (fil->y[1] - fil->y[0])*wz;
  hy = -wz*(fil->x[1] - fil->x[0]) + (fil->z[1] - fil->z[0])*wx;
  hz = -wx*(fil->y[1] - fil->y[0]) + (fil->x[1] - fil->x[0])*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;

  height[XX] = hx;
  height[YY] = hy;
  height[ZZ] = hz;
}


/* returns the minimum distance between an endpt of fil1 and one of fil2 */
double min_endpt_sep(fil1,fil2)
FILAMENT *fil1, *fil2;
{
  double min, tmp;
  int idx1, idx2;

  idx1 = 0;
  idx2 = 0;
  min = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);

  idx1 = 1;
  idx2 = 0;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  idx1 = 0;
  idx2 = 1;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  idx1 = 1;
  idx2 = 1;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  return min;
}

/* this finds the shortest distance between the line defined by
   r = s + tnew*D, t in [0,1] (which is along fil_line) and the endpoint
   on fil which is closer to fil_line (determined by the value of t) */
double dist_betw_pt_and_fil(fil_line, D, s, DD, fil,t)
FILAMENT *fil_line, *fil;
double *D, *s, t, DD;
{
  double e[3], sme[3], x,y,z;
  double tnew, Dsme;
  int idx;

  if (t < 0) {
    e[XX] = fil->x[0];
    e[YY] = fil->y[0];
    e[ZZ] = fil->z[0];
  }
  else if (t > 1) {
    e[XX] = fil->x[1];
    e[YY] = fil->y[1];
    e[ZZ] = fil->z[1];
  }
  else {
    /*fprintf(stderr, "Internal err: dist_bet_pt_and_fil: why is t = %lg?\n", t); */
      mexPrintf("Usage error!!");
      exit(1);
  }

  sme[XX] = s[XX] - e[XX];
  sme[YY] = s[YY] - e[YY];
  sme[ZZ] = s[ZZ] - e[ZZ];
  
  Dsme = vdotp(D,sme);
  
  tnew = -Dsme/DD;

  if (tnew <= 1 && tnew >= 0) {
    /* This will be the case when a small fil is near a big one (fil_line).*/
    /* Calculate r = (s - e) + tnew*D */
    getr(&x,&y,&z,sme,tnew,D);
    return sqrt(x*x + y*y + z*z);
  }
  else {
    /* just find the distance between the nearest endpt and e[] */
    if (tnew < 0)
      idx = 0;
    else
      idx = 1;

    return dist_between(e[XX], e[YY], e[ZZ], 
			fil_line->x[idx], fil_line->y[idx], fil_line->z[idx]);
  }
}



double dist_betw_fils(fil1, fil2, parallel)
FILAMENT *fil1, *fil2;
int *parallel;
{
  double x1,y1,z1,x2,y2,z2;
  double D1[3], D2[3], s1[3], s2[3], t1, t2, s1ms2[3];
  double D1D1, D1D2, D2D2, D1s1s2, D2s1s2;
  double tmp1;
  
  s1[XX] = fil1->x[0];
  s1[YY] = fil1->y[0];
  s1[ZZ] = fil1->z[0];
  s2[XX] = fil2->x[0];
  s2[YY] = fil2->y[0];
  s2[ZZ] = fil2->z[0];

  s1ms2[XX] = s1[XX] - s2[XX];
  s1ms2[YY] = s1[YY] - s2[YY];
  s1ms2[ZZ] = s1[ZZ] - s2[ZZ];

  getD(fil1, D1);
  getD(fil2, D2);

  D1D1 = vdotp(D1,D1);
  D1D2 = vdotp(D1,D2);
  D2D2 = vdotp(D2,D2);
  D1s1s2 = vdotp(D1,s1ms2);
  D2s1s2 = vdotp(D2,s1ms2);

  tmp1 = D1D2*D1D2/D1D1 - D2D2;

  if (fabs(tmp1/D2D2) < EPS) {
    /* fils are parallel */
    *parallel = 1;
    return min_endpt_sep(fil1,fil2);
  }
  else
    *parallel = 0;

  t2 = (D1D2*D1s1s2/D1D1 - D2s1s2)/tmp1;
  t1 = (t2*D1D2 - D1s1s2)/D1D1;

  if (t1 <= 1 && t1 >= 0) {
    if (t2 <= 1 && t2 >= 0) {
      getr(&x1,&y1,&z1,s1,t1,D1);
      getr(&x2,&y2,&z2,s2,t2,D2);
      return dist_between(x1,y1,z1,x2,y2,z2);
    }
    else
      /* nearest point along fil2 is outside line segment defining filament */
      return dist_betw_pt_and_fil(fil1, D1, s1, D1D1, fil2,t2);
  }
  else {
    if (t2 <= 1 && t2 >= 0) {
      /* nearest point along fil1 is outside line segment defining filament */
      return dist_betw_pt_and_fil(fil2, D2, s2, D2D2, fil1,t1);
    }    
    else 
      /* both point are out of range, just compare endpoints */
      return min_endpt_sep(fil1,fil2);
  }
  
}

/* This assumes the lengths of the fils are parallel and returns 1 if the 
   side faces are parallel also */
int edges_parallel(fil_j, fil_m, wid1, whperp)
FILAMENT *fil_j, *fil_m;
int *whperp;
double *wid1;
{
  double *wid_j = fil_j->widthdir;
  double *wid_m = fil_m->widthdir;
  double wid2[3];
  double prod,mj,mm;

  if (wid_j == NULL && wid_m == NULL) {
    /* both have unspecified width directions and their lengths are assumed
       parallel, so they are parallel */
    *whperp = 0;
    return TRUE;
  }
  else {
    if (wid_j == NULL) {
      /* get_wid(fil_j, wid1); */ 
      wid_j = wid1;
    }
    if (wid_m == NULL) {
      get_wid(fil_m, wid2);
      wid_m = wid2;
    }
    mj = mag(wid_j[XX],wid_j[YY],wid_j[ZZ]);
    mm = mag(wid_m[XX],wid_m[YY],wid_m[ZZ]);
    prod = dotp(wid_j[XX],wid_j[YY],wid_j[ZZ],wid_m[XX],wid_m[YY],wid_m[ZZ])
            / (mm*mj);
    if (fabs(prod) < EPS) {
      *whperp = 1;   /* width and height are perpend. to that on other fil*/
      return TRUE;
    }
    else if (fabs( fabs(prod) - 1 ) < EPS) {
      *whperp = 0;
      return TRUE;
    }
    else
      return FALSE;
  }
}



/*---------------------------------------------------------------------------------------- */
/*  Degenerate expressions */
/*---------------------------------------------------------------------------------------- */
#define LEN 4
#define WID 2
#define HEIGHT 1

double brick_to_brick(E,a,d,P,b,c,l3,l1,l2)
double E,a,d,P,b,c,l3,l1,l2;
{
  double q[4], r[4], s[4], totalM;
  int i,j,k, sign2;
  double eval_eq();

  fill_4(q, E,a,d);
  fill_4(r, P,b,c);
  fill_4(s, l3,l1,l2);
  
  totalM = 0;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      for(k = 0; k < 4; k++) {
	sign2 = ( (i+j+k)%2 == 0 ? 1 : -1);
	totalM += sign2*eval_eq(q[i],r[j],s[k], a);
      }

  return totalM;
}

double flat_to_flat_tape(E,a,d,P,l3,l1,l2)
double E,a,d,P,l3,l1,l2;
{
  double q[4], s[4], totalM;
  int i,k, sign2;
  double eval_eq_tape();

  fill_4(q, E,a,d);
  fill_4(s, l3,l1,l2);

  totalM = 0;
  
  for(i = 0; i < 4; i++)
      for(k = 0; k < 4; k++) {
	sign2 = ( (i+k)%2 == 0 ? 1 : -1);
	totalM += sign2*eval_eq_tape(q[i],P,s[k], a);
      }

  return totalM;
}

double eval_eq_tape(x,y,z,ref_len)
double x,y,z, ref_len;
{
  static double one_6 = 1.0/6.0;
  double retval;
  double len, xsq, ysq, zsq;
  double one_over_ref_len, one_over_ref_len_sq;
  int num_nearzero_sq;

  one_over_ref_len = 1.0/MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
  one_over_ref_len_sq = SQUARE(one_over_ref_len);

  xsq = x*x;
  ysq = y*y;
  zsq = z*z;

  len = sqrt(xsq + ysq + zsq);

  retval = -one_6*len*(xsq - 2*ysq + zsq);

  num_nearzero_sq = nearzero(xsq*one_over_ref_len_sq) 
                 + nearzero(ysq*one_over_ref_len_sq)
		 + nearzero(zsq*one_over_ref_len_sq);

  if (num_nearzero_sq < 2)
    retval += 0.5*( (xsq - ysq)*z*log(z + len) + (zsq - ysq)*x*log(x + len) ); 

  if (!nearzero(y*one_over_ref_len))
    retval -= x*y*z*atan(x*z/(y*len));

  return retval;
}

double flat_to_skinny_tape(E,a,P,c,l3,l1,l2)
double E,a,P,c,l3,l1,l2;
{
  double q[2], r[2], s[4], totalM;
  int i,j,k, sign2;
  double eval_eq_tape2();

  q[0] = E;
  q[1] = E - a;
  r[0] = P + c;
  r[1] = P;
  fill_4(s, l3,l1,l2);

  totalM = 0;
  
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      for(k = 0; k < 4; k++) {
	sign2 = ( (i+j+k)%2 == 0 ? 1 : -1);
	totalM += sign2*eval_eq_tape2(q[i],r[j],s[k], a);
      }

  return totalM;
}

double eval_eq_tape2(x,y,z,ref_len)
double x,y,z, ref_len;
{
  static double one_6 = 1.0/6.0;
  static double one_3 = 1.0/3.0;
  int num_nearzero;
  double retval;
  double len, xsq, ysq, zsq;
  double one_over_ref_len, one_over_ref_len_sq;
  double tan_tape();
  int nzxsq, nzysq, nzzsq;

  one_over_ref_len = 1.0/MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
  one_over_ref_len_sq = SQUARE(one_over_ref_len);

  xsq = x*x;
  ysq = y*y;
  zsq = z*z;

  len = sqrt(xsq + ysq + zsq);

  retval = -one_3*len*x*y;

  nzxsq = nearzero(xsq*one_over_ref_len_sq);
  nzysq = nearzero(ysq*one_over_ref_len_sq);
  nzzsq = nearzero(zsq*one_over_ref_len_sq);

  if (!(nzzsq && nzysq))
    retval += (0.5*zsq - one_6*ysq)*y*log(x + len);
  if (!(nzzsq && nzxsq))
    retval += (0.5*zsq - one_6*xsq)*x*log(y + len);
  if (!(nzzsq || nzysq || nzxsq))
    retval += x*y*z*log(z + len); 

  num_nearzero = nearzero(x*one_over_ref_len) 
                 + nearzero(y*one_over_ref_len)
		 + nearzero(z*one_over_ref_len);

  if (num_nearzero < 1)
    retval -= zsq*z*one_6*tan_tape(x,y,z,len) 
               + 0.5*z*(xsq*tan_tape(y,z,x,len) + ysq*tan_tape(x,z,y,len));

  return retval;
}

double tan_tape(x,y,z,len)
double x,y,z,len;
{
  return atan(x*y/(z*len));
}

enum degen_type find_deg_dims(fil)
FILAMENT *fil;
{
  double max;
  
  max = MAX(fil->length, fil->width);
  max = MAX(max, fil->height);

  return (fil->length/max < DEG_TOL)*LEN + (fil->width/max < DEG_TOL)*WID
         + (fil->height/max < DEG_TOL)*HEIGHT;
}

void setup_tape_to_tape(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m,
			  nfil_j, nfil_m, nx_j, ny_j)
FILAMENT *fil_j, *fil_m, *nfil_j, *nfil_m;
int whperp;
double *x_j, *y_j, **nx_j, **ny_j;  /* unit vectors in the fil coord sys */
enum degen_type deg_j, deg_m;
{

  if (deg_j == flat) {
    *nfil_j = *fil_j;
    *nfil_m = *fil_m;
    *nx_j = x_j;
    *ny_j = y_j;
  }
  else if (deg_j == skinny) {
    /* turn skinny into flat orientation */
    *nfil_j = *fil_j;
    *nfil_m = *fil_m;
    /* swap coord sys */
    *ny_j = x_j;
    *nx_j = y_j;
    /* swap height and width */
    nfil_j->width = fil_j->height;
    nfil_j->height = fil_j->width;
    nfil_m->width = fil_m->height;
    nfil_m->height = fil_m->width;
  }
}

double compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j, 
			      deg_j, deg_m, dist)
FILAMENT *fil_j, *fil_m;
int whperp;
double *x_j, *y_j;  /* unit vectors in the fil coord sys */
enum degen_type deg_j, deg_m;
double dist;
{

  FILAMENT nfil_j, nfil_m;   /* new temp fils */
  double *nx_j = NULL; 
  double *ny_j = NULL;

  if (deg_j == brick && deg_m == brick) {
    /* neither is degenerate, this shouldn't happen */
    fprintf(stderr,"Hey, compute_degenerate was called, impossible!\n");
    exit(1);
  }

  if ((deg_j == flat || deg_j == skinny)&&(deg_m == flat || deg_m == skinny)){
    setup_tape_to_tape(fil_j,fil_m,whperp,x_j,y_j,deg_j,deg_m,
		       &nfil_j,&nfil_m, &nx_j, &ny_j);
    return exact_mutual(&nfil_j, &nfil_m, whperp, nx_j, ny_j, deg_j, deg_m);
  }
  else if ( (deg_m == brick && (deg_j == flat || deg_j == skinny))
	   || (deg_j == brick && (deg_m == flat || deg_m == skinny)))
    return exact_mutual(&nfil_j, &nfil_m, whperp, nx_j, ny_j, deg_j, deg_m);
  else if ( deg_j == too_long && deg_m == too_long)
    return fourfil(fil_j, fil_m);
  else if (deg_j == too_long || deg_j == too_long)
    return fourfil(fil_j, fil_m);
  else
    return fourfil(fil_j, fil_m);

}

/*---------------------------------------------------------------------------------------- */
/*  EXACT_MUTUAL() */
/*---------------------------------------------------------------------------------------- */
/* exact mutual inductance based on C. Hoer and C.Love, 
   Journal of the National Bureau of Standards-C,  Vol. 69C, p 127-137, 1965.*/
   
#define log_term(x, xsq, ysq, zsq, len) \
           (((6*(zsq) - (ysq))*(ysq) - (zsq)*(zsq))*(x)*log( ((x) + (len))/sqrt((ysq) + (zsq))))

#define tan_term(x, y, z, zsq, len) \
           (x)*(y)*(z)*(zsq)*(atan((x)*(y)/((z)*(len))))

double eval_eq(double x, double y, double z, double ref_len)
{
  static double one_60 = 1.0/60.0;
  static double one_6 = 1.0/6.0;
  static double one_24 = 1.0/24.0;
  double retval;
  double len, xsq, ysq, zsq;
  int num_nearzero;
  int num_nearzero_sq;
  double one_over_ref_len;
  double one_over_ref_len_sq; 

  one_over_ref_len = 1.0/MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
  one_over_ref_len_sq = SQUARE(one_over_ref_len);

  xsq = x*x;
  ysq = y*y;
  zsq = z*z;

  len = sqrt(xsq + ysq + zsq);

  retval = one_60*len
              *(xsq*(xsq - 3*ysq) + ysq*(ysq - 3*zsq) + zsq*(zsq - 3*xsq));

  num_nearzero = nearzero(x*one_over_ref_len) 
                 + nearzero(y*one_over_ref_len)
		 + nearzero(z*one_over_ref_len);

  num_nearzero_sq = nearzero(xsq*one_over_ref_len_sq) 
                 + nearzero(ysq*one_over_ref_len_sq)
		 + nearzero(zsq*one_over_ref_len_sq);

  if (num_nearzero_sq < 2)
    retval += one_24*(log_term(x, xsq, ysq, zsq, len) 
		      + log_term(y, ysq, xsq, zsq, len) 
		      + log_term(z, zsq, xsq, ysq, len));

  if (num_nearzero < 1)
    retval -= one_6*(tan_term(x,y,z,zsq,len) + tan_term(x,z,y,ysq,len) 
		     + tan_term(z,y,x,xsq,len));

  return retval;
}

double exact_mutual(FILAMENT *fil_j, FILAMENT *fil_m,
int whperp,
double *x_j, double *y_j,  /* unit vectors in the fil coord sys */
int deg_j, int deg_m)
{
  double z_j[3];  /* unit vectors in the filament coord sys*/
  double origin[3];
  double ox,oy,oz, length;
  double a,b,c,d,l1,l2,l3,E,P, l3_1;
  double endx, endy, endz;
  int sign;
  double totalM=0; 
  int a_deg, b_deg, c_deg, d_deg;

  length = fil_j->length;
  z_j[XX] = fil_j->lenvect[XX]/length;
  z_j[YY] = fil_j->lenvect[YY]/length;
  z_j[ZZ] = fil_j->lenvect[ZZ]/length;

  a = fil_j->width;
  b = fil_j->height;
  if (whperp == 0) {
    c = fil_m->height;
    d = fil_m->width;
  }
  else {
    d = fil_m->height;
    c = fil_m->width;
  }
    
  ox = origin[XX] = fil_j->x[0] - x_j[XX]*a/2 - y_j[XX]*b/2;
  oy = origin[YY] = fil_j->y[0] - x_j[YY]*a/2 - y_j[YY]*b/2;
  oz = origin[ZZ] = fil_j->z[0] - x_j[ZZ]*a/2 - y_j[ZZ]*b/2;

  endx = fil_m->x[0] - ox;
  endy = fil_m->y[0] - oy;
  endz = fil_m->z[0] - oz;

  E = dotp(x_j[XX], x_j[YY], x_j[ZZ], endx, endy, endz) - d/2;

  P = dotp(y_j[XX], y_j[YY], y_j[ZZ], endx, endy, endz) - c/2;

  l3 = dotp(z_j[XX], z_j[YY], z_j[ZZ], endx, endy, endz);
  l3_1 = dotp(z_j[XX], z_j[YY], z_j[ZZ],fil_m->x[1] - ox, fil_m->y[1] - oy, 
	      fil_m->z[1] - oz);

  l1 = fil_j->length;
  l2 = fil_m->length;

  if ( fabs(fabs(l3 - l3_1) - l2)/l2 > EPS) {
    /*fprintf(stderr, "Huh?  filament not expected length\n"); */
    /*exit1); */
  }
  
  if (l3 <= l3_1)
    sign = 1;
  else {
    sign = -1;
    l3 = l3_1;
  }

  a_deg = a/MAX(l1,b) < DEG_TOL;
  b_deg = b/MAX(l1,a) < DEG_TOL;
  c_deg = c/MAX(l2,d) < DEG_TOL;
  d_deg = d/MAX(l2,c) < DEG_TOL;


  if (a_deg && b_deg && c_deg && d_deg) {
    /* two long filaments  */
    totalM = fourfil(fil_j, fil_m)/(sign*MUOVER4PI);
  }
  else if (!a_deg && b_deg && c_deg && d_deg) {
    /* one flat and one long  */
    /*totalM = tape_to_fil(E + d/2, a, P - b/2 + c/2,l3,l1,l2) / a; */
      mexPrintf("FastHenry never implemented tape to fil, but it should never have gotten here either.\n");
      exit(1);
  }
  else if (!a_deg && b_deg && c_deg && !d_deg) {
    /* two flat  */
    totalM = flat_to_flat_tape(E,a,d,P - b/2 + c/2 ,l3,l1,l2) / (a * d);
  }
  else if (!a_deg && b_deg && !c_deg && d_deg) {
    /* one flat and one skinny  */
    totalM = flat_to_skinny_tape(E + d/2,a,P - b/2 ,c,l3,l1,l2) / (a * c);
  }
  else if (!a_deg && !b_deg && c_deg && d_deg) {
    /* totalM = brick_to_fil(E + d/2,a,P + c/2,b,l3,l1,l2) / (a * b); */
    mexPrintf("FastHenry never implemented brick to fil, but it should never have gotten here either.\n");
      exit(1);
  }
  else if (deg_j != brick || deg_m != brick) {
    /*fprintf(stderr,"Internal Error: Bad Degenerate filament a/l1<tol\n"); */
    /*fprintf(stderr,"Using fourfil() instead...\n"); */
    totalM = fourfil(fil_j, fil_m)/(sign*MUOVER4PI);
  }
  else {
    /* all that are left are bricks  */
    /* this is the full 3-D filament calculation, no degeneracies */
    totalM = brick_to_brick(E,a,d,P,b,c,l3,l1,l2)/ (a * b * c * d);}
  return sign*MUOVER4PI*totalM;
}
    
int realcos_error = 0;

void print_infinity_warning(fil1, fil2)
     FILAMENT *fil1, *fil2;
{
  FILAMENT *fil;

  mexPrintf("Severe warning: mutual inductance = infinity for two filaments:\n");
  /*fprintf(stderr,"Severe warning: mutual inductance = infinity for two filaments:\n"); */

  fil = fil1;
  mexPrintf("  fil 1: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
  	  fil->x[0],fil->y[0],fil->z[0],fil->x[1],fil->y[1],fil->z[1],
   fil->width, fil->height);
  fil = fil2;
  mexPrintf("  fil 2: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
  	  fil->x[0],fil->y[0],fil->z[0],fil->x[1],fil->y[1],fil->z[1],
	  fil->width, fil->height);
  mexPrintf("Probably because there are overlapping but non-orthogonal segments in the input\n");
}

double parallel_fils(fil_j, fil_m, whperp, x_j, y_j, dist)
FILAMENT *fil_j, *fil_m;
int whperp;
double *x_j, *y_j;  /* unit vectors in the fil coord sys */
double dist;
{
  enum degen_type deg_j, deg_m;
  
  /* find degenerate dimensions */
  deg_j = find_deg_dims(fil_j);
  deg_m = find_deg_dims(fil_m);
  
  
  if (deg_j == brick && deg_m == brick) {
    /* no degenerate dimensions, both are bricks */
    

    return exact_mutual(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);

  }
  else {

  return compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j,
				  deg_j, deg_m, dist);

  }

}

/*---------------------------------------------------------------------------------------- */
/*  FOURFIL() */
/*---------------------------------------------------------------------------------------- */
/* calculates mid-distance interactions by the four filament quadrature approximation */

void findfourfils(fil, subfils)
FILAMENT *fil, subfils[4];
{
  double hx,hy,hz,mag,wx,wy,wz;
  int i;

  if (fil->widthdir != NULL) {
    wx = fil->widthdir[XX];
    wy = fil->widthdir[YY];
    wz = fil->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if ( fabs(wx/fil->length) < EPS && fabs(wy/fil->length) < EPS) {
      /* if all of x-y is perpendic to length, then choose x direction */
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }
  
  hx = -wy*(fil->z[1] - fil->z[0]) + (fil->y[1] - fil->y[0])*wz;
  hy = -wz*(fil->x[1] - fil->x[0]) + (fil->z[1] - fil->z[0])*wx;
  hz = -wx*(fil->y[1] - fil->y[0]) + (fil->x[1] - fil->x[0])*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;

  /* all mutualfil needs are the filament coordinates and length */
  for (i = 0; i < 2; i++) {
    subfils[0].x[i] = fil->x[i] + fil->width*wx/2; 
    subfils[0].y[i] = fil->y[i] + fil->width*wy/2; 
    subfils[0].z[i] = fil->z[i] + fil->width*wz/2; 

    subfils[1].x[i] = fil->x[i] - fil->width*wx/2; 
    subfils[1].y[i] = fil->y[i] - fil->width*wy/2; 
    subfils[1].z[i] = fil->z[i] - fil->width*wz/2; 

    subfils[2].x[i] = fil->x[i] + fil->height*hx/2; 
    subfils[2].y[i] = fil->y[i] + fil->height*hy/2; 
    subfils[2].z[i] = fil->z[i] + fil->height*hz/2; 

    subfils[3].x[i] = fil->x[i] - fil->height*hx/2; 
    subfils[3].y[i] = fil->y[i] - fil->height*hy/2; 
    subfils[3].z[i] = fil->z[i] - fil->height*hz/2; 
  }

  for(i = 0; i < 4; i++)
    subfils[i].length = fil->length;
}

double fourfil(fil_j, fil_m)
FILAMENT *fil_j, *fil_m;
{
  FILAMENT subfilj[4], subfilm[4];
  double totalM;
  int i;
  
  /* approximate 'filament' with width and length as a combination */
  /* of four filaments on the midpoints of the edges               */
  /* Known as the Rayleigh Quadrature formula. Grover p.11         */
      findfourfils(fil_j, subfilj);
      findfourfils(fil_m, subfilm);
      totalM = 0.0;
      for(i = 0; i < 4; i++)
        totalM += mutualfil(fil_j, &subfilm[i]);
      for(i = 0; i < 4; i++)
        totalM += mutualfil(fil_m, &subfilj[i]);
      totalM += -2.0*mutualfil(fil_j, fil_m);

      totalM = totalM/6.0;

  return totalM;
}


/*------------------------------------------------------------------------- */
/*  SPLITFIL() */
/*------------------------------------------------------------------------- */
/* Here, we will use the fourfil method but divide the filaments down along */
/* cross sections when their aspect ratios aren't great */
  
#define shiftFil(fil1, fil2, dist, dirj)  \
    (fil2)->x[0] = (fil1)->x[0] + dirj[0]*(dist); \
    (fil2)->y[0] = (fil1)->y[0] + dirj[1]*(dist); \
    (fil2)->z[0] = (fil1)->z[0] + dirj[2]*(dist); \
    (fil2)->x[1] = (fil1)->x[1] + dirj[0]*(dist); \
    (fil2)->y[1] = (fil1)->y[1] + dirj[1]*(dist); \
    (fil2)->z[1] = (fil1)->z[1] + dirj[2]*(dist)

#define shiftFac(N) (1/((double)N))
#define offsetFac(N) (1/(2*((double)N)) - 0.5)

    
double splitfil(fil_j, fil_m)
FILAMENT *fil_j, *fil_m;
{
  FILAMENT subfilj, subfilm;
  double totalM = 0.0;
  double widj[3], widm[3], heij[3], heim[3];
  int i,j,k,m;
  
  double N1=2; /* Width fil_j, iterator i */
  double M1=2; /* Height fil_j, iterator j */
  double N2=2; /* Width fil_m, iterator k */
  double M2=2; /* Height fil_m, iterator m */
  
  /* Determine the division fineness by examining the aspect ratio. */
  /* These are experimentally verified as being nice numbers */
  if (MIN(fil_j->length,fil_m->length) 
        < 10 * MAX(fil_j->width,fil_j->height)) 
  { 
    N1 = 4; M1 = 4;
    N2 = 3; M2 = 3;
  }

  get_wid(fil_j, widj);
  get_height(fil_j, widj, heij);
  get_wid(fil_m, widm);
  get_height(fil_m, widm, heim);
  
  /* Setup the virtual filaments */
  subfilj.length = fil_j->length;
  subfilm.length = fil_m->length;
  
  /* approximate 'filament' with width and length as a combination */
  /* of four filaments on the midpoints of the edges               */
  /* Known as the Rayleigh Quadrature formula. Grover p.11         */
  /* Cycle through the "to" filament */
  for (i = 0; i < N1; i++) {
    for (j = 0; j < M1; j++) {
        /* Set up "to" filament subfilj */
        shiftFil(fil_j, &subfilj, 
          (fil_j->width)*(offsetFac(N1)+i*shiftFac(N1)),
                widj);
        shiftFil(&subfilj, &subfilj, 
          (fil_j->height)*(offsetFac(M1)+j*shiftFac(M1)),
                heij);

        /* Cycle through filaments on horizontal edges */
        for (k=0; k<=N2; k++) {
            for (m=0; m<=M2; m++) {
                
                /* Cycle through filaments on horizontal edges */
                if (k != N2) {
                    shiftFil(fil_m, &subfilm, 
                      fil_m->width*(offsetFac(N2)+k*shiftFac(N2)),
                        widm);
                    shiftFil(&subfilm, &subfilm, 
                      fil_m->height*(m*shiftFac(M2)-0.5),
                        heim);
                    if (m==0 || m==M2)
                        totalM += mutualfil(&subfilj, &subfilm);
                    else
                        totalM += 2.0*mutualfil(&subfilj, &subfilm);
                }
        
                /* Cycle through filaments on vertical edges */
                if (m != M2) {
                    shiftFil(fil_m, &subfilm, 
                      fil_m->width*(k*shiftFac(N2)-0.5),
                        widm);
                    shiftFil(&subfilm, &subfilm, 
                      fil_m->height*(offsetFac(M2)+m*shiftFac(M2)),
                        heim);
                    if (k==0 || k==N2)
                        totalM += mutualfil(&subfilj, &subfilm);
                    else
                        totalM += 2.0*mutualfil(&subfilj, &subfilm);
                }
                    
                /* Cycle through filaments in the middle */
                if (k != N2 && m != M2) {
                    shiftFil(fil_m, &subfilm, 
                      fil_m->width*(offsetFac(N2)+k*shiftFac(N2)),
                        widm);
                    shiftFil(&subfilm, &subfilm, 
                      fil_m->height*(offsetFac(M2)+m*shiftFac(M2)),
                        heim);
                    totalM -= 2.0*mutualfil(&subfilj, &subfilm);
                }
            }
        }
    }
  }
  
      
  /* Now go around the other direction. */
  /* Cycle through the "to" filament */
  for (i = 0; i < N2; i++) {
    for (j = 0; j < M2; j++) {
        /* Set up "to" filament subfilm */
        shiftFil(fil_m, &subfilm, 
          (fil_m->width)*(offsetFac(N2)+i*shiftFac(N2)),
                widm);
        shiftFil(&subfilm, &subfilm, 
          (fil_m->height)*(offsetFac(M2)+j*shiftFac(M2)),
                heim);
        
        for (k=0; k<=N1; k++) {
            for (m=0; m<=M1; m++) {
                /* Cycle through filaments on horizontal edges */
                if (k != N1) {
                    shiftFil(fil_j, &subfilj, 
                      fil_j->width*(offsetFac(N1)+k*shiftFac(N1)),
                        widj);
                    shiftFil(&subfilj, &subfilj, 
                      fil_j->height*(m*shiftFac(M1)-0.5),
                        heij);
                    if (m==0 || m==M1)
                        totalM += mutualfil(&subfilm, &subfilj);
                    else
                        totalM += 2.0*mutualfil(&subfilm, &subfilj);
                }
                
                /* Cycle through filaments on vertical edges */
                if (m != M1) {
                    shiftFil(fil_j, &subfilj, 
                      fil_j->width*(k*shiftFac(N1)-0.5),
                        widj);
                    shiftFil(&subfilj, &subfilj, 
                      fil_j->height*(offsetFac(M1)+m*shiftFac(M1)),
                        heij);
                    if (k==0 || k==N1)
                        totalM += mutualfil(&subfilm, &subfilj);
                    else
                        totalM += 2.0*mutualfil(&subfilm, &subfilj);
                }
            }
        }
    }
  }

  totalM = totalM/(6.0*N1*M1*N2*M2);
  return totalM;
}    
    

    
/*---------------------------------------------------------------------------------------- */
/*  SELFTERM() */
/*---------------------------------------------------------------------------------------- */
/* calculates selfinductance of rectangular filament */
/* it uses an exact expression for the 6-fold integral from Ruehli */

double self(W,L,T)
double W,L,T; 
{

    double w,t,aw,at,ar,r, z;

    w = W/L; 
    t = T/L; 
    r = sqrt(w*w+t*t); 
    aw = sqrt(w*w+1.0); 
    at = sqrt(t*t+1.0); 
    ar = sqrt(w*w+t*t+1.0); 

    z = 0.25 * ((1/w) * asinh(w/at) + (1/t) * asinh(t/aw) + asinh(1/r)); 
    z += (1/24.0) * ((t*t/w) * asinh(w/(t*at*(r+ar))) + (w*w/t) * asinh(t/(w*aw*(r+ar))) + 
		     ((t*t)/(w*w)) * asinh(w*w/(t*r*(at+ar))) + ((w*w)/(t*t))*asinh(t*t/(w*r*(aw+ar))) + 
		     (1.0/(w*t*t)) * asinh(w*t*t/(at*(aw+ar))) + (1.0/(t*w*w))*asinh(t*w*w/(aw*(at+ar)))); 
    z -= (1.0/6.0) * ((1.0/(w*t)) * atan(w*t/ar) + (t/w) * atan(w/(t*ar)) + (w/t) * atan(t/(w*ar))); 
    z -= (1.0/60.0) * ( ((ar+r+t+at)*t*t)/((ar+r)*(r+t)*(t+at)*(at+ar)) 
		       + ((ar+r+w+aw)*(w*w)) / ((ar+r)*(r+w)*(w+aw)*(aw+ar)) 
		       + (ar+aw+1+at)/((ar+aw)*(aw+1)*(1+at)*(at+ar))); 
    z -= (1.0/20.0)*((1.0/(r+ar)) + (1.0/(aw+ar)) + (1.0/(at+ar))); 
    
    z *= (2.0/PI); 
    z *= L;  /* this is inductance */
    
    num_self += 1;
    
    return MU0*z;

}

/*---------------------------------------------------------------------------------------- */
/*  MUTUALFIL() */
/*---------------------------------------------------------------------------------------- */
/* calculates the mutual inductance between two filaments 
   from Grover, Chapter 7 

   this gives the mutual inductance of two filaments who represent
   opposite sides of a rectangle */
#define mut_rect(len, d) (sqrt((len)*(len) + (d)*(d)) - (len)*asinh((len)/(d)))
#define REPORTFIL(fil1,fil2) mexPrintf("fil1: x=%g y=%g z=%g, fil2:x=%g y=%g z=%g\n",\
    fil1.x[0],fil1.y[0],fil1.z[0],fil2.x[0],fil2.y[0],fil2.z[0]);

double mutualfil(fil1, fil2)
FILAMENT *fil1, *fil2;
{
  static double cose;

  /* R is the vector from fil1 start to fil2 start */
  double Rx = fil2->x[0] - fil1->x[0];  
  double Ry = fil2->y[0] - fil1->y[0];
  double Rz = fil2->z[0] - fil1->z[0];
  
  double R1sq = magdiff2(fil1,1,fil2,1);
  double R2sq = magdiff2(fil1,1,fil2,0);
  double R3sq = Rx*Rx + Ry*Ry + Rz*Rz;
  double R4sq = magdiff2(fil1,0,fil2,1);
  double R1 = sqrt(R1sq);
  double R2 = sqrt(R2sq);
  double R3 = sqrt(R3sq);
  double R4 = sqrt(R4sq);

  double minR = MIN(MIN(MIN(R1,R2),R3),R4);
  
  double l = fil1->length;
  double m = fil2->length;
  
  /* segments touching */
  if (nearzero(minR)) {
      double R;
      if (fabs(R1) < EPS)  R = R3;
      else if (fabs(R2) < EPS)  R = R4; 
      else if (fabs(R3) < EPS)  R = R1;
      else R = R2;

      return MUOVER4PI*2*(dotprod(fil1,fil2)/(l*m))
	   *(l*atanh(m/(l+R)) + m*atanh(l/(m+R)));
      /* note: dotprod should take care of signofM */
    }

  /* let's use the real cosine */
  cose = dotprod(fil1, fil2)/(l*m);

  /* Segments are perpendicular! */
  if (nearzero(cose))
    return 0.0;
  
  /* filaments parallel */
  if (nearzero(fabs(cose) - 1)) {

    /* Declare all our parallel-related variables. */
    double vx, vy, vz;
    double x2_0, x2_1;
    /* determine a vector in the direction of d with length d */
    /* (d is the distance between the lines made by the filament */
    /* u is the vector from fil1 start to fil1 end */
    double ux = (fil1->x[1] - fil1->x[0])/l;
    double uy = (fil1->y[1] - fil1->y[0])/l;
    double uz = (fil1->z[1] - fil1->z[0])/l;

    double dotp = ux*Rx + uy*Ry + uz*Rz;  /* component of R in direction of fil1 */

    /* d vector is R vector without its component in the direction of fils */
    double dx = Rx - dotp*ux;
    double dy = Ry - dotp*uy;
    double dz = Rz - dotp*uz;

    double x1_1 = l; /* Declare this to be a double for the precision it provides. */
    
    double d = sqrt(dx*dx + dy*dy + dz*dz);
    
    /* x2_0 = dotprod( fil2.node0 - (fil1.node0 + d), u ) */
    /* (dotproduct just gives it correct sign) */
    vx =  (fil2->x[0] - ( fil1->x[0] + dx));
    vy =  (fil2->y[0] - ( fil1->y[0] + dy));
    vz =  (fil2->z[0] - ( fil1->z[0] + dz));
    x2_0 = vx*ux + vy*uy + vz*uz;
    
    /* same thing for x2_1 */
    vx =  (fil2->x[1] - ( fil1->x[0] + dx));
    vy =  (fil2->y[1] - ( fil1->y[0] + dy));
    vz =  (fil2->z[1] - ( fil1->z[0] + dz));
    x2_1 = vx*ux + vy*uy + vz*uz;

    /* let fil1 be the x axis, with node 0 being origin and u be */
    /* its positive direction */
    if (nearzero(d))  { /* collinear! */
      return MUOVER4PI*(fabs(x2_1)*log(fabs(x2_1))
	             - fabs(x2_1 - x1_1)*log(fabs(x2_1 - x1_1))
		     - fabs(x2_0)*log(fabs(x2_0))
		     + fabs(x2_0 - x1_1)*log(fabs(x2_0 - x1_1)) );
    }  /* end collinear */
    return MUOVER4PI*(mut_rect(x2_1 - x1_1,d) - mut_rect(x2_1,d) 
		 - mut_rect(x2_0 - x1_1,d) + mut_rect(x2_0,d) );
  } else {  /* end if parallel filaments */
      /* Fill up the rest of the variables needed */
      double maxR = MAX(MAX(MAX(R1,R2),R3),R4);
      double maxlength = (l > m ? l : m);
      double alpha = R4sq - R3sq + R2sq - R1sq;
      
     /* the rest if for arbitrary filaments */
      double l2 = l*l;
      double m2 = m*m;
      double alpha2 = alpha*alpha;

      double u = l*(2*m2*(R2sq -R3sq - l2) + alpha*(R4sq - R3sq - m2))
            / (4*l2*m2 - alpha2);
      double v = m* (2*l2*(R4sq - R3sq - m2) + alpha*(R2sq - R3sq - l2))
            / (4*l2*m2 - alpha2);

      double u2 = u*u;
      double v2 = v*v;

      double sinesq = 1.0 - cose*cose;
      double sine = sqrt(sinesq);
      double omega, tmp1, tmp2;
      
      double d = (R3sq - u2 - v2 + 2*u*v*cose);
      if (nearzero(d/(R3sq + u2 + v2 + 1)*(maxlength*maxlength/(maxR*maxR))))
        d = 0.0;
      d = sqrt(d);
      tmp1 = d*d*cose;
      tmp2 = d*sine;
      /* sine will never be zero because we have already checked  */
      /* for (nearzero(fabs(cose) - 1)) above */
      if (fabs(d) < EPS) 
        omega = 0.0;   /* d is zero, so it doesn't matter */
      else
        omega = atan2( (tmp1 + (u+l)*(v + m)*sinesq),(tmp2*R1))
               - atan2( (tmp1 + (u + l)*v*sinesq),(tmp2*R2))
           + atan2( (tmp1 + u*v*sinesq),(tmp2*R3))
           - atan2( (tmp1 + u*(v + m)*sinesq),(tmp2*R4) );
      return MUOVER4PI*cose*(2*(  (u+l)*atanh( m/(R1 + R2)) 
          +(v+m)*atanh( l/(R1 + R4))
          -    u*atanh( m/(R3 + R4))
          -    v*atanh( l/(R2 + R3))  ) - omega*d/sine );
  }
}

/*---------------------------------------------------------------------------------------- */
/*  The main working function MUTUAL() */
/*---------------------------------------------------------------------------------------- */
/* calculates mutual inductance of "filaments" with width and height */
/* as a combination of filament approximations                       */
double mutual(fil_j, fil_m)
FILAMENT *fil_j, *fil_m;
{
  double totalM;
  double dist, rj, rm;
  int parallel, whperp;
  int edge_par = FALSE; /* Edge parallel checker */
                              
  double widj[3], heightj[3];
  dist = dist_betw_fils(fil_j, fil_m, &parallel); 
  rj = MAX(fil_j->width, fil_j->height)/2.0; 
  rm = MAX(fil_m->width, fil_m->height)/2.0;  
  
  /* MUTUALFIL if filaments are far apart. */
  if (MAX(rj,rm)*100 < dist) { 
    num_mutualfil++;
    totalM = mutualfil(fil_j, fil_m);
  } else {    /* close by */
        if (parallel == 1) {

          /* PARALLEL, check for edges parallel */
          get_wid(fil_j, widj);
          get_height(fil_j, widj, heightj);
          edge_par = edges_parallel(fil_j,fil_m,widj,&whperp);
        }
        else {

          /* PERPENDICULAR, set to zero. */
          if ( fabs(vdotp(fil_j->lenvect,fil_m->lenvect))
            /(fil_j->length*fil_m->length) < EPS  )  {
            num_perp++;
            totalM = 0.0;
          }
        }

        /* EXACT_MUTUAL for close, edge-parallel fils */
        if (edge_par && 2*MAX(rj,rm)*10 > dist){

            num_exact_mutual++;
            totalM = parallel_fils(fil_j, fil_m, whperp, widj, heightj, dist);
        } else {
        
      /* SPLITFIL for close-by filaments */
      if (6*MAX(rj,rm) > dist) {
        num_quadFil++;
        totalM = splitfil(fil_j, fil_m);
      } else {
        /* FOURFIL for all other filaments */
        num_fourfil++;
        totalM = fourfil(fil_j, fil_m);
      }
    }
  }
  
  if (!finite(totalM))
        print_infinity_warning(fil_j, fil_m);
  return totalM;
}

/*---------------------------------------------------------------------------------------- */
/* MATLAB MEX INTERFACE */
/*---------------------------------------------------------------------------------------- */
#define XX0 0
#define YY0 1
#define ZZ0 2
#define XX1 3
#define YY1 4
#define ZZ1 5
#define DEBUG 1

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{   
	double *O, *L, *W, *H;
    double *Oin, *Lin, *Win, *Hin;
    double *out, *stats;
    double P[6], Q[6];
    int ind, ii;
    int n;
    FILAMENT fil1, fil2;
    int chkslf;

    int dispflag = DEBUG;
    
    if (nlhs>2)
    mexErrMsgTxt("Wrong number of output parameters, usage:  M = mutual(O, L, W, H)");
    if (nrhs!=4)
    mexErrMsgTxt("Wrong number of input parameters, usage:  M = mutual(O, L, W, H)");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
    mexErrMsgTxt("mutual: Input arguments must be double.");
 

    /* Prepare input */
    Oin = mxGetPr(prhs[0]); /* O - origin */
    Lin = mxGetPr(prhs[1]); /* L - length */
    Win = mxGetPr(prhs[2]); /* W - width */
    Hin = mxGetPr(prhs[3]); /* H - height */
    
    /* number of elements */
    n = (int) mxGetN(prhs[0]);
    
    /* Error check for m == 6 */
    if(!(6 == mxGetM(prhs[0]) &&
    		 6 == mxGetM(prhs[1]) &&
    		 6 == mxGetM(prhs[2]) &&
    		 6 == mxGetM(prhs[3])))
    		 mexErrMsgTxt("mutual: Must have six rows for each of the four inputs");
    if(n != mxGetN(prhs[1]) ||
       n != mxGetN(prhs[2]) ||
       n != mxGetN(prhs[3]))
    		mexErrMsgTxt("mutual: Must have same number of columns");   		   
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    
    /* Display output? */
    if (nlhs == 2) {
        plhs[1] = mxCreateDoubleMatrix(6, 1, mxREAL);
        stats = mxGetPr(plhs[1]);
    } else {
        dispflag = TRUE;
    }
            
    
    /* Reset counters */
    num_exact_mutual=0;
    num_fourfil=0;
    num_mutualfil=0;
    num_quadFil=0;
    num_found=0;
    num_perp=0;
    num_self=0;
    
    #ifdef _OPENMP
        omp_set_num_threads(omp_get_num_procs());  
    #endif
    
    #pragma omp parallel for default(none) shared(out,n,Oin,Lin,Win,Hin) private(O,L,W,H,fil1,fil2,chkslf,P,Q,ind,ii)
    for (ind=0;ind<n;ind++) {
          /* Increment pointers */
          O = Oin + 6*ind; L = Lin + 6*ind; W = Win + 6*ind; H = Hin + 6*ind;
          
          chkslf = 1 == 1;
          for (ii=0; ii<3; ii++) {
              chkslf = chkslf && nearzero(O[ii]-O[ii+3]);
              chkslf = chkslf && nearzero(L[ii]-L[ii+3]);
              chkslf = chkslf && nearzero(W[ii]-W[ii+3]);
              chkslf = chkslf && nearzero(H[ii]-H[ii+3]);
          }
          
          /* Convert to the FILAMENT format         */
          /* Do fil1 first. */
          /* length, width, height */
          fil1.length = sqrt(SQUARE(L[XX0])+SQUARE(L[YY0])+SQUARE(L[ZZ0]));
          fil1.width = sqrt(SQUARE(W[XX0])+SQUARE(W[YY0])+SQUARE(W[ZZ0]));
          fil1.height = sqrt(SQUARE(H[XX0])+SQUARE(H[YY0])+SQUARE(H[ZZ0]));

          if(chkslf) {
              out[ind] = self(fil1.width, fil1.length, fil1.height);
          } else {
          
          /* length, width, height */
          fil2.length = sqrt(SQUARE(L[XX1])+SQUARE(L[YY1])+SQUARE(L[ZZ1]));
          fil2.width = sqrt(SQUARE(W[XX1])+SQUARE(W[YY1])+SQUARE(W[ZZ1]));
          fil2.height = sqrt(SQUARE(H[XX1])+SQUARE(H[YY1])+SQUARE(H[ZZ1]));
          
          /* Make end points.  */
          /* Note that FastHenry defines the origin at the center of the  */
          /* square, but the OLWH notation defines the origin at the corner. */
          /* We define P to be the FastHenry origin and Q to be the end point. */
          for (ii = 0; ii < 6; ii++) {
              P[ii] = O[ii] + W[ii]/2 + H[ii]/2;  
              Q[ii] = O[ii] + L[ii] + W[ii]/2 + H[ii]/2;
          }

          fil1.x[0] = P[XX0]; fil1.y[0] = P[YY0]; fil1.z[0] = P[ZZ0];
          fil1.x[1] = Q[XX0]; fil1.y[1] = Q[YY0]; fil1.z[1] = Q[ZZ0];
          fil2.x[0] = P[XX1]; fil2.y[0] = P[YY1]; fil2.z[0] = P[ZZ1];
          fil2.x[1] = Q[XX1]; fil2.y[1] = Q[YY1]; fil2.z[1] = Q[ZZ1];

          /* lenvect */
          fil1.lenvect[XX] = L[XX0]; fil2.lenvect[XX] = L[XX1];
          fil1.lenvect[YY] = L[YY0]; fil2.lenvect[YY] = L[YY1];
          fil1.lenvect[ZZ] = L[ZZ0]; fil2.lenvect[ZZ] = L[ZZ1];

          /* Widthdir */
          fil1.widthdir[XX] = W[XX0] / fil1.width;
          fil1.widthdir[YY] = W[YY0] / fil1.width;
          fil1.widthdir[ZZ] = W[ZZ0] / fil1.width;

          fil2.widthdir[XX] = W[XX1] / fil2.width;
          fil2.widthdir[YY] = W[YY1] / fil2.width;
          fil2.widthdir[ZZ] = W[ZZ1] / fil2.width;
          /* Output */
          out[ind] = mutual(&fil1,&fil2);   
          }
          
            
    }
    /*More debug tools */
/*     if (dispflag){ */
/*         mexPrintf("Calls to exact_mutual: %15d\n",num_exact_mutual); */
/*         mexPrintf("         quadFil:      %15d\n",num_quadFil); */
/*         mexPrintf("         fourfils:     %15d\n",num_fourfil); */
/*         mexPrintf("         mutualfil:    %15d\n",num_mutualfil); */
/*         mexPrintf("Number found in table: %15d\n",num_found); */
/*         mexPrintf("Number perpendicular:  %15d\n",num_perp); */
/*         mexPrintf("Number self terms:     %15d\n",num_self); */
/*         mexPrintf("\n"); */
/*     } else { */
/*         stats[0] = num_exact_mutual; */
/*         stats[1] = num_fourfil; */
/*         stats[2] = num_mutualfil; */
/*         stats[3] = num_found; */
/*         stats[4] = num_perp; */
/*         stats[5] = num_self; */
/*     } */
    return;
}


