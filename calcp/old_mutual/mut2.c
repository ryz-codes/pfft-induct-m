#include <math.h>
#include <stdlib.h>
#include "mex.h"

#define FILS_PER_MESH 2         // filaments per mesh for breaking big loops.*/
#define PI 3.14159265358979323846
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define K 0.2235
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13
#define QUICK_EVAL 1          // do one filament approximations if == 1 */
#define MAXsubfils 6          // maximum points for Gaussian quadrature */
#define NUM_DIMS 10           //Number of distinguishing parms. table lookup*/

/* for deg_mutual, how small a dimension must be to be degenerate */
#define DEG_TOL 1.1e-4

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

#define nearzero(x) (fabs(x) < EPS)

#define TRUE (1==1)
#define FALSE (1==0)

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#define atanh(x) 0.5*log((1+x)/(1-x))
#define asinh(x) log(x+sqrt(x*x+1)) 
#define finite(x) (x != HUGE_VAL && x != -HUGE_VAL)

// COUNTERS
int num_mutualfil = 0;
int num_perp = 0;
int num_found = 0;
int num_fourfil = 0;
int num_exact_mutual = 0;
int num_self = 0;


typedef struct _alloc_list {
  char *ptr;
  struct _alloc_list *next;
} AllocList;

typedef struct _alloc_info {
  int size;
  int blocksize;      /* how many elements to allocate at once */
  int elems_left;
  char *next_elem;
  AllocList *head;
} AllocInfo;

/* stuff for mutual terms lookup table */
typedef struct _table {
  double val;
  struct _table *next_val;
  union {
    struct _table *next_dim;
    double *mut_term;
  } u;
} Table;

typedef struct Node {
   char *name;
   int number;
   int index;           /* internal number, for placement in A matrix */
   struct Node *equiv;  /* electically equivalence */
   double x, y, z;
   struct spath *to_end;
   int num_to_end;

   struct seglist *connected_segs;

   int type;                   /* NORMAL or (GPTYPE or GPHOLE) */

   struct Groundplane *gp;      /* CMS 6/7/92 ---- pointer to a groundplane */
   int s1, s2;                  /* indices into groundplane node array */
   struct tree *treeptr;
   char examined;      /* 1 = examined or never to be examined */
   int level;          /* number of nodes away from root of tree */
   //struct sseg_ptr pred;      /* predecessor, really just the branch of tree */

   struct g_nodes *gp_node;  /* node of nonuni gp that i really correspond to*/

   struct Node *next;
} NODES;

typedef struct Filament {
  double x[2], y[2], z[2];  // endpoints
  double length, area, width, height;
  double lenvect[3];        // vector along the length of filament
  int filnumber;
  struct Segment *segm;
  struct charge *pchg;      // 'charge' to send to multipole routines
} FILAMENT;

typedef struct Segment {
   char *name;
   double *widthdir;   //if width is not || to x-y plane and perpendicular to
                       // the length, then this is 3 element vector in      
                       // in the direction of width
   int number;         // an arbitrary number for the segment
   int type;    // CMS 8/21/92 -- type of structure the segment is in
   double length;      
   double area;        // area of cross section
   double width, height;  //width and height to cross section
   int hinc, winc;             // number of filament divisions in each dir
   NODES *node[2];                // nodes at the ends
   double sigma;              // conductivity
   double r_width, r_height;  //ratio of adjacent fil widths(see assignFil())
   int num_fils;               // hinc*winc
   FILAMENT *filaments;        // this segment's filaments
//   struct npath *conds;      // linked list of conductors which this seg is in
   struct pathlist *loops;   // loops in which this segment is a member
   int is_deleted;           // has this segment been used already

   //struct g_nodes *gp_node[2];/ nonuni_gp nodes that are really the ends

   struct _table *table;          // lookup table for mutual terms
   struct Segment *next;      // next segment in list
 } SEGMENT;

enum degen_type {brick = 0, flat = 1, skinny = 2, too_long = 3, too_short = 4, 
		   short_flat = 5, short_skinny = 6, impossible = 7};
//  
 
/* Calculates the minimum distance between two filaments */
/* Among other things, like Gaussian Quadrature weights */
//#include "induct.h"

/* this is where the Gaussian Quadrature are defined */
double **Gweight, **Gpoint;

double dist_between();
double min_endpt_sep();
double dist_betw_pt_and_fil();
double mutualfil(FILAMENT *fil1, FILAMENT *fil2);

double fourfil(FILAMENT *fil_j, FILAMENT *fil_m);
double parallel_fils(FILAMENT *fil_j, FILAMENT *fil_m,
    int whperp, double *x_j, double *y_j, double dist);

void savemat_mod();
double **MatrixAlloc();
char *MattAlloc();
void fillA();
void fillM();
void fillZ();
double resistance();
double selfterm();
double mutual();
double mutualfil();
double exact_mutual();
/*double magdiff();*/
double magdiff2();
double mut_rect();
extern int gmres();
extern int matvec();
extern int directmatvec();
extern int SetupComputePsi();
double fourfil();
double selfterm(FILAMENT *fil);
char *AllocAnEntry(AllocInfo *allocptr);

double brick_to_brick();
double parallel_fils();
double flat_to_flat_tape();
double flat_to_skinny_tape();
double do_tape_to_brick();
double tape_to_fil();  /* not implemented */
double brick_to_fil(); /* not implemented */


getD(fil, D)
FILAMENT *fil;
double *D;
{
  D[XX] = fil->x[1] - fil->x[0];
  D[YY] = fil->y[1] - fil->y[0];
  D[ZZ] = fil->z[1] - fil->z[0];
}

getr(x,y,z,s,t,D)
double *x,*y,*z;
double *s,t,*D;
{
  *x = s[XX] + t*D[XX];
  *y = s[YY] + t*D[YY];
  *z = s[ZZ] + t*D[ZZ];
}

double vdotp(v1, v2)
double *v1,*v2;
{
  return v1[XX]*v2[XX] + v1[YY]*v2[YY] + v1[ZZ]*v2[ZZ];
}


double magdiff2(fil1, node1, fil2, node2)
FILAMENT *fil1, *fil2;
int node1, node2;
{
   return ( SQUARE(fil1->x[node1] - fil2->x[node2])
	       +SQUARE(fil1->y[node1] - fil2->y[node2])
	       +SQUARE(fil1->z[node1] - fil2->z[node2])
	       );
}

/* this gives the mutual inductance of two filaments who represent */
/* opposite sides of a rectangle */

double mut_rect(len, d)
double len,d;
{
  double temp,temp1;

  temp = sqrt(len*len + d*d);
  temp1 = len*asinh(len/d) ;
  return temp - temp1;
}

/* returns the dotproduct of the vector from node0 to node1 of fil1 */
/* with that of fil2 */
double dotprod(fil1, fil2)
FILAMENT *fil1, *fil2;
{
  return(  (fil1->x[1] - fil1->x[0])*(fil2->x[1] - fil2->x[0])
	 + (fil1->y[1] - fil1->y[0])*(fil2->y[1] - fil2->y[0])
	 + (fil1->z[1] - fil1->z[0])*(fil2->z[1] - fil2->z[0]) );
}

double mag(x1,y1,z1)
double x1,y1,z1;
{
  return sqrt( x1*x1 + y1*y1 + z1*z1 );
}

double magsq(x1,y1,z1)
double x1,y1,z1;
{
  return ( x1*x1 + y1*y1 + z1*z1 );
}

double dotp(x1,y1,z1,x2,y2,z2)
double x1,y1,z1,x2,y2,z2;
{
  return x1*x2 + y1*y2 + z1*z2;
}

void fill_4(double *vec, double E, double a, double d)
{
  vec[0] = E - a;
  vec[1] = E + d - a;
  vec[2] = E + d;
  vec[3] = E;
}

double eval_eq(double x, double y, double z, double ref_len)
{
  static double one_60 = 1.0/60.0;
  static double one_6 = 1.0/6.0;
  static double one_24 = 1.0/24.0;
  double retval;
  double len, xsq, ysq, zsq;
  int num_nearzero;
  int num_nearzero_sq;
  double log_term(), tan_term();
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

double log_term(double x, double xsq, double ysq, double zsq, double len)
{
  double retval;
  retval = ((6*zsq - ysq)*ysq - zsq*zsq)*x*log( (x + len)/sqrt(ysq + zsq));
  return retval;
}

double tan_term(double x,double y,double z,double zsq,double len)
{
  double retval;
  retval =  x*y*z*zsq*atan(x*y/(z*len));
  return retval;
}

double dist_between(x1,y1,z1,x2,y2,z2)
double x1,y1,z1,x2,y2,z2;
{
  return sqrt(SQUARE(x1 - x2) + SQUARE(y1 - y2) + SQUARE(z1 - z2));
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
    //fprintf(stderr, "Internal err: dist_bet_pt_and_fil: why is t = %lg?\n", t);
    //exit1);
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

double aspectratio(fil)
FILAMENT *fil;
{
  double rat;

  rat = fil->height/fil->width;
  if (rat >= 1.0)
    return rat;
  else
    return 1.0/rat;
}

// fill_Gquad()
// {
//   int i,j;
//   
//   
//   Gweight = (double **)MattAlloc(MAXsubfils+1, sizeof(double *));
//   Gpoint = (double **)MattAlloc(MAXsubfils+1, sizeof(double *));
//   
//   for(i = 1; i <= MAXsubfils; i++) {
//     Gweight[i] = (double *)MattAlloc(i, sizeof(double));
//     Gpoint[i] = (double *)MattAlloc(i, sizeof(double));
//   }
// 
//   Gweight[1][0] = 2.0;
//   Gpoint[1][0] = 0.0;
// 
//   for(i = 2; i <= MAXsubfils; i++)
//     /* subtract 1 from pointers since this function starts with p[1] */
//     gquad_weights(i,Gpoint[i] - 1, Gweight[i] - 1);
// }

findnfils(fil, subfils, nfils)
FILAMENT *fil, *subfils;
int nfils;
{
  double hx,hy,hz,mag,wx,wy,wz,dx,dy,dz;
  int i,j;

  if (fil->segm->widthdir != NULL) {
    wx = fil->segm->widthdir[XX];
    wy = fil->segm->widthdir[YY];
    wz = fil->segm->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if ( fabs(wx/fil->segm->length) < EPS && fabs(wy/fil->segm->length) < EPS) {
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

  if (fil->width > fil->height) {
    dx = fil->width*wx/2;
    dy = fil->width*wy/2;
    dz = fil->width*wz/2;
  }
  else {
    dx = fil->height*hx/2;
    dy = fil->height*hy/2;
    dz = fil->height*hz/2;
  }

  /* all mutualfil needs are the filament coordinates and length */
  for (j = 0; j < nfils; j++) {
    for(i = 0; i < 2; i++) {
      subfils[j].x[i] = fil->x[i] + dx*Gpoint[nfils][j];
      subfils[j].y[i] = fil->y[i] + dy*Gpoint[nfils][j];
      subfils[j].z[i] = fil->z[i] + dz*Gpoint[nfils][j];
    }
    subfils[j].length = fil->length;
  }
}

double dist_betw_fils(fil1, fil2, parallel)
FILAMENT *fil1, *fil2;
int *parallel;
{
  double k1,k2,k3,k4,c1,c2,c3,c4,a1,a2,a3,b1,b2,b3;
  double x1,y1,z1,x2,y2,z2;
  double D1[3], D2[3], s1[3], s2[3], e[3], *D, t1, t2, s1ms2[3], s1me[3];
  double D1D1, D1D2, D2D2, D1s1s2, D2s1s2;
  double tmp1, tmp2;
  
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


#if 1==0
/* main for testing if these work */
main()
{
  FILAMENT fil1, fil2;

  while(1) {
    printf("fil1 1: ");
    scanf("%lf %lf %lf", &(fil1.x[0]), &(fil1.y[0]), &(fil1.z[0]));
    printf("fil1 2: ");
    scanf("%lf %lf %lf", &(fil1.x[1]), &(fil1.y[1]), &(fil1.z[1]));
    
    printf("fil2 1: ");
    scanf("%lf %lf %lf", &(fil2.x[0]), &(fil2.y[0]), &(fil2.z[0]));
    printf("fil2 2: ");
    scanf("%lf %lf %lf", &(fil2.x[1]), &(fil2.y[1]), &(fil2.z[1]));

    printf("dist = %lg\n",dist_betw_fils(&fil1, &fil2));
  }
}
#endif

/* this file contains the functions for exact calculation of the
   self-inductance of a rectangular bar and the mutual inductance of
   two parallel rectangular bars.

   Also, it contains the code for the lookup table for these inductances */

/* this counts on all arguments being positive! */
#define compare(x,y,eps) (  (((x)==0 && (y)==0) || (fabs((x) - (y)) < eps*((x) + (y)) )) \
  ? 0 : (  ((x) > (y)) ? 1 : -1 )  )

#define nearzero(x) (fabs(x) < EPS)
  
/* self inductance */
double self(W,L,T)
double W,L,T; 
{

    double w,t,aw,at,ar,r, z;
    //double asinh(), atan(), sqrt(); 

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
    
    return z; 

}

/* This assumes the lengths of the fils are parallel and returns 1 if the 
   side faces are parallel also */
edges_parallel(fil_j, fil_m, wid1, whperp)
FILAMENT *fil_j, *fil_m;
int *whperp;
double *wid1;
{
  double *wid_j = fil_j->segm->widthdir;
  double *wid_m = fil_m->segm->widthdir;
  double wid2[3];
  double prod,mj,mm;

  if (wid_j == NULL && wid_m == NULL) {
    /* both have unspecified width directions and their lengths are assumed
       parallel, so they are parallel */
    *whperp = 0;
    return 1==1;
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
      return 1==1;
    }
    else if (fabs( fabs(prod) - 1 ) < EPS) {
      *whperp = 0;
      return 1==1;
    }
    else
      return 1==0;
  }
}

/* calculates direction of width if not specified */
get_wid(fil, wid)
FILAMENT *fil;
double *wid;
{

  double wx,wy,wz;
  double mag;

  if (fil->segm->widthdir != NULL) {
    wx = fil->segm->widthdir[XX];
    wy = fil->segm->widthdir[YY];
    wz = fil->segm->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if (fabs(wx/fil->segm->length) < EPS && fabs(wy/fil->segm->length) < EPS) {
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
get_height(fil, wid, height)
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
    
/* the following is code for the lookup table for common mutual terms */

/* The lookup table for gp's will be a global variable. Yuck! */

/* The lookup table head */
Table *table = NULL;
AllocInfo table_alloc;
AllocInfo double_alloc;
int table_init = 0;

#define TABLE_SMALL OFF


/* lookup mutual term in table */

int lookup (FILAMENT *fil_j, FILAMENT *fil_m,
int whperp,
double *widj, double *heightj,   /* width and height vectors */
double *retval,
double *dims,
Table ***lastptr,
int *dim_count, int *p_num_dims)
{
  int num_dims = NUM_DIMS;
  Table **s_table;

  if (whperp == 1) {
    *lastptr = NULL;
    return 0;
  }

//#if TABLE_SMALL == ON
//
//  /* this was an attempted shorter table but works worse than original! */
//
//  /* for ground planes, look in global table */
//
//  if (fil_j->segm->node[0]->gp != NULL && fil_m->segm->node[0] != NULL) {
//    num_dims = *p_num_dims = 9;
//    fill_dims(fil_j, fil_m, widj, heightj, dims,num_dims);
//    s_table = &table;  /* global table */
//  }
//  else if (fil_j->segm == fil_m->segm) {
//    /* if fils on same seg, look up individual segment tables */
//    num_dims = *p_num_dims = 6;
//    fill_dims_seg(fil_j, fil_m, widj, heightj, dims, num_dims);
//    s_table = &fil_j->segm->table;
//  }
//  else {
//    /* otherwise don't look for it */
//    *lastptr = NULL;
//    return 0;
//  }
//#else
    /* this is the original lookup table which works better! */
  num_dims = *p_num_dims = 10;
  fill_dims(fil_j, fil_m, widj, heightj, dims,num_dims);
  s_table = &table;  /* global table */

//#endif

  if (find_dims(dims, num_dims, s_table, retval,
		dim_count, lastptr) == 0)
    return 0;
  else {
    if (dotprod(fil_j,fil_m) < 0)
      *retval *= -1;
    return 1;
  }
    
}

/* this fills the vector dims with the dimension information for fil_j
   and fil_m to later determine if the pair has been previous computed and
   is in the lookup table.  All dims should be positive!! */
fill_dims(fil_j, fil_m, widthj, heightj, dims,num_dims)
FILAMENT *fil_j, *fil_m;
double *dims, *widthj, *heightj;
{
  int j_first;
  int is_same;
  int i = 0;
  double x,y,z;
  double z_j[3],length;

#if TABLE_SMALL != ON
  if (fil_j->segm->node[0]->gp != NULL && fil_m->segm->node[0]->gp != NULL)
    dims[i++] = 1.0;
  else
    dims[i++] = 0.0;
#endif

  length = fil_j->length;
  /* vector along length */
  z_j[XX] = fil_j->lenvect[XX]/length;
  z_j[YY] = fil_j->lenvect[YY]/length;
  z_j[ZZ] = fil_j->lenvect[ZZ]/length;
  
  x = (fil_j->x[0]+fil_j->x[1]) - (fil_m->x[0]+fil_m->x[1]);
  y = (fil_j->y[0]+fil_j->y[1]) - (fil_m->y[0]+fil_m->y[1]);
  z = (fil_j->z[0]+fil_j->z[1]) - (fil_m->z[0]+fil_m->z[1]);

  /* compute center to center distances in fil_j coordinates */
  /* (putting height component first should sort the information by */
  /* plane (for multiple planes only)) */
/*
  dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ],x,y,z));
  dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ],x,y,z));
  dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ],x,y,z));
*/

  if ( (is_same = compare(fil_j->height, fil_m->height, EPS)) == -1)
    j_first = 1;
  else if (is_same == 1)
    j_first = 0;
  else if ( (is_same = compare(fil_j->width, fil_m->width, EPS)) == -1)
    j_first = 1;
  else if (is_same == 1)
    j_first = 0;
  else if ( (is_same = compare(fil_j->length, fil_m->length, EPS)) == -1)
    j_first = 1;
  else if (is_same == 1 || is_same == 0)
    j_first = 0;
    
  if (j_first == 1) {
    dims[i++] = fil_j->height;
    dims[i++] = fil_j->width;
    dims[i++] = fil_j->length;
    dims[i++] = fil_m->height;
    dims[i++] = fil_m->width;
    dims[i++] = fil_m->length;
  }
  else {
    dims[i++] = fil_m->height;
    dims[i++] = fil_m->width;
    dims[i++] = fil_m->length;
    dims[i++] = fil_j->height;
    dims[i++] = fil_j->width;
    dims[i++] = fil_j->length;
  }

  dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ],x,y,z));
  dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ],x,y,z));
  dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ],x,y,z));

  if (i != num_dims) {
    //fprintf(stderr, "bad number for num_dims %d\n",num_dims);
    //exit1);
  }
}
/* this fills the vector dims with the dimension information for fil_j
   and fil_m to later determine if the pair has been previous computed and
   is in the lookup table.  It is for fils on the same seg only. */
fill_dims_seg(fil_j, fil_m, widthj, heightj, dims,num_dims)
FILAMENT *fil_j, *fil_m;
double *dims, *widthj, *heightj;
{
  int j_first;
  int is_same;
  int i = 0;
  double x,y,z;
  double z_j[3],length;

  x = (fil_j->x[0]+fil_j->x[1]) - (fil_m->x[0]+fil_m->x[1]);
  y = (fil_j->y[0]+fil_j->y[1]) - (fil_m->y[0]+fil_m->y[1]);
  z = (fil_j->z[0]+fil_j->z[1]) - (fil_m->z[0]+fil_m->z[1]);
  dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ],x,y,z));
  dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ],x,y,z));

  /* should always be zero for fils on same seg */
  /*  dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ],x,y,z)); */

  /* choose which fil to put first based on height, then width */
  if ( (is_same = compare(fil_j->height, fil_m->height, EPS)) == -1)
    j_first = 1;
  else if (is_same == 1)
    j_first = 0;
  else if ( (is_same = compare(fil_j->width, fil_m->width, EPS)) == -1)
    j_first = 1;
  else if (is_same == 1 || is_same == 0)
    j_first = 0;
    
  if (j_first == 1) {
    dims[i++] = fil_j->height;
    dims[i++] = fil_j->width;
    dims[i++] = fil_m->height;
    dims[i++] = fil_m->width;
  }
  else {
    dims[i++] = fil_m->height;
    dims[i++] = fil_m->width;
    dims[i++] = fil_j->height;
    dims[i++] = fil_j->width;
  }

  if (i != num_dims) {
    //fprintf(stderr, "bad number for num_dims %d\n",num_dims);
    //exit1);
  }
}
  
/*
int compare(x, y, eps)
double x,y,eps;
{
  if ( (x==0 && y==0) || (fabs(x - y)/fabs(x+y) < eps))
    return 0;
  else if (x > y)
    return 1;
  else
    return -1;
}
*/

find_dims(dims, num_dims, a_table, retval, ret_dim_count, ret_lastptr)
double *dims;
int num_dims;
double *retval;
int *ret_dim_count;
Table ***ret_lastptr, **a_table;
{
  Table *entry, **lastptr;
  int is_same, maybe_its_there;
  int dim_count, i;

  maybe_its_there = TRUE;
  dim_count = 0;
  /*  entry = table;  lastptr = &table;  the old code */
  entry = *a_table;
  lastptr = a_table;
  while(entry!=NULL && maybe_its_there == TRUE) {
    if ( (is_same = compare(dims[dim_count], entry->val, EPS)) == 1)
      /* not found */
      maybe_its_there = FALSE;
    else if (is_same == -1) {
      lastptr = &(entry->next_val);
      entry = entry->next_val;
    }      
    else {
      if (dim_count < num_dims - 1) {
	dim_count++;
	lastptr = &(entry->u.next_dim);
	entry = entry->u.next_dim;
      }
      else {
	*retval = *(entry->u.mut_term);
	   /*
	    printf("Found!    ");
	    for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
	    printf("%13.6lg\n",*retval);
	   */
	return 1;
      }
    }
  }

  /* info for put_in_table to quickly insert in table */
  (*ret_lastptr) = lastptr;
  (*ret_dim_count) = dim_count;
      /*
	printf("Not Found!");
	for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
        printf("\n");
      */
  return 0;
}
      
put_in_table(fil_j, fil_m, whperp, mutterm, dims, dim_count, lastptr, num_dims)
FILAMENT *fil_j, *fil_m;
int whperp;
double mutterm;
double *dims;
int dim_count, num_dims;
Table **lastptr;
{
  Table *entry;
  int i;

  if (lastptr == NULL)
    return;

  entry = (Table *)AllocAnEntry(&table_alloc);
  entry->next_val = (*lastptr);
  (*lastptr) = entry;
  entry->val = dims[dim_count++];
  lastptr = &(entry->u.next_dim);
  
  for(i = dim_count; i < num_dims; i++) {
    entry = (Table *)AllocAnEntry(&table_alloc);
    (*lastptr) = entry;
    entry->val = dims[i];
    entry->next_val = NULL;
    lastptr = &(entry->u.next_dim);
  }

  entry->u.mut_term = (double *)AllocAnEntry(&double_alloc);
  *(entry->u.mut_term) = fabs(mutterm);

  
  printf("put in:   ");
  for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
  printf("%13.6lg\n",fabs(mutterm));
  
}

/* no logic, just a nice number */
#define ALLOCBLOCK 256 

/* initialize info for allocating table so we can free it later */
init_table()
{
  table_alloc.size = sizeof(Table);
  table_alloc.blocksize = ALLOCBLOCK;
  table_alloc.elems_left = 0;
  table_alloc.head = NULL;
  double_alloc.size = sizeof(double);
  double_alloc.blocksize = ALLOCBLOCK;
  double_alloc.elems_left = 0;
  double_alloc.head = NULL;
  table_init = 1;
  mexPrintf("Table initiated.\n");
}

get_table_mem()
{
  return MemoryForEntries(&table_alloc) + MemoryForEntries(&double_alloc);
}

destroy_table()
{
  DestroyEntries(table_alloc);
  DestroyEntries(double_alloc);
}


/* allocates table entries in blocks for more efficient memory */
char *AllocAnEntry(allocptr)
AllocInfo *allocptr;
{
  int blocksize, size;
  AllocList *elem;

  if (allocptr->elems_left > 0) {
    allocptr->elems_left--;
    return (allocptr->next_elem += allocptr->size);
  }
  else {
    blocksize = allocptr->blocksize;
    elem = (AllocList *) mxMalloc(sizeof(AllocList));
    if (elem == NULL) {
        mexPrintf("Fatal error: AllocList failed");
        exit;
    }
    elem->next = allocptr->head;
    allocptr->head = elem;
    elem->ptr = allocptr->next_elem = 
                        (char *) mxMalloc(allocptr->size*blocksize);
    allocptr->elems_left = blocksize - 1;
    if (allocptr->next_elem == NULL) {
      //fprintf(stderr, "Out of memory in AllocAnEntry. size = %d, block = %d\n",
	  //    allocptr->size, blocksize);
      //exit1);
    }
    return allocptr->next_elem;
  }
}

DestroyEntries(allocinfo)
AllocInfo allocinfo;
{
  AllocList *lastelem, *head;

  head = allocinfo.head;

  while(head != NULL) {
    lastelem = head;
    free(head->ptr);
    head = head->next;
    free(lastelem);
  }

  allocinfo.head = NULL;
  allocinfo.elems_left = 0;
}

MemoryForEntries(allocptr)
AllocInfo *allocptr;
{
  AllocList *entry;
  int count = 0;

  for(entry = allocptr->head; entry != NULL; entry = entry->next)
    count++;

  return count*allocptr->blocksize*allocptr->size;
}

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
  static double one_3 = 1.0/3.0;
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
  int num_nearzero_sq;
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

/*  retval += (0.5*zsq - one_6*ysq)*y*log(x + len)
                + (0.5*zsq - one_6*xsq)*x*log(y + len)
		  + x*y*z*log(z + len); 
*/

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

double tape_to_fil(E,a,P,l3,l1,l2)
double E,a,P,l3,l1,l2;
{
  /* I have not implemented this degenerate case.  It should be done
     by fourfil() and never get to this point */
  //fprintf(stderr,"Hey, tape_to_fil should not have been called!\n");
  //exit1);
}

double brick_to_fil(E,a,P,b,l3,l1,l2)
double E,a,P,b,l3,l1,l2;
{
  /* I have not implemented this degenerate case.  It should be done
     by fourfil() and never get to this point */
  //fprintf(stderr,"Hey, brick_to_fil should not have been called!\n");
  //exit1);
}

#define LEN 4
#define WID 2
#define HEIGHT 1

enum degen_type find_deg_dims(fil)
FILAMENT *fil;
{
  double max;
  
  max = MAX(fil->length, fil->width);
  max = MAX(max, fil->height);

  return (fil->length/max < DEG_TOL)*LEN + (fil->width/max < DEG_TOL)*WID
         + (fil->height/max < DEG_TOL)*HEIGHT;
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
  double *nx_j, *ny_j;

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
  else if ( deg_m == brick && (deg_j == flat || deg_j == skinny)
	   || deg_j == brick && (deg_m == flat || deg_m == skinny))
    return do_tape_to_brick(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);
  else if ( deg_j == too_long && deg_m == too_long)
    return fourfil(fil_j, fil_m);
  else if (deg_j == too_long || deg_j == too_long)
    return fourfil(fil_j, fil_m);
  else
    return fourfil(fil_j, fil_m);

}

setup_tape_to_tape(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m,
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

double do_tape_to_brick(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m)
FILAMENT *fil_j, *fil_m;
int whperp;
double *x_j, *y_j;  /* unit vectors in the fil coord sys */
enum degen_type deg_j, deg_m;
{

  FILAMENT nfil_j, nfil_m;
  double *nx_j, *ny_j, *dR;
  double wid_brick[3], hei_brick[3], orig_x[2], orig_y[2], orig_z[2];
  double x_flat[3], y_flat[3];
  double small_dim, sum;
  int i,j,gpoints;
  extern double **Gweight, **Gpoint;    /* gaussian quad weights. */
  enum degen_type ndeg_j, ndeg_m;

/*
  if ( deg_m == brick && (deg_j == flat || deg_j == skinny)
	   || deg_j == brick && (deg_m == flat || deg_m == skinny))
    return do_tape_to_brick(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);
*/

  if (deg_j == flat) {
    nfil_j = *fil_j;
    nfil_m = *fil_m;
    nx_j = x_j;
    ny_j = y_j;
    get_wid(fil_m,wid_brick);
    get_height(fil_m,wid_brick,hei_brick);
  }
  else if (deg_j == skinny) {
    /* turn skinny into flat orientation */
    nfil_j = *fil_j;
    nfil_m = *fil_m;
    /* swap coord sys */
    ny_j = x_j;
    nx_j = y_j;
    /* swap height and width */
    nfil_j.width = fil_j->height;
    nfil_j.height = fil_j->width;
    nfil_m.width = fil_m->height;
    nfil_m.height = fil_m->width;
    /* get them swapped */
    get_wid(fil_m,hei_brick);
    get_height(fil_m,hei_brick,wid_brick);
  }
  else if (deg_j == brick) {
    /* swap j and m */
    nfil_j = *fil_m;
    nfil_m = *fil_j;
    get_wid(fil_m,x_flat);
    get_height(fil_m,x_flat,y_flat);

    if (deg_m == flat) {
      nx_j = x_flat;
      ny_j = y_flat;
      for(i = 0; i < 3; i++) {
        wid_brick[i] = x_j[i];
        hei_brick[i] = y_j[i];
      }
    
    }
    else {
      nx_j = y_flat;
      ny_j = x_flat;
      nfil_j.width = fil_m->height;
      nfil_j.height = fil_m->width;
      nfil_m.width = fil_j->height;
      nfil_m.height = fil_j->width;
      for(i = 0; i < 3; i++) {
        wid_brick[i] = y_j[i];
        hei_brick[i] = x_j[i];
      }
    
    }
  }

  /* store original brick position */
  for(i = 0; i < 2; i++) {
    orig_x[i] = nfil_m.x[i];
    orig_y[i] = nfil_m.y[i];
    orig_z[i] = nfil_m.z[i];
  }
  
  if (nfil_m.width > nfil_m.height) {
    /* the height direction will be done discretely */
    small_dim = nfil_m.height/2;
    nfil_m.height = 0;
    dR = hei_brick;
    if (whperp == 0)   /* useful for testing only. if forced == 1 */
      ndeg_m = flat;
    else 
      ndeg_m = skinny;
  }
  else {
    /* the width direction will be done discretely */
    small_dim = nfil_m.width/2;
    nfil_m.width = 0;
    dR = wid_brick;
    if (whperp == 0)
      ndeg_m = skinny;
    else
      ndeg_m = flat;
  }

  /* if forced == 1, then setting ndeg_j matters */
  ndeg_j = flat;
  nfil_j.height = 0.0;   /* insure we use the middle of filament x-section*/
    
  gpoints = 3;
  /* now do gaussian quadrature of tape_to_tape to approximate */
  sum = 0;
  for(i = 0; i < gpoints; i++) {
    for(j = 0; j < 2; j++) {
      nfil_m.x[j] = orig_x[j] + dR[XX]*small_dim*Gpoint[gpoints][i];
      nfil_m.y[j] = orig_y[j] + dR[YY]*small_dim*Gpoint[gpoints][i];
      nfil_m.z[j] = orig_z[j] + dR[ZZ]*small_dim*Gpoint[gpoints][i];
    }
    sum += Gweight[gpoints][i]*exact_mutual(&nfil_j, &nfil_m, whperp, 
					    nx_j, ny_j, ndeg_j, ndeg_m);
  }

  return sum/2.0;
}


/* exact mutual inductance based on C. Hoer and C.Love, 
   Journal of the National Bureau of Standards-C,  Vol. 69C, p 127-137, 1965.*/
double exact_mutual(FILAMENT *fil_j, FILAMENT *fil_m,
int whperp,
double *x_j, double *y_j,  /* unit vectors in the fil coord sys */
int deg_j, int deg_m)
{
  double z_j[3];  /* unit vectors in the filament coord sys*/
  double origin[3];
  double dumb, ox,oy,oz, length;
  double a,b,c,d,l1,l2,l3,E,P, l3_1;
  double endx, endy, endz;
  int sign, sign2;
  double q[4], r[4], s[4];
  double totalM, eval_eq;
  int i,j,k;
  int a_deg, b_deg, c_deg, d_deg;
  int forced;

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

  l1 = fil_j->segm->length;
  l2 = fil_m->segm->length;

  if ( fabs(fabs(l3 - l3_1) - l2)/l2 > EPS) {
    //fprintf(stderr, "Huh?  filament not expected length\n");
    //exit1);
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

  if (forced) {
    if (deg_j == brick && deg_m == brick) {a_deg=d_deg=b_deg=c_deg=0;}
       //fprintf(stderr,"fdbb");}
    if (deg_j == flat && deg_m == flat) {a_deg=d_deg=0; b_deg=c_deg=1;}
       //fprintf(stderr,"fdff ");}
    if (deg_j == flat && deg_m == skinny) {a_deg=c_deg=0; b_deg=d_deg=1;}
       //fprintf(stderr,"fdfs ");}
    if (deg_j == skinny && deg_m == flat) {a_deg=c_deg=0; b_deg=d_deg=1;}
       //fprintf(stderr,"fdsf ");}
    if (deg_j == skinny && deg_m == skinny) {a_deg=d_deg=0; b_deg=c_deg=1;}
       //fprintf(stderr,"fdss ");}
  }


  if (a_deg && b_deg && c_deg && d_deg) {
    // two long filaments 
    totalM = fourfil(fil_j, fil_m)/(sign*MUOVER4PI);
  }
  else if (!a_deg && b_deg && c_deg && d_deg) {
    // one flat and one long 
    totalM = tape_to_fil(E + d/2, a, P - b/2 + c/2,l3,l1,l2) / a;
  }
  else if (!a_deg && b_deg && c_deg && !d_deg) {
    // two flat 
    totalM = flat_to_flat_tape(E,a,d,P - b/2 + c/2 ,l3,l1,l2) / (a * d);
  }
  else if (!a_deg && b_deg && !c_deg && d_deg) {
    // one flat and one skinny 
    totalM = flat_to_skinny_tape(E + d/2,a,P - b/2 ,c,l3,l1,l2) / (a * c);
  }
  else if (!a_deg && !b_deg && c_deg && d_deg) {
    totalM = brick_to_fil(E + d/2,a,P + c/2,b,l3,l1,l2) / (a * b);
  }
  else if (deg_j != brick || deg_m != brick) {
    //fprintf(stderr,"Internal Error: Bad Degenerate filament a/l1<tol\n");
    //fprintf(stderr,"Using fourfil() instead...\n");
    totalM = fourfil(fil_j, fil_m)/(sign*MUOVER4PI);
  }
  else {
    // all that are left are bricks 
    // this is the full 3-D filament calculation, no degeneracies
    totalM = brick_to_brick(E,a,d,P,b,c,l3,l1,l2)/ (a * b * c * d);}
  return sign*MUOVER4PI*totalM;
}
    


/* these are missing in some math.h files */
//extern double asinh();
//extern double atanh();

int realcos_error = 0;



print_infinity_warning(fil1, fil2)
     FILAMENT *fil1, *fil2;
{
  FILAMENT *fil;

  mexPrintf("Severe warning: mutual inductance = infinity for two filaments:\n");
  //fprintf(stderr,"Severe warning: mutual inductance = infinity for two filaments:\n");

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

/* calculates selfinductance of rectangular filament */
/* it uses an exact expression for the 6-fold integral from Ruehli */
/* (the actual function comes is in file joelself.c */
double selfterm(fil)
FILAMENT *fil;
{
   double self();
   double approx, joelself;

/*   approx = fil->length*MUOVER4PI*2
               *(log(2*fil->length/(K*(fil->width + fil->height))) - 1); */
   joelself = MU0*self(fil->width, fil->length, fil->height);
/*   printf("Joel's function: %lg,  my approx: %lg\n",joelself, approx); */
   num_self += 1;
   
   return joelself; 
}



double fourfil(fil_j, fil_m)
FILAMENT *fil_j, *fil_m;
{
  FILAMENT subfilj[MAXsubfils], subfilm[MAXsubfils];
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
  
  /* printf("4: %14.8le ",totalM); */

  return totalM;
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

    /* //fprintf(stderr,"  degenerate: fil %d to fil %d: %13.6lg\n",
           fil_j->filnumber,
	   fil_m->filnumber, 
	    compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j,
				  deg_j, deg_m, dist));
     */

  }

}

findfourfils(fil, subfils)
FILAMENT *fil, subfils[4];
{
  double hx,hy,hz,mag,wx,wy,wz;
  int i;

  if (fil->segm->widthdir != NULL) {
    wx = fil->segm->widthdir[XX];
    wy = fil->segm->widthdir[YY];
    wz = fil->segm->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if ( fabs(wx/fil->segm->length) < EPS && fabs(wy/fil->segm->length) < EPS) {
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



/* calculates the mutual inductance between two filaments */
/* from Grover, Chapter 7 */
double mutualfil(fil1, fil2)
FILAMENT *fil1, *fil2;
{
  double R, R1, R2, R3, R4, l, m;
  double cose, sine, u,v,d;
  double alpha, tmp1, tmp2, tmp3, tmp4, tmp5;
  double maxR, maxlength, sinsq, blah, minR;
  double scaleEPS, realcos;
  double vtemp;

  double omega,M;
  int signofM;

  double junk;

  /* for parallel filaments */
  double ux, uy, uz, Rx, Ry, Rz, x1_0, x1_1, x2_0, x2_1, magu;
  double dx, dy, dz, dotp, vx, vy, vz;
  double R1sq, R2sq, R3sq, R4sq, m2, l2, u2, v2, alpha2;

  R1sq = magdiff2(fil1,1,fil2,1);
  R2sq = magdiff2(fil1,1,fil2,0);
  R3sq = magdiff2(fil1,0,fil2,0);
  R4sq = magdiff2(fil1,0,fil2,1);
  R1 = sqrt(R1sq);
  R2 = sqrt(R2sq);
  R3 = sqrt(R3sq);
  R4 = sqrt(R4sq);

  maxR = minR = R1;
  if (R2 > maxR) maxR = R2;
  if (R3 > maxR) maxR = R3;
  if (R4 > maxR) maxR = R4;

  if (R2 < minR) minR = R2;
  if (R3 < minR) minR = R3;
  if (R4 < minR) minR = R4;

  l = fil1->length;
  m = fil2->length;
  maxlength = (l > m ? l : m);

  scaleEPS = minR/maxlength*10;
  if (scaleEPS < 1) scaleEPS = 1;
  if (scaleEPS > 100) scaleEPS = 100;

  alpha = R4sq - R3sq + R2sq - R1sq;
  signofM = 1;

  /* segments touching */
  if ( (fabs(R1) < EPS)||(fabs(R2) < EPS)||(fabs(R3) < EPS)||(fabs(R4) < EPS) )
    {
      if (fabs(R1) < EPS)  R = R3;
      else if (fabs(R2) < EPS)  R = R4; 
      else if (fabs(R3) < EPS)  R = R1;
      else R = R2;

      M = MUOVER4PI*2*(dotprod(fil1,fil2)/(l*m))
	   *(l*atanh(m/(l+R)) + m*atanh(l/(m+R)));
      /* note: dotprod should take care of signofM */

      return M;
    }

  cose = alpha/(2*l*m);
  if (fabs(cose) > 1) cose = (cose < 0 ? -1.0 : 1.0);
  blah = 1.0 - fabs(cose);

  /* let's use the real cosine */
  realcos = dotprod(fil1, fil2)/(l*m);
  /*realcos = dotprod(fil1, fil2)/(fil1->length*fil2->length);*/

  /* Segments are perpendicular! */
  if (fabs(realcos) < EPS)
    return 0.0;

  if (fabs((realcos - cose)/cose) > 0.1) 
    if (realcos_error == 0) {
      //fprintf(stderr, "Internal Warning: realcos = %lg,  cose = %lg\n",realcos, cose);
      //fprintf(stderr,"  This may be due to two filaments that are separated \n\
by a distance 1e10 times their length\n");
      realcos_error = 1;
    }

  cose = realcos;

  /* filaments parallel */
  tmp1 = fabs( fabs(cose) - 1);
/*  if ( fabs( fabs(cose) - 1) < scaleEPS*EPS*10.0) { */
  if ( fabs( fabs(cose) - 1) < EPS) {
    /* determine a vector in the direction of d with length d */
    /* (d is the distance between the lines made by the filament */
    Rx = fil2->x[0] - fil1->x[0];  /* vector from fil1 to fil2 */
    Ry = fil2->y[0] - fil1->y[0];
    Rz = fil2->z[0] - fil1->z[0];
    ux = fil1->x[1] - fil1->x[0];
    uy = fil1->y[1] - fil1->y[0];
    uz = fil1->z[1] - fil1->z[0];
    magu = sqrt(ux*ux + uy*uy + uz*uz);
    ux = ux/magu;    /* unit vector in direction of fil1 */
    uy = uy/magu;    
    uz = uz/magu;

    dotp = ux*Rx + uy*Ry + uz*Rz;  /* component of R in direction of fil1 */

    /* d vector is R vector without its component in the direction of fils */
    dx = Rx - dotp*ux;
    dy = Ry - dotp*uy;
    dz = Rz - dotp*uz;
    d = sqrt(dx*dx + dy*dy + dz*dz);

    /* let fil1 be the x axis, with node 0 being origin and u be */
    /* its positive direction */
    x1_0 = 0;
    x1_1 = l;
    
    /* x2_0 = dotprod( fil2.node0 - (fil1.node0 + d), u ) */
    /* (dotproduct just gives it correct sign) */
    vx =  (fil2->x[0] - ( fil1->x[0] + dx));
    vy =  (fil2->y[0] - ( fil1->y[0] + dy));
    vz =  (fil2->z[0] - ( fil1->z[0] + dz));
    x2_0 = vx*ux + vy*uy + vz*uz;
    vtemp = sqrt(vx*vx + vy*vy + vz*vz);

    /* same thing for x2_1 */
    vx =  (fil2->x[1] - ( fil1->x[0] + dx));
    vy =  (fil2->y[1] - ( fil1->y[0] + dy));
    vz =  (fil2->z[1] - ( fil1->z[0] + dz));
    x2_1 = vx*ux + vy*uy + vz*uz;

    if ( fabs( (sqrt(vx*vx + vy*vy + vz*vz) - fabs(x2_1))
	      /(MAX(fabs(x2_0)+d,fabs(x2_1)+d))) 
	> EPS) {
      //printf("uh oh, segs don't seem parallel %lg\n",(sqrt(vx*vx + vy*vy * vz*vz) - fabs(x2_1)));
    }

    if ( fabs( (vtemp - fabs(x2_0))/(MAX(fabs(x2_0)+d,fabs(x2_1)+d))) 
	> EPS) {
      //printf("uh oh, segs don't seem parallel\n");
    }

    if ( fabs(d) < EPS )  { /* collinear! */
      M = MUOVER4PI*(fabs(x2_1 - x1_0)*log(fabs(x2_1 - x1_0))
	             - fabs(x2_1 - x1_1)*log(fabs(x2_1 - x1_1))
		     - fabs(x2_0 - x1_0)*log(fabs(x2_0 - x1_0))
		     + fabs(x2_0 - x1_1)*log(fabs(x2_0 - x1_1)) );
      return M;
    }  /* end collinear */

    M = MUOVER4PI*(mut_rect(x2_1 - x1_1,d) - mut_rect(x2_1 - x1_0,d) 
		 - mut_rect(x2_0 - x1_1,d) + mut_rect(x2_0 - x1_0,d) );

    return M;
  }  /* end if parallel filaments */

 /* the rest if for arbitrary filaments */

  l2 = l*l;
  m2 = m*m;
  alpha2 = alpha*alpha;

  u = l*(2*m2*(R2sq -R3sq - l2) + alpha*(R4sq - R3sq - m2))
        / (4*l2*m2 - alpha2);
  v = m* (2*l2*(R4sq - R3sq - m2) + alpha*(R2sq - R3sq - l2))
        / (4*l2*m2 - alpha2);

  u2 = u*u;
  v2 = v*v;

  d = (R3sq - u2 - v2 + 2*u*v*cose);
  if (fabs(d/(R3sq + u2 + v2 + 1)*(maxlength*maxlength/(maxR*maxR))) < EPS)
    d = 0.0;
  d = sqrt(d);

  sinsq = 1.0 - cose*cose;
  if (fabs(sinsq) < EPS) sinsq = 0.0;
  sine = sqrt(sinsq);
  tmp1 = d*d*cose;
  tmp2 = d*sine;
  tmp3 = sine*sine;

  if (fabs(d) < EPS) 
    omega = 0.0;   /* d is zero, so it doesn't matter */
  else
    omega = atan2( (tmp1 + (u+l)*(v + m)*tmp3),(tmp2*R1))
           - atan2( (tmp1 + (u + l)*v*tmp3),(tmp2*R2))
	   + atan2( (tmp1 + u*v*tmp3),(tmp2*R3))
	   - atan2( (tmp1 + u*(v + m)*tmp3),(tmp2*R4) );

  tmp4 =(  (u+l)*atanh( m/(R1 + R2)) 
	  +(v+m)*atanh( l/(R1 + R4))
	  -    u*atanh( m/(R3 + R4))
	  -    v*atanh( l/(R2 + R3))  );

  if ( fabs(sine) < 1e-150) tmp5 = 0.0;
    else tmp5 = omega*d/sine;

  M = MUOVER4PI*cose*(2*tmp4 - tmp5 );

  return M;
}




/* calculates mutual inductance of "filaments" with width and height */
/* as a combination of filament approximations                       */
/* some functions it uses are in dist_betw_fils.c */
double mutual(fil_j, fil_m)
FILAMENT *fil_j, *fil_m;
{
  double totalM;
  int i,j,ij;
  double aspect_j, aspect_m;
  static double cutoff = 3;
  int nfilsj, nfilsm;
  double dist, rj, rm, sum1, sum2;
  int parallel, whperp;
  int edge_par, num_dims;
  extern double **Gweight;    /* gaussian quad weights. */
                              /*  Filled in dist_betw_fils.c*/
  double widj[3], heightj[3];
  double dims[NUM_DIMS];
  Table **lastptr;
  int dim_count, test_self;

  dist = dist_betw_fils(fil_j, fil_m, &parallel); 
  rj = MAX(fil_j->width, fil_j->height)/2.0; 
  rm = MAX(fil_m->width, fil_m->height)/2.0;  

  if (MAX(rj,rm)*100 < dist) { 
    /* fils are far apart */ 
    num_mutualfil++;

    totalM = mutualfil(fil_j, fil_m);

    if (!finite(totalM))
      print_infinity_warning(fil_j, fil_m);

    return totalM;
  }
  else {    
    if (parallel == 1) {
      get_wid(fil_j, widj);
      get_height(fil_j, widj, heightj);
    }
    else
      if ( fabs(vdotp(fil_j->lenvect,fil_m->lenvect))
	    /(fil_j->length*fil_m->length) < EPS  )  {
        num_perp++;
	/* fils are perpendicular */
	return 0.0;
      }

    edge_par = parallel == 1
                 && edges_parallel(fil_j,fil_m,widj,&whperp);
    if (edge_par)
      if (lookup(fil_j, fil_m, whperp, widj, heightj, 
		 &totalM, dims, &lastptr, &dim_count,&num_dims)==1)
	{
	  num_found++;
	  return totalM;
	}
    if (edge_par && 2*MAX(rj,rm)*10 > dist){
      /* fils are close enough to use exact integrals  (6/94) */
        
        // This is actually a self term
        if (nearzero(dist)) {
            totalM = selfterm(fil_j);
            return totalM;
        }
        
        num_exact_mutual++;
        totalM = parallel_fils(fil_j, fil_m, whperp, widj, heightj, dist);
      put_in_table(fil_j, fil_m, whperp, totalM, dims, dim_count, lastptr,
		   num_dims);
      /* printf("ex: %14.8lg \n", totalM);*/

      if (!finite(totalM))
	print_infinity_warning(fil_j, fil_m);
      
      return totalM;
    }
    else {
      /* do 5 filament approximation to the filament */
      num_fourfil++;
      totalM = fourfil(fil_j, fil_m);
      
      if (edge_par) put_in_table(fil_j, fil_m, whperp, totalM,
				 dims, dim_count, lastptr, num_dims);
      /* printf("4fil: %14.8lg\n ", totalM); */
      
      if (!finite(totalM))
	print_infinity_warning(fil_j, fil_m);
      
      return totalM;
    }
  }
  
}



#define XX0 0
#define YY0 1
#define ZZ0 2
#define XX1 3
#define YY1 4
#define ZZ1 5
#define DEBUG FALSE

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{   
	double *O, *L, *W, *H;
    double *Oin, *Lin, *Win, *Hin;
    double *out, *stats, totalM;
    double P[6], Q[6];
    int m, n, ind, ii, dispflag;
    FILAMENT fil1, fil2;
    FILAMENT *fil;
    SEGMENT seg1, seg2;

    if (nlhs>2)
    mexErrMsgTxt("Wrong number of output parameters, usage:  M = mutual(O, L, W, H)");
    if (nrhs!=4)
    mexErrMsgTxt("Wrong number of input parameters, usage:  M = mutual(O, L, W, H)");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
    mexErrMsgTxt("mutual: Input arguments must be double.");
 

    // Prepare input
    Oin = mxGetPr(prhs[0]); // O - origin
    Lin = mxGetPr(prhs[1]); // L - length
    Win = mxGetPr(prhs[2]); // W - width
    Hin = mxGetPr(prhs[3]); // H - height
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    // Error check
    if(!(6 == mxGetM(prhs[0]) &&
    		 6 == mxGetM(prhs[1]) &&
    		 6 == mxGetM(prhs[2]) &&
    		 6 == mxGetM(prhs[3])))
    		 mexErrMsgTxt("mutual: Must have six rows for each of the four inputs");
    if(n != mxGetN(prhs[1]) ||
       n != mxGetN(prhs[2]) ||
       n != mxGetN(prhs[3]))
    		mexErrMsgTxt("mutual: Must have same number of columns");   		   
    
    // Output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    
    // Display output?
    if (nlhs == 2) {
        plhs[1] = mxCreateDoubleMatrix(6, 1, mxREAL);
        stats = mxGetPr(plhs[1]);
    } else
        dispflag == TRUE;
            
    
    // Reset counters
    num_exact_mutual=0;
    num_fourfil=0;
    num_mutualfil=0;
    num_found=0;
    num_perp=0;
    num_self=0;
    
    //if (table_init == 0) {
        init_table();
        table_init = 1;
    //}
    
    for (ind=0;ind<n;ind++) {
          // Increment pointers
          O = Oin + 6*ind; L = Lin + 6*ind; W = Win + 6*ind; H = Hin + 6*ind;
          
          // Convert to the FILAMENT format        
          // Do fil1 first.
          // length, width, height
          fil1.length = sqrt(SQUARE(L[XX0])+SQUARE(L[YY0])+SQUARE(L[ZZ0]));
          fil1.width = sqrt(SQUARE(W[XX0])+SQUARE(W[YY0])+SQUARE(W[ZZ0]));
          fil1.height = sqrt(SQUARE(H[XX0])+SQUARE(H[YY0])+SQUARE(H[ZZ0]));
          
          // If this is a self term, we can skip the rest of the stuff early.
          

          // length, width, height
          fil2.length = sqrt(SQUARE(L[XX1])+SQUARE(L[YY1])+SQUARE(L[ZZ1]));
          fil2.width = sqrt(SQUARE(W[XX1])+SQUARE(W[YY1])+SQUARE(W[ZZ1]));
          fil2.height = sqrt(SQUARE(H[XX1])+SQUARE(H[YY1])+SQUARE(H[ZZ1]));

          // Make end points. 
          // Note that FastHenry defines the origin at the center of the 
          // square, but the OLWH notation defines the origin at the corner.
          // We define P to be the FastHenry origin and Q to be the end point.
          for (ii = 0; ii < 6; ii++) {
              P[ii] = O[ii] + W[ii]/2;  
              Q[ii] = O[ii] + L[ii] + W[ii]/2;
          }

          fil1.x[0] = P[XX0]; fil1.y[0] = P[YY0]; fil1.z[0] = P[ZZ0];
          fil1.x[1] = Q[XX0]; fil1.y[1] = Q[YY0]; fil1.z[1] = Q[ZZ0];
          fil2.x[0] = P[XX1]; fil2.y[0] = P[YY1]; fil2.z[0] = P[ZZ1];
          fil2.x[1] = Q[XX1]; fil2.y[1] = Q[YY1]; fil2.z[1] = Q[ZZ1];

          // lenvect
          fil1.lenvect[XX] = L[XX0]; fil2.lenvect[XX] = L[XX1];
          fil1.lenvect[YY] = L[YY0]; fil2.lenvect[YY] = L[YY1];
          fil1.lenvect[ZZ] = L[ZZ0]; fil2.lenvect[ZZ] = L[ZZ1];

          // Make segments
          // SEGMENT.widthdir is a unit vector in the direction of the width
          seg1.widthdir = (double *) mxMalloc(3* sizeof(double));
          seg1.widthdir[XX] = W[XX0] / fil1.width;
          seg1.widthdir[YY] = W[YY0] / fil1.width;
          seg1.widthdir[ZZ] = W[ZZ0] / fil1.width;
          seg1.length = fil1.length;

          seg2.widthdir = (double *) mxMalloc(3* sizeof(double));
          seg2.widthdir[XX] = W[XX1] / fil2.width;
          seg2.widthdir[YY] = W[YY1] / fil2.width;
          seg2.widthdir[ZZ] = W[ZZ1] / fil2.width;
          seg2.length = fil2.length;

          // Point the segments
          fil1.segm = &seg1;
          fil2.segm = &seg2;
          
          // Output
          out[ind] = mutual(&fil1,&fil2);          
          
          // Debug code if needed
          if (DEBUG) {
              fil = &fil1;
              mexPrintf("  fil 1: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
              fil->x[0],fil->y[0],fil->z[0],fil->x[1],fil->y[1],fil->z[1],
              fil->width, fil->height);
              mexPrintf("  fil 1: widthdir:%lg %lg %lg\n", 
              fil->segm->widthdir[0],fil->segm->widthdir[1],fil->segm->widthdir[2]);
              fil = &fil2;
              mexPrintf("  fil 2: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
              fil->x[0],fil->y[0],fil->z[0],fil->x[1],fil->y[1],fil->z[1],
              fil->width, fil->height);
              mexPrintf("  fil 2: widthdir:%lg %lg %lg\n", 
              fil->segm->widthdir[0],fil->segm->widthdir[1],fil->segm->widthdir[2]);
          }
    }

    mxFree(seg1.widthdir);
    mxFree(seg2.widthdir);
    //destroy_table();
    // More debug tools
    if (dispflag){
        mexPrintf("Calls to exact_mutual: %15d\n",num_exact_mutual);
        mexPrintf("         fourfils:     %15d\n",num_fourfil);
        mexPrintf("         mutualfil:    %15d\n",num_mutualfil);
        mexPrintf("Number found in table: %15d\n",num_found);
        mexPrintf("Number perpendicular:  %15d\n",num_perp);
        mexPrintf("Number self terms:     %15d\n",num_self);
        mexPrintf("\n");
    } else {
        stats[0] = num_exact_mutual;
        stats[1] = num_fourfil;
        stats[2] = num_mutualfil;
        stats[3] = num_found;
        stats[4] = num_perp;
        stats[5] = num_self;
    }
    return;
}



