#include <math.h>

#define FILS_PER_MESH 2         /* filaments per mesh for breaking big loops.*/
#define PI 3.14159265358979323846
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define K 0.2235
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13
#define QUICK_EVAL 1          /* do one filament approximations if == 1 */
#define MAXsubfils 6          /* maximum points for Gaussian quadrature */
#define NUM_DIMS 10           /*Number of distinguishing parms. table lookup*/

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

typedef struct Filament {
  double x[2], y[2], z[2];  /* endpoints */
  double length, area, width, height;
  double lenvect[3];        /* vector along the length of filament */
  int filnumber;
  struct Segment *segm;
  struct charge *pchg;      /* 'charge' to send to multipole routines */
} FILAMENT;

typedef struct Segment {
   char *name;
   double *widthdir;   /*if width is not || to x-y plane and perpendicular to*/
                       /* the length, then this is 3 element vector in       */
                       /* in the direction of width*/
   int number;         /* an arbitrary number for the segment */
   int type;    /* CMS 8/21/92 -- type of structure the segment is in */
   double length;      
   double area;        /* area of cross section */
   double width, height;  /*width and height to cross section */
   int hinc, winc;             /* number of filament divisions in each dir */
//   NODES *node[2];                /* nodes at the ends */
   double sigma;              /* conductivity */
   double r_width, r_height;  /*ratio of adjacent fil widths(see assignFil())*/
   int num_fils;               /* hinc*winc */
   FILAMENT *filaments;        /* this segment's filaments */
/*   struct npath *conds;  */    /* linked list of conductors which this seg is in */
   struct pathlist *loops;   /* loops in which this segment is a member */
   int is_deleted;           /* has this segment been used already */

   /*struct g_nodes *gp_node[2];*//* nonuni_gp nodes that are really the ends*/

   /*struct _table *table;*/          /* lookup table for mutual terms */
   struct Segment *next;      /* next segment in list */
 } SEGMENT;
 
 
 //extern double magdiff2(FILAMENT *fil1, int node1, FILAMENT *fil2, int node2);
 //extern double dotprod(FILAMENT *fil1, FILAMENT *fil2);
 extern double mut_rect(double len, double d);