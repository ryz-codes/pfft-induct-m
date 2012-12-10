#include <math.h>
#include "mex.h"

#define PI 3.14159265358979323846
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13

// my own OLWH format
#define XX0 0
#define YY0 1
#define ZZ0 2
#define XX1 3
#define YY1 4
#define ZZ1 5

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

/* these are missing in some math.h files */
#define atanh(x) 0.5*log((1+x)/(1-x))
#define asinh(x) log(x+sqrt(x*x+1)) 

/* Turn these commonly used funcitons into macros */
#define magdiff2(x1,y1,z1,x2,y2,z2,node1,node2) ( SQUARE(x1[node1] - x2[node2]) \
	       +SQUARE(y1[node1] - y2[node2]) \
	       +SQUARE(z1[node1] - z2[node2]) \
	       )
#define dotprod(x1,y1,z1,x2,y2,z2)  ((x1[1] - x1[0])*(x2[1] - x2[0]) \
	 + (y1[1] - y1[0])*(y2[1] - y2[0]) \
	 + (z1[1] - z1[0])*(z2[1] - z2[0]) )
	 

     
// Function prototypes
double mutual(double* O, double* L, double* W, double* H);
double self(double W, double L, double T);
double mutualfil(double *O1, double *L1, double *O2, double *L2);
double mut_rect(double len, double d);
double getDist(double* O, double* L, double* W, double* H);
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]);

int realcos_error = 0;


/* calculates mutual inductance of "filaments" with width and height */
/* as a combination of filament approximations                       */
// Missing the analytical solutions, the single filament approximation
double mutual(
        double *O,
        double *L,
        double *W,
        double *H
        )
{
  	double *x1,*x2,*y1,*y2,*z1,*z2;
  	double len[2], wid[2], hei[2];
  	double P0[6], P1[6], P2[6], P3[6], P4[6];
    int ii;
  	double totalM = 0;

  	// This part is dedicated to parsing the OLWH input
	len[0] = sqrt(SQUARE(L[XX0])+SQUARE(L[YY0])+SQUARE(L[ZZ0]));
	len[1] = sqrt(SQUARE(L[XX1])+SQUARE(L[YY1])+SQUARE(L[ZZ1]));
	
	wid[0] = sqrt(SQUARE(W[XX0])+SQUARE(W[YY0])+SQUARE(W[ZZ0]));
	hei[0] = sqrt(SQUARE(H[XX0])+SQUARE(H[YY0])+SQUARE(H[ZZ0]));

  	// Perpendicular
	if ( fabs(L[XX0]*L[XX1]+L[YY0]*L[YY1]+L[ZZ0]*L[ZZ1]) // dot product of the L vector
	    /(len[0] * len[1]) < EPS  )  {
		return 0.0;   
	    mexPrintf("Perp");
   	}

	//"if they are close enough, then just calculate self term"   
   if (checkSelf(O,L,W,H)) {
	    mexPrintf("Self");
   		totalM = self(wid[0],len[0],hei[0]);
   		return totalM;
   	}

    /* do 5 filament approximation to the filament */
    // Generate the ten filaments
    for (ii = 0; ii < 6; ii++) {
    	// central filament
	    P0[ii] = O[ii] + W[ii]/2 + H[ii]/2; 
	    
	    // side filament
	    P1[ii] = O[ii] + W[ii]/2;
	    P2[ii] = O[ii] + W[ii]/2 + H[ii];
	    P3[ii] = O[ii] + H[ii]/2;
	    P4[ii] = O[ii] + H[ii]/2 + W[ii];
	}	

    totalM += mutualfil(P0+3,L+3,P1,L);
    totalM += mutualfil(P0+3,L+3,P2,L);
    totalM += mutualfil(P0+3,L+3,P3,L);
    totalM += mutualfil(P0+3,L+3,P4,L);
            
    totalM += mutualfil(P0,L,P1+3,L+3);
    totalM += mutualfil(P0,L,P2+3,L+3);
    totalM += mutualfil(P0,L,P3+3,L+3);
    totalM += mutualfil(P0,L,P4+3,L+3);
    
    totalM -= 2*mutualfil(P0,L,P0+3,L+3);
      
    return totalM/6.0;
    
}

int checkSelf(double* O, double* L, double* W, double* H) 
{
	int ii;
	double dsum;
	for (ii=0;ii<2;ii++) {
		dsum += fabs(O[ii] - O[ii+3]);
		dsum += fabs(L[ii] - L[ii+3]);
		dsum += fabs(W[ii] - W[ii+3]);
		dsum += fabs(H[ii] - H[ii+3]);
	}
	return dsum < (EPS*12);
}

/* self inductance */
double self(double W, double L, double T)
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
    z *= MU0*L;  /* this is inductance */
    
    return z; 

}

/* calculates the mutual inductance between two filaments */
/* from Grover, Chapter 7 */
double mutualfil(double *O1, double *L1, double *O2, double *L2)
{
  double x1[2],x2[2],y1[2],y2[2],z1[2],z2[2];
  double P1[3], P2[3];
  double R, R1, R2, R3, R4, l, m;
  double cose, sine, u,v,d;
  double alpha, tmp1, tmp2, tmp3, tmp4, tmp5;
  double maxR, maxlength, sinsq, blah, minR;
  double scaleEPS, realcos;
  double vtemp;

  double omega,M;
  int signofM, ii;

  double junk;

  /* Load OLWH format */
  // Starting points
  x1[0] = O1[XX]; y1[0] = O1[YY]; z1[0] = O1[ZZ];
  x2[0] = O2[XX]; y2[0] = O2[YY]; z2[0] = O2[ZZ];
  
  // Ending points
  for (ii = 0; ii < 3; ii++) {
  		P1[ii] = O1[ii] + L1[ii];
  		P2[ii] = O2[ii] + L2[ii];
  }
  x1[1] = P1[XX]; y1[1] = P1[YY]; z1[1] = P1[ZZ];
  x2[1] = P2[XX]; y2[1] = P2[YY]; z2[1] = P2[ZZ];  
  
  /* for parallel filaments */
  double ux, uy, uz, Rx, Ry, Rz, x1_0, x1_1, x2_0, x2_1, magu;
  double dx, dy, dz, dotp, vx, vy, vz;
  double R1sq, R2sq, R3sq, R4sq, m2, l2, u2, v2, alpha2;

  R1sq = magdiff2(x1,y1,z1,x2,y2,z2,1,1);
  R2sq = magdiff2(x1,y1,z1,x2,y2,z2,1,0);
  R3sq = magdiff2(x1,y1,z1,x2,y2,z2,0,0);
  R4sq = magdiff2(x1,y1,z1,x2,y2,z2,0,1);
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

  l = sqrt(dotprod(x1,y1,z1,x1,y1,z1));
  m = sqrt(dotprod(x2,y2,z2,x2,y2,z2));
 
  if (l<EPS || m<EPS)
  	mexPrintf("Warning: zero-length filament detected\n"); 
 
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

      M = MUOVER4PI*2*(dotprod(x1,y1,z1,x2,y2,z2)/(l*m))
	   *(l*atanh(m/(l+R)) + m*atanh(l/(m+R)));
      /* note: dotprod should take care of signofM */

      return M;
    }

  cose = alpha/(2*l*m);
  if (fabs(cose) > 1) cose = (cose < 0 ? -1.0 : 1.0);
  blah = 1.0 - fabs(cose);

  /* let's use the real cosine */
  realcos = dotprod(x1,y1,z1,x2,y2,z2)/(l*m);
  
  /* Segments are perpendicular! */
  if (fabs(realcos) < EPS)
    return 0.0;

  if (fabs((realcos - cose)/cose) > 0.1) 
    if (realcos_error == 0) {
      fprintf(stderr, "Internal Warning: realcos = %lg,  cose = %lg\n",realcos, cose);
      fprintf(stderr,"  This may be due to two filaments that are separated \n\
by a distance 1e10 times their length\n");
      realcos_error = 1;
    }

  cose = realcos;

  /* filaments parallel */
  tmp1 = fabs( fabs(cose) - 1);
  if ( fabs( fabs(cose) - 1) < EPS) {
    /* determine a vector in the direction of d with length d */
    /* (d is the distance between the lines made by the filament */
    Rx = x2[0] - x1[0];  /* vector from fil1 to fil2 */
    Ry = y2[0] - y1[0];
    Rz = z2[0] - z1[0];
    ux = x1[1] - x1[0];
    uy = y1[1] - y1[0];
    uz = z1[1] - z1[0];
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
    vx =  (x2[0] - ( x1[0] + dx));
    vy =  (y2[0] - ( y1[0] + dy));
    vz =  (z2[0] - ( z1[0] + dz));
    x2_0 = vx*ux + vy*uy + vz*uz;
    vtemp = sqrt(vx*vx + vy*vy + vz*vz);

    /* same thing for x2_1 */
    vx =  (x2[1] - ( x1[0] + dx));
    vy =  (y2[1] - ( y1[0] + dy));
    vz =  (z2[1] - ( z1[0] + dz));
    x2_1 = vx*ux + vy*uy + vz*uz;
    
    if ( fabs( (sqrt(vx*vx + vy*vy + vz*vz) - fabs(x2_1))
	      /(MAX(fabs(x2_0)+d,fabs(x2_1)+d))) 
	> EPS) {
      printf("uh oh, segs don't seem parallel %lg\n",(sqrt(vx*vx + vy*vy * vz*vz) - fabs(x2_1)));
    }

    if ( fabs( (vtemp - fabs(x2_0))/(MAX(fabs(x2_0)+d,fabs(x2_1)+d))) 
	> EPS) {
      printf("uh oh, segs don't seem parallel\n");
    }


    if ( fabs(d) < EPS )  { /* collinear! */
      M = MUOVER4PI*(fabs(x2_1 - x1_0)*log(fabs(x2_1 - x1_0))
	             - fabs(x2_1 - x1_1)*log(fabs(x2_1 - x1_1))
		     - fabs(x2_0 - x1_0)*log(fabs(x2_0 - x1_0))
		     + fabs(x2_0 - x1_1)*log(fabs(x2_0 - x1_1)) );
      return M;
    }  /* end collinear */

    M = MUOVER4PI*(mut_rect((x2_1 - x1_1),d) - mut_rect((x2_1 - x1_0),d) 
		 - mut_rect((x2_0 - x1_1),d) + mut_rect((x2_0 - x1_0),d) );
    
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

double mut_rect(double len, double d)
{
  double temp,temp1;

  temp = sqrt(len*len + d*d);
  temp1 = len*asinh(len/d) ;
  return temp - temp1;
}


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{   
	double *O, *L, *W, *H;
	double *this_O, *this_L, *this_W, *this_H;
    double *out;
    int m, n, ind;

    if (nlhs>1)
    mexErrMsgTxt("Wrong number of output parameters, usage:  M = mutual(O, L, W, H)");
    if (nrhs!=4)
    mexErrMsgTxt("Wrong number of input parameters, usage:  M = mutual(O, L, W, H)");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
    mexErrMsgTxt("mutual: Input arguments must be double.");
 

    // Prepare input
    O = mxGetPr(prhs[0]); // O - origin
    L = mxGetPr(prhs[1]); // L - length
    W = mxGetPr(prhs[2]); // W - width
    H = mxGetPr(prhs[3]); // H - height
    
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
    
    for (ind=0;ind<n;ind++) {
        out[ind] = mutual(O+ind*6, L+ind*6, W+ind*6, H+ind*6);
    }
    return;
}