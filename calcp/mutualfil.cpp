#include <math.h>
#include "mex.h"
#include "mutualfil.h"
//#define dotp(x1,y1,z1,x2,y2,z2) x1*x2+y1*y2+z1*z2 //fast dot product macro

/* these are missing in some math.h files */
#define atanh(x) 0.5*log((1+x)/(1-x))
#define asinh(x) log(x+sqrt(x*x+1)); 

/* Turn these commonly used funcitons into macros */
#define magdiff2(x1,y1,z1,x2,y2,z2,node1,node2) ( SQUARE(x1[node1] - x2[node2]) \
	       +SQUARE(y1[node1] - y2[node2]) \
	       +SQUARE(z1[node1] - z2[node2]) \
	       )
#define dotprod(x1,y1,z1,x2,y2,z2)  ((x1[1] - x1[0])*(x2[1] - x2[0]) \
	 + (y1[1] - y1[0])*(y2[1] - y2[0]) \
	 + (z1[1] - z1[0])*(z2[1] - z2[0]) )
     

int realcos_error = 0;

/* calculates the mutual inductance between two filaments */
/* from Grover, Chapter 7 */
double mutualfil_joel(
        double *x1,
        double *y1,
        double *z1,
        double *x2,
        double *y2,
        double *z2)
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

// Four input form M = mutualfil(from1,to1,from2,to2)
#define IN_from1     prhs[0]
#define IN_to1       prhs[1]
#define IN_from2     prhs[2]
#define IN_to2       prhs[3]

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
    FILAMENT fil1[1], fil2[1];
    double *from1_in, *to1_in, *from2_in, *to2_in, x1[2], y1[2], z1[2], x2[2], y2[2], z2[2];
    double *from1, *to1, *from2, *to2, *out;
    int ind, m, n;
    
    // Parse inputs
    if (nrhs == 4) {
        from1_in = mxGetPr(IN_from1);
        to1_in = mxGetPr(IN_to1);
        from2_in = mxGetPr(IN_from2);
        to2_in = mxGetPr(IN_to2);
    } else {
        mexErrMsgTxt("mutualfil: Incorrect use of function.\n Inputs are M = mutualfil(from1,to1,from2,to2)");
    }
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    // Error check
    if(!(3 == mxGetM(prhs[0]) &&
    		 3 == mxGetM(prhs[1]) &&
    		 3 == mxGetM(prhs[2]) &&
    		 3 == mxGetM(prhs[3])))
    		 mexErrMsgTxt("mutualfil: Must have three rows for each of the four inputs");
    if(n != mxGetN(prhs[1]) ||
       n != mxGetN(prhs[2]) ||
       n != mxGetN(prhs[3]))
    		mexErrMsgTxt("mutualfil: Must have same number of columns");   		   
    
    // Output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    
    for (ind=0;ind<n;ind++) {
    
        from1 = from1_in + 3*ind;
        to1 = to1_in + 3*ind;
        from2 = from2_in + 3*ind;
        to2 = to2_in + 3*ind;
        
        // Set up the c structs
        x1[0] = from1[XX];
        y1[0] = from1[YY];
        z1[0] = from1[ZZ];
        x1[1] = to1[XX];
        y1[1] = to1[YY];
        z1[1] = to1[ZZ];

        x2[0] = from2[XX];
        y2[0] = from2[YY];
        z2[0] = from2[ZZ];
        x2[1] = to2[XX];
        y2[1] = to2[YY];
        z2[1] = to2[ZZ];

        // Prepare for output.
        out[ind] = mutualfil_joel(x1,y1,z1,x2,y2,z2);
    
    }
    return;
}
    
/* MATLAB CODE
USEDBLQUAD=0; % debug with dblquad
u1 = fil1.u;
u2 = fil2.u;
L1 = fil1.lenvec;
L2 = fil2.lenvec;

% Obtain the x directional unit vector in fil1's frame of reference
ux = [u1(1:2) 0];
ux = ux./norm(ux);

% z directional pointer is just z
uz = [0 0 1];

% Rise angle calculations
cos_t1 = dot(u1,ux); % yes
sin_t1 = dot(u1,uz); % yes
cos_t2 = dot(u2,ux); % yes
sin_t2 = dot(u2,uz); % yes

% define distance vectors
rx1 = @(lu) fil1.from(1) + lu*L1(1); % yes
ry1 = @(lu) fil1.from(2) + lu*L1(2); % yes
rz1 = @(lu) fil1.from(3) + lu*L1(3); % yes
rx2 = @(lu) fil2.from(1) + lu*L2(1); % yes
ry2 = @(lu) fil2.from(2) + lu*L2(2); % yes
rz2 = @(lu) fil2.from(3) + lu*L2(3); % yes

rval = @(lu1,lu2) ((rx1(lu1)-rx2(lu2)).^2+(ry1(lu1)-ry2(lu2)).^2).^0.5; % yes
% Dyadic, anisotropic Green's function. Integrate Gxx and Gzz
% separately. Assumes that all cross terms are zero.
if abs(cos_t1)>hObj.EPS && abs(cos_t2)>hObj.EPS  
    if USEDBLQUAD==1
        Zx_int = @(l1,l2) hObj.g.get_Gr(rval(l1,l2),rz1(l1),rz2(l2),f*ones(size(l1)));
        Zx = dblquad(Zx_int,0,1,0,1,1e-12)...
            .*1j.*2*pi*f*cos_t1*cos_t2*fil1.len*fil2.len;
    else
        Zx_int = @(l) hObj.g.get_Gr(rval(l(:,1),l(:,2)), ...
                       rz1(l(:,1)),rz2(l(:,2)), ...
                       f*ones(length(l),1));
        Zx = hObj.sp2.spgr(Zx_int);
        Zx = Zx.*1j.*2*pi*f*cos_t1*cos_t2*fil1.len*fil2.len;
    end
else % one of the filaments is perp to the xy plane
    Zx=0;
end

if abs(sin_t1)>hObj.EPS && abs(sin_t2)>hObj.EPS
    if USEDBLQUAD==1
        Zz_int = @(l1,l2) hObj.g.get_Gz(rval(l1,l2),rz1(l1),rz2(l2),f*ones(size(l1)));
        Zz = dblquad(Zz_int,0,1,0,1,1e-12)...
           .*1j.*2*pi*f*sin_t1*sin_t2*fil1.len*fil2.len;
    else
        Zz_int = @(l) hObj.g.get_Gz(rval(l(:,1),l(:,2)), ... Radial values
                        rz1(l(:,1)),rz2(l(:,2)), ... z and d
                        f*ones(length(l),1)); % frequency
        Zz = hObj.sp2.spgr(Zz_int);
        Zz = Zz.*1j.*2*pi*f*sin_t1*sin_t2*fil1.len*fil2.len;
    end
else % one of the filaments is parallel to the xy plane
    Zz = 0;
end

% Combine
Zm_mf = Zx + Zz;

*/

/* Alternate method to populate the output, one element at a time */
/* for (j = 0; j < 3; j++) */
/* { */
/*   output[j] = data[j]; */
/* } */

