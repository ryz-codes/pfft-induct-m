#include "mex.h"
#include "lru.h"

class key_box
{
    public:
     double a;
     double b;

     key_box(double a = 0, double b = 0) : a(a), b(b) {}

     bool operator<( const key_box & n ) const {
         //mexPrintf("Compared %g, %g with %g, %g\n",this->a,this->b,n.a,n.b);
       return (this->a < n.a) || (this->b < n.b);   // for example
     }
};

double dbl_fun (double in1, double in2) {
    return in1*2 + in2*2;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
  const mxArray *prhs[]) {
    double *in1, *in2, *out, *my_in1, *my_in2, *my_out;
    int t, t2, ii;
    LRUCache<key_box,double> mycache(10);
    key_box mykey;
    
    if (nrhs != 2)
        mexErrMsgTxt("Must have two inputs\n");
    
    
    // Set up input
    in1 = mxGetPr(prhs[0]);
    t = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    
    in2 = mxGetPr(prhs[1]);
    t2 = mxGetM(prhs[1]) * mxGetN(prhs[1]);
    
    if (t != t2)
        mexErrMsgTxt("numel of both inputs must be equal!");
    
    // Setup output
    plhs[0] = mxCreateDoubleMatrix(t, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    
    // iterate through the first column
    for (ii=0; ii<t; ii++) {
        my_in1 = in1+ii;
        my_in2 = in2+ii;
        my_out = out+ii;
        
        mykey = key_box(*my_in1,*my_in2);
        *my_out = mycache.get(mykey);
        
        // Handle not finding the key-value pair.
        if (*my_out == -1) {
            *my_out = dbl_fun(*my_in1, *my_in2);
            mycache.put(mykey,*my_out);
            mexPrintf("Put %g, %g\n",*my_in1, *my_in2);
        }
    }
}