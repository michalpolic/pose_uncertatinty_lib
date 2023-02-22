// Copyright Richard Hartley, 2010
static const char *copyright = "Copyright Richard Hartley, 2010";

//--------------------------------------------------------------------------
// LICENSE INFORMATION
//
// 1.  For academic/research users:
//
// This program is free for academic/research purpose:   you can redistribute
// it and/or modify  it under the terms of the GNU General Public License as 
// published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// Under this academic/research condition,  this program is distributed in 
// the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along 
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
// 2.  For commercial OEMs, ISVs and VARs:
// 
// For OEMs, ISVs, and VARs who distribute/modify/use this software 
// (binaries or source code) with their products, and do not license and 
// distribute their source code under the GPL, please contact NICTA 
// (www.nicta.com.au), and NICTA will provide a flexible OEM Commercial 
// License. 
//
//---------------------------------------------------------------------------

//
// Matlab usage:
//
//    E = solveE_nister (q1, q2)
//
// where
//    q1, q2 are matched points (each of dimension 3 x 6)
//    E is a list of essential matrices returned
//
// method used is the one described by Li and Hartley
//

#include "string.h"
#include "mex.h"
#include "hidden6.h"

typedef double Ematrix[3][3];

void compute_E_matrices_LAF(
	Observations q, Observations qp, Affinity a, Affinity ap, Ematrix Ematrices[10], int &nroots,
	bool optimized = true);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray *prhs[])
   {
   // Fast way to compute the Essential matrix from 6 points

   // -------------------------------------------------------------------------
   // Checking the arguments.

   // Check right number of arguments
   if (nrhs  != 4)
      {
      mexPrintf ("solveE: Wrong number of arguments.\n");
      mexPrintf ("Usage: [E,f] = solveE_nister_LAF(q1, q2)\n");
      return; 
      }

   // Check the input
   mxArray *q1 = prhs[0];
   mxArray *q2 = prhs[1];
   mxArray *a1 = prhs[2];
   mxArray *a2 = prhs[3];

   // Check the dimensions of the arguments
   size_t ndimensions1 = mxGetNumberOfDimensions (q1);
   const mwSize *q1dim = mxGetDimensions (q1);
   size_t ndimensions2 = mxGetNumberOfDimensions (q2);
   const mwSize *q2dim = mxGetDimensions (q2);
   size_t ndimensions3 = mxGetNumberOfDimensions(a1);
   const mwSize *q3dim = mxGetDimensions(a1);
   size_t ndimensions4 = mxGetNumberOfDimensions(a2);
   const mwSize *q4dim = mxGetDimensions(a2);

   // Now check them
   if (ndimensions1 != 2 || q1dim[0] != 2 || q1dim[1] != 2 ||
       ndimensions2 != 2 || q2dim[0] != 2 || q2dim[1] != 2 ||
	   ndimensions3 != 2 || q3dim[0] != 2 || q3dim[1] != 2 ||
	   ndimensions4 != 2 || q4dim[0] != 2 || q4dim[1] != 2 )
      {
      mexPrintf ("Bad input to mex function solveE\n");
      mexPrintf ("Inputs q1, q2, a1, a2 must have dimensions [2,2]\n");
      return; 
      }

   if (q1dim[1] != q2dim[1])
      {
      mexPrintf ("Bad input to mex function solveE\n");
      mexPrintf ( "Inputs q1 and q2 must have same dimensions.\n");
      return; 
      }

   // -------------------------------------------------------------------------
   // Read and reformat the input
   double q[2][2], qp[2][2], a[2][2], ap[2][2];

   double *p1 = (double *) mxGetData(q1);
   memcpy (&(q[0][0]),  p1, 4*sizeof(double));
   double *p2 = (double *) mxGetData(q2);
   memcpy (&(qp[0][0]),  p2, 4*sizeof(double));

   double *x1 = (double *)mxGetData(a1);
   memcpy(&(a[0][0]), x1, 4 * sizeof(double));
   double *x2 = (double *)mxGetData(a2);
   memcpy(&(ap[0][0]), x2, 4 * sizeof(double));

   /*mexPrintf("q[0][0] = %f\n", q[0][0]);
   mexPrintf("q[0][1] = %f\n", q[0][1]);
   mexPrintf("q[1][0] = %f\n", q[1][0]);
   mexPrintf("q[1][1] = %f\n", q[1][1]);

   mexPrintf("a[0][0] = %f\n", a[0][0]);
   mexPrintf("a[0][1] = %f\n", a[0][1]);
   mexPrintf("a[1][0] = %f\n", a[1][0]);
   mexPrintf("a[1][1] = %f\n", a[1][1]);*/

   // -------------------------------------------------------------------------
   // Do the computation
   Ematrix Ematrices[Maxdegree];
   double flengths[Maxdegree];
   int nroots;
   compute_E_matrices_LAF (q, qp, a, ap, Ematrices, nroots);


   // -------------------------------------------------------------------------
   // Return the results
   // Return the cameras MM
   if (nlhs > 0) 
      {
      // Create an array of the right size and fill with the two matrices
	  size_t dims[3];
      dims[0] = 3;
      dims[1] = 3; 
      dims[2] = nroots;
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[0]), Ematrices, nroots*sizeof(Ematrix));
      }

   // Return the cameras MM
   if (nlhs > 1) 
      {
      // Create an array of the right size and fill with the two matrices
	  size_t dims[2];
      dims[0] = nroots;
      dims[1] = 1; 
      plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[1]), flengths, nroots*sizeof(double));
      }
   }
