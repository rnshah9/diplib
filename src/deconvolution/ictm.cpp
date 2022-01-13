/*
 * DIPlib 3.0
 * This file contains functions for ICTM deconvolution
 *
 * (c)2022, Cris Luengo.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "diplib.h"
#include "diplib/deconvolution.h"
#include "diplib/transform.h"
#include "diplib/generation.h"

#include "get_otf.h"

namespace dip {

namespace {

dfloat MeanRelativeError( Image const& a, Image const& b ) {
   Image error = a - b;
   SafeDivide( error, a, error );
   MeanAbs( error, {}, error );
   std::cout << "MRE = " << error.As< dfloat >() << '\n';
   return error.As< dfloat >();
}

} // namespace

void IterativeConstrainedTikhonovMiller(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations,
      dfloat regularization,
      dfloat stepSize,
      StringSet const& options
) {
   DIP_THROW_IF( !in.IsForged() || !psf.IsForged(), E::IMAGE_NOT_FORGED );
   DIP_THROW_IF( !in.IsScalar() || !psf.IsScalar(), E::IMAGE_NOT_SCALAR );
   DIP_THROW_IF( !in.DataType().IsReal(), E::DATA_TYPE_NOT_SUPPORTED );
   bool isOtf = false;
   bool simple = false;
   for( auto const& opt: options ) {
      if( opt == "OTF" ) {
         isOtf = true;
      } else if( opt == "simple" ) {
         simple = true;
      } else {
         DIP_THROW_INVALID_FLAG( opt );
      }
   }
   // Fourier transform of inputs
   Image H;
   DIP_STACK_TRACE_THIS( H = GetOTF( psf, in.Sizes(), isOtf ));
   Image G = FourierTransform( in );
   // Regularization (an ideal Laplacian)
   dip::uint nD = G.Dimensionality();
   Image CtC( nD );
   for( dip::uint ii = 0; ii < nD; ++ii ) {
      Image ramp;
      CreateRamp( ramp, G.Sizes(), ii, { S::FREQUENCY } );
      ramp.UnexpandSingletonDimensions();
      Power( ramp, 2, ramp );
      if( ii == 0 ) {
         CtC = ramp;
      } else {
         CtC += ramp;
      }
   }
   CtC *= pi * pi * regularization;
   SquareModulus( CtC, CtC );
   // HtH
   Image A = SquareModulus( H );
   // A = HtH + regularization CtC
   A += CtC;
   CtC.Strip();
   // H^T g
   Image HtG = MultiplyConjugate( G, H );
   H.Strip();
   /*
   // Our first guess for the output is the unconstrained Tikhonov solution
   Image F = G * HtG;
   G.Strip();
   SafeDivide( F, A, F, F.DataType() );
   */
   // Our first guess for the output is the input
   Image F = G;
   G.Strip();
   // Initialize the remaining intermediate images used in the iterative process
   Image d{ 0.0f };
   Image r{ 1.0f };
   Image outPrev = in;
   Image Tf( true );
   Image rPrevSqMod( 1.0f );
   while( true ) {
      // r = A * f - H^T * g
      r = A * F;
      r -= HtG;
      // d = r + |r|^2 / |rPrev|^2 * dPrev
      if( simple ) {
         d = r;
      } else {
         // d = r + |r|^2 / |rPrev|^2 * dPrev
         /*
         Image rSqMod = SquareModulus( r );
         d *= rSqMod;
         SafeDivide( d, rPrevSqMod, d );
         rPrevSqMod = rSqMod;
         d += r;
          */
         d = r;
      }
      dfloat beta = -stepSize;
      if( !simple ) {
         // beta = - ( d^T T(f) r ) / ( d^T T(f) A T(f) d )
         Image dtTf = FourierTransform( d, { S::INVERSE, S::REAL } );
         dtTf *= Tf; // this is both "d^T T(f)" and "T(f) d" (with "d" in the spatial domain)
         Image ATfd = FourierTransform( dtTf );
         ATfd *= A;
         FourierTransform( ATfd, ATfd, { S::INVERSE, S::REAL } ); // ATfd is "A T(f) d"
         Image rSpatial = FourierTransform( r, { S::INVERSE, S::REAL } );
         beta = -InProduct( dtTf, rSpatial ) // dot product of "d^T T(f)" and "r" (in spatial domain)
               / InProduct( dtTf, ATfd );    // dot product of "d^T T(f)" and "A T(f) d"
      }
      // f = P( fPrev + beta * d )
      std::cout << "beta = " << beta << '\n';
      F += beta * d;
      DIP_STACK_TRACE_THIS( FourierTransform( F, out, { S::INVERSE, S::REAL } ));
      Tf = out >= 0; // save this for next iteration
      std::cout << "Number of negative values corrected: " << (Tf.NumberOfPixels() - Count( Tf )) << '\n';
      out *= Tf; // P(.) sets negative pixels to 0.
      // Do we stop iterating?
      if( --maxIterations == 0 ) {
         std::cout << "Terminating because maxIterations\n";
         return;
      }
      if( MeanRelativeError( out, outPrev ) < 1e-12 ) { // TODO: what is a good threshold here?
         std::cout << "Terminating because MRE < 1e-12\n";
         return;
      }
      if( !simple ) {
         outPrev = out.Copy(); // save for next iteration
      }
      DIP_STACK_TRACE_THIS( FourierTransform( out, F ));
   }
}

} // namespace dip
