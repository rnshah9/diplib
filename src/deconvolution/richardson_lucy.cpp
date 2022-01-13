/*
 * DIPlib 3.0
 * This file contains functions for RL deconvolution
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

#include "get_otf.h"

namespace dip {

namespace {

bool ParseRichardsonLucyOptions( StringSet const& options ) {
   bool isOtf = false;
   for( auto const& opt: options ) {
      if( opt == "OTF" ) {
         isOtf = true;
      } else {
         DIP_THROW_INVALID_FLAG( opt );
      }
   }
   return isOtf;
}

} // namespace

void RichardsonLucy(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations,
      StringSet const& options
) {
   DIP_THROW_IF( !in.IsForged() || !psf.IsForged(), E::IMAGE_NOT_FORGED );
   DIP_THROW_IF( !in.IsScalar() || !psf.IsScalar(), E::IMAGE_NOT_SCALAR );
   DIP_THROW_IF( !in.DataType().IsReal(), E::DATA_TYPE_NOT_SUPPORTED );
   bool isOtf;
   DIP_STACK_TRACE_THIS( isOtf = ParseRichardsonLucyOptions( options ));
   // Fourier transforms etc.
   Image H;
   DIP_STACK_TRACE_THIS( H = GetOTF( psf, in.Sizes(), isOtf ));
   Image G = FourierTransform( in );

   DIP_THROW( E::NOT_IMPLEMENTED );
   ( void )maxIterations;

   // Inverse Fourier transform
   DIP_STACK_TRACE_THIS( FourierTransform( G, out, { S::INVERSE, S::REAL } ));
}

void RichardsonLucyTotalVariation(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations,
      dfloat regularization,
      StringSet const& options
) {
   DIP_THROW_IF( !in.IsForged() || !psf.IsForged(), E::IMAGE_NOT_FORGED );
   DIP_THROW_IF( !in.IsScalar() || !psf.IsScalar(), E::IMAGE_NOT_SCALAR );
   DIP_THROW_IF( !in.DataType().IsReal(), E::DATA_TYPE_NOT_SUPPORTED );
   bool isOtf;
   DIP_STACK_TRACE_THIS( isOtf = ParseRichardsonLucyOptions( options ));
   // Fourier transforms etc.
   Image H;
   DIP_STACK_TRACE_THIS( H = GetOTF( psf, in.Sizes(), isOtf ));
   Image G = FourierTransform( in );

   DIP_THROW( E::NOT_IMPLEMENTED );
   ( void )maxIterations;
   ( void )regularization;

   // Inverse Fourier transform
   DIP_STACK_TRACE_THIS( FourierTransform( G, out, { S::INVERSE, S::REAL } ));
}

} // namespace dip
