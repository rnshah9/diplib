/*
 * (c)2016-2017, Cris Luengo.
 * Based on original DIPlib code: (c)1995-2014, Delft University of Technology.
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
#include "diplib/framework.h"
#include "diplib/overload.h"
#include "diplib/saturated_arithmetic.h"

namespace dip {

//
void Add(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt
) {
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_ALL( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return dip::saturated_add( *its[ 0 ], *its[ 1 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt, *scanLineFilter ));
}

//
void Subtract(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt
) {
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_ALL( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return dip::saturated_sub( *its[ 0 ], *its[ 1 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt, *scanLineFilter ));
}

//
namespace {

template< typename TPI >
class MultiplyLineFilter : public Framework::ScanLineFilter {
   public:
      virtual dip::uint GetNumberOfOperations( dip::uint, dip::uint, dip::uint ) override {
         return nRows_ * nColumns_ * nInner_;
      }
      virtual void Filter( Framework::ScanLineFilterParameters const& params ) override {
         // This function is only called for two non-scalar images.
         DIP_ASSERT( params.inBuffer.size() == 2 );
         DIP_ASSERT( params.outBuffer.size() == 1 );
         TPI const* lhs = static_cast< TPI const* >( params.inBuffer[ 0 ].buffer );
         TPI const* rhs = static_cast< TPI const* >( params.inBuffer[ 1 ].buffer );
         TPI* out = static_cast< TPI* >( params.outBuffer[ 0 ].buffer );
         dip::sint const lhsStride = params.inBuffer[ 0 ].stride;
         dip::sint const rhsStride = params.inBuffer[ 1 ].stride;
         dip::sint const outStride = params.outBuffer[ 0 ].stride;
         dip::sint const lhsTensorStride = params.inBuffer[ 0 ].tensorStride;
         dip::sint const rhsTensorStride = params.inBuffer[ 1 ].tensorStride;
         dip::sint const outTensorStride = params.outBuffer[ 0 ].tensorStride;
         DIP_ASSERT( params.inBuffer[ 0 ].tensorLength == nRows_ * nInner_ );
         DIP_ASSERT( params.inBuffer[ 1 ].tensorLength == nInner_ * nColumns_ );
         DIP_ASSERT( params.outBuffer[ 0 ].tensorLength == nRows_ * nColumns_ );
         dip::uint const bufferLength = params.bufferLength;
         for( dip::uint ii = 0; ii < bufferLength; ++ii ) {
            TPI const* rhsT = rhs;
            TPI* outT = out;
            for( dip::uint col = 0; col < nColumns_; ++col ) {
               TPI const* lhsT = lhs;
               for( dip::uint row = 0; row < nRows_; ++row ) {
                  TPI const* lhsTT = lhsT;
                  TPI const* rhsTT = rhsT;
                  FlexType< TPI > v = 0;
                  for( dip::uint jj = 0; jj < nInner_; ++jj ) {
                     v += static_cast< FlexType< TPI >>( *lhsTT ) * static_cast< FlexType< TPI >>( *rhsTT );
                     lhsTT += static_cast< dip::sint >( nRows_ ) * lhsTensorStride;
                     rhsTT += rhsTensorStride;
                  }
                  *outT = clamp_cast< TPI >( v );
                  lhsT += lhsTensorStride;
                  outT += outTensorStride;
               }
               rhsT += static_cast< dip::sint >( nInner_ ) * rhsTensorStride;
            }
            lhs += lhsStride;
            rhs += rhsStride;
            out += outStride;
         }
      }
      MultiplyLineFilter( dip::uint nRows, dip::uint nColumns, dip::uint nInner )
            : nRows_( nRows ), nColumns_( nColumns ), nInner_( nInner ) {}
   private:
      dip::uint nRows_;     // == lhs.TensorRows
      dip::uint nColumns_;  // == rhs.TensorColumns
      dip::uint nInner_;    // == lhs.TensorColumns == rhs.TensorRows
};

template< typename TPI >
class MultiplySymmetricLineFilter : public Framework::ScanLineFilter {
   public:
      virtual dip::uint GetNumberOfOperations( dip::uint, dip::uint, dip::uint ) override {
         return nOuter_ * ( nOuter_ + 1 ) * nInner_ / 2;
      }
      virtual void Filter( Framework::ScanLineFilterParameters const& params ) override {
         // This function is only called for one non-scalar image.
         DIP_ASSERT( params.inBuffer.size() == 1 ); // RHS matrix, meaning the inner dimension is the columns
         DIP_ASSERT( params.outBuffer.size() == 1 );
         TPI const* in = static_cast< TPI const* >( params.inBuffer[ 0 ].buffer );
         TPI* out = static_cast< TPI* >( params.outBuffer[ 0 ].buffer );
         dip::sint const inStride = params.inBuffer[ 0 ].stride;
         dip::sint const outStride = params.outBuffer[ 0 ].stride;
         dip::sint const inTensorStride = params.inBuffer[ 0 ].tensorStride;
         dip::sint const outTensorStride = params.outBuffer[ 0 ].tensorStride;
         DIP_ASSERT( params.inBuffer[ 0 ].tensorLength == nOuter_ * nInner_ );
         DIP_ASSERT( params.outBuffer[ 0 ].tensorLength == ( nOuter_ * ( nOuter_ + 1 )) / 2 );
         dip::uint const bufferLength = params.bufferLength;
         for( dip::uint ii = 0; ii < bufferLength; ++ii ) {
            TPI const* inT = in;
            TPI* outT = out;
            // Compute diagonal elements first
            for( dip::uint col = 0; col < nOuter_; ++col ) {
               FlexType< TPI > v = 0;
               for( dip::uint jj = 0; jj < nInner_; ++jj ) {
                  v += static_cast< FlexType< TPI >>( *inT ) * static_cast< FlexType< TPI >>( *inT );
                  inT += inTensorStride;
               }
               *outT = clamp_cast< TPI >( v );
               outT += outTensorStride;
            }
            // Elements above diagonal are stored column-wise
            dip::sint colSkip = static_cast< dip::sint >( nInner_ ) * inTensorStride;
            TPI const* rhsT = in + colSkip;
            for( dip::uint col = 1; col < nOuter_; ++col ) { // Elements above diagonal stored column-wise
               TPI const* lhsT = in;
               for( dip::uint row = 0; row < col; ++row ) {
                  TPI const* lhsTT = lhsT;
                  TPI const* rhsTT = rhsT;
                  FlexType< TPI > v = 0;
                  for( dip::uint jj = 0; jj < nInner_; ++jj ) {
                     v += static_cast< FlexType< TPI >>( *lhsTT ) * static_cast< FlexType< TPI >>( *rhsTT );
                     lhsTT += inTensorStride;
                     rhsTT += inTensorStride;
                  }
                  *outT = clamp_cast< TPI >( v );
                  lhsT += colSkip;
                  outT += outTensorStride;
               }
               rhsT += colSkip;
            }
            in += inStride;
            out += outStride;
         }
      }
      MultiplySymmetricLineFilter( dip::uint nOuter, dip::uint nInner )
            : nOuter_( nOuter ), nInner_( nInner ) {}
   private:
      dip::uint nOuter_;    // == lhs.TensorRows == rhs.TensorColumns
      dip::uint nInner_;    // == lhs.TensorColumns == rhs.TensorRows
};

template< typename TPI >
class MultiplyDiagonalLineFilter : public Framework::ScanLineFilter {
   public:
      virtual dip::uint GetNumberOfOperations( dip::uint, dip::uint, dip::uint ) override {
         return lhsTensor_.Elements();
      }
      virtual void Filter( Framework::ScanLineFilterParameters const& params ) override {
         DIP_ASSERT( params.inBuffer.size() == 2 );
         DIP_ASSERT( params.outBuffer.size() == 1 );
         TPI const* lhs = static_cast< TPI const* >( params.inBuffer[ 0 ].buffer );
         TPI const* rhs = static_cast< TPI const* >( params.inBuffer[ 1 ].buffer );
         TPI* out = static_cast< TPI* >( params.outBuffer[ 0 ].buffer );
         dip::sint const lhsStride = params.inBuffer[ 0 ].stride;
         dip::sint const rhsStride = params.inBuffer[ 1 ].stride;
         dip::sint const outStride = params.outBuffer[ 0 ].stride;
         dip::sint const lhsTensorStride = params.inBuffer[ 0 ].tensorStride;
         dip::sint const rhsTensorStride = params.inBuffer[ 1 ].tensorStride;
         dip::sint const outTensorStride = params.outBuffer[ 0 ].tensorStride;
         dip::uint const N = lhsTensor_.Columns();
         DIP_ASSERT( params.inBuffer[ 0 ].tensorLength == lhsTensor_.Elements() );
         DIP_ASSERT( params.inBuffer[ 1 ].tensorLength == N );
         DIP_ASSERT( params.outBuffer[ 0 ].tensorLength == lhsTensor_.Elements() );
         dip::uint const bufferLength = params.bufferLength;
         if( lhsTensor_.IsSymmetric() ) {
            // The symmetric case: symm * diag
            for( dip::uint ii = 0; ii < bufferLength; ++ii ) {
               TPI const* lhsT = lhs;
               TPI const* rhsT = rhs;
               TPI* outT = out;
               for( dip::uint m = 0; m < N; ++m ) {
                  *outT = static_cast< TPI >( *lhsT * *rhsT ); // integer promotion causes compiler warnings
                  lhsT += lhsTensorStride;
                  rhsT += rhsTensorStride;
                  outT += outTensorStride;
               }
               rhsT = rhs;
               for( dip::uint n = 1; n < N; ++n ) {
                  for( dip::uint m = 0; m < n; ++m ) {
                     *outT = static_cast< TPI >( *lhsT * *rhsT ); // integer promotion causes compiler warnings
                     lhsT += lhsTensorStride;
                     outT += outTensorStride;
                  }
                  rhsT += rhsTensorStride;
               }
               lhs += lhsStride;
               rhs += rhsStride;
               out += outStride;
            }
         } else {
            // The full case: full * diag (full can be column-major or row-major)
            dip::uint const M = lhsTensor_.Rows();
            dip::sint lhsRowStride = lhsTensorStride;
            dip::sint lhsColStride = lhsTensorStride * static_cast< dip::sint >( M );
            if( !lhsTensor_.HasNormalOrder() ) {
               lhsRowStride = lhsTensorStride * static_cast< dip::sint >( N );
               lhsColStride = lhsTensorStride;
            }
            dip::sint outRowStride = outTensorStride;
            dip::sint outColStride = outTensorStride * static_cast< dip::sint >( M );
            if( transposeOutput_ ) {
               outRowStride = outTensorStride * static_cast< dip::sint >( N );
               outColStride = outTensorStride;
            }
            for( dip::uint ii = 0; ii < bufferLength; ++ii ) {
               TPI const* lhsC = lhs;
               TPI const* rhsT = rhs;
               TPI* outC = out;
               for( dip::uint n = 0; n < N; ++n ) {
                  TPI const* lhsR = lhsC;
                  TPI* outR = outC;
                  for( dip::uint m = 0; m < M; ++m ) {
                     *outR = static_cast< TPI >( *lhsR * *rhsT ); // integer promotion causes compiler warnings
                     lhsR += lhsRowStride;
                     outR += outRowStride;
                  }
                  lhsC += lhsColStride;
                  rhsT += rhsTensorStride;
                  outC += outColStride;
               }
               lhs += lhsStride;
               rhs += rhsStride;
               out += outStride;
            }
         }
      }
      MultiplyDiagonalLineFilter( Tensor lhsTensor, bool transposeOutput )
            : lhsTensor_( lhsTensor ), transposeOutput_( transposeOutput ) {}
   private:
      Tensor lhsTensor_; // lhs is either symmetric or full; rhs is diagonal with same sizes
      bool transposeOutput_;
};

} // namespace

void Multiply(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   if( lhs.IsScalar() || rhs.IsScalar() ) {
      DIP_STACK_TRACE_THIS( MultiplySampleWise( lhs, rhs, out, dt_out ));
      return;
   }
   if( lhs.TensorColumns() != rhs.TensorRows() ) {
      DIP_THROW( "Inner tensor dimensions must match in multiplication" );
   }
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   Tensor lhsTensorTransposed = lhs.Tensor();
   lhsTensorTransposed.Transpose();
   if(( lhsTensorTransposed == rhs.Tensor() ) && lhs.IsIdenticalView( rhs )) {
      // a' * a  or  a * a' : produces a symmetric matrix
      dip::uint nOuter = lhs.TensorRows();
      dip::uint nInner = lhs.TensorColumns();
      Tensor outTensor = Tensor( Tensor::Shape::SYMMETRIC_MATRIX, nOuter, nOuter );
      std::unique_ptr< Framework::ScanLineFilter > scanLineFilter;
      DIP_OVL_NEW_FLEXBIN( scanLineFilter, MultiplySymmetricLineFilter, ( nOuter, nInner ), dt );
      ImageRefArray outar{ out };
      DIP_STACK_TRACE_THIS( Framework::Scan( { rhs }, outar, { dt }, { dt }, { dt_out }, { outTensor.Elements() }, *scanLineFilter,
                                             Framework::ScanOption::ExpandTensorInBuffer + Framework::ScanOption::NotInPlace ));
      out.ReshapeTensor( outTensor );
   } else {
      bool oneIsDiagonal = lhs.Tensor().IsDiagonal() || rhs.Tensor().IsDiagonal();
      bool oneIsTriangular = lhs.Tensor().IsTriangular() || rhs.Tensor().IsTriangular();
      if( oneIsDiagonal && !oneIsTriangular ) {
         Image lhs_copy = lhs;
         Image rhs_copy = rhs;
         if( lhs_copy.IsVector() ) {
            rhs_copy.ReshapeTensor( lhs_copy.Tensor() ); // convert the other into a same shape vector
         } else if( rhs_copy.IsVector() ) {
            lhs_copy.ReshapeTensor( rhs_copy.Tensor() ); // convert the other into a same shape vector
         }
         if(( lhs_copy.Tensor().IsDiagonal() && rhs_copy.Tensor().IsDiagonal() ) ||
            ( lhs_copy.IsVector() && rhs_copy.IsVector() )) {
            // Here we deal with: diag*diag, vector*diag, diag*vector
            MultiplySampleWise( lhs_copy, rhs_copy, out, dt );
            if( rhs.IsVector() ) {
               out.ReshapeTensor( rhs.Tensor() );
            } else {
               out.ReshapeTensor( lhs.Tensor() );
            }
         } else {
            // Here we deal with:
            //  - full*diag or symm*diag: multiply each lhs column with corresponding diag element
            //  - diag*full or diag*symm: multiply each rhs row with corresponding diag element
            bool transposeOutput = false;
            Tensor outTensor = lhs_copy.Tensor();
            if( lhs_copy.Tensor().IsDiagonal() ) {
               // compute the reverse operation, then transpose the result
               lhs_copy.swap( rhs_copy );
               outTensor = lhs_copy.Tensor();
               lhs_copy.Transpose();
               transposeOutput = true;
            }
            if( outTensor.TensorShape() == Tensor::Shape::ROW_MAJOR_MATRIX ) {
               outTensor.ChangeShape( outTensor.Rows() ); // Force column-major matrix
            }
            DIP_ASSERT( lhs_copy.TensorColumns() == rhs_copy.TensorElements() );
            // We've transformed the problem to one of the two cases: full*diag or symm*diag
            std::unique_ptr< Framework::ScanLineFilter > scanLineFilter;
            DIP_OVL_NEW_FLEXBIN( scanLineFilter, MultiplyDiagonalLineFilter, ( lhs_copy.Tensor(), transposeOutput ), dt );
            ImageRefArray outar{ out };
            DIP_STACK_TRACE_THIS( Framework::Scan( { lhs_copy, rhs_copy }, outar, { dt, dt }, { dt }, { dt_out }, { outTensor.Elements() }, *scanLineFilter,
                                                   Framework::ScanOption::NotInPlace ));
            out.ReshapeTensor( outTensor );
         }
      } else {
         // General case: tri*diag, diag*tri, or anything not involving a diagonal matrix
         Tensor outTensor = Tensor( lhs.TensorRows(), rhs.TensorColumns() );
         std::unique_ptr< Framework::ScanLineFilter > scanLineFilter;
         DIP_OVL_NEW_FLEXBIN( scanLineFilter, MultiplyLineFilter, ( lhs.TensorRows(), rhs.TensorColumns(), lhs.TensorColumns() ), dt );
         ImageRefArray outar{ out };
         DIP_STACK_TRACE_THIS( Framework::Scan( { lhs, rhs }, outar, { dt, dt }, { dt }, { dt_out }, { outTensor.Elements() }, *scanLineFilter,
                                                Framework::ScanOption::ExpandTensorInBuffer + Framework::ScanOption::NotInPlace ));
         out.ReshapeTensor( outTensor );
      }
   }
}

void MultiplySampleWise(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_ALL( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return dip::saturated_mul( *its[ 0 ], *its[ 1 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt_out, *scanLineFilter ));
}

void MultiplyConjugate(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   if( rhs.DataType().IsComplex() && dt.IsComplex() ) {
      std::unique_ptr< Framework::ScanLineFilter > scanLineFilter;
      DIP_OVL_CALL_ASSIGN_COMPLEX( scanLineFilter, Framework::NewDyadicScanLineFilter, (
            []( auto its ) { return dip::saturated_mul( *its[ 0 ], std::conj( *its[ 1 ] )); }, 4
      ), dt );
      DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt_out, *scanLineFilter ));
   } else {
      DIP_STACK_TRACE_THIS( MultiplySampleWise( lhs, rhs, out, dt_out ));
   }
}

//
void Divide(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_FLEXBIN( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return dip::saturated_div( *its[ 0 ], *its[ 1 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt_out, *scanLineFilter ));
}

//
void SafeDivide(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   if( dt.IsBinary() ) {
      Divide( lhs, rhs, out, dt_out );
      return;
   }
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_FLEX( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return dip::saturated_safediv( *its[ 0 ], *its[ 1 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt_out, *scanLineFilter ));
}

//
void Modulo(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt
) {
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   if( dt.IsFloat() ) {
      DIP_OVL_CALL_ASSIGN_FLOAT( scanLineFilter, Framework::NewDyadicScanLineFilter, (
            []( auto its ) { return std::fmod( *its[ 0 ], *its[ 1 ] ); }
      ), dt );
   } else {
      DIP_OVL_CALL_ASSIGN_INTEGER( scanLineFilter, Framework::NewDyadicScanLineFilter, (
            []( auto its ) { return static_cast< decltype( *its[ 0 ] ) >( *its[ 0 ] % *its[ 1 ] ); }
      ), dt );
   }
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt, *scanLineFilter ));
}

//
void Power(
      Image const& lhs,
      Image const& rhs,
      Image& out,
      DataType dt_out
) {
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DataType dt = DataType::SuggestArithmetic( lhs.DataType(), rhs.DataType() );
   if( dt.IsBinary() ) {
      dt = DataType::SuggestFlex( dt_out );
   }
   DIP_OVL_CALL_ASSIGN_FLEX( scanLineFilter, Framework::NewDyadicScanLineFilter, (
         []( auto its ) { return std::pow( *its[ 0 ], *its[ 1 ] ); }, 20 // Rough guess at the cost
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanDyadic( lhs, rhs, out, dt, dt, dt_out, *scanLineFilter ));
}

//
template< typename TPI >
class InvertLineFilter : public Framework::ScanLineFilter {
   public:
      virtual void Filter( Framework::ScanLineFilterParameters const& params ) override {
         TPI const* in = static_cast< TPI const* >( params.inBuffer[ 0 ].buffer );
         TPI* out = static_cast< TPI* >( params.outBuffer[ 0 ].buffer );
         dip::sint const inStride = params.inBuffer[ 0 ].stride;
         dip::sint const outStride = params.outBuffer[ 0 ].stride;
         dip::uint const bufferLength = params.bufferLength;
         for( dip::uint ii = 0; ii < bufferLength; ++ii ) {
            // Tensor dimension is 1 because we request `ScanOption::TensorAsSpatialDim`
            *out = saturated_inv( *in );
            in += inStride;
            out += outStride;
         }
      }
};

void Invert(
      Image const& in,
      Image& out
) {
   DataType dt = in.DataType();
   std::unique_ptr< Framework::ScanLineFilter >scanLineFilter;
   DIP_OVL_CALL_ASSIGN_ALL( scanLineFilter, Framework::NewMonadicScanLineFilter, (
         []( auto its ) { return saturated_inv( *its[ 0 ] ); }
   ), dt );
   DIP_STACK_TRACE_THIS( Framework::ScanMonadic( in, out, dt, dt, 1, *scanLineFilter, Framework::ScanOption::TensorAsSpatialDim ));
}


} // namespace dip

#ifdef DIP_CONFIG_ENABLE_DOCTEST
#include "doctest.h"
#include "diplib/math.h"
#include "diplib/testing.h"

DOCTEST_TEST_CASE("[DIPlib] testing the matrix multiplication operation") {
   // Case 1: general case:
   dip::Image matrix2x3( { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
   matrix2x3.ReshapeTensor( 2, 3 );
   dip::Image tmp;
   DOCTEST_CHECK_THROWS( tmp = matrix2x3 * matrix2x3 );
   dip::Image matrix3x2( { 0.5, 1.0, 2.0, 5.0, 10.0, 20.0 } );
   matrix3x2.ReshapeTensor( 3, 2 );
   dip::Image out = matrix2x3 * matrix3x2;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 2, 2 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 13.5, 17.0, 135.0, 170.0 } ));
   out = dip::Transpose( matrix3x2 ) * dip::Transpose( matrix2x3 );
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 2, 2 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 13.5, 135.0, 17.0, 170.0 } ));
   out = matrix3x2 * matrix2x3;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 3 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 10.5, 21.0, 42.0, 21.5, 43.0, 86.0, 32.5, 65.0, 130.0 } ));
   dip::Image vector3( { 1.0, 10.0, 100.0 } );
   out = matrix2x3 * vector3;
   DOCTEST_CHECK( out.TensorElements() == 2 );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_VECTOR );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 531.0, 642.0 } ));
   out = dip::Transpose( vector3 ) * matrix3x2;
   DOCTEST_CHECK( out.TensorElements() == 2 );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::ROW_VECTOR );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 210.5, 2105.0 } ));
   dip::Image otherVector3( { 1.0, 2.0, 3.0 } );
   otherVector3.Transpose();
   out = otherVector3 * vector3;
   DOCTEST_CHECK( out.TensorElements() == 1 );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 321.0 } ));
   out = vector3 * otherVector3;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 3 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0 } ));
   dip::Image symmetric( { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
   symmetric.ReshapeTensor( dip::Tensor( dip::Tensor::Shape::SYMMETRIC_MATRIX, 3, 3 ));
   out = symmetric * matrix3x2;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 2 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 14.5, 16.0, 14.5, 145.0, 160.0, 145.0 } ));

   // Case 2: diagonal matrix case:
   dip::Image diag( { 1.0, 10.0, 100.0 } );
   diag.ReshapeTensorAsDiagonal();
   out = matrix2x3 * diag;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 2, 3 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 1.0, 2.0, 30.0, 40.0, 500.0, 600.0 } ));
   out = diag * dip::Transpose( matrix2x3 );
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 2 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 1.0, 30.0, 500.0, 2.0, 40.0, 600.0 } ));
   out = diag * matrix3x2;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 2 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_MAJOR_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 0.5, 10.0, 200.0, 5.0, 100.0, 2000.0 } ));
   dip::Image diag2( { 1.0, 2.0, 3.0 } );
   diag2.ReshapeTensorAsDiagonal();
   out = diag * diag2;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3, 3 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::DIAGONAL_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 1.0, 20.0, 300.0 } ));
   out = diag2 * vector3;
   DOCTEST_CHECK( out.TensorSizes() == dip::UnsignedArray{ 3 } );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::COL_VECTOR );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 1.0, 20.0, 300.0 } ));

   // Case 3: creating symmetric matrices when multiplying by transposed self
   out = matrix2x3 * Transpose( matrix2x3 );
   DOCTEST_CHECK( out.TensorElements() == 3 );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::SYMMETRIC_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 35.0, 56.0, 44.0 } ));
   out = Transpose( matrix2x3 ) * matrix2x3;
   DOCTEST_CHECK( out.TensorElements() == 6 );
   DOCTEST_CHECK( out.TensorShape() == dip::Tensor::Shape::SYMMETRIC_MATRIX );
   DOCTEST_CHECK( out.At( 0 ) == dip::Image::Pixel( { 5.0, 25.0, 61.0, 11.0, 17.0, 39.0 } ));
}

DOCTEST_TEST_CASE("[DIPlib] testing the operator overloading") {
   // We want all these things to compile:
   dip::Image lhs( { 1.0, 2.0, 3.0 } );
   dip::Image rhs( { 1.0, 10.0, 100.0 } );
   dip::Image out;
   Add( lhs, rhs, out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs + rhs, dip::Option::CompareImagesMode::FULL ));
   Add( lhs[ 0 ], rhs, out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs[ 0 ] + rhs, dip::Option::CompareImagesMode::FULL ));
   Add( lhs, rhs[ 0 ], out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs + rhs[ 0 ], dip::Option::CompareImagesMode::FULL ));
   Add( lhs[ 0 ], rhs[ 0 ], out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs[ 0 ] + rhs[ 0 ], dip::Option::CompareImagesMode::FULL ));
   Add( lhs, 1, out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs + 1, dip::Option::CompareImagesMode::FULL ));
   Add( 1, rhs, out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, 1 + rhs, dip::Option::CompareImagesMode::FULL ));
   Add( lhs[ 0 ], 1, out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, lhs[ 0 ] + 1, dip::Option::CompareImagesMode::FULL ));
   Add( 1, rhs[ 0 ], out );
   DOCTEST_CHECK( dip::testing::CompareImages( out, 1 + rhs[ 0 ], dip::Option::CompareImagesMode::FULL ));
}

DOCTEST_TEST_CASE("[DIPlib] inplace arithmetic produces correct values") {
   dip::Image floatim( { 1.0 } );
   dip::Image intim = dip::Convert( floatim, dip::DT_SINT32 );
   DOCTEST_CHECK( dip::testing::CompareImages( intim * 2, floatim * 2 ));
   DOCTEST_CHECK( dip::testing::CompareImages( intim * 0.5, floatim * 0.5 ));
   floatim *= 2;
   intim *= 2;
   DOCTEST_CHECK( dip::testing::CompareImages( intim, floatim ));
   floatim *= 0.5;
   intim *= 0.5; // There's a bug where this value is converted to int before multiplying
   DOCTEST_CHECK( dip::testing::CompareImages( intim, floatim ));
}

DOCTEST_TEST_CASE( "[DIPlib] testing inplace arithmetic operators on dip::Image::View" ) {
   dip::Image im( { 5 }, 3 );
   im.Fill( 0.7 );
   im[ 1 ].Fill( 0.2 );
   dip::Image res = im.Copy();

   // Regular view
   dip::Image tmp = res[ 0 ];
   dip::Add( tmp, res[ 1 ], tmp, tmp.DataType() );
   im[ 0 ] += im[ 1 ]; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));

   dip::MultiplySampleWise( tmp, res[ 2 ], tmp, tmp.DataType() );
   im[ 0 ] *= im[ 2 ]; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));

   // Irregular view
   res.At( 4 ) += 1;
   im.At( dip::CoordinateArray{{ 4 }} ) += 1; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));

   res.At( 3 ) *= 3;
   res.At( 4 ) *= 3;
   im.At( dip::CoordinateArray{{ 3 }, { 4 }} ) *= 3; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));

   // Idem with bitwise operators
   im = im > 1;
   res = res > 1;
   tmp = res[ 0 ];
   dip::And( tmp, res[ 1 ], tmp );
   im[ 0 ] &= im[ 1 ]; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));

   res.At( 4 ) &= 0;
   im.At( dip::CoordinateArray{{ 4 }} ) &= 0; // operation we're testing
   DOCTEST_CHECK( dip::testing::CompareImages( im, res ));
}

DOCTEST_TEST_CASE("[DIPlib] testing the scan framework with singleton-expanded dimensions") {

   // Generate test image with a singleton-expanded dimension
   dip::Image img{ dip::UnsignedArray{ 256, 1 }, 1, dip::DT_SFLOAT };
   img.Fill( 0 );
   img.ExpandSingletonDimension( 1, 256 );
   DOCTEST_CHECK( img.Sizes() == dip::UnsignedArray{ 256, 256 } );

   // Arithmetic operation with a scalar -- both dimensions will get singleton expanded
   dip::Image out = img + 1;
   DOCTEST_CHECK( out.Sizes() == dip::UnsignedArray{ 256, 256 } );
   DOCTEST_CHECK( out.IsSingletonExpanded() );
   DOCTEST_CHECK( out.Stride( 0 ) == 1 );
   DOCTEST_CHECK( out.Stride( 1 ) == 0 );
   DOCTEST_CHECK( out.At( 0, 0 ) == 1 ); // ensure we're producing the right values
   DOCTEST_CHECK( out.At( 50, 200 ) == 1 ); // ensure we're producing the right values
}

#endif // DIP_CONFIG_ENABLE_DOCTEST
