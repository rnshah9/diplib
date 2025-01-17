/*
 * (c)2017, Cris Luengo.
 *
 * Encapsulates code Taken from OpenCV 3.1: (c)2000, Intel Corporation.
 * (see src/transform/opencv_dxt.cpp)
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

#ifndef DIP_DFT_H
#define DIP_DFT_H

#include <vector>
#include <complex>
#include <limits>

#include "diplib/library/export.h"

/// \file
/// \brief Declares an interface to a DFT function.
/// See \ref transform.


namespace dip {


/// \addtogroup transform


/// \brief An object that encapsulates the Discrete Fourier Transform (DFT).
///
/// Usage:
///
/// ```cpp
/// DFT dft( size, inverse );               // creates the object with all the data ready to start running DFTs.
/// std::vector< std::complex< T >> buf( opts.BufferSize() ); // creates a buffer
/// dft.Apply( in, out, buf.data() );                         // computes a DFT, repeat as necessary
/// dft.Initialize( size2, inverse );                         // changes the options for the new size / direction
/// buf.resize( opts.BufferSize() );                          // resizes the buffer
/// dft.Apply( in, out, buf.data() );                         // computes a different DFT, repeat as necessary
/// ```
///
/// Note that this code uses `int` for sizes, rather than \ref dip::uint. \ref maximumDFTSize is the largest length
/// of the transform.
///
/// The template can be instantiated for `T = float` or `T = double`. Linker errors will result for other types.
///
/// The DFT is computed using an FFT algorithm that is optimized for lengths that are a multiple of 2, 3 and 5.
/// The larger the factors above 5, the longer the algorithm will take.
template< typename T >
class DFT {
   public:

      /// \brief A default-initialized `DFT` object is useless. Call `Initialize` to make it useful.
      DFT() = default;

      /// \brief Construct a `DFT` object by specifying the size and direction of the transform.
      /// Note that this is not a trivial operation.
      DFT( std::size_t size, bool inverse ) {
         Initialize( size, inverse );
      }

      /// \brief Re-configure a `DFT` object to the given transform size and direction.
      /// Note that this is not a trivial operation.
      DIP_EXPORT void Initialize( std::size_t size, bool inverse );

      /// \brief Apply the transform that the `DFT` object is configured for.
      ///
      /// `source` and `destination` are pointers to contiguous buffers with \ref TransformSize elements.
      /// This is the value of the `size` parameter of the constructor or \ref Initialize. `buffer` is a pointer
      /// to a contiguous buffer used for intermediate data. It should have \ref BufferSize elements.
      ///
      /// `scale` is a real scalar that the output values are multiplied by. It is typically set to `1/size` for
      /// the inverse transform, and 1 for the forward transform.
      ///
      /// `source` and `destination` can only point to the same buffer if all factors of \ref TransformSize are
      /// the same. One should avoid this in general situations.
      DIP_EXPORT void Apply(
            const std::complex< T >* source,
            std::complex< T >* destination,
            std::complex< T >* buffer,
            T scale
      ) const;

      /// \brief Returns true if this represents an inverse transform, false for a forward transform.
      bool IsInverse() const { return inverse_; }

      /// \brief Returns the size that the transform is configured for.
      std::size_t TransformSize() const { return static_cast< std::size_t >( nfft_ ); }

      /// \brief Returns the size of the buffer expected by `Apply`.
      std::size_t BufferSize() const { return static_cast< std::size_t >( sz_ ); }

   private:
      int nfft_ = 0;
      bool inverse_ = false;
      std::vector< int > factors_;
      std::vector< int > itab_;
      std::vector< std::complex< T >> wave_;
      int sz_ = 0; // Size of the buffer to be passed to DFT.
};

/// \brief Returns a size equal or larger to `size0` that is efficient for our DFT implementation.
///
/// Set `larger` to false to return a size equal or smaller instead.
///
/// Returns 0 if `size0` is too large for our DFT implementation.
///
/// Prefer to use \ref dip::OptimalFourierTransformSize in your applications, it will throw an error if
/// the transform size is too large.
DIP_EXPORT std::size_t GetOptimalDFTSize( std::size_t size0, bool larger = true );

/// \brief The largest size supported by \ref DFT and \ref FourierTransform, equal to 2^31^-1.
/// Both the built-in DFT and the FFTW library use `int` for array sizes.
constexpr std::size_t maximumDFTSize = static_cast< std::size_t >( std::numeric_limits< int >::max() );


/// \endgroup

} // namespace dip

#endif // DIP_DFT_H
