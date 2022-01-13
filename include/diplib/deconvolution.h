/*
 * DIPlib 3.0
 * This file contains declarations for deconvolution functions
 *
 * (c)2017-2022, Cris Luengo.
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

#ifndef DIP_DECONVOLUTION_H
#define DIP_DECONVOLUTION_H

#include "diplib.h"


/// \file
/// \brief Microscopy-related functionality.
/// See \ref microscopy.


namespace dip {


/// \group deconvolution Deconvolution
/// \brief Deconvolution algorithms (inverse filtering).
/// \addtogroup


/// \brief Wiener Deconvolution using estimates of signal and noise power
///
/// If $G$ is the Fourier transform of `in`, $H$ is the Fourier transform of `psf`,
/// and $F$ is the Fourier transform of `out`, then this function estimates the $F$ that optimally
/// (in the least squares sense) satisfies $G = FH$ (that is, `in` is the result of the convolution of
/// `out` with `psf`).
///
/// Finding `out` requires knowledge of the power spectrum of the signal and the noise. The Wiener deconvolution
/// filter is defined in the frequency domain as
///
/// $$ H_\text{inv} = \frac{H^* S}{ H^* H S + N } \; , $$
///
/// where $S$ is `signalPower`, and $N$ is `noisePower`. These functions are typically not known, but:
///
/// - `signalPower` can be estimated as the Fourier transform of the autocorrelation of `in`. If a raw image
///   is passed for this argument (`dip::Image{}`), then it will be computed as such.
///
/// - `noisePower` can be estimated as a flat function. A 0D image can be given here, it will be expanded to
///   the size of the other images. `noisePower` should not be zero anywhere, as that might lead to division
///   by zero and consequently meaningless results.
///
/// The other syntax for \ref dip::WienerDeconvolution takes an estimate of the noise-to-signal
/// ratio instead of the signal and noise power spectra. Note that $H_\text{inv}$ can be rewritten as
///
/// $$ H_\text{inv} = \frac{H^*}{ H^* H  + \frac{N}{S} } = \frac{H^*}{ H^* H + K } \; , $$
///
/// where $K$ is the noise-to-signal ratio. If $K$ is a constant, then the Wiener deconvolution filter is
/// equivalent to the Tikhonov regularized inverse filter.
///
/// `psf` is given in the spatial domain, and will be zero-padded to the size of `in` and Fourier transformed.
/// The PSF (point spread function) should sum to one in order to preserve the mean image intensity.
/// If the OTF (optical transfer function, the Fourier transform of the PSF) is known, it is possible to pass
/// that as `psf`; add the string `"OTF"` to `options`.
///
/// All input images must be real-valued and scalar, except if the OFT is given instead of the PSF, in which
/// case `psf` could be complex-valued.
DIP_EXPORT void WienerDeconvolution(
      Image const& in,
      Image const& psf,
      Image const& signalPower,
      Image const& noisePower,
      Image& out,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image WienerDeconvolution(
      Image const& in,
      Image const& psf,
      Image const& signalPower,
      Image const& noisePower,
      StringSet const& options = {}
) {
   Image out;
   WienerDeconvolution( in, psf, signalPower, noisePower, out, options );
   return out;
}

/// \brief Wiener Deconvolution using an estimate of noise-to-signal ratio
///
/// See the description of the function with the same name above. The difference here is that a single number,
/// `regularization`, is given instead of the signal and noise power spectra. We then set $K$ (the
/// noise-to-signal ratio) to `regularization * dip::Maximum(P)`, with `P` equal to $H^* H$.
///
/// This formulation of the Wiener deconvolution filter is equivalent to the Tikhonov regularized inverse filter.
DIP_EXPORT void WienerDeconvolution(
      Image const& in,
      Image const& psf,
      Image& out,
      dfloat regularization = 1e-4,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image WienerDeconvolution(
      Image const& in,
      Image const& psf,
      dfloat regularization = 1e-4,
      StringSet const& options = {}
) {
   Image out;
   WienerDeconvolution( in, psf, out, regularization, options );
   return out;
}

//   G.M.P. van Kempen, Image Restoration in Fluorescence Microscopy,
//         PhD Thesis, Delft University of Technology, Delft, The Netherlands, 1998.
//   P.J. Verveer and T.M. Jovin, Acceleration of the ICTM image restoration algorithm,
//         Journal of Microscopy 188(3):191-195, 1997.
DIP_EXPORT void IterativeConstrainedTikhonovMiller(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations = 10,
      dfloat regularization = 0.1,
      dfloat stepSize = 1.0,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image IterativeConstrainedTikhonovMiller(
      Image const& in,
      Image const& psf,
      dip::uint maxIterations = 10,
      dfloat regularization = 0.1,
      dfloat stepSize = 1.0,
      StringSet const& options = {}
) {
   Image out;
   IterativeConstrainedTikhonovMiller( in, psf, out, maxIterations, regularization, stepSize, options );
   return out;
}

//   G.M.P. van Kempen, Image Restoration in Fluorescence Microscopy,
//         PhD Thesis, Delft University of Technology, Delft, The Netherlands, 1998.
//   W.H. Richardson, Bayesian-based iterative method of image restoration,
//        Journal of the Optical Society of America 62(1):55–59, 1972.
//   L.B. Lucy, An iterative technique for the rectification of observed distributions,
//        Astronomical Journal 79(6):745–754, 1974.
//
// also sometimes called the expectation maximization (EM) method.
DIP_EXPORT void RichardsonLucy(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations = 10,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image RichardsonLucy(
      Image const& in,
      Image const& psf,
      dip::uint maxIterations = 10,
      StringSet const& options = {}
) {
   Image out;
   RichardsonLucy( in, psf, out, maxIterations, options );
   return out;
}

//   N. Dey, L. Blanc-Féraud, C. Zimmer, P. Roux, Z. Kam, J. Olivo-Marin, J. Zerubia,
//        Richardson–Lucy algorithm with total variation regularization for 3D confocal microscope deconvolution,
//        Microscopy Research & Technique 69(4):260–266, 2006.
DIP_EXPORT void RichardsonLucyTotalVariation(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations = 10,
      dfloat regularization = 0.1,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image RichardsonLucyTotalVariation(
      Image const& in,
      Image const& psf,
      dip::uint maxIterations = 10,
      dfloat regularization = 0.1,
      StringSet const& options = {}
) {
   Image out;
   RichardsonLucyTotalVariation( in, psf, out, maxIterations, regularization, options );
   return out;
}

//   A. Beck, M. Teboulle, A fast iterative shrinkage-thresholding algorithm for linear inverse problems,
//        SIAM Journal on Imaging Sciences 2(1):183–202, 2009.
DIP_EXPORT void FastIterativeSoftThresholding(
      Image const& in,
      Image const& psf,
      Image& out,
      dip::uint maxIterations = 10,
      dfloat stepSize = 1.0,
      dfloat regularization = 0.1,
      dip::uint scale = 3,
      StringSet const& options = {}
);
DIP_NODISCARD inline Image FastIterativeSoftThresholding(
      Image const& in,
      Image const& psf,
      dip::uint maxIterations = 10,
      dfloat stepSize = 1.0,
      dfloat regularization = 0.1,
      dip::uint scale = 3,
      StringSet const& options = {}
) {
   Image out;
   FastIterativeSoftThresholding( in, psf, out, maxIterations, stepSize, regularization, scale, options );
   return out;
}

/// \endgroup

} // namespace dip

#endif // DIP_DECONVOLUTION_H
