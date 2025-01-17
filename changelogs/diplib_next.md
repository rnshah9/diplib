---
layout: post
title: "Changes DIPlib 3.x.x"
---

## Changes to *DIPlib*

### New functionality

- Added `dip::Kernel::IgnoreZeros()`, which causes subsequent calls to `dip::Kernel::PixelTable()` for that
  object to also exclude zero values in the gray-scale kernel image.

- Added `dip::ExtendImageToSize()`, which is like `dip::ExtendImage()` but allows instead specifying the
  size of the output image, like in `dip::Image::Pad()`.

- Added a static function `dip::Image::CropWindow()`, similar to the existing non-static one, which takes
  the image size as input.

- Added functions `dip::MeanRelativeError()` and `dip::MaximumRelativeError()`.

- Added `dip::Option::CompareImagesMode::APPROX_REL` as an option to `dip::testing::CompareImages()`, to
  compare pixel values with a relative tolerance.

- The function `dip::ConvolveFT()` has a new parameter `boundaryCondition`, if used the input image will be
  extended before computing the FFT, and the result cropped after the inverse FFT.

- Added `dip::Convolution()`, a gateway to the three implementations of the convolution.

- `dip::GeneralConvolution()` accepts a complex-valued filter kernel.

- `dip::Kernel`, `dip::PixelTable` and `dip::PixelTableOffsets` can now represent complex-valued weights. 
  Added `dip::Kernel::HasComplexWeights()`, `dip::PixelTable::WeightsAreComplex()` and
  `dip::PixelTableOffsets::WeightsAreComplex()`.

### Changed functionality

- `dip::testing::Timer::CpuResolution()` and `WallResolution()` are static members.

- Stream output operator for `dip::testing::Timer` decouples the two values, choosing appropriate units
  for them independently.

- Added some more `[[nodiscard]]` annotations.

- Updated included *libics*, which fixes some potential bugs, and adds support for metadata related
  to multi-detector microscopes.

### Bug fixes

- `dip::ResampleAt(in, map)` didn't copy the color space information from the input image to the output image.

- `dip::GeneralConvolution()` skips zero pixels in the kernel image, as was described in the documentation.
  This makes the operation significantly faster if the kernel has many zero pixels.

- `dip::ExpandTensor()` didn't write to an existing data segment in its output image, even if it was allocated
  to the right dimensions.

- `dip::BoundaryCondition::ASYMMETRIC_MIRROR` didn't work in some functions because it had the same value as
  `dip::BoundaryCondition::ALREADY_EXPANDED`.

- `dip::Histogram::Configuration`, when given an upper and lower bounds and a bin size, computed the number of
  bins in such a way that, for small difference between the bounds, the bin size was not respected.

- `dip::PixelTable::Mirror()` mirrored the footprint of the kernel, but not the weights.

- `dip::SeparateFilter()` yielded the complex conjugate of the correct 1D filter for all dimensions except the first.




## Changes to *DIPimage*

### New functionality

### Changed functionality

- Added a new parameter to `convolve`, which allows the selection of the computation method. The default, automatic
  selection of the method is now handled by *DIPlib* though the `dip::Convolution()` function. When computing
  through the Fourier domain, the periodic boundary condition is no longer the default; the operation should now
  produce the same result (up to rounding error) no matter which computation method is chosen.

(See also changes to *DIPlib*.)

### Bug fixes

(See also bugfixes to *DIPlib*.)




## Changes to *PyDIP*

### New functionality

- Added bindings for `dip::testing::Timer`.

- Added functions `dip.MeanRelativeError()` and `dip.MaximumRelativeError()`.

- Added `dip.Convolution()`, a gateway to the three implementations of the convolution.

### Changed functionality

(See also changes to *DIPlib*.)

### Bug fixes

(See also bugfixes to *DIPlib*.)




## Changes to *DIPviewer*

### New functionality

### Changed functionality

### Bug fixes




## Changes to *DIPjavaio*

### New functionality

### Changed functionality

### Bug fixes
