%GABOR   Gabor filter
%
% SYNOPSIS:
%  image_out = gabor(image_in,sigma,frequency,direction,boundaryCondition,process,truncation)
%  image_out = gabor(image_in,sigma,frequencies,boundaryCondition,process,truncation)
%
% PARAMETERS:
%  sigma:       Sigma in the spatial domain
%  frequency:   Magnitude of the frequency in pixel, [0,0.5]
%  direction:   Direction in the fourier domain where to filter, [0,2pi]
%  frequencies: Frequencies in cartesian coordinates, [-0.5,0.5]
%  boundary_condition: Defines how the boundary of the image is handled.
%                      See HELP BOUNDARY_CONDITION
%  process:     Dimensions to process, set to [] to process all dimensions
%  truncation:  Determines the size of the Gaussian filters.
%
% DEFAULTS:
%  bounary_condition = 'mirror'
%  process = []
%  truncation = 3
%
% NOTES:
%  The first form is specific to 2D images, and specifies the frequencies
%  in polar coordinates.
%
%  The second form is for arbitrary dimensionality, and specifies the
%  frequencies as a vector with coordinates in the frequency domain.
%
%  The coordinate system in the Fourier domain is as follows:
%  The origin is at IMSIZE(IMAGE_IN)/2. Positive frequencies are to the
%  right and down. The angle is computed in mathematical sense, with
%  0 degrees along the positive x-axis (to the right), and 90 degrees
%  along the positive y-axis (downwards).
%
% DEFAULTS:
%  sigma = 5
%  frequency = 0.15
%  direction = pi
%
%
% DIPlib:
%  This function calls the DIPlib functions dip::GaborIIR.

% (c)2018, Cris Luengo.
% Based on original DIPlib code: (c)1995-2014, Delft University of Technology.
% Based on original DIPimage code: (c)1999-2014, Delft University of Technology.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
