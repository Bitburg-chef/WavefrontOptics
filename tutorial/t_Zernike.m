% t_Zernike
% Tutorial on representing wavefront aberrations using Zernike polynomials
%
% Class:     Psych 221/EE 362
% Tutorial:  Wavefront Toolbox - Zernike Polynomials
% Author:    BW, MDL, KP
% Purpose:   Introduce the concept of Zernike polynomials; show the pupil
% function and how it is formed using Zernkie polynomials; show the
% associated point-spread functions for given pupil functions; look at
% measured human data and show how Zernike polynomial representation of
% aberrations allow us to correct certain parts
%
%
% Date:      3/12/2012	
% Duration:  --
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Things to include: 
% - a pathdef file which includes the wavefront toolbox?
% - SCE explanation in tutorial
% - anything else??
%
%
% This tutorial is meant to explain a method of representing the wavefront
% aberration function using a set of functions known as Zernike polynomials.  
% The Zernike polynomials form an orthogonal basis set over a unit disk.
% They are useful because they can isolate aberrations into separate
% components, each of which is given a weight and has potential for being
% corrected. For example, rather than seeing an entire aberrated wavefront,
% we can instead look at the amount of astigmatism in the 45 degree
% direction and how it contributes to the PSF on its own by knowing the
% measured Zernike coefficient associated with it.
% The Zernike polynomials allow us to represent the pupil function, which 
% contains the phase deviation information of a wavefront.
% The Fourier Transform of this pupil function gives us the PSF.
%
%%
%________________________________________________________________________
% Some details about Zernike polynomials that will be helpful for
% understanding the rest of the tutorial:
%
% Zernike polynomials consist of:
% -a weighting coefficient (Zernike coeff),
% -a normalization factor,
% -a radially-dependent polynomial, and
% -an azimuthally-dependent sinusoid
% For example, one such polynomial (known as Astigmatism along 45degrees) 
% is given by:
% Z3 * sqrt(6) * ro^2 * cos(2*theta)
% where ro and theta are natural polar coordinates representing radial norm
% and angle on a disk. These can be easily converted to rectangular
% coordinates, making the Zernike polynomial representation useful for
% computing the wavefront aberrations. 
%
%
% Zernike polynomials can be expressed using two indices, one representing
% the highest order of the radial polynomial term (n), and the other 
% representing the frequency of the azimuthal sinusoid (m). 
% (See wiki for more details.)
% The polynomials can also be represented in a single-indexing scheme (j) 
% using OSA standards, which is easier to manage in vector form for Matlab,
% so we will use it here.
%
% Each order has order+1 number of polynomial terms (and therefore, order+1
% number of Zernike coefficients to represent them). Thus, 0th order has 1
% term, 1st order has 2 terms, etc. 
% In this tutorial we will be working with up to 10 orders of Zernike
% polynomials. Counted out, this represents 1+2+3+...+11 = 66 terms for
% orders 0 through 10. We will actually treat this as 65 coefficients since
% the 0th order term (j=0) is constant. Additionally, the 1st order terms 
% (coeffs 1 and 2, known as tip and tilt) only serve to shift the PSF along
% the x or y axis, and don't represent true aberrations. They will be left
% at 0 for this tutorial, although they still make up the 65 coefficients,
% unlike the j=0 term, which is removed from the index. 
% (Which is also convenient because Matlab indexes starting from 1!)
%_________________________________________________________________________
%%
% Looking at how each Zernike polynomial specifies a pupil function and a
% point spread function


clear all; close all;

% We start by specifying a column vector of Zernike polynomial weighting
% coefficents:

Zcoeffs = zeros(65,1);

% This represents an unaberrated wavefront, one whose PSF is just the
% diffraction-limited PSF.
% We can create a wavefront function to see this.

wvf0 = wvfCreate %creates the default wavefront variable and its params
% This wavefront already by default has the 0's for all zernike coeffs
% Notice that the calcpupilMM is by default 3, meaning we are simulating
% the wavefront PSF for a pupil of 3MM diameter

% We now compute the PSF for this wavefront
wvf0 = wvfComputePSF(wvf0); %computes the PSF and stores it back in

% We can look at the plot of the normalized PSF within 1 mm of the center
vcNewGraphWin;
maxMM = 1; %MM from the center of the PSF
wvfPlot(wvf0,'2dpsf space normalized','mm',maxMM);

%
% perhaps include 1d comparison to diffraction limited psf from
% airydisk 
%

% The plot shows an airy disk, representing the diffraction-limited PSF
% when the Zernike coefficients are all zero. 

% Now we can look at the first non-zero Zernike coefficient and see how it
% independently contributes to the PSF.

Zcoeffs(3) = 0.75; %just a non-zero weight

% This 3rd term (or 4th term when counting the 0th order constant as in the
% j indexing scheme) is known as astigmatism with axis at 45 degrees.

% We create a new wvf with this coefficient weighting

wvf3 = wvfSet(wvf0,'zcoeffs',Zcoeff);
% Here we have started with the default, wvf0, and chosen to set the
% zcoeff column vector to be our new non-zero vector. 

% Before we look at the PSF, let's look at the pupil function for 
% astigmatism with axis at 45 degrees. 

wvf3 = wvfComputePupilFunction(wvf3);
% We have used this function separately here, but it is actually also
% contained within wvfComputePSF, which we will use from now on.

% Now we plot the pupil function, which captures phase information about
% the wavefront aberrations

pupilfuncrange = 2;

vcNewGraphWin;
wvfPlot(wvf3,'2d pupil function space','mm',pupilfuncrange);

% We can see that the phase changes seem to be aligned with the + and - 45
% degree axes. Unfortunately, while these pupil functions are well
% specified by Zernike polynomials, it's hard to get meaning from them.
% We'd much prefer to look at the PSF, which gives us an idea of how the
% pupil will blur an image. 
% This is essentially done by applying a Fourier Transform to the pupil
% function.

wvf3 = wvfComputePSF(wvf3); 

% Now we can plot the normalized PSF for a pupil only whose only aberration
% is the 45 degree astigmatism.

vcNewGraphWin;
maxMM = 2;
wvfPlot(wvf3, '2d psf space normalized','mm',maxMM);

% As you can see, this no longer looks like the narrower
% diffraction-limited PSF. It has also lost its radial symmetry. We will
% see that the higher the order of Zernike polynomial, the more complex the
% associated PSF will be.

% We can quickly look at the 5th coeff (j index 6), which is astigmatism
% along the 0 or 90 degree axis.

Zcoeffs = zeros(65,1); %reset our coeffs
Zcoeffs(5) = 0.75; %some non-zero weight

wvf5 = wvfSet(wvf0,'zcoeffs',Zcoeffs);
wvf5 = wvfComputePSF(wvf5);

vcNewGraphWin;
wvfPlot(wvf5,'2d pupil function space','mm',3);
% We can see that unlike the 3rd coefficient, this coefficient for
% astigmatism is aligned to the x and y axes

vcNewGraphWin;
wvfPlot(wvf5,'2d psf space normalized','mm',maxMM);

%_______________________________________________________________________
%%
% 
% Some plots of pupil functions and their respective point-spread functions
% for different Zernike polynomials of 2nd and 3rd orders
% (j index 3 through 9)

clear all; close all;

wvf0 = wvfCreate;

jindices = 3:9;
% figure;
PupilFuncRange = 2; 
maxMM = 3; 
for ii = jindices
    vcNewGraphWin;
    zcoeffs = zeros(65,1);
    zcoeffs(ii) = 0.75;
    wvf = wvfSet(wvf0,'zcoeffs',zcoeffs);
    wvf = wvfComputePSF(wvf);
%     subplot(2,7,ii-2);
    wvfPlot(wvf,'2d pupil function space','mm',PupilFuncRange);
%     subplot(2,7,ii+5);
%     wvfPlot(wvf,'2d psf space','mm',maxMM);
end

% this subplotting is really ugly. need to consider another way.

%_______________________________________________________________________
%
% How chromatic aberration affects the PSF
%

    







