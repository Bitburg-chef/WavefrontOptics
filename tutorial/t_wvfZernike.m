% t_Zernike
%
% Representing wavefront aberrations using Zernike polynomials.
%
% This tutorial explains a method of representing the wavefront aberration
% function using a set of functions known as Zernike polynomials. The
% wavefront aberration function models the effect of the human cornea,
% lens, and pupil on the optical wavefront propagating through them.
% Absorption is modeled by an amplitude < 1, and phase aberrations are
% modeled by a complex phasor of the form exp(i*2*pi*[summation of Zernike
% polynomials]/wavelength). From Fourier optics, the eye's point spread
% function (PSF) can be computed from the wavefront aberration function, or
% pupil function, by taking the Fourier transform: PSF = fft2(pupil
% function).
%
% The Zernike polynomials form an orthogonal basis set over a unit disk.
% They are useful because they can isolate aberrations into separate
% components, each of which is given a weight and has potential for being
% corrected. For example, rather than seeing an entire aberrated wavefront,
% we can instead look at the amount of astigmatism in the 45 degree
% direction and how it contributes to the PSF on its own by knowing the
% measured Zernike coefficient associated with it.
%
% The tutorial introduces the concept of Zernike polynomials; shows the pupil
% function and how it is formed using Zernkie polynomials; shows the
% associated point-spread functions for given pupil functions; demonstrates
% and explains longitudinal chromatic aberration; 
% demonstrates and explains Stiles-Crawford effect; looks at measured human
% data and shows how eyeglasses only allow us to correct certain wavefront
% aberrations.
%
% See:
%  http://white.stanford.edu/teach/index.php/Wavefront_optics_toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class:     Psych 221/EE 362
% Tutorial:  Wavefront Toolbox - Zernike Polynomials
% Author:    BW, MDL, KP
% Purpose:   
%
% History:
%   3/12/2012	  Created.
%   4/29/12  dhb  Tried to tighten vertical spacing
%                 and apply uniform convention.  Some editing of
%                 text for clarity.
%	
% (c) Wavefront Toolbox Team 2011, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTES
%   a) Add and call a pathdef file which includes the wavefront toolbox?
%   See function AiryPattern in PTB, or perhaps function airy that comes
%   with Matlab.
%
%   b) Introduce notion of a figure of merit for quality of PSF, and
%   compute some explicitly for various cases considered.
%
%   c) The fact that for an aberrated eye, the best optical quality does
%   not occur when nominal defocus wl matches the calculated wavelength is
%   not considered here, but can be quite important when thinking about
%   real optical quality.
%

%% Zernike polynomials
%
% Zernike polynomials consist of:
%   -a weighting coefficient (Zernike coeff),
%   -a normalization factor,
%   -a radially-dependent polynomial,
%   -an azimuthally-dependent sinusoid
% 
% For example, one such polynomial (known as Astigmatism along 45 degrees) 
% is given by:
%
%          Z3 * sqrt(6) * rho^2 * cos(2*theta)
%
% where rho and theta are natural polar coordinates representing radial
% norm and angle on a disk. These can be easily converted to rectangular
% coordinates, making the Zernike polynomial representation useful for
% computing the wavefront aberrations.
%
% Zernike polynomials can be expressed using two indices, one representing
% the highest order of the radial polynomial term (n), and the other
% representing the frequency of the azimuthal sinusoid (m). (See wiki for
% more details.)
%
% The polynomials can also be represented in a single-indexing scheme (j)
% using OSA standards, which is easier to manage in vector form for Matlab,
% so we will use it here.
%
% Each radial order has (order+1) number of polynomial terms (and
% therefore, order+1 number of Zernike coefficients to represent them).
% Thus, 0th order has 1 term, 1st order has 2 terms, etc.
%
% In this tutorial we will be working with up to 10 orders of Zernike
% polynomials. Counted out, this represents 1+2+3+...+11 = 66 terms for
% orders 0 through 10. We will actually treat this as 65 coefficients since
% the 0th order term (j=0) is constant. Additionally, the 1st order terms
% (coeffs 1 and 2, known as tip and tilt) only serve to shift the PSF along
% the x or y axis, and don't represent true aberrations. They will be left
% at 0 for this tutorial, although they still make up the 65 coefficients,
% unlike the j=0 term, which is removed from the index (which is also
% convenient because Matlab indexes starting from 1!)

%% Clear
s_initISET

% The tutorial only uses 1 wavelength at a time. So, for plotting, we use
% this index.
waveIdx = 1;
maxMM = 2;
maxUM = 20;      
pupilfuncrangeMM = 5;

%% Use Zernike polynomials to specify a diffraction limited PSF.
%
% Use wvfCreate to create a wavefront variable to explore with.
%
% This wavefront by default has the 0's for all zernike coeffs
% Notice that the calcpupilMM is by default 3, meaning we are simulating
% the wavefront PSF for a pupil of 3MM diameter.  This code dumps
% out the structure so you can get a sense of what is in it.
% 
% The validation script v_wvfDiffractionPSF compares the diffraction
% limited PSFs obtained in this manner with those obtained by
% computing an Airy disk and shows that they match.
wvf0 = wvfCreate                    

% Look at the plot of the normalized PSF within 1 mm of the center.
% Variable maxUM is used to specify size of plot from center of the PSF.
%
% The plot shows an airy disk computed from the Zernike polynomials; that
% is representing the diffraction-limited PSF obtained when the Zernike
% coefficients are all zero.
%
% You might wonder, where does the psf get computed.  The answer is that
% wvfPlot calls wvfGet, which computes what it needs to get desired
% quantities, as necessary.
vcNewGraphWin;
wvfPlot(wvf0,'2dpsf space normalized','um',waveIdx,maxUM);

%% Examine how the first non-zero Zernike coefficient contributes to the PSF.

% The 3rd term (or 4th term when counting the 0th order constant as in the
% j indexing scheme) is known as astigmatism with axis at 45 degrees.
%
% We start with the default structure , wvf0 created above, and set the
% zcoeff column vector to be a new non-zero vector.
zcoeffs = zeros(65,1);
zcoeffs(3) = 0.75;                              % Just a non-zero weight
wvf3 = wvfSet(wvf0,'zcoeffs',zcoeffs);

% Look at the pupil function for astigmatism with axis at 45 degrees.
%
% We have used wvfComputePupilFunction separately here, but it is actually also
% contained within wvfComputePSF, which we will use from now on.
wvf3 = wvfGet(wvf3,'pupil function');

% Now we plot the pupil function, which captures phase information about
% the wavefront aberrations.
%
% We can see that the phase changes seem to be aligned with the + and - 45
% degree axes. 
vcNewGraphWin;
wvfPlot(wvf3,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);

%% Plot a PSF

% While the pupil functions are well specified by Zernike polynomials, it's
% hard to get meaning from them. We'd much prefer to look at the PSF, which
% gives us an idea of how the pupil will blur an image. 
% This is essentially done by applying a Fourier Transform to the pupil
% function.
wvf3 = wvfComputePSF(wvf3); 

% Now we can plot the normalized PSF for a pupil only whose only aberration
% is the 45 degree astigmatism.
%
% As you can see, this no longer looks like the narrower
% diffraction-limited PSF. It has also lost its radial symmetry. We will
% see that the higher the order of Zernike polynomial, the more complex the
% associated PSF will be.
vcNewGraphWin;
wvfPlot(wvf3, '2d psf space normalized','um',waveIdx,maxUM);

%% Examine effect of  the 5th coeff (j index 6), which is astigmatism
% along the 0 or 90 degree axis.
%
% We can see that unlike the 3rd coefficient, this coefficient for
% astigmatism is aligned to the x and y axes
zcoeffs = zeros(65,1);                            
zcoeffs(5) = 0.75;                              % Just a non-zero weight
wvf5 = wvfSet(wvf0,'zcoeffs',zcoeffs);
wvf5 = wvfComputePSF(wvf5);

vcNewGraphWin;
wvfPlot(wvf5,'2d pupil phase space','mm',waveIdx,maxMM);
vcNewGraphWin;
wvfPlot(wvf5,'2d psf space normalized','um',waveIdx,maxUM);

%% Go wild and make plots of various pupil functions and their respective
% point-spread functions for different Zernike polynomials of 2nd and 3rd orders
% (j index 3 through 9)
wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'calculated pupil',wvfGet(wvf0,'measured pupil','mm'));
pupilfuncrangeMM = 4;
jindices = 3:9;
maxMM = 4; 
for ii = jindices
    vcNewGraphWin;
    zcoeffs = zeros(65,1);
    zcoeffs(ii) = 0.75;
    wvf = wvfSet(wvf0,'zcoeffs',zcoeffs);
    wvf = wvfComputePSF(wvf);

    subplot(2,1,1);
    wvfPlot(wvf,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);

    subplot(2,1,2);
    wvfPlot(wvf,'2d psf space','mm',waveIdx,maxMM);
end

%% How longitudinal chromatic aberration (LCA) affects the PSF / "Defocus"
%
% What happens if we want to know how the PSF looks for different
% wavelengths? You may have learned that optical systems can have chromatic
% aberration, where one wavelength is brought into focus but others may be
% blurry because they are refracted closer or farther from the imaging
% plane. In this case, the PSF is dependent on wavelength.

% We can set this using the  "in-focus wavelength" of our wvf.
% This code indicates that the data is given for a nominal focus of 550 nm,
% which is also the default in wvfCreate.  We also now explictly set the
% wavelength for which the PSF is calculated to 550 nm (this is also the
% default.  
wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'measured wl',550); 
wvf0 = wvfSet(wvf0,'wavelength',550); 

% It turns out that all aberrations other than "Defocus" are known
% to vary only slightly with wavelength. As a result, the Zernike
% coefficients don't have to be modified, apart from one.
% zcoeff(4) is the "Defocus" of a pupil function. It is what typical eyeglasses
% correct for using + or - diopters lenses. The wavefront toolbox combines
% the longitudinal chromatic aberration (LCA) into this coefficient.
wvf0 = wvfComputePSF(wvf0);
vcNewGraphWin;
maxMM = 3; 
wvfPlot(wvf0,'1dpsfspacenormalized','mm',waveIdx,maxMM);
hold on;

% Keep the calculated wavelength at default 550 nm but change
% the nominal in-focus wavelength, then add to the plot.
%
% The new psf is wider due to longitudinal chromatic aberration, even though
% it's still just the diffraction-limited wavefront function (the Zernike
% coefficients are still 0).%
%
% To put it another way, the code produces the
% PSF at the specified wavelengths (set explicitly to just 550 above) given
% that the wavelength of nominal focus is as specified.  Since the two differ
% here, we see the effect of LCA.
wvf1 = wvfCreate;
lcaDiopters = wvfLCAFromWavelengthDifference(wvfGet(wvf1,'measured wl'),600);
wvf1 = wvfSet(wvf1,'calc observer focus correction',-lcaDiopters);
wvf1 = wvfComputePSF(wvf1);
wvfPlot(wvf1,'1dpsfspacenormalized','mm',waveIdx,maxMM);

%% Verify that The LCA effect is contained solely within the Defocus
% coefficient.
%
% First we explicitly computethe defocus implied by the wavelength difference above.
% Function wvfGetDefocusFromWavelengthDifference takes in a wavefront which
% contains information about the specified nominal wavelength and calculated
% wavelength. It computes the defocus in diopters from using unmatched wavelengths, 
% and returns the defocus converted into microns. This is important because Zernike
% coefficients are assumed to be given in microns.
defocusMicrons = wvfGetDefocusFromWavelengthDifference(wvf1);

% Make a new wavefront which does not have the mismatched wavelengths.
% Because these both default to 550 nm, they match here.
wvf2 = wvfCreate;
nominalWl = wvfGet(wvf2,'nominalfocuswl')
calcWl = wvfGet(wvf2,'wavelength')

% Make our adjustment purely to the j=4 Zernike
% coefficient (you may remember from our earlier plotting that this term on
% its own has a radially symmetric PSF which widens the diffraction limited
% PSF). We'll plot this PSF with a thinner blue line and overlay it.
%
% % The two aberrated plots are identical. The defocus of a pupil can be
% measured separately, whether using Zernike coeffs or in diopters, but any
% chromatic aberration is added solely into this coefficient.
zcoeffs = zeros(65,1);
zcoeffs(4) = defocusMicrons;
wvf2 = wvfSet(wvf2,'zcoeffs',zcoeffs);
wvf2 = wvfComputePSF(wvf2);
[udataS, pData] = wvfPlot(wvf2,'1dpsfspacenormalized','mm',waveIdx,maxMM);
set(pData,'color','b','linewidth',2);

%%  How cone geometry affects the PSF: the Stiles-Crawford effect (SCE)
%
% The cones that line our retinas are tall rod-shaped cells. They act like
% waveguides, such that rays parallel to their long axis excite the photoreceptors
% more readily than rays that are travelling at skewed angles. This has the
% benefit of reducing the chance that scattered or aberrated rays will
% excite the cone cells. Although this effect physically comes from the
% retina, it can be modeled using the pupil function discussed above. The
% amplitude of the pupil function is altered, such that it decays in the
% form exp(-alpha*((x-x0)^2+(y-y0)^2)). This physically attenuates rays
% that enter near the edges of the pupil and lens. Since the phase aberration at
% the edges of the pupil is usually most severe, SCE can actually improve vision. 
% Note that generally the position of the pupil with highest transmission 
% efficiency does not lie in the exact center of the pupil (x0,y0 are nonzero).

% Begin with an unaberrated pupil and see what its amplitude and phase look
% like. We'll also plot the associated diffraction-limited PSF.
wvf0 = wvfCreate;
wvf0 = wvfComputePSF(wvf0);
maxMM = 2; %MM from the center of the PSF
pupilfuncrangeMM = 5;
vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf0,'2d pupil amplitude space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,2);
wvfPlot(wvf0,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,3:4);
wvfPlot(wvf0,'2d psf space','mm',waveIdx,maxMM);

% To this unaberrated pupil function, we add the Stiles-Crawford
% parameters, as measured by Berendshot et al. (see sceCreate for
% details). This adds a decaying exponential amplitude to the pupil
% function, causing less light to be transmitted to the retina.
%
% Compare the diffraction-limited PSF without SCE to the one with SCE. What
% are the differences? Is the amplitude different? Why? Is the width of the
% PSF different? Why?
wvf0SCE = wvfSet(wvf0,'sceParams',sceCreate(wvfGet(wvf0,'wave'),'berendshot'));
wvf0SCE = wvfComputePSF(wvf0SCE);
vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf0SCE,'2d pupil amplitude space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,2);
wvfPlot(wvf0SCE,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,3:4);
wvfPlot(wvf0SCE,'2d psf space','mm',waveIdx,maxMM);

% Compare the above with how the SCE affects an aberrated PSF. Let's create a
% PSF with moderate astigmatism along the xy axes.
zcoeffs = zeros(65,1);
zcoeffs(5) = 0.75;
wvf5 = wvfSet(wvf0,'zcoeffs',zcoeffs);
wvf5 = wvfComputePSF(wvf5);
vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf5,'2d pupil amplitude space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,2);
wvfPlot(wvf5,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,3:4);
wvfPlot(wvf5,'2d psf space','mm',waveIdx,maxMM);

% Add SCE to the aberrated pupil function.
%
% Compare the two aberrated PSFs. How do their peak amplitudes compare?
% How do their widths compare? How did the symmetry of the PSF change?
% Which PSF would create a "better image" on the retina?
wvf5SCE = wvfSet(wvf5,'sceParams',sceCreate(wvfGet(wvf5,'wave'),'berendshot'));
wvf5SCE = wvfComputePSF(wvf5SCE);
vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf5SCE,'2d pupil amplitude space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,2);
wvfPlot(wvf5SCE,'2d pupil phase space','mm',waveIdx,pupilfuncrangeMM);
subplot(2,2,3:4);
wvfPlot(wvf5SCE,'2d psf space','mm',waveIdx,maxMM);

%% Wavefront measurements of human eyes and the effects of single-vision
% corrective eyeglasses 
%
% We have access to measurements of the pupil function of real human eyes. The
% optics of these eyes are not perfect, so they have interesting pupil functions
% and PSF shapes.

% Set up the wvf structure
measMM = 6;
calcMM = 3;
maxMM = 3;
theWavelengthNM = 550;
wvfHuman0 = wvfCreate('measured pupil',measMM,'calculated pupil',calcMM);
wvfHuman0 = wvfSet(wvfHuman0,'wavelength',theWavelengthNM);

% Load in some measured data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = load(sDataFile);
whichSubjects = [3 7];
theZernikeCoeffs = theZernikeCoeffs(:,whichSubjects);
nSubjects = size(theZernikeCoeffs,2);
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Stiles Crawford
wvfHuman0 = wvfSet(wvfHuman0,'sceParams',sceCreate(wvfGet(wvfHuman0,'wave'),'berendshot'));

% Plot subject PSFs, one by one
for ii = 1:nSubjects
    fprintf('** Subject %d\n',whichSubjects(ii))

    wvfHuman = wvfSet(wvfHuman0,'zcoeffs',theZernikeCoeffs(:,ii));
    wvfHuman = wvfComputePSF(wvfHuman);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfHuman,'2d pupil amplitude space','mm',waveIdx,calcMM);
    subplot(2,2,2);
    wvfPlot(wvfHuman,'2d pupil phase space','mm',waveIdx,calcMM);
    subplot(2,2,3:4);
    wvfPlot(wvfHuman,'2d psf space','mm',waveIdx,maxMM);
end

%% Single-vision eyewear generally corrects only the lowest-order
% Zernike aberrations (defocus given in diopters) and astigmatism (cylinder
% correction also given in diopters). The Zernike coefficients give us an
% easy and convenient way to simulate corrective lenses; we can simply set
% those Zernike coefficients to zero and see what the PSFs look like!
%
% Plot their corrected PSFs, one by one, How do the corrected PSFs compare
% to the uncorrected ones? their peaks? their widths?
%
% Try changing the whichSubjects array above to look at other sample data. Do
% eyeglasses help correct the aberrations in those subjects?
%
% If you were to spend thousands of dollars on laser eye surgery, would you
% want them to only correct the first order of wavefront aberrations, like
% eyeglasses, or do a full wavefront measurement?
% 
% Suppose you knew that such surgery would correct some of the lower order aberrations but some of the
% higher order aberrations worse.  How would you compute the net effect of
% something like that?
for ii = 1:nSubjects
    fprintf('** Subject %d corrected\n',whichSubjects(ii))
    
    % Correct defocus and astigmatism
    zCoeffs = theZernikeCoeffs(:,ii);
    zCoeffs(3:5) = 0;
    wvfHuman = wvfSet(wvfHuman0,'zcoeffs',zCoeffs);
    wvfHuman = wvfComputePSF(wvfHuman);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfHuman,'2d pupil amplitude space','mm',waveIdx,calcMM);
    subplot(2,2,2);
    wvfPlot(wvfHuman,'2d pupil phase space','mm',waveIdx,calcMM);
    subplot(2,2,3:4);
    wvfPlot(wvfHuman,'2d psf space','mm',waveIdx,maxMM);
end

%% End

