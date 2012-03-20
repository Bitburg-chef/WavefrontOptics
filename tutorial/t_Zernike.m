% t_Zernike
% Tutorial on representing wavefront aberrations using Zernike polynomials
%
% Class:     Psych 221/EE 362
% Tutorial:  Wavefront Toolbox - Zernike Polynomials
% Author:    BW, MDL, KP
% Purpose:   Introduce the concept of Zernike polynomials; show the pupil
% function and how it is formed using Zernkie polynomials; show the
% associated point-spread functions for given pupil functions; demonstrate and 
% explain Stiles-Crawford effect; look at measured human data and show how 
% eyeglasses only allow us to correct certain wavefront aberrations
%
%
% Date:      3/12/2012	
% Duration:  --
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Things to include: 
% - a pathdef file which includes the wavefront toolbox? - eh, probably
% will be just included in ISET later
% - anything else??
%
%
% This tutorial is meant to explain a method of representing the wavefront
% aberration function using a set of functions known as Zernike polynomials.  
% The wavefront aberration function models the effect of the human cornea,
% lens, and pupil on the optical wavefront propagating through them.
% Absorption is modeled by an amplitude < 1, and phase aberrations are
% modeled by a complex phasor of the form exp(i*2*pi*[summation of Zernike polynomials]/wavelength).
% From Fourier optics, the eye's point spread function (PSF) can be computed 
% from the wavefront aberration function, or pupil function, by taking the
% Fourier transform: PSF = fft2(pupil function).
%
% The Zernike polynomials form an orthogonal basis set over a unit disk.
% They are useful because they can isolate aberrations into separate
% components, each of which is given a weight and has potential for being
% corrected. For example, rather than seeing an entire aberrated wavefront,
% we can instead look at the amount of astigmatism in the 45 degree
% direction and how it contributes to the PSF on its own by knowing the
% measured Zernike coefficient associated with it.
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
% Z3 * sqrt(6) * rho^2 * cos(2*theta)
% where rho and theta are natural polar coordinates representing radial norm
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
% Each radial order has (order+1) number of polynomial terms (and therefore, order+1
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

wvf3 = wvfSet(wvf0,'zcoeffs',Zcoeffs);
% Here we have started with the default, wvf0, and chosen to set the
% zcoeff column vector to be our new non-zero vector. 

% Before we look at the PSF, let's look at the pupil function for 
% astigmatism with axis at 45 degrees. 

wvf3 = wvfComputePupilFunction(wvf3);
% We have used this function separately here, but it is actually also
% contained within wvfComputePSF, which we will use from now on.

% Now we plot the pupil function, which captures phase information about
% the wavefront aberrations

pupilfuncrange = 5;

vcNewGraphWin;
wvfPlot(wvf3,'2d pupil phase space','mm',pupilfuncrange);

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
wvfPlot(wvf5,'2d pupil phase space','mm',3);
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
wvf0 = wvfSet(wvf0,'calculated pupil',wvfGet(wvf0,'measured pupil','mm'));

jindices = 3:9;  %3:14;
% figure;
PupilFuncRange = 4;  %2;
maxMM = 4; 
for ii = jindices
    vcNewGraphWin;
    zcoeffs = zeros(65,1);
    zcoeffs(ii) = 0.75;
    wvf = wvfSet(wvf0,'zcoeffs',zcoeffs);
    wvf = wvfComputePSF(wvf);

    subplot(2,1,1);
    wvfPlot(wvf,'2d pupil phase space','mm',PupilFuncRange);

    subplot(2,1,2);
    wvfPlot(wvf,'2d psf space','mm',maxMM);
%     psf = wvfGet(wvf,'psf');
%     imagesc(psf);
end


%_______________________________________________________________________
%%
% How chromatic aberration affects the PSF / "Defocus"
%

clear all; close all;

% What happens if we want to know how the PSF looks for different wavelengths?
% You may have learned that optical systems can have chromatic aberration,
% where one wavelength is brought into focus but others may be blurry
% because they are refracted closer or farther from the imaging plane. In
% this case, the PSF is dependent on wavelength. 

% We can always set this into "in-focus" wavelength into our wvf:

wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'nominalfocuswl',550); 
% This indicates that the data is given for a nominal focus of 550nm, which
% is also the default in wvfCreate. 

% It turns out that all aberrations other than "Defocus" have variations
% with wavelength that are known to be small. As a result, the Zernike
% coefficients don't have to be modified, apart from one.
% zcoeff(4) is the "Defocus" of a pupil. It is what typical eyeglasses
% correct for using + or - diopters lenses. The wavefront toolbox combines
% the longitudinal chromatic aberration (LCA) into this coefficient.

% We can first make a diffraction-limited PSF from wvf0:
wvf0 = wvfComputePSF(wvf0);
vcNewGraphWin;
maxMM = 3; 
wvfPlot(wvf0,'1dpsfspacenormalized','mm',maxMM);
hold on;

% Let's keep the calculated wavelength at default 550nm but change
% the in focus wavelength:

wvf1 = wvfCreate;
wvf1 = wvfSet(wvf1,'nominalfocuswl',600); %sets nominal wavelength to 600nm
wvf1 = wvfComputePSF(wvf1);
wvfPlot(wvf1,'1dpsfspacenormalized','mm',maxMM);

% The new plot looks wider due to the chromatic aberration, even though
% it's still just the diffraction-limited wavefront function (the Zernike
% coefficients are still 0).

% Let's see how this LCA effect is contained solely within the Defocus
% coefficient.

defocusMicrons = wvfGetDefocusFromWavelengthDifference(wvf1);
% This function takes in a wavefront which contains information about the
% specified nominal wavelength and calculated wavelength. It computes the
% defocus in diopters from using unmatched wavelengths, and returns the
% defocus converted into microns. This is important because Zernike
% coefficients are assumed to be given in microns.

% Now that we can make a new wavefront which does not have the mismatched
% wavelengths:
wvf2 = wvfCreate;
nominalWl = wvfGet(wvf2,'nominalfocuswl')
calcWl = wvfGet(wvf2,'wave')
% We can see that unlike before, the two wavelengths are identical.

% Instead, we will make our adjustment purely to the j=4 Zernike
% coefficient (you may remember from our earlier plotting that this term on
% its own has a radially symmetric PSF which widens the diffraction limited
% PSF). We'll plot this PSF with a thinner blue line and overlay it.

Zcoeffs = zeros(65,1);
Zcoeffs(4) = defocusMicrons;

wvf2 = wvfSet(wvf2,'zcoeffs',Zcoeffs);
wvf2 = wvfComputePSF(wvf2);
[udataS, pData] = wvfPlot(wvf2,'1dpsfspacenormalized','mm',maxMM);
set(pData,'color','b','linewidth',2);

% The two aberrated plots are identical. The defocus of a pupil can be
% measured separately, whether using Zernike coeffs or in diopters, but any
% chromatic aberration is added solely into this coefficient.


%_______________________________________________________________________
%%
% How cone geometry affects the PSF: the Stiles-Crawford effect (SCE)
%

clear all; close all;

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

% Let's began with an unaberrated pupil and see what its amplitude and phase look
% like. We'll also plot the associated diffraction-limited PSF.

wvf0 = wvfCreate;
wvf0 = wvfComputePSF(wvf0);
maxMM = 2; %MM from the center of the PSF
PupilFuncRange = 5;

vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf0,'2d pupil amplitude space','mm',PupilFuncRange);
subplot(2,2,2);
wvfPlot(wvf0,'2d pupil phase space','mm',PupilFuncRange);
subplot(2,2,3:4);
wvfPlot(wvf0,'2d psf space','mm',maxMM);

% To this unaberrated pupil function, we'll add the Stiles-Crawford
% parameters, as measured by Berendshot et al. (see sceCreate for
% details). This adds a decaying exponential amplitude to the pupil
% function, causing less light to be transmitted to the retina.

wvf0SCE = wvfSet(wvf0,'sceParams',sceCreate(wvfGet(wvf0,'wave'),'berendshot'));
wvf0SCE = wvfComputePSF(wvf0SCE);

vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf0SCE,'2d pupil amplitude space','mm',PupilFuncRange);
subplot(2,2,2);
wvfPlot(wvf0SCE,'2d pupil phase space','mm',PupilFuncRange);
subplot(2,2,3:4);
wvfPlot(wvf0SCE,'2d psf space','mm',maxMM);

% Compare the diffraction-limited PSF without SCE to the one with SCE. What
% are the differences? Is the amplitude different? Why? Is the width of the
% PSF different? Why?



% Now, let's compare how the SCE affects an aberrated PSF. Let's create a
% PSF with moderate astigmatism along the xy axes.

Zcoeffs = zeros(65,1); %reset our coeffs
Zcoeffs(5) = 0.75; %some non-zero weight

wvf5 = wvfSet(wvf0,'zcoeffs',Zcoeffs);
wvf5 = wvfComputePSF(wvf5);

vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf5,'2d pupil amplitude space','mm',PupilFuncRange);
subplot(2,2,2);
wvfPlot(wvf5,'2d pupil phase space','mm',PupilFuncRange);
subplot(2,2,3:4);
wvfPlot(wvf5,'2d psf space','mm',maxMM);

% We now add SCE to the aberrated pupil function.

wvf5SCE = wvfSet(wvf5,'sceParams',sceCreate(wvfGet(wvf5,'wave'),'berendshot'));
wvf5SCE = wvfComputePSF(wvf5SCE);

vcNewGraphWin;
subplot(2,2,1);
wvfPlot(wvf5SCE,'2d pupil amplitude space','mm',PupilFuncRange);
subplot(2,2,2);
wvfPlot(wvf5SCE,'2d pupil phase space','mm',PupilFuncRange);
subplot(2,2,3:4);
wvfPlot(wvf5SCE,'2d psf space','mm',maxMM);

% Now compare the two aberrated PSFs. How do their peak amplitudes compare?
% How do their widths compare? How did the symmetry of the PSF change?
% Which PSF would create a "better image" on the retina?


%_______________________________________________________________________
%%
% Studying wavefront measurements of human eyes and the effects of
% single-vision corrective eyeglasses
%

clear all; close all;

% In this wavefront toolbox, we have access to measurements of the pupil
% function of real human eyes. Many of these eyes are not perfect, so they
% have interesting pupil functions and PSF shapes.

measMM = 6;
calcMM = 3;
maxMM = 3;
% Set values in millimeters
wvfParams0 = wvfCreate('measured pupil',measMM,'calculated pupil',calcMM);

% load the patient data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = load(sDataFile);

whichSubjects = [3 7];
theZernikeCoeffs = theZernikeCoeffs(:,whichSubjects);
nSubjects = size(theZernikeCoeffs,2);
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Stiles Crawford
theWavelength = 550;
wvfParams0.sceParams = sceCreate(theWavelength,'berendshot');

% Now, let's plot their PSFs, one by one
for ii = 1:nSubjects
    fprintf('** Subject %d\n',ii)

    wvfParams = wvfSet(wvfParams0,'zcoeffs',theZernikeCoeffs(:,ii));
    wvfParams = wvfComputePSF(wvfParams);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfParams,'2d pupil amplitude space','mm',calcMM);
    subplot(2,2,2);
    wvfPlot(wvfParams,'2d pupil phase space','mm',calcMM);
    subplot(2,2,3:4);
    wvfPlot(wvfParams,'2d psf space','mm',maxMM);
end

% Single-vision eyewear generally corrects only the lowest-order
% Zernike aberrations (defocus given in diopters) and astigmatism (cylinder
% correction also given in diopters). The Zernike coefficients give us an
% easy and convenient way to simulate corrective lenses; we can simply set
% those Zernike coefficients to zero and see what the PSFs look like!

% Let's plot their corrected PSFs, one by one
for ii = 1:nSubjects
    fprintf('** Subject %d\n',ii)
    
    % correct defocus and astigmatism
    zCoeffs = theZernikeCoeffs(:,ii);
    zCoeffs(3:5) = 0;

    wvfParams = wvfSet(wvfParams0,'zcoeffs',zCoeffs);
    wvfParams = wvfComputePSF(wvfParams);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfParams,'2d pupil amplitude space','mm',calcMM);
    subplot(2,2,2);
    wvfPlot(wvfParams,'2d pupil phase space','mm',calcMM);
    subplot(2,2,3:4);
    wvfPlot(wvfParams,'2d psf space','mm',maxMM);
end

% How do the corrected PSFs compare to the uncorrected ones? their peaks?
% their widths?

% Try changing the whichSubjects array to look at other sample data. Do
% eyeglasses help correct the aberrations in those subjects?

% If you were to spend thousands of dollars on laser eye surgery, would you
% want them to only correct the first order of wavefront aberrations, like
% eyeglasses, or do a full wavefront measurement?

