function vesselness = vesselness3D(I, sigmas, spacing, tau, brightondark)
% calculates vesselness probability map (local tubularity) of a 3D input 
% image
% 
% vesselness = vesselness3D(V, sigmas, spacing, tau, brightondark)
% 
% inputs,
%   I : 3D image
%   sigmas : vector of scales on which the vesselness is computed
%   spacing : input image spacing resolution - during hessian matrix 
%       computation, the gaussian filter kernel size in each dimension can 
%       be adjusted to account for different image spacing for different
%       dimensions 
%   tau : (between 0.5 and 1) : parameter that controls response uniformity
%       - lower tau -> more intense output response            
%   brightondark: (true/false) : are vessels (tubular structures) bright on 
%       dark background or dark on bright (default for 3D is true)
%
% outputs,
%   vesselness: maximum vesselness response over scales sigmas
%
% example:
%   V = vesselness3D(I, 1:5, [1;1;1], 1, true);
%
% Function was written by T. Jerman, University of Ljubljana (October 2014)
% Based on code by D. Kroon, University of Twente (May 2009)

verbose = 1;

if nargin<5
    brightondark = true; % default
end

I(~isfinite(I)) = 0;
I = single(I);

for j = 1:length(sigmas)
    
    if verbose   
        disp(['Current Filter Sigma: ' num2str(sigmas(j)) ]);
    end
    
    [~, Lambda2, Lambda3] = volumeEigenvalues(I,sigmas(j),spacing,brightondark);
    if brightondark == true
        Lambda2 = -Lambda2;
        Lambda3 = -Lambda3;
    end
    
    % proposed filter      
    Lambda_rho = Lambda3;
    Lambda_rho(Lambda3 > 0 & Lambda3 <= tau .* max(Lambda3(:))) = tau .* max(Lambda3(:));
    Lambda_rho(Lambda3 <= 0) = 0;
    response = Lambda2.*Lambda2.*(Lambda_rho-Lambda2).* 27 ./ (Lambda2 + Lambda_rho).^3;    
    
    response(Lambda2 >= Lambda_rho./2 & Lambda_rho > 0) = 1;    
    response(Lambda2 <= 0 | Lambda_rho <= 0) = 0;
    response(~isfinite(response)) = 0;

    %keep max response
    if(j==1)
        vesselness = response;
    else        
        vesselness = max(vesselness,response);
    end
        
    clear response Lambda2 Lambda3 Lambda3M    
     
end

vesselness = vesselness ./ max(vesselness(:));
vesselness(vesselness < 1e-2) = 0;    


function [Lambda1, Lambda2, Lambda3] = volumeEigenvalues(V,sigma,spacing,brightondark)
% calculates the three eigenvalues for each voxel in a volume

% Calculate 3D hessian
[Hxx, Hyy, Hzz, Hxy, Hxz, Hyz] = Hessian3D(V,sigma,spacing);

% Correct for scaling
c=sigma.^2;
Hxx = c*Hxx; Hxy = c*Hxy;
Hxz = c*Hxz; Hyy = c*Hyy;
Hyz = c*Hyz; Hzz = c*Hzz;

% reduce computation by computing vesselness only where needed
% S.-F. Yang and C.-H. Cheng, “Fast computation of Hessian-based
% enhancement filters for medical images,” Comput. Meth. Prog. Bio., vol.
% 116, no. 3, pp. 215–225, 2014.
B1 = - (Hxx + Hyy + Hzz);
B2 = Hxx .* Hyy + Hxx .* Hzz + Hyy .* Hzz - Hxy .* Hxy - Hxz .* Hxz - Hyz .* Hyz;
B3 = Hxx .* Hyz .* Hyz + Hxy .* Hxy .* Hzz + Hxz .* Hyy .* Hxz - Hxx .* Hyy .* Hzz - Hxy .* Hyz .* Hxz - Hxz .* Hxy .* Hyz;

T = ones(size(B1));

if brightondark == true
    T(B1<=0) = 0;
    T(B2<=0 & B3 == 0) = 0;
    T(B1>0 & B2>0 & B1 .* B2 < B3) = 0;
else
    T(B1>=0) = 0;
    T(B2>=0 & B3 == 0) = 0;
    T(B1<0 & B2<0 & (-B1) .* (-B2) < (-B3)) = 0;
end
clear B1 B2 B3;

indeces = find(T==1);

Hxx = Hxx(indeces);
Hyy = Hyy(indeces);
Hzz = Hzz(indeces);
Hxz = Hxz(indeces);
Hyz = Hyz(indeces);
Hxy = Hxy(indeces);

% Calculate eigen values
[Lambda1i,Lambda2i,Lambda3i]=eig3volume(Hxx,Hxy,Hxz,Hyy,Hyz,Hzz);

% Free memory
clear Hxx Hyy Hzz Hxy Hxz Hyz;

Lambda1 = zeros(size(T));
Lambda2 = zeros(size(T));
Lambda3 = zeros(size(T));

Lambda1(indeces) = Lambda1i;
Lambda2(indeces) = Lambda2i;
Lambda3(indeces) = Lambda3i;

% some noise removal
Lambda1(~isfinite(Lambda1)) = 0;
Lambda2(~isfinite(Lambda2)) = 0;
Lambda3(~isfinite(Lambda3)) = 0;

Lambda1(abs(Lambda1) < 1e-4) = 0;
Lambda2(abs(Lambda2) < 1e-4) = 0;
Lambda3(abs(Lambda3) < 1e-4) = 0;


function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma,spacing)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma,spacing)
% 
% inputs,
%   I : The image volume, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%   spacing : input image spacing
%
% outputs,
%   Dxx, Dyy, Dzz, Dxy, Dxz, Dyz: The 2nd derivatives
%
% Function is written by D.Kroon University of Twente (June 2009)
% defaults
if nargin < 2, Sigma = 1; end

if(Sigma>0)
    %F=imbigaussian(Volume,Sigma,0.5);
    F=imgaussian(Volume,Sigma,spacing);
else
    F=Volume;
end

% Create first and second order diferentiations
Dz=gradient3(F,'z');
Dzz=(gradient3(Dz,'z'));
clear Dz;

Dy=gradient3(F,'y');
Dyy=(gradient3(Dy,'y'));
Dyz=(gradient3(Dy,'z'));
clear Dy;

Dx=gradient3(F,'x');
Dxx=(gradient3(Dx,'x'));
Dxy=(gradient3(Dx,'y'));
Dxz=(gradient3(Dx,'z'));
clear Dx;

function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');

[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
        
function I=imgaussian(I,sigma,spacing,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with 
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
% 
% Function is written by D.Kroon University of Twente (September 2009)

if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)

    % Filter each dimension with the 1D Gaussian kernels\
    x=-ceil(siz/spacing(1)/2):ceil(siz/spacing(1)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(1))^2)));
    H = H/sum(H(:));    
    Hx=reshape(H,[length(H) 1 1]);
    
    x=-ceil(siz/spacing(2)/2):ceil(siz/spacing(2)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(2))^2)));
    H = H/sum(H(:));    
    Hy=reshape(H,[1 length(H) 1]);

    x=-ceil(siz/spacing(3)/2):ceil(siz/spacing(3)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(3))^2)));
    H = H/sum(H(:));    
    Hz=reshape(H,[1 1 length(H)]);
    
    I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
end

