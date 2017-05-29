%% Enhancement of the 3D cerebral vasculature

% load input image
load('volume.mat');

%rotate volume for visualization
I = permute(I,[3,1,2]);
I = I(end:-1:1,:,:);

%normalize input a little bit
I = I - min(I(:));
I = I / prctile(I(I(:) > 0.5 * max(I(:))),90);
I(I>1) = 1;

% compute enhancement for two different tau values
V = vesselness3D(I, 1:4, [1;1;1], 0.75, true);

% display result
figure; 
subplot(1,2,1)
% maximum intensity projection
imshow(max(I,[],3))
title('MIP of the input image')
axis image

subplot(1,2,2)
% maximum intensity projection
imshow(max(V,[],3))
title('MIP of the filter enhancement (tau=0.75)')
axis image