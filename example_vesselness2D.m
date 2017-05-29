%% Enhancement of the 2D retinal vasculature

% load input image
I = imread('fundus2D.png');

% compute enhancement for two different tau values
V1 = vesselness2D(I, 0.5:0.5:2.5, [1;1], 1, false);
V2 = vesselness2D(I, 0.5:0.5:2.5, [1;1], 0.5, false);

% display result
figure; 
subplot(2,2,1)
imshow(I)
title('input image')
axis image

subplot(2,2,3)
imshow(V1)
title('Filter enhancement (tau=1)')
axis image

subplot(2,2,4)
imshow(V2)
title('Filter enhancement (tau=0.5)')
axis image