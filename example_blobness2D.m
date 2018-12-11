%% Enhancement of a synthetic 2D blob

% load input image
I = rgb2gray(imread('blob2d.png'));

% preprocess the input a little bit
Ip = single(I);
thr = prctile(Ip(Ip(:)>0),1) * 0.9;
Ip(Ip<=thr) = thr;
Ip = Ip - min(Ip(:));
Ip = Ip ./ max(Ip(:));    

% compute enhancement for two different tau values
B1 = blobness2D(Ip, 2.5:2.5:7.5, [1;1], 1, true);
B2 = blobness2D(Ip, 2.5:2.5:7.5, [1;1], 0.5, true);

% display result
figure; 
subplot(2,2,1)
imshow(I)
title('input image')
axis image

subplot(2,2,3)
imshow(B1)
title('Filter enhancement (tau=1)')
axis image

subplot(2,2,4)
imshow(B2)
title('Filter enhancement (tau=0.5)')
axis image