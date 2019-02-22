clear all
close all
clc

% params of filter mentioned in the paper
sj=0.03;
q=2;


% create random image for test n size of 401x401
I=randn(401,401);

% calculate centered frequcy of given image
F=fftshift(fft2(I));

figure;imagesc(I);
figure;mesh(abs(F));

% spatial domain filter (filter size nxn)
n=51;
[x y]=meshgrid(-(n-1)/2:(n-1)/2);
r=sqrt(x.^2+y.^2);
% calculate first part of filter
p1=2*pi*sj^2*(q*exp(-2*pi^2*q^2*sj^2*r.^2));
% then second part
p2=2*pi*sj^2*(exp(-2*pi^2*sj^2*r.^2));
% make sure both parts are normalized
p1=p1/sum(p1(:));
p2=p2/sum(p2(:));
% and set the spatial filter
sDoG=p1-p2;



% frequency domain filter
% create meshgrid in the same resolution of given image
[U V]=meshgrid(-fix(size(F,2)/2):fix(size(F,2)/2),-fix(size(F,1)/2):fix(size(F,1)/2));
% note the last index must be -0.5 to 0.5 hence we must use normalized
% frequency
U=0.5*U/max(U(:));
V=0.5*V/max(V(:));
% calculate ro value
p=sqrt(U.^2+V.^2);

fDoG=exp(-0.5*(p/(q*sj)).^2)-exp(-0.5*(p/sj).^2);


% filter image in frequency domain
FD=F.*fDoG;
% transform the filtered image in to spatial domain
fI=ifft2(ifftshift(FD));

% filter image in spatial domain
oI=conv2(I,sDoG,'same');


figure;
subplot(1,3,1);imagesc(fI);title('frequency filter result');
subplot(1,3,2);imagesc(oI);title('spatial filter result');
subplot(1,3,3);imagesc(fI-oI);title('differences');

figure;
subplot(1,3,1);mesh(abs(fftshift(fft2(fI))));title('frequency of frequency filter result');
subplot(1,3,2);mesh(abs(fftshift(fft2(oI))));title('frequency spatial filter result');
subplot(1,3,3);imagesc(abs(fftshift(fft2(fI)))-abs(fftshift(fft2(oI))));title('differences of remaining freq');


% impulse response
[H,f1,f2] = freqz2(sDoG);
figure;
subplot(1,3,1);mesh(imresize(fDoG,[size(H,1) size(H,2)]));
title('frequency filter''s impulse response');
subplot(1,3,2);mesh(H);title('spatial filter''s impulse response');
subplot(1,3,3);mesh(H-imresize(fDoG,[size(H,1) size(H,2)]));
title('impulse responses differences');




