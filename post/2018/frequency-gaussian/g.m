I = im2double(imread('input.jpg'));
[H, W, ~] = size(I);
N = 2*H;
M = 2*W;
Ip = zeros(N,M, 3);
Ip(1:H,1:W,:) = I;
[x,y]=meshgrid(1:M,1:N);
Ic = Ip.*repmat((-1).^(x+y),1,1,3);
F=fft2(Ic);

Mag=abs(F);
Pha = angle(F);
Vmag = log(Mag)/8;
Vpha=(Pha/pi+1)/2;


[u,v] = meshgrid(1:M,1:N);
cu = (1+M)/2;
cv = (1+N)/2;
sigma = 10;
h = exp(-2*pi^2*sigma^2*((u-cu).^2 + (v-cv).^2)/(M*N));

G = F.*repmat(h,1,1,3);
Is=real(ifft2(G));
t=Is.*repmat((-1).^(x+y),1,1,3);
Iout = zeros(H, W, 3);
Iout(1:H, 1:W, :) = t(1:H, 1:W, :);

groundtruth=imgaussfilt(I,sigma);

imwrite(Vmag, 'mag.jpg');
imwrite(Vpha, 'phase.jpg');
imwrite(Ip, 'Ip.jpg');

imwrite(h, 'hfilter.jpg');
imwrite(Is, 'Is.jpg');
imwrite(Iout, 'Iout.jpg');

imwrite(groundtruth, 'groundtruth.jpg');

dif = abs(groundtruth-Iout);
imshow(dif);

DftOut = Iout;