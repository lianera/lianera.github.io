I = im2double(imread('input.jpg'));
[H, W, ~] = size(I);
[u,v]=meshgrid(1:W,1:H);
sigma = 10;
h = exp(-2*pi^2*sigma^2*(u.^2 + v.^2)/(4*H*W));
Iout = zeros(H,W,3);
for i = 1:3
    s = I(:,:,i);
    C = dct2(s);
    Cf = C.*h;
    Iout(:,:,i) = idct2(Cf);
end

DctOut = Iout;