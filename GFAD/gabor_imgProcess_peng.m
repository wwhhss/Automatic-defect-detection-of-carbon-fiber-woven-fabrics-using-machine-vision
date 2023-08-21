function Ig=gabor_imgProcess_peng(I,ksize,lambda,theta,phase,sigma,ratio)
% input
%   I: input gray image
%   ksize: kernel size
%   lambda: wavelength
%   theta: orientation
%   phase: pahse angle
%   sigma: variation
%   ratio: spatial aspect ratio
% output
%   g: gabor filter
 
[m,n] = size(I);
d = ksize/2;
 
% pad image
Ip = zeros(m+ksize, n+ksize);
Ip(d+1:d+m, d+1:d+n)=I;
 
g = gabor_func_peng(ksize,lambda,theta,phase,sigma,ratio);
g = real(g); % only use the real part of gabor filter
Ig = zeros(m,n);
% disp('get gabor');
for x = 1:m
    for y = 1:n
        Ig(x,y) = sum(sum(Ip(x:x+ksize-1,y:y+ksize-1).*g));
    end
end
 
Ig = uint8(Ig);
Ig = min(255, max(0, Ig));
end

function g=gabor_func_peng(ksize,lambda,theta,phase,sigma,ratio)
% input
%   ksize: kernel size
%   lambda: wavelength
%   theta: orientation
%   phase: pahse angle
%   sigma: variation
%   ratio: spatial aspect ratio
% output
%   g: gabor filter
 
d = ksize/2;
g = zeros(ksize, ksize);
for x = 1:ksize
    xd = x - d;  % distance from the center
    for y = 1:ksize
        yd = y - d;
        xn = xd*cos(theta) + yd*sin(theta);
        yn = -xd*sin(theta) + yd*cos(theta);
        g(x,y) = exp(-(xn^2+ratio^2*yn^2)/(2*sigma^2))*exp(1i*(2*pi*xn/lambda+phase));
    end
end
 
end