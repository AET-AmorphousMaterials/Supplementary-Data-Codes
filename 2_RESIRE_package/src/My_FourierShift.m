function img2 = My_FourierShift(img,dy,dx)

ny = size(img,1); nx = size(img,2);
[X Y] = meshgrid(-ceil((nx-1)/2):floor((nx-1)/2),-ceil((ny-1)/2):floor((ny-1)/2));
F = my_ifft(img);
Pfactor = exp(2*pi*i*(dx*X/nx + dy*Y/ny));
img2 = my_fft(F.*Pfactor);

end