function kout = my_fft(img)
kout = fftshift(fftn((ifftshift(img))));
end