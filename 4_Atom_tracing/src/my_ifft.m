function realout = my_ifft(k)
realout =fftshift((ifftn(ifftshift(k))));
%realout =ifftshift((ifftn(fftshift(k))));

end