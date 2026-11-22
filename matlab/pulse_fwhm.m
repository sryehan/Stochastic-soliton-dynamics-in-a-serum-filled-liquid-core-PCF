function w = pulse_fwhm(x, I)
%PULSE_FWHM  Full width at half maximum of a single-peaked profile I(x), by linear interpolation.
[Imax, j] = max(I); h = Imax/2;
a = find(I(1:j) < h, 1, 'last');  b = j - 1 + find(I(j:end) < h, 1, 'first');
xl = x(a) + (h - I(a))*(x(a+1) - x(a))/(I(a+1) - I(a));
xr = x(b-1) + (h - I(b-1))*(x(b) - x(b-1))/(I(b) - I(b-1));
w = xr - xl;
end
