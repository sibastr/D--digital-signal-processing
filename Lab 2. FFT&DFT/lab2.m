% Main function
function main()

% Input parameters
T = 2.0;
sigma = 0.5;

% Borders of calculation
delta = 5;
t = -delta:0.05:delta;

% Calculation of pulse functions
x1 = zeros(size(t));
x1(abs(t) - T < 0) = 1;
x1(abs(t) == T) = 0.5;
x2 = exp(-(t/sigma).^2);

% FFT of functions (with "twin" effect)
y_rec = fft(x1);
y_gauss = fft(x2);

% FFT of functions (without "twin" effect)
y_rec_without_twin = fftshift(y_rec);
y_gauss_without_twin = fftshift(y_gauss);

% DFT of functions (with "twin" effect)
z_rec = dft(x1);
z_gauss = dft(x2);

% DFT of functions (without "twin" effect)
z_rec_without_twin = fftshift(z_rec);
z_gauss_without_twin = fftshift(z_gauss);

%
% PLOTTING
%

xs = 0:length(t)-1;

figure (1);
subplot(2,1,1);
plot(xs, y_rec,xs,abs(y_rec_without_twin)/length(xs));
title('FFT: amplitude spectrum of rectangle pulse');
legend('With "twin" effect','Without "twin" effect');
subplot(2,1,2);
plot(xs,angle(y_rec)/length(xs),xs,angle(y_rec_without_twin)/length(xs));
title('FFT: phase spectrum of rectangle pulse');
legend('With "twin" effect','Without "twin" effect');
print -dpng plot2_1.png;

figure (2);
subplot(2,1,1);
plot(xs,abs(y_gauss)/length(xs),xs,abs(y_gauss_without_twin)/length(xs));
title('FFT: amplitude spectrum of gaussian pulse');
legend('With "twin" effect','Without "twin" effect');
subplot(2,1,2);
plot(xs,angle(y_gauss)/length(xs),xs,angle(y_gauss_without_twin)/length(xs));
title('FFT: phase spectrum of gaussian pulse');
legend('With "twin" effect','Without "twin" effect');
print -dpng plot2_2.png;

figure (3);
subplot(2,1,1);
plot(xs,abs(z_rec)/length(xs),xs,abs(z_rec_without_twin)/length(xs));
title('DFT: amplitude spectrum of rectangle pulse');
legend('With "twin" effect','Without "twin" effect');
subplot(2,1,2);
plot(xs,angle(z_rec)/length(xs),xs,angle(z_rec_without_twin)/length(xs));
title('DFT: phase spectrum of rectangle pulse');
legend('With "twin" effect','Without "twin" effect');
print -dpng plot2_3.png;

figure (4);
subplot(2,1,1);
plot(xs,abs(z_gauss)/length(xs),xs,abs(zg2)/length(xs));
title('DFT: amplitude spectrum of gaussian pulse');
legend('With "twin" effect','Without "twin" effect');
subplot(2,1,2);
plot(xs,angle(z_gauss)/length(xs),xs,angle(zg2)/length(xs));
title('DFT: phase spectrum of gaussian pulse');
legend('With "twin" effect','Without "twin" effect');
print -dpng plot2_4.png;

end

% Discrete Fourier Transform function
function y = dft(x)
a = 0:length(x)-1;
b = -2 * pi * sqrt(-1) * a / length(x);
for i = 1:length(a)
a(i) = 0;
for j = 1:length(x)
a(i) = a(i) + x(j) * exp(b(i) * j);
end
end
y = a;
end