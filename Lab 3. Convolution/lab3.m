function main

% Input parameters
T = 2.0;
sigma = 1.0;

% Borders of calculation
mult = 5;
step = 0.05;
t = -mult:step:mult;

% Pulse generation
x1 = [rectpls(t,T) zeros(1,length(t))];
x2 = [gauspls(t, sigma) zeros(1,length(t))];
x3 = [rectpls(t,T/2) zeros(1,length(t))];
x4 = [gauspls(t, sigma/2) zeros(1,length(t))];

% Convolution
% Фурье-образ взаимной свертки равен произведению Фурье-образов свертываемых функций.
y1 = ifft(fft(x1).*fft(x2))*step;
y2 = ifft(fft(x1).*fft(x3))*step;
y3 = ifft(fft(x2).*fft(x4))*step;

% Normalize convolution
start = fix((length(y1)-length(t))/2);
y1 = y1(start+1:start+length(t));
y2 = y2(start+1:start+length(t));
y3 = y3(start+1:start+length(t));

%
% PLOTTING
%

figure(1)
plot(t,x1(1:201),t,x2(1:201),t,y1);
title('Rectangular and Gaussian convolution');
legend('Pulse 1','Pulse 2','Convolution');
print -dpng plot3_1.png;

figure(2)
plot(t,x1(1:201),t,x3(1:201),t,y2);
title('Two Rectangular convolution');
legend('Pulse 1','Pulse 2','Convolution');
print -dpng plot3_2.png;

figure(3)
plot(t,x2(1:201),t,x4(1:201),t,y3);
title('Two Gaussian convolution');
legend('Pulse 1','Pulse 2','Convolution');
print -dpng plot3_3.png;

end

% Rectangular pulse generation
function y = rectpls(x,T)
y = zeros(size(x));
y(abs(x) - T < 0) = 1;
y(abs(x) == T) = 1/2;
end

% Gaussian pulse generation
function y = gauspls(x,s)
y = exp(-(x/s).^2);
end
