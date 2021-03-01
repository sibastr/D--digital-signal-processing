% Source signal param
sigma = 5;
tt = 5;

% Discrete signal
n = input('Input number of samples: ');
dt = input('Input sample step: ');
t_max = dt*(n-1)/2;
t = -t_max:dt:t_max; 

gauss_discrete = exp(-(t/sigma).^2);
rect_discrete = zeros(size(t));
rect_discrete(abs(t) - tt < 0) = 1;

% Reference signal
x = -t_max:0.005:t_max;
gauss_ref = exp(-(x/sigma).^2);
rect_ref = zeros(size(x));
rect_ref(abs(x) - tt < 0) = 1;



% Signal restore (
gauss_restored = zeros(1, length(x));
rect_restored = zeros(1, length(x));
for i=1:length(x)
   for j = 1:n
       gauss_restored(i) = gauss_restored(i) + gauss_discrete(j)*sinc((x(i)-t(j))/dt);
       rect_restored(i) = rect_restored(i) + rect_discrete(j)*sinc((x(i)-t(j))/dt);
   end
end

figure;
subplot(2,1,1);
title('Gaussian filter');
hold on;
grid on;
plot(x, gauss_ref, 'b');
plot(x, gauss_restored, 'k');
plot(t, gauss_discrete, '.m');
legend('Reference', 'Restored', 'Discrete');

subplot(2,1,2);
title('Rectangular function');
hold on;
grid on;
plot(x, rect_ref, 'b');
plot(x, rect_restored, 'k');
plot(t, rect_discrete, '.m');
legend('Reference', 'Restored', 'Discrete');

print -dpng plot1.png;