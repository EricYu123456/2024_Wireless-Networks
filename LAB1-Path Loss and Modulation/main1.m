d = 600;
sub1(d);


% 提取实部和虚部
real_part = real(bpsk_received ./ h);
imag_part = imag(bpsk_received ./ h);

% 绘制在复数平面上
plot(real_part, imag_part, '.');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Complex Plane Plot');
axis equal;
grid on;

% BPSK with 600m
