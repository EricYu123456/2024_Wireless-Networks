% This is a script, not a function
% Task 8: Plot Figures
% Constellation Diagrams
figure;
subplot(2,2,1);
plot(bpsk_received(1:500), 'x'); 
title('BPSK, 40m');
axis([-2 2 -2 2]);
grid on;

subplot(2,2,2);  
plot(qpsk_received(1:500), 'x');
title('QPSK, 20m');
axis([-2 2 -2 2]);
grid on;

subplot(2,2,3);
plot(qam16_received(1:500), 'x');  
title('16QAM, 10m');
axis([-4 4 -4 4]);
grid on;

subplot(2,2,4);
plot(qam64_received(1:500), 'x');
title('64QAM, 5m');  
axis([-8 8 -8 8]);
grid on;

% Throughput vs. Distance
distances = 50:50:600;
bpsk_throughputs = zeros(size(distances));
qpsk_throughputs = zeros(size(distances));
qam16_throughputs = zeros(size(distances));
qam64_throughputs = zeros(size(distances));

for i = 1:length(distances)
    d = distances(i);
    
    % Calculations from Task 1 to Task 7
    % ... (omitted for brevity)
    script LAB1_111550016.m

    bpsk_throughputs(i) = bpsk_throughput_empirical;
    qpsk_throughputs(i) = qpsk_throughput_empirical;
    qam16_throughputs(i) = qam16_throughput_empirical;
    qam64_throughputs(i) = qam64_throughput_empirical;
end

figure;
plot(distances, bpsk_throughputs, '-o', distances, qpsk_throughputs, '-s', ...
     distances, qam16_throughputs, '-d', distances, qam64_throughputs, '-v');
legend('BPSK', 'QPSK', '16QAM', '64QAM', 'Location', 'best');
xlabel('Distance (m)');
ylabel('Throughput (bit/s)');
title('Throughput vs. Distance');
grid on;