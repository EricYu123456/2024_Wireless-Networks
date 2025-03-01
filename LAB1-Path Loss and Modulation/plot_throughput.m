% Throughput vs. Distance
distances = 50:50:600;
bpsk_throughputs = zeros(size(distances));
qpsk_throughputs = zeros(size(distances));
qam16_throughputs = zeros(size(distances));
qam64_throughputs = zeros(size(distances));

for i = 1:length(distances)
    d = distances(i);
    
    [bpsk_throughput_empirical, qpsk_throughput_empirical, qam16_throughput_empirical, qam64_throughput_empirical] = sub1(d);
    
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
