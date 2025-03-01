close all; clear; clc;
addpath ./tasks;

origin = [0, 0];
tx_location = origin;
P_tx_dBm = 10;          % Transmission power of Tx (dBm)
N0_dBm = -95;           % Assume noise power is -90 dBm

% q1
distances = 50:50:500;
antenna_nums = [8, 16];
avg_snr = zeros(length(distances), length(antenna_nums));
avg_inr = zeros(length(distances), length(antenna_nums));

for d_idx = 1:length(distances)
    d = distances(d_idx);
    for ant_idx = 1:length(antenna_nums)
        ant_num = antenna_nums(ant_idx);
        snr_arr = zeros(1,10);
        inr_arr = zeros(1,10);
        for topo_idx = 1:10
            numbers = 0:10:180;
            random_index1 = randi(length(numbers));
            random_index2 = randi(length(numbers));
            random_number1 = numbers(random_index1);
            random_number2 = numbers(random_index2);

            offset = -5 + 10*rand();
            rx1_x = d * cosd(random_number1 + offset);
            rx1_y = d * sind(random_number1 + offset);
            rx1_location = [rx1_x, rx1_y];
            
            offset = -5 + 10*rand();
            rx2_x = d * cosd(random_number2 + offset);
            rx2_y = d * sind(random_number2 + offset);
            rx2_location = [rx2_x, rx2_y];
            
            [snr, inr, ~] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, ant_num, 0, 19);
            snr_arr(topo_idx) = snr;
            inr_arr(topo_idx) = inr;
        end
        avg_snr(d_idx, ant_idx) = mean(snr_arr);
        avg_inr(d_idx, ant_idx) = mean(inr_arr);
    end
end
fprintf('Average SNR (dBm) for different distances and antenna numbers:\n');
disp(avg_snr);

fprintf('Average INR (dBm) for different distances and antenna numbers:\n');
disp(avg_inr);


% q2
d = 200;
ant_num = 16;
snr_arr = zeros(1,10);
inr_arr = zeros(1,10);
for topo_idx = 1:10
    numbers = 0:10:180;
    random_index1 = randi(length(numbers));
    random_index2 = randi(length(numbers));
    random_number1 = numbers(random_index1);
    random_number2 = numbers(random_index2);

    offset = -5 + 10*rand();
    rx1_x = d * cosd(random_number1 + offset);
    rx1_y = d * sind(random_number1 + offset); 
    rx1_location = [rx1_x, rx1_y];
    
    offset = -5 + 10*rand();
    rx2_x = d * cosd(random_number2 + offset);
    rx2_y = d * sind(random_number2 + offset);
    rx2_location = [rx2_x, rx2_y];
    
    [snr, inr, ~] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, ant_num, 0, 19);
    snr_arr(topo_idx) = snr;
    inr_arr(topo_idx) = inr;
end

figure;
subplot(1,2,1);
plot(1:10, snr_arr, '-o');
xlabel('Topology Index');
ylabel('SNR (dBm)');
title('SNR for 10 Topologies, d=200m, 16 Antennas');

subplot(1,2,2); 
plot(1:10, inr_arr, '-o');
xlabel('Topology Index');
ylabel('INR (dBm)');
title('INR for 10 Topologies, d=200m, 16 Antennas');

% q3
d = 200;
ant_num = 16;
codebook_sizes = [19, 37, 73];
prx_array = zeros(10, length(codebook_sizes));

for topo_idx = 1:10
    numbers = 0:10:180;
    random_index1 = randi(length(numbers));
    random_number1 = numbers(random_index1);

    offset = -5 + 10*rand();
    rx1_x = d * cosd(random_number1 + offset);
    rx1_y = d * sind(random_number1 + offset);
    rx1_location = [rx1_x, rx1_y];
    
    for cb_idx = 1:length(codebook_sizes)
        cb_size = codebook_sizes(cb_idx);
        [~, ~, prx] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, ant_num, 0, cb_size);
        prx_array(topo_idx, cb_idx) = prx;
    end
end

figure;
plot(1:10, prx_array(:,1), '-o', 1:10, prx_array(:,2), '-x', 1:10, prx_array(:,3), '-s');
legend('Codebook Size 19', 'Codebook Size 37', 'Codebook Size 73');
xlabel('Topology Index');
ylabel('Prx,1 (dBm)');
title('Prx,1 for Different Codebook Sizes, d=200m, 16 Antennas');