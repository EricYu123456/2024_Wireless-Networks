function [rx1_SNR_dbm, rx2_INR_dbm, Prx1_dBm] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_antenna_number, detail, cb_size)
addpath ./ewa_function;

% Generate codebook for analog beamforming
tx_beam_direction = 0:(180/(cb_size - 1)):180; % degree
d = 0.5;                      % Distance between antennas (multiple number of wavelength)
% Generate gain table: calculate power gain of beam at each direction(beam direction * sector angle)
% gain = gain(beam, sector angle)
resolution = 360;                                           % Number of angles dividing 180 degrees
a = zeros(numel(tx_beam_direction), tx_antenna_number);     % Antennas' array coefficient(weights), direction * number of antenna
A = zeros(numel(tx_beam_direction), resolution);            % Antennas' array factor
gain = zeros(numel(tx_beam_direction), resolution);         % Antennas' power gain
phi = (1 : resolution) * pi / resolution;                   % Equally-spaced over [0, pi] into resolution angle
phi_deg = phi * 180 / pi;                                   % Convert phi to degree
psi = 2 * pi * d * cos(phi);                                % Antenna phase shift of angle phi

for tx_beam_idx = 1:numel(tx_beam_direction)
    ph0 = tx_beam_direction(tx_beam_idx);
    a(tx_beam_idx,:) = uniform(d, ph0, tx_antenna_number);  % beam steering weights
    A(tx_beam_idx,:) = dtft(a(tx_beam_idx,:), -psi);        % Array factor, note dtft(a,psi)=dtft(a,-psi)
    gain(tx_beam_idx,:) = abs(A(tx_beam_idx,:)).^2;         % Power gain
end

% TODO1: Calculate the closest_sector_index and beam_index using gain table
% Hint: rx1_sector_index: Find the closest sector angle based on rx1_theta and rx2_theta
% Hint: rx1_beam_index: Find the optimal beam with the highest power gain for rx1
rx1_distance = norm(rx1_location - tx_location); % Distance between receiver 1 and transmitter (m)
rx2_distance = norm(rx2_location - tx_location); % Distance between receiver 2 and transmitter (m)
rx1_theta = atan2(rx1_location(2), rx1_location(1)); % Direction of receiver 1 (rad)
rx2_theta = atan2(rx2_location(2), rx2_location(1)); % Direction of receiver 2 (rad)
rx1_sector_index = 0;
rx1_beam_index = 0;

% 將rx1_theta和rx2_theta從弧度轉換為度
rx1_theta_deg = rx1_theta * 180 / pi;
rx2_theta_deg = rx2_theta * 180 / pi;

% 找到rx1和rx2的最近扇區角度索引
[~, rx1_sector_index] = min(abs(phi_deg - rx1_theta_deg));
[~, rx2_sector_index] = min(abs(phi_deg - rx2_theta_deg));

% 找到接收端1的最佳波束索引（最大增益）
[~, rx1_beam_index] = max(gain(:, rx1_sector_index));

% Uncomment to plot the beam pattern and show the results of optimal beam
if detail ~= 0
    figure(1); 
    polarplot(phi, gain(rx1_beam_index,:));
    title('Power Gain in different directions');
    
    figure(2); 
    plot(phi_deg, gain(rx1_beam_index, :)); 
    xlabel('Angle (degree)');
    xlim([0 180]);
    ylabel('Power Gain');
    title('Power Gain vs. Angle');
    
    fprintf('Beam Scanning Results\n');
    fprintf('\tReceiver 1 Distance: %f m\tAngle of departure: %2f degree\n', rx1_distance, phi_deg(rx1_sector_index));
    fprintf('\tReceiver 2 Distance: %f m\tAngle of departure: %2f degree\n', rx2_distance, phi_deg(rx2_sector_index));
    fprintf('\tClosest beam direction is %.2f (index: %d)\n', tx_beam_direction(rx1_beam_index), rx1_beam_index);
    fprintf('\t%d-th beam is the best beam for rx1 (with power gain %f)\n', rx1_beam_index, gain(rx1_beam_index, rx1_sector_index));
    fprintf('\t%d-th beam side lobe for rx2 (with power gain %f)\n', rx1_beam_index, gain(rx1_beam_index, rx2_sector_index));
end

% TODO2: Calculate receiving power and SNR (Prx-noise)
% Hint: Use friis equation with gain
freq = 24e9;
tx_gain_1 = gain(rx1_beam_index, rx1_sector_index); % 接收端1的增益
tx_gain_2 = gain(rx1_beam_index, rx2_sector_index); % 接收端2的增益 (主波束旁瓣)

% 计算接收功率
Prx1_dBm = P_tx_dBm + friis_equation(freq, tx_gain_1, 1, rx1_distance);
Prx2_dBm = P_tx_dBm + friis_equation(freq, tx_gain_2, 1, rx2_distance);

% 计算SNR和INR
rx1_SNR_dbm = Prx1_dBm - N0_dBm; % 接收端1的信噪比 (dB)
rx2_INR_dbm = Prx2_dBm - N0_dBm; % 接收端2的干扰噪声比 (dB)

% 输出结果
% fprintf('Receiver 1 SNR: %.2f dB\n', rx1_SNR_dbm);
% fprintf('Receiver 2 INR: %.2f dB\n', rx2_INR_dbm);

