function [y] = filterMyData(y,fs,orderFilt,fc_lowpass,fc_highPass)
% [y] = filterMyData(y,fs,orderFilt,fcUp,fc_highPass) extract the background
% component of the bridge vertical acceleration response. Since it requires
% low vibrational frequencies, a Zero-Pole-Gain design is preferred to the 
% Transfer Function design (see
% https://se.mathworks.com/help/signal/ref/butter.html).
% 
% Input
% y: [1xN] or [Nx1] time series of the bridge displacement response
% fs: [1x1] double: sampling frequency
% orderFilt [1x1] integer: order of the high-pas and low-pass filters 
% fc_lowpass [1x1] double: cut-off frequency of the low-pass filter
% fc_highPass [1x1] double: cut-off frequency of the high-pass filter
% 
% Output
% y: [1xN] filtered bridge displacement response
% 
% Author: E. Cheynet - UiB -03-08-2020


y = y(:)';

WnU=fc_lowpass/(fs/2); % Normalized cutoff frequency for the low-pass filter


% remove the frequencies above fc_lowpass
[z,p1,k] = butter(orderFilt,WnU,'low'); % low pass filter design
[sos, g0] = zp2sos(z, p1, k); % Zero-pole-gain to second-order sections model conversion.
y  =filtfilt(sos, g0, y);  % Zero-phase forward and reverse digital IIR filtering.

% remove the frequencies below fc_highPass
WnB=fc_highPass/(fs/2); % Normalized cutoff frequency for the high-pass filter
[z,p1,k] = butter(orderFilt,WnB,'high');
[sos, g0] = zp2sos(z, p1, k);
y  =filtfilt(sos, g0, y);


end