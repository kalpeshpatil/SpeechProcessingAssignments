clear all

//Function to calculate filter response using difference equation
function [y]= time_response (x,num,den,n_samples)
    y = zeros(n_samples,1)
    //numerator is constant (all pole filter)
    y(1) = num(1)*x(1)
    //response by taking coefficients for denominator
    for ii =2:n_samples
        temp = k*x(ii)
        for jj = 1:min(ii-1,length(den)-1)
            temp = temp - den(jj+1)*y(ii - jj)
        end
        y(ii) = temp
    end
endfunction

//specifying parameters
F1 = 1000
B1 = 200
Fs = 16000
//numerator of filter assumed to be 1
k = 1

//number of samples
n_samples = 200
//number of points for fft
n_fft = n_samples*10

//computing poles 
r = exp(-B1*%pi/Fs);
theta = 2*%pi*F1/Fs

//specifying coefficients of discrete time filter
num = k
den = [1 -2*r*cos(theta) r^2]

//impulse input 
delta_n = zeros(n_samples,1)
delta_n(1) = 1

//impulse response
h = time_response(delta_n,num,den,n_samples)
h_padded = zeros(n_fft,1)
h_padded(1:length(h)) = h
n_time_samples = 100
time_array = [0:n_time_samples-1]*1000/Fs
fig = scf()
plot(time_array, h(1:n_time_samples))
xtitle('impulse response','time (ms)','h(n)')
xs2jpg(gcf(), '../plots/Q1/impulse response.jpg');

//Frequency response of the filter
fig = scf()
H = fftshift(fft(h_padded))
H_mag = abs(H(n_fft/2 +1:n_fft))
freq_array = linspace(0,Fs/2,n_fft/2)
plot(freq_array,20*log(H_mag))
xtitle('Single formant resonator (Q.1)','Frequency (Hz)','Magnitude in dB')
xs2jpg(gcf(), '../plots/Q1/single_formant_resonator_magnitude_plot.jpg');
