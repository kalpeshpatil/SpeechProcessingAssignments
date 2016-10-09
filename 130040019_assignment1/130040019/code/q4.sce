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

//function to find discrete coefficients of cascade of multiple filters
function [num,den] = find_cascade_filter(F_list,B_list,Fs)
    num = 1
    den = [1]
    for iter = 1:length(F_list)
        F = F_list(iter)
        B = B_list(iter)
        //finding poles for current iteration
        r = exp(-B*%pi/Fs)
        theta = 2*%pi*F/Fs
        //current filter
        num_curr = 1
        den_curr = [1 -2*r*cos(theta) r^2]
        //multiplication of polynomials can be computed using their 
        //convolution. Numerator is constant(1). Hence only 
        //denominator needs to be multiplied
        den = conv(den,den_curr)
    end
endfunction

//plots frequency response given impulse response 
function plot_frequency_response(h,n_fft,vowel)
    h_padded = zeros(n_fft,1)
    h_padded(1:length(h)) = h
    //Frequency response of the filter
    fig = scf()
    H = fftshift(fft(h_padded))
    H_mag = abs(H(n_fft/2 +1:n_fft))
    freq_array = linspace(0,Fs/2,n_fft/2)
    plot(freq_array,20*log(H_mag))
    plot_title = strcat(['Frequency response of vowel ',vowel])
    xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
    xs2jpg(gcf(), strcat(['../plots/Q4/',plot_title,'.jpg']));    
endfunction

//Plots output signal, plays and store the sound
function plot_y(y,F0,Fs,vowel)
    //plotting time domain for response
    t_obs = 25 //we will observe signal for 25ms
    n_time_samples = floor(t_obs*Fs/1000)
    time_array = linspace(0,t_obs,n_time_samples)
    temp = find(y~=0)
    init = temp(1) //capture n_time_samples here onwards
    fig = scf()
    plot(time_array, y(init:init+n_time_samples-1))
    plot_title = strcat(['Output signal of vowel ',vowel,' for F0 = ',string(F0)])
    xtitle(plot_title,'time (ms)','y(n)')
    xs2jpg(gcf(), strcat(['../plots/Q4/',plot_title,'.jpg']));
    
    //converting to y to sound by limiting amplitude in [-1,1]
    //as required by the wavewrite 
    y_snd = y'/max(y')
    playsnd(y_snd,Fs);    
    wavfile = strcat(['../wav_files/Q4/',plot_title,'.wav'])
    wavwrite(y_snd, Fs, wavfile);
endfunction

F0_list = [120,220]
vowel_frequency_matrix = [730, 1090, 2440; //F1,F2,F3 for vowell a
                          270, 2290, 3010; //F1,F2,F3 for vowell i
                          300, 870, 2240]  //F1,F2,F3 for vowell u

vowel_list = ['a','i','u']

Fs = 16000



for i_f0 = 1:length(F0_list)
    F0 = F0_list(i_f0)
    t_duration = 0.5
    //numerator of filter assumed to be 1
    k = 1

    //specifying input signal
    n_samples = t_duration*Fs + 1
    n_fft = n_samples*10;
    x = zeros(n_samples,1);
    //impulse is approximated by narrow triangular pulse 
    for i = 2:t_duration*F0
        x(floor((i-1)*Fs/F0))= 1
        x(floor((i-1)*Fs/F0 + 1))= 0.75
        x(floor((i-1)*Fs/F0 + 2))= 0.5
        x(floor((i-1)*Fs/F0 + 3))= 0.25
        x(floor((i-1)*Fs/F0 - 1))= 0.75
        x(floor((i-1)*Fs/F0 - 2))= 0.5    
        x(floor((i-1)*Fs/F0 - 3))= 0.25    
    end

    //Create delta input for finding impulse response
    delta_n = zeros(n_samples,1)
    delta_n(1) = 1
    
    for j = 1:length(vowel_frequency_matrix(1,:))
        F_list = vowel_frequency_matrix(j,:)
        //Bandwidth is constant for all formants
        B_list = [100,100,100]
        
        //find output of cascade filter
        [num,den] = find_cascade_filter(F_list,B_list,Fs)
        h = time_response(delta_n,num,den)
        y = time_response(x,num,den)
        
        //plot responses and save sound files
        plot_frequency_response(h,n_fft,vowel_list(j))
        plot_y(y,F0,Fs,vowel_list(j))        
   end
end



