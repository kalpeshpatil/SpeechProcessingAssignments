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


F0_list = [120 120 180]
F1_list = [300 1200 300]
B1_list = [100 200 100]

for iter = 1:length(F0_list)
    //specifying parameters
    F1 = F1_list(iter)
    B1 = B1_list(iter)
    Fs = 16000
    F0 = F0_list(iter)
    t_duration = 0.5
    //numerator of filter assumed to be 1
    k = 1
    //specifying input signal
    n_samples = t_duration*Fs + 1
    n_fft = n_samples*10
    
    x = zeros(n_samples,1)
    
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
    
    //computing poles 
    r = exp(-B1*%pi/Fs);
    theta = 2*%pi*F1/Fs
    
    //specifying coefficients of discrete time filter
    num = k
    den = [1 -2*r*cos(theta) r^2]

    //output from filter
    y = time_response(x,num,den,n_samples)
    
    //plotting time domain for response
    t_obs = 25 //we will observe signal for 25ms
    n_time_samples = floor(t_obs*Fs/1000)
    time_array = linspace(0,t_obs,n_time_samples)
    temp = find(y~=0)
    init = temp(1) //capture n_time_samples here onwards
    fig = scf()
    plot(time_array, y(init:init+n_time_samples-1))
    plot_title = strcat(['Output signal for F0 = ',string(F0),'Hz,F1 = ',string(F1),'Hz,B1 = ',string(B1),'Hz'])
    xtitle(plot_title,'time (ms)','y(n)')
    xs2jpg(gcf(), strcat(['../plots/Q3/',plot_title,'.jpg']));
    
    //converting to y to sound by limiting amplitude in [-1,1]
    //as required by the wavewrite 
    y_snd = y'/max(y')
    playsnd(y_snd,Fs);
    
    wavfile = strcat(['../wav_files/Q3/',plot_title,'.wav']);
    wavwrite(y_snd, Fs, wavfile);

end
