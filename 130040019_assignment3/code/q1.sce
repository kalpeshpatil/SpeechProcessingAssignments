
clear all


//Function to calculate filter response using difference equation
function [y]= time_response (x,num,den)
    n_samples = length(x);  
    y = zeros(n_samples,1)
    //numerator is constant (all pole filter)
    y(1) = num(1)*x(1)
    k = num(1)
    //response by taking coefficients for denominator
    
    for ii =2:n_samples
        temp = k*x(ii)
        for jj = 1:min(ii-1,length(den)-1)
            temp = temp - den(jj+1)*y(ii - jj)
        end
        y(ii) = temp
    end
endfunction

//function to read previously generated wavefiles
function y = read_file(wavefile)
    y = loadwave(strcat(['../wav_files/input/',wavefile]));
endfunction

//pre_emphasis filter
function y = pre_emphasis_filter(x,alpha)
    y = zeros(length(x))
    y(1) = x(1)
    for i = 2:length(x)
        y(i) = x(i)-alpha*x(i-1)
    end
endfunction

//de_emphasis filter
function y = de_emphasis_filter(x,alpha)
    y = zeros(length(x))
    y(1) = alpha*x(1)
    for i = 2:length(x)
        y(i) = alpha*y(i-1) + x(i)
    end
endfunction


//find windowed fft
function [X_mag] = find_windowed_FFT (x,window_type,N,n_fft)
    //we can put this window anywhere on the signal. We will put it
    //at the centre
    x_init = (length(x)+1)/2 - (N-1)/2  
    select window_type
    case 'hamming' then
        win_hamming = window('hm', N);
        windowed_x = x(x_init: x_init + N - 1).*win_hamming;
    case 'rect' then
        windowed_x = x(x_init: x_init + N - 1)       
    end
    //zero padding
    pad = zeros(1,ceil((n_fft - N)/2));
    padded_x = [pad windowed_x pad]    
    //finding fft
    freq_x = fftshift(fft(padded_x))
    X_mag = abs(freq_x(n_fft/2 +1:n_fft))    
endfunction


// function to find output of LP aproximation
function y = find_output_of_LP (num, den, input_type,Fs,F0,t_duration)
        n_samples = t_duration*Fs + 1
        x = zeros(n_samples,1)
        //impulse train generation
        if(input_type == "voiced")
            for i = 2:t_duration*F0
                x(floor((i-1)*Fs/F0))= 1  
            end
        else
            x = rand(x,'normal')
        end
        y = zeros(n_samples,1)
        y = time_response(x,num,den)
        
        if(input_type == 'voiced')
            y = de_emphasis_filter(y,0.95)
        end
        
endfunction


files_list = ["a_male_8k","i_male_8k","n_male_8k","s_male_16k"]
Fs_list = [8000,8000,8000,16000]


do_emphasis_list = [1,1,1,0]
input_type = ["voiced","voiced","voiced","unvoiced"]
selected_order_list = [10,10,10,20]

for iii = 1:length(Fs_list)
        Fs = Fs_list(iii)
        x = read_file(files_list(iii))'
        if(do_emphasis_list(iii)==1)
            x = pre_emphasis_filter(x,0.95)
        end
        x = x'
        window_time = 30
        N = window_time*Fs/1000 + 1
        n_fft = 10*N
        x_init = (length(x)+1)/2 - (N-1)/2  
        win_hamming = window('hm', N);
        y = x(x_init: x_init + N - 1).*win_hamming;

        r = zeros(length(y),1)  
        for ki=0:length(x) 
             for ni=ki:N-1
                r(ki+1)=r(ki+1)+y(ni+1)*y(ni-ki+1);
             end
        end
        
        delta_n = zeros(N,1)
        delta_n(1) = 1

        h_total = []
        //LP recursion
        p_max = selected_order_list(iii)
        e_vec = [r(1)]
        i = 1
        k = r(2)/e_vec(1)
        a_vec = [1,k]
        e_vec = [e_vec,(1-(k)^2)*e_vec(1)]
        
        for i = 2:p_max
            k = r(i+1)
            for j = 1:i-1
                k = k - a_vec(j+1)*r(i-j+1)
            end
            k = k/e_vec(i)
            a_vec_old = a_vec
            a_new = k
            for j = 1:i-1
                a_vec(j+1) = a_vec_old(j+1)-k*a_vec_old(i-j+1)
            end
            a_vec = [a_vec , a_new]
            e_vec(i+1) = (1-k^2)*e_vec(i)
                   
        end
        

        
               fig = scf()
               h_padded = zeros(n_fft,1)               
               h_padded(1:length(y)) = y'
               H = fftshift(fft(h_padded))    
               H_mag = abs(H(n_fft/2 +1:n_fft))
               
               freq_array = linspace(0,Fs/2,n_fft/2)
               plot(freq_array,20*log(H_mag))
               num = sqrt(e_vec(i+1))
               den = -a_vec(2:length(a_vec))
               den = [1,den]
               [h,w] = frmag(num,den,n_fft)
               plot(w*Fs,20*log(abs(h)),'red')
               plot_title = strcat(['LP approximation of plots ',files_list(iii),' for p = ',string(selected_order_list(iii))])
               xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
               xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']));  
               
               
               
                t_duration = 0.3
                //calculated from Assignment 2
                F0 = 125
                y_LP = find_output_of_LP (num, den, input_type(iii),Fs,F0,t_duration) 
                //converting to y to sound by limiting amplitude in [-1,1]
                //as required by the wavewrite 
                
                t_obs = 25 //we will observe signal for 25ms
                n_time_samples = floor(t_obs*Fs/1000)
                time_array = linspace(0,t_obs,n_time_samples)
                temp = find(y_LP~=0)
                init = temp(1) //capture n_time_samples here onwards
                fig = scf()
                plot(time_array, y_LP(init:init+n_time_samples-1))
                plot_title = strcat(['LP approximation of signal ',files_list(iii),' for p = ',string(selected_order_list(iii))])
                xtitle(plot_title,'time (ms)','y(n)')
                xs2jpg(gcf(), strcat(['../plots/Q1/',plot_title,'.jpg']));
                
                y_snd = y_LP'/max(y_LP')
                playsnd(y_snd,Fs);    
                wavfile = strcat(['../wav_files/Q1/',plot_title,'.wav'])
                wavwrite(y_snd, Fs, wavfile);           

end


xdel(winsid())








