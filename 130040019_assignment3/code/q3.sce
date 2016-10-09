

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
        
endfunction


function [num, den] = find_LP_coefs_using_levinson (y,p_max)
        r = zeros(length(y),1)  
        for ki=0:length(x) 
             for ni=ki:N-1
                r(ki+1)=r(ki+1)+y(ni+1)*y(ni-ki+1);
             end
        end
        


        //LP recursion

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
        
        num = sqrt(e_vec(i+1))
        den = -a_vec(2:length(a_vec))
        den = [1,den]


endfunction


function ceptral_coefs = find_ceptral_coefs(x,n_fft)
    N = length(x)
    //zero padding
    x_padded = zeros(n_fft,1)    
    x_padded(1:length(x)) = x'
    //finding fft
    freq_x = fftshift(fft(x_padded))
    //X_mag = abs(freq_x(n_fft/2 +1:n_fft))    
    X_mag = (abs(freq_x))
    log_mag = log(X_mag)
    log_mag_padded = zeros(n_fft,1)
    log_mag_padded(1:length(log_mag)) = log_mag
   // ceptral_coefs = ifft(log_mag_padded)
   ceptral_coefs = ifft(log(abs(fft(x_padded))))    
    
endfunction
// Compute and plot the narrowband spectrum using a Hamming window 
// of duration = 30ms before and after pre-emphasis.

files_list = ["a_male_8k","i_male_8k","n_male_8k","s_male_16k"]
Fs_list = [8000,8000,8000,16000]

// Using a 30 ms Hamming window centered in the segment of the 
//  waveform preemphasised for the voiced sounds)

do_emphasis_list = [1,1,1,0]
p_list = [4,6,8,10,12,20]
input_type = ["voiced","voiced","voiced","unvoiced"]
selected_order_list = [10,10,10,18]
ceptral_lengths = [15,15,15,15]
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
        x = x(x_init: x_init + N - 1).*win_hamming;
        
        p_max = selected_order_list(iii)
        [num, den] = find_LP_coefs_using_levinson (x,p_max)
        
        n_fft = 10*length(x)
        ceptral_length = ceptral_lengths(iii)
        ceptral_coefs = find_ceptral_coefs(x,n_fft)
        padded_ceptral_coefs = zeros(n_fft,1)
        padded_ceptral_coefs(1:ceptral_length) = ceptral_coefs(1:ceptral_length)
        padded_ceptral_coefs(n_fft - ceptral_length+1:n_fft) = ceptral_coefs(n_fft - ceptral_length+1:n_fft)        
                // LP approximation
                t_duration = 0.3
                F0 = 125
                y = find_output_of_LP (num, den, input_type(iii),Fs,F0,t_duration) 
             
               fig = scf()
               plot([1:99],ceptral_coefs(2:100))
               plot_title = strcat(['Cepstral coefficients (zoomed) ',files_list(iii),' for p = ',string(selected_order_list(iii))])
               xtitle(plot_title,'n','C[n]')
               xs2jpg(gcf(), strcat(['../plots/Q3/',plot_title,'.jpg']));  
        
               
               fig = scf()
               h_padded = zeros(n_fft,1)               
               h_padded(1:length(x)) = x'
               H = fftshift(fft(h_padded))    
               H_mag = abs(H(n_fft/2 +1:n_fft))          
               //H_mag = abs(H(1:n_fft/2))          
               freq_array = linspace(0,Fs/2,n_fft/2)
               plot(freq_array,20*log(H_mag))
               [h,w] = frmag(num,den,n_fft)
               plot(w*Fs,20*log(abs(h)),'red')
               
               H = fftshift(fft(padded_ceptral_coefs))    
               H_mag = abs(H(n_fft/2 +1:n_fft))          
               freq_array = linspace(0,Fs/2,n_fft/2)
               plot(freq_array,20*log(abs(exp(H(n_fft/2 +1:n_fft)))),'k')
               
               plot_title = strcat(['Cepstral estimation and LP approximation for natural sound',files_list(iii),' for p = ',string(selected_order_list(iii))])
               xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')               
               legend (['Signal Spectrum','LP approximation','Cepstral Estimate'],3)
               xs2jpg(gcf(), strcat(['../plots/Q3/',plot_title,'.jpg']));  
        
end
//xdel(winsid())








