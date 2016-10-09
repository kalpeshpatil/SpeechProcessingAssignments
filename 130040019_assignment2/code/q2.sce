clear all

//function to read previously generated wavefiles
function y = read_file(wavefile)
    y = loadwave(strcat(['../wav_files/Q2/',wavefile]));
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
        y(i) = y(i-1) + alpha*x(i)
    end
endfunction


//inverse filter
function y = inverse_filter(x,num,den)
    y = zeros(length(x),1)
    for i = 1:length(x)
        for j = 1:min(length(den),i)
            y(i) = y(i) + den(j)*x(i-j+1)
        end
        y(i) = y(i)/num
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



//1. Compute and plot the narrowband spectrum using a Hamming window 
// of duration = 30ms before and after pre-emphasis.

files_list = ["a_male_8k","i_male_8k","n_male_8k","s_male_16k"]
Fs_list = [8000,8000,8000,16000]

for i = 1:length(Fs_list)
    x = read_file(files_list(i))
    window_time = 30
    Fs = Fs_list(i)
    N = window_time*Fs/1000 + 1
    n_fft = 10*N
    
   [X_mag] = find_windowed_FFT (x,'hamming',N,n_fft)
   freq_array = linspace(0,Fs/2,n_fft/2)
   fig = scf()
   x_preemph = pre_emphasis_filter(x,0.95)
   [X_preemph_mag] = find_windowed_FFT (x_preemph','hamming',N,n_fft)
   plot2d(freq_array,[20*log(X_preemph_mag)',20*log(X_mag)'],[2,5])
   plot_title = strcat(['Spectrum before and after preemphasis of ',files_list(i)])
   xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
   legends(['after preemphasis','before preemphasis'],[2,5],opt = 6)
   xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']));    
   
end 

//2. Using a 30 ms Hamming window centered in the segment of the 
//  waveform preemphasised for the voiced sounds)

do_emphasis_list = [1,1,1,0]
p_list = [4,6,8,10,12,20]
pole_zero_list = [6,10]
do_inverse_filtering = [1,0,0,1]
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
             for ni=ki:L-1
                r(ki+1)=r(ki+1)+y(ni+1)*y(ni-ki+1);
             end
        end
        
        delta_n = zeros(N,1)
        delta_n(1) = 1

        h_total = []
        //LP recursion
        p_max = 20
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
            
            //overlapping LP approximation plot on windowed fft       
            if (find(i==p_list) ~= [])
               fig = scf()
               h_padded = zeros(n_fft,1)               
               h_padded(1:length(y)) = y'
               H = fftshift(fft(h_padded))    
               H_mag = abs(H(n_fft/2 +1:n_fft))
               if(h_total == [])
                   h_total = [h_total,H_mag]
               end
               freq_array = linspace(0,Fs/2,n_fft/2)
               plot(freq_array,20*log(H_mag))
               num = sqrt(e_vec(i+1))
               den = -a_vec(2:length(a_vec))
               den = [1,den]
               [h,w] = frmag(num,den,n_fft)
               [h2,w2] = frmag(num,den,n_fft/2)
               h_total = [h_total,h2']
               plot(w*Fs,20*log(abs(h)),'red')
               plot_title = strcat(['LP approximation of ',files_list(iii),' for p = ',string(i)])
               xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
               xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']));    
               
            end
            
            //plotting pole zero plots
            if(find(i == pole_zero_list)~=[])
                fig = scf()
                den = -a_vec(2:length(a_vec))
                den = [1,den]
                poles = roots(den)
                plot2d([real(poles),zeros(length(poles),1)],[imag(poles),zeros(length(poles),1)],[-2,-3])
                legends(['poles','zeros'],[-2,-3],opt = 'lr')
                plot_title = strcat(['pole zero plot of ',files_list(iii),' for p = ',string(i)])
                xtitle(plot_title,'real','imaginary')
                ax=get("current_axes")
                ax.data_bounds = [-1,-1;1,1]
                xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg'])); 
            end
        
            //inverse filtering 
            if((do_inverse_filtering(iii)) & (i == 10))
                den = -a_vec(2:length(a_vec))
                den = [1,den]
                num = sqrt(e_vec(i+1))
                residual_signal = inverse_filter(x,num,den)
                fig = scf()
                plot(residual_signal)
                plot_title = strcat(['Residual signal ',files_list(iii),' for p = ',string(i)])
                xtitle(plot_title,'n','s[n]')
                xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']))
                fig = scf()
                h_padded = zeros(n_fft,1)               
                h_padded(1:length(residual_signal)) = residual_signal
                H = fftshift(fft(h_padded))    
                H_mag = abs(H(n_fft/2 +1:n_fft))
                freq_array = linspace(0,Fs/2,n_fft/2)
                plot(freq_array,20*log(H_mag))
                plot_title = strcat(['Residual signal spectrum ',files_list(iii),' for p = ',string(i)])
                xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
                xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg'])); 
   
                //3. finding pitch period for voiced signals using ACF of residual signal 
                // and plotting magnitude spectrum of residual signal
                if(do_emphasis_list(iii)==1)
                    temp_res = length(residual_signal)*xcorr(residual_signal,["biased"])
                    acf_residual = temp_res(length(residual_signal):length(temp_res))  
                    fig = scf()
                    plot(acf_residual)
                    plot_title = strcat(['Residual signal acf of ',files_list(iii),' for p = ',string(i)])
                    xtitle(plot_title,'lag(k)','ACF')
                    xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg'])); 
                    //63.5 is period
                    acf_zoomed = acf_residual(1:200)
                    fig = scf()
                    plot(acf_zoomed)
                    plot_title = strcat(['Residual signal acf zoommed of ',files_list(iii),' for p = ',string(i)])
                    xtitle(plot_title,'lag(k)','ACF')
                    xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg'])); 
                    pitch_period = Fs_list(iii)*63.5     
                end       
            end           
        end
        
        //plotting error vs p graph
        fig = scf()
        plot([0:p_max],e_vec)
        plot_title = strcat(['Error signal energy of ',files_list(iii)])
        xtitle(plot_title,'p','Error signal energy')
        xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']));       
        
        //plotting all LP together
        fig = scf()
        plot2d(w2,20*log(h_total),[1:(length(p_list)+1)])
        plot_title = strcat(['All lp plots together for ',files_list(iii)])
        xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
        legends(['original','p = 4','p = 6','p = 8','p = 10','p = 12','p = 20'],[1:(length(p_list)+1)],opt = 3)
        xs2jpg(gcf(), strcat(['../plots/Q2/',plot_title,'.jpg']));    

end
xdel(winsid())








