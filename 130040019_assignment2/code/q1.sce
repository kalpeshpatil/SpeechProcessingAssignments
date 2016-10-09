
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

//produce vowel /a/ for given F0
function a = produce_vowel(F0)
    vowel_formant_list = [730, 1090, 2440]; //F1,F2,F3 for vowell a
    Fs = 8000
    t_duration = 0.5
    //numerator of filter assumed to be 1
    k = 1
    //specifying input signal
    n_samples = t_duration*Fs + 1
    n_fft = n_samples*10;
    x = zeros(n_samples,1);
    //impulse train generation
    for i = 2:t_duration*F0
        x(floor((i-1)*Fs/F0))= 1  
    end
    //Create delta input for finding impulse response
    delta_n = zeros(n_samples,1)
    delta_n(1) = 1
    F_list = vowel_formant_list
    //Bandwidth is constant for all formants
    B_list = [50,50,50]
    //find output of cascade filter
    [num,den] = find_cascade_filter(F_list,B_list,Fs)
    h = time_response(delta_n,num,den)
    y = time_response(x,num,den)
    a = y
    //plot responses and save sound files
//    plot_frequency_response(y,n_fft,F0,Fs)
//    plot_y(y,F0,Fs)  
endfunction



//produce vowels
Fs = 8000
F0_list = [120, 300]

for iii = 1:length(F0_list)
        //window signal near centre
        window_time = 30
        F0 = F0_list(iii)
        x = produce_vowel(F0)'
        N = window_time*Fs/1000 + 1
        n_fft = 10*N
        //applying hamming window
        x_init = (length(x)+1)/2 - (N-1)/2  
        win_hamming = window('hm', N);
        y = x(x_init: x_init + N - 1).*win_hamming;
        //finding autocorrelation
        r = zeros(length(y),1)  
        for ki=0:length(x) 
             for ni=ki:N-1
                r(ki+1)=r(ki+1)+y(ni+1)*y(ni-ki+1);
             end
        end
        
        delta_n = zeros(N,1)
        delta_n(1) = 1
        h_total = []
        //Levinson Durbin recursion
        p_max = 10     
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
            
                
            //Plotting LP approximation overlapping with the original spectrum
            if (modulo(i,2)==0)
               fig = scf()
               h_padded = zeros(n_fft,1)               
               h_padded(1:length(y)) = y'
               H = fftshift(fft(h_padded))    
               H_mag = abs(H(n_fft/2 +1:n_fft))
               freq_array = linspace(0,Fs/2,n_fft/2)
               
               plot(freq_array,20*log(H_mag))
               if(h_total == [])
                   h_total = H_mag
               end
               ax = gca()
               original_frequencies = [730, 1090, 2440]
               temp_d = (n_fft*original_frequencies)/Fs
               temp_d2 = zeros(length(freq_array),1)
               temp_d2(floor(temp_d)+1) = max(ax.y_ticks.locations)
               
               plot(freq_array,temp_d2,'k-.')
               num = sqrt(e_vec(i+1))
               den = -a_vec(2:length(a_vec))
               den = [1,den]
               [h,w] = frmag(num,den,n_fft)
               [h2,w2] = frmag(num,den,n_fft/2)
               h_total = [h_total,h2']
               plot(w*Fs,20*log(abs(h)),'red')
               plot_title = strcat(['LP approximation of vowel a for F0 = ',string(F0),' p = ',string(i)])
               xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
               xs2jpg(gcf(), strcat(['../plots/Q1/',plot_title,'.jpg'])); 
            end
        end
        
        //plotting all LP together
        fig = scf()
        plot2d(w2,20*log(h_total),[1:(p_max/2+1)])
        plot_title = strcat(['All lp plots together for F0 ',string(F0)])
        xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
        legends(['original','p = 2','p = 4','p = 6','p = 8','p = 10'],[1:(p_max/2+1)],opt = 3)
        xs2jpg(gcf(), strcat(['../plots/Q1/',plot_title,'.jpg']));    
end

xdel(winsid())
