clear all

//function to read previously generated wavefiles
function y = read_vowel(vowel,F0)
    wavefile = strcat(['Output signal of vowel ',vowel,' for F0 = ',string(F0),'.wav'])
    y = loadwave(strcat(['../wav_files/Q4/',wavefile]));
endfunction


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


window_times = [5, 10, 20, 40] //time in ms
Fs = 16000
//finding window lengths
window_lengths = window_times*Fs/1000 + 1
n_fft = max(window_lengths)*10


//vowel_list = ['a','i','u','a','i','u']
//F0_list = [120 120 120 220 220 220]
vowel_list = ['u']
F0_list = [220]
for p = 1:length(F0_list)
    vowel = vowel_list(p)
    F0 = F0_list(p)
    for window_type = ['hamming','rect']
        for i = 1:length(window_times)
           N = window_lengths(i)
           x = read_vowel(vowel,F0)
           [X_mag] = find_windowed_FFT (x,window_type,N,n_fft)
           freq_array = linspace(0,Fs/2,n_fft/2)
           fig = scf()
           plot(freq_array,20*log(X_mag))
           plot_title = strcat(['FFT for vowel ',vowel,' using ',window_type,' window of ',string(window_times(i)),' ms',' for F0 = ',string(F0)])
           xtitle(plot_title,'Frequency (Hz)','Magnitude in dB')
           xs2jpg(gcf(), strcat(['../plots/Q5/',plot_title,'.jpg']));
        end
    end
end
