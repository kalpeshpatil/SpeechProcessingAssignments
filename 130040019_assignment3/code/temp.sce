 function [startPoint,endPoint] = endpointer(Energy,Thres1,Thres2,Thres3)
 // startPoint - index of starting point of word
 // endPoint - index of ending point of word
 // Energy - energy of the signal
 // Thres1,Thres2,Thres3 - Energy Thresholds

 startPoint = zeros(20,1);
 endPoint = zeros(20,1);
 count = 0;
 decrease = %f;
 point1 = 1;

 for i = 2:length(Energy)

 if(Energy(i) > Thres1 & Energy(i-1) <= Thres1)
 point1 = i;
 decrease = %f;
 end

 if(Energy(i) > Thres2 & Energy(i-1) <= Thres2)
 decrease = %f;
 end

 if(Energy(i) > Thres3 & Energy(i-1) <= Thres3)
 decrease = %f;
 end

 if(Energy(i) < Thres3 & Energy(i-1) >= Thres3)
 decrease = %t;
 end

 if(Energy(i) < Thres2 & Energy(i-1) >= Thres2)
 decrease = %t;
 end

 if(Energy(i) < Thres1 & Energy(i-1) >= Thres1 & decrease)
 point6 = i;
 count = count + 1;
 startPoint(count) = point1;
 endPoint(count) = point6;
 end
 end

endfunction 


Fs = 8000
window_length_ms = 0.025
window_hop_ms = 0.010

name = "divyansh"
file_to_read = "../data/Digits male 8Khz Updated/zero_to_nine_Hitesh1.wav"
[y,Fs,bits] = wavread(file_to_read);
y = y(2000:length(y))



y = y/max(abs(y))

window_len = window_length_ms*Fs
window_hop = window_hop_ms*Fs

t = 1/Fs:1/Fs:(length(y)/Fs)

subplot(3,1,1)
plot(t,y)

sum1 = 0;
energy = 0;
w = window('hm',window_len)

jj = 1
for i = 1:(floor((length(y))/window_hop)-ceil(window_len/window_hop))
    for j = (((i-1)*window_hop + 1):(((i-1)*window_hop)+window_len))
        y(j)=y(j)*w(jj)
        jj = jj + 1
        yy = y(j)*y(j)
        sum1 = sum1 + yy
    end
    energy(i) = sum1
    sum1 = 0
    jj = 1
end

w = 0
c = energy

tt = 1/Fs:(window_hop_ms):(length(energy)*window_hop_ms)

subplot(3,1,2)
plot(tt,energy)

energy_interp = interp1(tt,energy,t,'linear')

[startPoint,endPoint] = endpointer(energy_interp,25e-4,0.05,0.15)



p = zeros(1,length(y))
q = zeros(1,length(y))
p(startPoint) = 1
q(endPoint) = 1


for i = 2:length(startPoint)
    if(endPoint(i-1)-startPoint(i) < 1200)
        p(i) = 0
        q(i-1) = 0 
    end
end


subplot(3,1,3)
plot([p',q',y'])






