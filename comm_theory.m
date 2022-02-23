clc 
clear all
close all

a = 9 ;
b = 1 ;
Ns = 2 ;
L = 4 ; %L = 4   L = 8

%--------------------------------------------------------------------------
%Sampling
%xs ---> Sampled signal 

t = 0:0.01:30 ;
x = a*sin(0.5*b*t) ;
ts = t(1:Ns:end) ;
xs = x(1:Ns:end) ;
xq = zeros(1, length(xs)) ;
##plot(ts,xs);
Vmax = a ;
Vmin = -a ; 
delta = (Vmax-Vmin)/L ;

%--------------------------------------------------------------------------
%Quantizing 
%xq ---> Quantized signal 
%qz ---> the levels of the quantizer

function [qz,xq] = Quantize(xs, L, delta)
  
  count = 0 ;
  for i=1:L/2 
    q(i) = delta/2 + count ;
    count = count + delta ;
  end
  
  maxi = -max(q) ;
  
  for i=1:L 
    qz(i) = maxi ;
    maxi = maxi + delta ;
  end   
  
  % we round down if xs is exactly in between the two levels!
  for i=1:length(xs)
    min = 100;
    for j=1:L
      if (abs(qz(j)-xs(i)) < min)
        min = abs(qz(j)-xs(i)) ;
        current = qz(j) ;
      end
      xq(i) = current ;   
    end    
  end
  
endfunction  

[qz, xq] = Quantize(xs, L, delta) ; 

%the mean absolute quantization error 
mean_abs_qe = mean(abs(xq-xs)) ;
nl = [2 4 8 16 32 64] ;
mean = [2.7364 1.2858 0.6226 0.302 0.148 0.072995] ;
##figure()
##plot(nl, mean) ;
##xlabel("The number of levels") ;
##ylabel("The mean absolute Q. Error") ;

%the practical variance of the quantization error 
var_pr_qe = var(abs(xq-xs)) ;
var_pr = [1.8958 0.4781 0.1161 0.028506 7.0072e-03 1.7303e-03] ;

%the theoretical variance of the quantization error
var_theo_qe = (delta)^2 / 12 ;
var_theo = [6.75 1.6875 0.4219 0.1055 0.026367 6.5918e-03] ;

##figure()
##plot(nl, var_pr) ;
##hold on 
##plot(nl, var_theo) ;
##xlabel("The number of levels") ;
##ylabel("The variance of Q. Error") ;
##legend("Practical", "Theoretical") ;

%the signal to quantization noise ratio
sqnr = 3 * (L)^2 
sq_theo = [12 48 192 768 3072 12288] ;

%the signal to quantization noise ratio practical
sqnr_p = (a)^2 / var_pr_qe 
sq_p = [42.726 169.4 697.49 2481.5 1.1560e+04 4.6812e+04] ;

##figure()
##plot(nl, sq_theo) ;
##hold on 
##plot(nl, sq_p) ;
##xlabel("The number of levels") ;
##ylabel("The SQNR") ;
##legend("Theoretical", "Practical") ;

%--------------------------------------------------------------------------
%Encoding 
%xe ---> Encoded signal
%s---> the levels(symbols) of the quantized signal

for i=1:length(qz)
  s(i)=i ;  % if on octave, please make it s(i) = i 
end

for i=1:length(xq)
  xe(i) = s(find(qz==xq(i))) ;
end  

%--------------------------------------------------------------------------
%Huffmann Source encoding
%xhe ---> The huffmann encoded signal
%s ---> The symbols 
%p ---> Their respective probabilities

for i=1:length(s)
  p(i) = length(find(xe==s(i))) / length(xe) ;
end

dict = huffmandict(s, p) ;
xhe = huffmanenco(xe, dict) ;

%--------------------------------------------------------------------------
%Huffmann Source Decoding
%xhd ---> The huffmann decoded signal

xhd = huffmandeco(xhe,dict);
isequal(xe,xhd) ; % the logical AND between input to source encoder and output; equals 1

%--------------------------------------------------------------------------
%Decoding
%xd ---> The decoded signal 

for i= 1:length(xhd)
     xd(i)=qz(find(s==xhd(i))) ;
end

isequal(xd, xq) ;

%--------------------------------------------------------------------------
%The input/output Graph 

figure()
plot(t, x)
hold on 
plot(ts, xd)
xlabel("time (s)")
ylabel("Signal")
legend("The input signal", "The output signal")

%--------------------------------------------------------------------------
%Calculating efficiency and compression rate

SI = -1*log2(p) ;
H = sum(p.* SI) ;
Vk = ceil(SI) ;
L = sum(Vk.* p) ;
n = max(max(Vk)) ;

E = H / L ;
CR = n / L ;
%--------------------------------------------------------------------------





 






