%% PPT_16   Design a radio propagation simulator that is capable to find the median path loss 
% using Okumura’s model in a suburban area for a given T-R distance- d, base station
% antenna height-hte, mobile antenna height-hre, transmitted RF power and carrier 
% frequency. Also find the RF power at the receiver. Show the variation of path loss 
% obtained using Okumura’s model with different carrier frequencies and analyze the graph. 
% Take proper assumptions wherever needed. Put values that are needed from the 
% Okumura curves in a spread sheet and read from the spread sheet.  (6)

clc;
clear all;
close all;

f= 150:1:1920;   %frequency 150Mhz<f<1920Mhz
Hte=50;         % transmitter height
Hre=5;          % Mobile Antenna Height
d = 30;         % Distance 30 Km
    
c=3*10^8;
lamda=(c)./(f.*10^6);
Lf = 10*log((lamda.^2)./((4*pi)^2)*d^2); %   Free Space Propagation Loss
Amu = 35;       % Median Attenuation Relative to Free Space (900 MHz and 30 Km)
Garea = 9;      % Gain due to the Type of Environment (Suburban Area)
Ghte = 20*log(Hte/200); % Base Station Antenna Height Gain Factor
if(Hre>3)
Ghre = 20*log(Hre/3);
else
Ghre = 10*log(Hre/3);
end
%   Propagation Path Loss
L50 = Lf+Amu-Ghte-Ghre-Garea;
display('Propagation pathloss is : ');
disp(L50);
plot(f,L50,'LineWidth',1.5);
title('Okumura Model');
xlabel('Frequency(MHz)');
ylabel('Propagation Path loss (dB)');
grid on;


%% PPT_16  Design a radio propagation simulator that is capable to find the median path loss 
% using Okumura’s model in a suburban area for a given T-R distance- d, base station 
% antenna height-hte, mobile antenna height-hre, transmitted RF power and carrier 
% frequency. Also find the RF power at the receiver. Show the variation of path loss 
% obtained using Okumura’s model with different base station antenna height and analyze
% the graph. Take proper assumptions wherever needed. Put values that are needed from 
% the Okumura curves in a spread sheet and read from the spread sheet.  (7)



clc;
clear all;
close all;

Hte=30:1:100;    % Base Station Antenna Height 
Hre=5;          % Mobile Antenna Height
d = 30;         % Distance 30 Km
f=900;       %frequency 150Mhz<f<1920Mhz
c=3*10^8;
lamda=(c)/(f*10^6);
Lf = 10*log((lamda^2)/((4*pi)^2)*d^2); %   Free Space Propagation Loss
Amu = 35;       % Median Attenuation Relative to Free Space
Garea = 9;      % Gain due to the Type of Environment (Suburban Area)
Ghte = 20*log(Hte/200); % Base Station Antenna Height Gain Factor
if(Hre>3)
Ghre = 20*log(Hre/3);
else
Ghre = 10*log(Hre/3);
end
%   Propagation Path Loss
L50 = Lf+Amu-Ghte-Ghre-Garea;
display('Propagation pathloss is : ');
disp(L50);
plot(Hte,L50,'LineWidth',1.5);
title('Okumura Model');
xlabel('Transmitter antenna Height (Km)');
ylabel('Propagation Path loss (dB)');
grid on;




%% PPT_16  Design a radio propagation simulator that is capable to determine median path loss for
% urban areas using Hata model for a given T-R distance- d, base station antenna 
% height-hte, mobile antenna height-hre, transmitted RF power and carrier frequency-fc. 
% The carrier frequency will be less than 300MHz. Plot graphs between the path loss 
% obtained using Hata model with different carrier frequencies less than 300MHz and 
% analyze the graph.  (8) 


clc;
clear all;
close all;
                    

Hte=150;      
Hre=10;            
d =20000;         
f=1:1:300;  %Frequency                  
cm=0;
city_type= 1; % 1 for small cities and 2 for large cities   
for i=1:length(f)
if(city_type==1)                                                                                               % for small city
CH = 0.8 +((1.1.*log(f(i)))-0.7).*(Hre) - 1.56.*log(f(i));
else                                                                                                           % for large city
if(f(i)>=150 && f(i)<=200)
CH=8.29.*(log(1.54*Hre)).*2-1.1;
elseif(f(i)>200 && f(i)<=1500)
CH=3.2.*(log(11.75*Hre)).*2-4.97;
end
end
LU1(i)=69.55+26.16.*log(f(i))-13.82.*log(Hte)-CH+(44.9-6.55.*log(Hte)).*log(d);                                     % loss for hata model
% LU2(i)=46.3+33.9.*log(f(i))-13.82.*log(Hte)-CH+(44.9-6.55.*log(Hte)).*log(d)+cm;                                    % loss for extended hata model
end

figure(1)
plot(f,LU1,'Linewidth',3)
title('Hata Model Frequency (MHz) vs Loss (dB) ');
xlabel('Frequency (MHz)');
ylabel('Propagation Path loss(dB)');
grid on;



%% PPT_16 Design a radio propagation simulator that is capable to determine median path loss for 
% suburban areas using Hata model for a given T-R distance- d, base station antenna 
% height-hte, mobile antenna height-hre, transmitted RF power and carrier frequency-fc.
% Plot graphs between the path loss obtained using Hata model in suburban area with 
% different carrier frequencies and for different T-R separation and
% analyze the graph. (9)



clc;
clear all;
close all;


Hte=150;      
Hre=10;
dis=1000:20000;    
d = 20000;     
f=1:1:300;  %Frequency                  
cm=0;
city_type= 1; % 1 for small cities and 2 for large cities   
for i=1:length(f)
if(city_type==1)                                                                                         
CH = 0.8 +((1.1.*log(f(i)))-0.7).*(Hre) - 1.56.*log(f(i));
else                                                                                                        
if(f(i)>=150 && f(i)<=200)
CH=8.29.*(log(1.54*Hre)).*2-1.1;
elseif(f(i)>200 && f(i)<=1500)
CH=3.2.*(log(11.75*Hre)).*2-4.97;
end
end

LU1(i)=69.55+26.16.*log(f(i))-13.82.*log(Hte)-CH+(44.9-6.55.*log(Hte)).*log(d)-2.*(log(f(i)./28))^2 - 5.4;                              
% LU2(i)=46.3+33.9.*log(f(i))-13.82.*log(Hte)-CH+(44.9-6.55.*log(Hte)).*log(d)+cm;                           
end

CH2=47;
freq=300;

for i=1000:length(dis)
   LU2(i)=69.55+26.16*log(freq)-13.82*log(Hte)-CH2+(44.9-6.55*log(Hte)).*log(dis(i)) -2*(log(freq/28))^2 - 5.4;   
end



figure(1)
subplot(2,1,1)
plot(f,LU1,'Linewidth',3)
title('Hata Model Frequency (MHz) vs Loss (dB) ');
xlabel('Frequency (MHz)');
ylabel('Propagation Path loss(dB)');
grid on;

subplot(2,1,2)
plot(dis,LU2,'Linewidth',3)
title('Distance Separation (m) vs Loss (dB) ');
xlabel('Distance Separation (m)');
ylabel('Propagation Path loss(dB)');
grid on;



%% PPT_18  Implement even parity check coder and decoder using Matlab and demonstrate
% the results.  (10)


clc;
clear all;
close all;


x = input('Enter the bit sequence to test for Even parity: ');
t = 0;
for i = 1:length(x) 
%also do by 'for' loop just by t=sum(x)        
if(x(i))
        t = t + 1; %increment by one if a bit is one
        end
end

if(mod(t,2)~=0) %check if not even then attach another '1' to make the parity even
   y = [x 1]; disp('Parity bit generated : 1');
else %check if even then attach another '0' to let the parity be even
  y = [x 0]; disp('Parity bit generated : 0');
end
disp('Input bit sequence:');
disp(x); %display the input bit sequence
disp('Bit sequence with parity (Even) bit : ');
disp(y); %display the resultant bit sequence after parity bit addition

% Decoder
disp("Decoder Result");
% Y=sum(y(:) == 1);
parity_bit=y(end);
y(end)=[];
no_of_ones= mod(sum(y(:) == 1),2);


if(length(y)==length(x))
 if (parity_bit==1 &&  no_of_ones ==1)
     disp("Even parity with parity bit 1")
  elseif (parity_bit==0 &&  no_of_ones ==0)
     disp("Even Parity with parity bit 0")
 else
     disp("Error detected")
 end
else
    disp("Error detected")
end


% Enter the bit stream : [1 0 1 1 0 1 1]


%% PPT_18  Implement odd parity check coder and decoder using Matlab and demonstrate 
% the results.   (11)

clc;
clear all;
close all;

x = input('Enter the bit sequence to test for Odd parity: ');
t = 0;
for i = 1:length(x)       
if(x(i))
        t = t + 1; %increment by one if a bit is one
        end
end

if(mod(t,2)~=1) %check if not odd then attach another '1' to make the parity odd
   y = [x 1]; disp('Parity bit generated : 1');
else %check if odd then attach another '0' to let the parity be odd
  y = [x 0]; disp('Parity bit generated : 0');
end
disp('Input bit sequence:');
disp(x); %display the input bit sequence
disp('Bit sequence with parity (Odd) bit : ');
disp(y); 


% Decoder
disp("Decoder Result");
% Y=sum(y(:) == 1);
parity_bit=y(end);
y(end)=[];
no_of_ones= mod(sum(y(:) == 1),2);

if(length(y)==length(x))

 if (parity_bit==1 &&  no_of_ones ==0)
     disp("Odd Parity with parity bit 1")
 elseif (parity_bit==0 &&  no_of_ones ==1)
     disp("Odd Parity with parity bit 0")
 else
     disp("Error detected")
 end
else
     disp("Error detected")
 end



%  Enter the bit stream :[1 0 1 1 0 1 0]

%% PPT_18  Implement a CRC encoder and decoder. Demonstrate the results of both. (12)

clc;
clear all;
close all;

% msg = randi([0 1],6,2);
% disp(msg)

%msg=[1 1 1 0 0 0 1 1 ]
msg=  [1 1 0 1 0 1 1 0 1 1]; %input('Input Message sequence :');
%poly=[1 1 0 0 1 1]
poly= [1 0 0 1 1]; %input('Input Generator Polynomial :');
[M N]=size(poly);
mseg=[msg zeros(1,N-1)];
[q r]=deconv(mseg,poly);
r=abs(r);

for i=1:length(r)
    a=r(i);
    if ( mod(a,2)== 0 )
        r(i)=0;
    else
        r(i)=1;
    end
end
crc=r(length(msg)+1:end);
Tx = bitor(mseg,r);  % Transmitter


% Decoder

[quo rem]= deconv([Tx],poly);  % Receiver

rem=abs(rem);

for i=1:length(rem)
    a=rem(i);
    if ( mod(a,2)== 0 )
        rem(i)=0;
    else
        rem(i)=1;
    end
end
crc=rem(length(msg)+1:end);
Rx = bitor(Tx,rem);



disp('Transmitted Frame')
disp(Tx)

disp('Received Frame')  % Here I have not allowed user to give input
disp(Tx)

disp("CRC Detection Result")
    if ( Tx == Rx )
        disp("No Error");
    else
        disp("Error Detected")
    end




%% PPT_18  Write a Matlab code to determine the hamming distance of the coding 
% scheme with a set of code words.  (13)

clc;
clear all;
close all;

x=input('Enter the first string:');
y=input('Enter the second string:');
x1=double(x);
y1=double(y);
c=0;
for i=1:length(x1)
    if(x(i)~=y(i))
        c=c+1;
    end
end
disp("Hamming Distance");
disp(c);

% Enter the first string:'hi'
% Enter the second string:'hello'


%% PPT_18   Implement Hamming code encoder and decoder using Matlab/C/C++. 
% Take appropriate assumptions on the length of code words.  (14)
clc;
clear all;
close all;

M = [1 0 0 1]; 
C = encode(M,7,4);  % Encode
disp('Input Message')
disp(M);
E = [ 0 0 0 1 0 0 0];  % single bit error
CC = bitxor(C,E);  % Add the single bit error to the codeword
de=decode(CC,7,4);  % Decode
disp('Decoded message')
disp(de);



%% PPT_18  Implement block interleaving using MATLAb or using Simulink and 
% demonstrate the results. (15)

clc;
clear all;
close all;

st1 = 27221; st2 = 4831; % States for random number generator
n = 7; k = 4; % Parameters for Hamming code
msg = randi([0 1],k*500,1); % Data to encode
code = encode(msg,n,k,'hamming/binary'); % Encoded data
% Create a burst error that will corrupt two adjacent codewords.
errors = zeros(size(code)); errors(n-2:n+3) = [1 1 1 1 1 1];

% With Interleaving
%------------------
inter = randintrlv(code,st2); % Interleave.
inter_err = bitxor(inter,errors); % Include burst error.
deinter = randdeintrlv(inter_err,st2); % Deinterleave.
decoded = decode(deinter,n,k,'hamming/binary'); % Decode.
disp('Number of errors and error rate, with interleaving:');
[number_with,rate_with] = biterr(msg,decoded) % Error statistics

% Without Interleaving
%---------------------
code_err = bitxor(code,errors); % Include burst error.
decoded = decode(code_err,n,k,'hamming/binary'); % Decode.
disp('Number of errors and error rate, without interleaving:');
[number_without,rate_without] = biterr(msg,decoded) % Error statistics