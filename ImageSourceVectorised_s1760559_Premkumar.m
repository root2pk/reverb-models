% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%        VECTORISED IMAGE SOURCE METHOD FOR REGULAR
%        ROOM REVERBERATION MODELLING
% 
%        USES THE NDGRID() FUNCTION TO CREATE A 3D ARRAY AND PERFORMS 
%        COMPUTATIONS ACCORDINGLY. ALSO COMBINES EVEN AND ODD TERMS
%        INTO ONE EXPRESSION. 
%
%            RUTHU PREM KUMAR
%            FEBRUARY 2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all;

%% Parameters %%

Fs = 48000;                             % Initialising Sample rate
Cair = 343;                             % Speed of sound in air (m/s)


%% Room dimensions (in meters) %%


% Lengths along x,y and z directions
Lx = 5;
Ly = 4;
Lz = 3;

V = Lx*Ly*Lz;                           % Volume of room

A1 = Lx*Ly;                             % Floor and ceiling area
A2 = Lx*Lz ;                            % Area of wall 1 and wall 3
A3 = Ly*Lz;                             % Area of wall 2 and wall 4


%% Absorption and Reflection Coefficients %% (Change as necessary)

% Walls
alpha(1) = 0.1; 
alpha(2) = 0.1;
alpha(3) = 0.1;
alpha(4) = 0.1;

% Ceiling
alpha(5) = 0.1;

% Floor
alpha(6) = 0.1;

% Reflection Coefficients
n = (1:6);
R(n) = sqrt(1-alpha(n));


%% Calculation of T60 of room %%

% Formula for alpha_bar
alpha_bar = alpha(1)*A2 + alpha(2)*A2 + alpha(3)*A3 + alpha(4)*A3 + alpha(5)*A1 + alpha(6)*A1;
T_sixty = (24*log(10)*V)/(Cair*alpha_bar);     % Formula for T60

% Number of reflections Nx,Ny,Nz 

Nx = Cair*T_sixty/Lx;
Ny = Cair*T_sixty/Ly;
Nz = Cair*T_sixty/Lz;

%Initialising vector for impulse respone output (stereo)
IR1 = zeros(ceil(T_sixty*Fs),1);
IR2 = IR1;

%% Source and Receiver Positions

% Position of source (p,q,r)
p = 2;
q = 2;
r = 2;

% Position of Receiver 1 (For stereo output)
a1 = Lx/sqrt(2);
b1 = Ly/sqrt(2);
c1 = Lz/sqrt(2);          

% Position of Receiver 2 (For stereo output)
a2 = Lx/sqrt(2);
b2 = Ly/sqrt(2);
c2 = Lz/sqrt(2) + 0.1;              % 10 cm away

%% Computation %%
%   This method stores every possible combination of values for d,e and f
%   in a 3D Matrix using ndgrid, and computes the corresponding
%   values for Ad,Be and Cf, taking into consideration odd and even terms.

% Creating a 3D Matrix
[d,e,f] = ndgrid(-Nx:Nx,-Ny:Ny,-Nz:Nz);

% Values of Ad,Be and Cf 
% This formula differentiates odd and even values of d,e and f using the
% mod() function and properties of (-1)^
Ad1 = (d + mod(d,2))*Lx + (-1.^d) - a1;
Ad2 = (d + mod(d,2))*Lx + (-1.^d) - a2;
Be1 = (e + mod(e,2))*Lx + (-1.^e) - b1;
Be2 = (e + mod(e,2))*Lx + (-1.^e) - b2;
Cf1 = (f + mod(f,2))*Lx + (-1.^f) - c1;
Cf2 = (f + mod(f,2))*Lx + (-1.^f) - c2;

% Calculation of distances
l1 = sqrt((Ad1).^2 + (Be1).^2 + (Cf1).^2);
l2 = sqrt((Ad2).^2 + (Be2).^2 + (Cf2).^2);

% Number of reflections 
w = abs(d) + abs(e) + abs(f);
                    
% Calculation of modified reflection constant contribution to g

R_total = 1;
for n = 1:6
    R_total = R_total.*(power(R(n),w/6));
end

% Calculation of reflected impulse magnitude
g1 = R_total./l1;
g2 = R_total./l2;
            
% Calculation of arrival times
t1 = l1/Cair;
t2 = l2/Cair;
                        
% Calculation of sample bins 
bin1 = round(t1*Fs);
bin2 = round(t2*Fs);

% Impulse response and its magnitude

% Rearranging the bin and impulse magnitude arrays into vectors and
% rearranging them.
bin1 = bin1(:); [bin1,ind1] = sort(bin1,'ascend');
bin2 = bin2(:); [bin2,ind2] = sort(bin2,'ascend');
g1 = g1(:); g1 = g1(ind1);
g2 = g2(:); g2 = g2(ind2);

% Using accumarray to accumulate magnitude values corresponding to the same
% bin values
accum1 = accumarray(bin1,g1); 
accum2 = accumarray(bin2,g2);

% Putting the values into the corresponding impulse response vectors
for n = 1:length(accum1)
    
    if accum1(n)~=0
        
        if bin1(n)<=length(IR1)        
            IR1(bin1(n)) = accum1(n);
        end
    end
end
for n = 1:length(accum2)
    
    if accum2(n)~=0
        
        if bin2(n)<=length(IR2)
            IR2(bin2(n)) = accum2(n);
        end
    end
end



 
%%  Convolution with input audio %%

% Input Audio 
[input_audio,Fs] = audioread('speechdirectsound_48.wav');
% In case of stereo,convert to mono
input_audio = 0.5*sum(input_audio,2);   

% Output Audio
output_audio(:,1) = conv(input_audio,IR1);
output_audio(:,2) = conv(input_audio,IR2);
            
% Play impulse response
soundsc(output_audio, Fs);

% Plot impulse responses

subplot(2,1,1);
plot((0:length(IR1)-1)/Fs, IR1);
title('Plot of Impluse Response vs. time (Channel 1)');
xlabel('Time(s)');ylabel('Magnitude');

subplot(2,1,2);
plot((0:length(IR2)-1)/Fs, IR2);
title('Plot of Impluse Response vs. time (Channel 2)');
xlabel('Time(s)');ylabel('Magnitude');

