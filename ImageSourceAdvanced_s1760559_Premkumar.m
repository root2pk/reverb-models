% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%       ADVANCED IMAGE SOURCE METHOD FOR REGULAR 
%       ROOM REVERBERATION MODELLING
% 
%       USES DIFFERENT VALUES OF ABSORPTION COEFFICIENTS
%       FOR VARIOUS SURFACES AND PRODUCES A STEREO OUTPUT 
%       FOR TWO RECEIVERS.
%
%       ALSO TAKES WEIGHTED IMPULSE MAGNITUDE AND SPREADS ACROSS TWO 
%       NEAREST INTEGER BINS  
% 
%            RUTHU PREM KUMAR
%            FEBRUARY 2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all;

%% Parameters %%

Fs = 48000;                                % Initialising Sample rate
Cair = 343;                                % Speed of sound in air (m/s)


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
alpha(5) = 0.3;

% Floor
alpha(6) = 0.3;

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

% Position of source in cartesian (p,q,r)
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
c2 = Lz/sqrt(2) + 0.1; % 10 cm away

%% Computation %%

for d = -Nx:Nx
    if rem(d,2)~=0                       % Odd case 
        Ad1 = (d+1)*Lx-p-a1;
        Ad2 = (d+1)*Lx-p-a2;
    else                                 % Even case
        Ad1 = d*Lx+p-a1;
        Ad2 = d*Lx+p-a2;
    end
    
    for e = -Ny:Ny
        if rem(e,2)~=0                   % Odd case
            Be1 = (e+1)*Ly-q-b1;
            Be2 = (e+1)*Ly-q-b2;
        else                             % Even case
            Be1 = e*Lx+q-b1;             
            Be2 = e*Lx+q-b2;
        end
    
        
        for f = -Nz:Nz
            if rem(f,2)~=0               % Odd case
                Cf1 = (f+1)*Lz-r-c1;
                Cf2 = (f+1)*Lz-r-c2;
            else                         % Even case
                Cf1 = f*Lz+r-c1;
                Cf2 = f*Lz+r-c2;
            end
            
            %% Calculation of impulse response %%
            
            % Calculation of distances
            l1 = sqrt(Ad1*Ad1 + Be1*Be1 + Cf1*Cf1);
            l2 = sqrt(Ad2*Ad2 + Be2*Be2 + Cf2*Cf2);
            
            % Number of reflections 
            w = abs(d) + abs(e) + abs(f);
            
            % Calculation of modified reflection constant contribution to g
            R_total = 1;
            
            for n = 1:6
                R_total = R_total*(power(R(n),w/6));
            end
            
            % Calculation of reflected impulse magnitude
            g1 = R_total/l1;
            g2 = R_total/l2;
            
            % Calculation of arrival times
            t1 = l1/Cair;
            t2 = l2/Cair;
            
            % Calculation of sample bin 
            bin1 = t1*Fs;          
            bin2 = t2*Fs;
            
            % Impulse response and its magnitude, using 2-point spreading
            
           if bin1 <= length(IR1)
               
               % If bin falls on integer sample
               if mod(bin1,1)==0                    
                   IR1(bin1) = IR(bin1) + g1;
                   
               % If bin falls on non integer sample
               % Take weighted impulse magnitude and add to neighbouring
               % integer bins
               else                                
                   IR1(floor(bin1)) = IR1(floor(bin1)) + abs(ceil(bin1)-bin1)*g1;
                   IR1(ceil(bin1)) = IR1(ceil(bin1)) + abs(floor(bin1)-bin1)*g1;  
               end
           end
           
           if bin2 <= length(IR2)
               
               % If bin falls on integer sample
               if mod(bin2,1)==0                  
                   IR1(bin2) = IR(bin2) + g2;
                  
               % If bin falls on non integer sample   
               % Take weighted impulse magnitude and add to neighbouring
               % integer bins
               else
                   IR2(floor(bin2)) = IR2(floor(bin2)) + abs(ceil(bin2)-bin2)*g2;
                   IR2(ceil(bin2)) = IR2(ceil(bin2)) + abs(floor(bin2)-bin2)*g2;  
               end
           end
                            
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
            
% Play Output Audio
soundsc(output_audio, Fs);

% Plot impulse responses

subplot(2,1,1);
plot((0:length(IR1)-1)/Fs, IR1);
title('Plot of Impluse Response vs. time (Channel 1)');
xlabel('Time(s)');ylabel('Magnitude');

subplot(2,1,2);
plot((0:length(IR2)-1)/Fs, IR2);
title('Plot of Impluse Response vs. time (Channel 1)');
xlabel('Time(s)');ylabel('Magnitude');


