% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%          BASIC IMAGE SOURCE METHOD FOR REGULAR SHAPE
%          ROOM REVERBERATION MODELLING
% 
%          USES CONSTANT ABSORPTION COEFFICIENT FOR ALL SURFACES
%          AND PRODUCES A MONO AUDIO OUTPUT
% 
%            RUTHU PREM KUMAR
%            FEBRUARY 2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all;

%% Parameters %%

Fs = 48000;                             % Initialising Sample rate
Cair = 343;                                % Speed of sound in air (m/s)



alpha = 0.1;                             % Acoustic Absorption
R = sqrt(1 - alpha);                     % Reflection coefficient
%% Room dimensions (in meters) %%

% Lengths along x,y and z directions
Lx = 5;                                
Ly = 4;
Lz = 3;

V = Lx*Ly*Lz;                            % Volume of room

A1 = Lx*Ly ; A2 = Lx*Lz ; A3 = Ly*Lz;    % Area of room walls



%% Calculation of T60 of room %%

alpha_bar = 2*alpha*(A1 + A2 + A3);            % Formula for alpha_bar
T_sixty = (24*log(10)*V)/(Cair*alpha_bar);     % Formula for T60
N = (Cair*T_sixty)/min([Lx,Ly,Lz]);            % Number of reflections 

%Initialising vector for impulse respone output
IR = zeros(ceil(T_sixty*Fs),1);


%% Source and Receiver Positions %%

% Position of source in cartesian coordinates(p,q,r)
p = 2;
q = 2;
r = 2;

% Position of observer in cartesian coordinates(a,b,c)
a = Lx/sqrt(2);
b = Ly/sqrt(2);
c = Lz/sqrt(2); 

%% Computation %%

for d = -N:N
    if rem(d,2)~=0                          % Odd case
        Ad = (d+1)*Lx-p-a;
    else                                    % Even case
        Ad = d*Lx+p-a;
    end
    
    for e = -N:N
        if rem(e,2)~=0                      % Odd case
            Be = (e+1)*Ly-q-b;
        else                                % Even case
            Be = e*Lx+q-b;
        end
    
        
        for f = -N:N
            if rem(f,2)~=0                  % Odd case
                Cf = (f+1)*Lz-r-c;
            else                            % Even case
                Cf = f*Lz+r-c;
            end
            
            %% Calculation of impulse response %%
            
            % Calculation of distances
            l = sqrt(Ad*Ad + Be*Be + Cf*Cf); 
            
            % Number of reflections
            w = abs(d) + abs(e) + abs(f);
            
            % Calculation of reflected impulse magnitude
            g = power(R,w)/l;
            
            % Calculation of arrival times
            t = l/Cair;
                                                         
            bin = round(t*Fs);           % Calculation of sample bin 
            if bin <= length(IR)
                IR(bin) = IR(bin) + g;   % Impulse response and its magnitude
            end                            
        end
    end
end

%%  Convolution with input audio %%

%Input Audio
[input_audio,Fs] = audioread('speechdirectsound_48.wav');
% In case of stereo,convert to mono
input_audio = 0.5*sum(input_audio,2);   

% Output Audio
output_audio = conv(input_audio,IR);
            
% Play Output Audio
soundsc(output_audio, Fs);

% Plot impulse response
plot([0:length(IR)-1]/Fs, IR);
title('Plot of Impluse Response vs. time');
xlabel('Time(s)');ylabel('Magnitude');

