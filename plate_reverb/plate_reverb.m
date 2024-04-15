%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&
%  PBMMI ASSIGNMENT 6
%  MODAL PLATE REVERBERATION 
%  MAIN SCRIPT
%  
%  PROGRAM THAT SIMULATES THE MODAL PLATE
%  REVERBERATION MODEL USING A 
%  LINEAR PLATE MODEL AND FINITE 
%  DIFFERENCE SCHEME
%
%  CAN BE COMPUTED EITHER USING 
%  THE REGULAR METHOD OR EXACT METHOD 
% 
% 
%          RUTHU PREM KUMAR
%            APRIL 2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&

clear all; close all;

%% User defined Parameters

Fs = 44100;           % Sample Rate (kHz)
update_method = 1;    % 1 - Regular method, 2 - 'Exact' method
input_type = 2;       % 1 - Audio Clip, 2 - Kronecker Delta IR
material = 2;         % 1 - Steel, 2 - Aluminium

% Plate Parameters

Lx = 2;               % Plate width(m)
Ly = 1;               % Plate height(m)
a = 5e-4;             % Plate Thickness (m)
T = 700;              % Tension per unit length of the plate (N/m)

if material==1   
    rho = 8000;           % Plate Density (kg/m^3)    
    E = 2e11;             % Young's Modulus of plate (N/m^2)
    v = 0.3;              % Poisson's ratio
else
    rho = 2710;           % Plate Density (kg/m^3)
    E = 6.9e10;           % Young's Modulus of plate (N/m^2)
    v = 0.334;            % Poisson's ratio
end
    
% T60 values (seconds)

T60min = 1;             
T60max = 4;

% Normalized Input and Output Coordinates (between 0 and 1)

input_param = [0.2,0.6];        % Input forcing signal [x,y]
outputL = [0.2,0.8];            % Left Output signal [x,y]
outputR = [0.7,0.3];            % Right Output signal [x,y]

%% Input Signal

if input_type==1
    % Audio Input
    [F,Fs]=audioread('piano_short.wav');
    
    % In case of stereo, to mono
    F = 0.5*sum(F,2);
    
else
    % Kronecker Delta impulse input
    F = [1;0];
end
    
%% Derived Parameters

% Simulation Parameteres

k = 1/Fs;                        % Sample Period (s)
wmax = 2/k;                      % Maximmum value of w
fmax = wmax/(2*pi);              % Corresponding maximum value of f

% Exact input and output positions based on Lx and Ly

xi = input_param(1)*Lx; yi = input_param(2)*Ly;   % Actual position of input signal(m)
xoL = outputL(1)*Lx; yoL = outputL(2)*Ly;         % Actual position of left output signal(m)
xoR = outputR(1)*Lx; yoR = outputR(2)*Ly;         % Actual position of right output signal(m)

% Plate Parameters

c = sqrt(T/(rho*a));
K = sqrt((E*a^2)/(12*rho*(1-v^2)));    % Stiffness Factor

% Loss Parameters
sigma_max = 6*log(10)/T60min;          % Max value of sigma
sigma_min = 6*log(10)/T60max;          % Min value of sigma

beta_max_sq = (sqrt(c^4+4*K^2*wmax^2)-c^2)/(2*K^2);  % Max value of beta^2 (at Qx,Qy)
beta_min_sq = (pi/max([Lx,Ly]))^2;                   % Min value of beta^2 ((0,1)or(1,0))

% Loss coefficients
epsilon1 = (sigma_max - sigma_min)/(beta_max_sq - beta_min_sq);
epsilon0 = sigma_max - epsilon1*beta_max_sq;

% Assert that they are non negative
assert(epsilon0 >= 0);
assert(epsilon1 >= 0);

%% Finding number of modes to simulate

% Maximum horizontal limit is when qy = 1 and beta=beta_max
Qx = floor(sqrt(beta_max_sq-(pi/Ly)^2)*Lx/pi);
% Maximum vertical limit is when qx = 1 and beta=beta_max
Qy = floor(sqrt(beta_max_sq-(pi/Lx)^2)*Ly/pi);

% 2D grid
[qx,qy]=meshgrid(1:Qx,1:Qy);

%% Calculations and loop

% Vector of all possible values from q = 1...Q

beta_sq = (qx*pi/Lx).^2+(qy*pi/Ly).^2;         % beta^2
wq_sq = c^2*beta_sq + K^2*(beta_sq.^2);        % wq^2

% Using a logic mask to limit values where wq>w_max

beta_sq = beta_sq(wq_sq<wmax^2);
wq_sq = wq_sq(wq_sq<=wmax^2);
sigma_q = epsilon0+epsilon1*beta_sq;

% Phi values at input and output locations

phi_input=(2/sqrt(Lx*Ly))*sin(qx*pi*xi/Lx).*sin(qy*pi*yi/Ly);
phi_input=phi_input(wq_sq<wmax^2);

phi_outputL=(2/sqrt(Lx*Ly))*sin(qx*pi*xoL/Lx).*sin(qy*pi*yoL/Ly);
phi_outputL=phi_outputL(wq_sq<wmax^2);
  
phi_outputR=(2/sqrt(Lx*Ly))*sin(qx*pi*xoR/Lx).*sin(qy*pi*yoR/Ly);
phi_outputR=phi_outputR(wq_sq<wmax^2);

%% Calculations

% Finding length of output signal
dur = length(F) + T60max*Fs;
outputL = zeros(dur,1);
outputR = outputL;

% Extending force input with zeros to fit into loop
F = vertcat(F,zeros(T60max*Fs,1));

v0=0;          % Initial Velocity
x0=0;          % Initial Displacement      
p2=x0;         % Value at time step n=1
p1=x0+(k*v0);  % Value at time step n=2

% Parameters for Exact method

pow = sqrt((sigma_q.^2)-wq_sq)*k;
coeff1 = exp(-sigma_q*k).*(exp(pow)+ exp(-pow));
coeff2 = exp(-2*sigma_q*k);

 for n=1:dur
    % Accurate Scheme
    if update_method==1
        p0 = ((2-k^2*wq_sq).*p1 + (k*sigma_q-1).*p2 + k^2*phi_input*F(n))./(1 + k*sigma_q);  
    
    % Exact Scheme
    elseif update_method==2
        p0 = coeff1.*p1 - coeff2.*p2 + k^2*phi_input*F(n);
      
    end
          
     % Left Output
     outputL(n)=sum(p0.*phi_outputL,'all');
     % Right Output
     outputR(n)=sum(p0.*phi_outputR,'all');
     
     % Value updates fo next loop    
     p2=p1;
     p1=p0; 
       
end
 

% Normalize Output channels 
outputL = outputL/max(abs(outputL));
outputR = outputR/max(abs(outputR));

% Stereo Signal 
output_signal=[outputL outputR];

% Play Output Signal
soundsc(output_signal,Fs);

%% Write Output Signal

if material == 1 && input_type == 1
    audiowrite('steel_audio.wav',output_signal,Fs);
elseif material == 2 && input_type == 1
    audiowrite('aluminium_audio.wav',output_signal,Fs);
elseif material == 1 && input_type == 2
    audiowrite('steel_delta.wav',output_signal,Fs);
else
    audiowrite('aluminium_delta.wav',output_signal,Fs);
end

%% Spectrogram
if input_type == 2
    
    % Stereo to mono
    signal = outputL + outputR;
    NFFT=2*1024;
    
    % Window length
    window_length=floor(2*round(0.031*Fs));

    % Hann window 
    n = 1:window_length;
    window(n)= 0.5*(1-cos(2*pi*n/window_length)); 

    % Number of windows samples without overlapping (0.45 overlap factor)  
    overlap=floor(0.45*window_length); 
    
    % Spectrogram
    [S,F,T] = spectrogram(signal,window,window_length-overlap,NFFT,Fs); 
    
    % Properties
    [Nf,Nw]=size(S);
    figtitle1 = 'Impulse response Spectrogram';
    figure('name',figtitle1)
    factor_colour=0.000001;
    S_one_sided=max(S(1:fix(length(F)/2),:),factor_colour); 
    pcolor(T,F(1:fix(Nf/2)),10*log10(abs(S_one_sided)));
    shading interp;
    colormap('jet');
    colorbar;
    title('Impulse response');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    if material == 1
        image_filename='steel_IRspec.png';
        saveas(gcf,image_filename);
    else
        image_filename='aluminium_IRspec.png';
        saveas(gcf,image_filename);
    end
end
