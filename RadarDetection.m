clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77 GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Definition
v_max = 100;
rangeRes = 1;
rangeMax = 200; 
c = 3e8; 


%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
% Target initial position and velocity
targetPosInit = 100; 
targetVel = 50; 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B = c/(2*rangeRes); 
T_chirp = 5.5*2*rangeMax/c;
alpha = B/T_chirp; 

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                      
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps 

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*T_chirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    targetPos = targetPosInit + t(i)*targetVel; 
    r_t(i) = targetPos;
    tau = (2*targetPos) / c; % trip time for the received signal
    td(i) = tau;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + (alpha*t(i)*t(i))/2)); 
    Rx(i)  = cos(2*pi*(fc*(t(i) - tau) + (alpha*((t(i) - tau)*(t(i) - tau)))/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i); 
end


%% RANGE MEASUREMENT


% *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr,Nd]); 

% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
dim = 1; % for operating along the Nr dimension
%signal_fft = fft(Mix,[],dim); 
signal_fft= fft(Mix,Nr);


% *%TODO* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft/Nr); 

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:Nr/2+1,:); 

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% *%TODO* :
% plot FFT output 
r = linspace(1, (Nr/2)+1, (Nr/2)+1);  % Range resolution is 1 m, so the last point in the FFT signal would correspond to Nr m.   
bin_1 = reshape(signal_fft(:,1),[1,(Nr/2)+1]); 
plot(r, bin_1); 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 8; % Range
Td = 8;  % Doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4; % Range
Gd = 4; % Doppler

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 10; 

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR

% Create dedicated matrix for CFAR outcome
CFAR_Map = zeros(Nr/2,Nd);

for i = Tr+Gr+1:(Nr/2)-(Gr+Tr) % Iterate over RANGE dimension but cut don't consider the training and guard cell boundaries 
    for j = Td+Gd+1:Nd-(Gd+Td) % Iterate over DOPPLER dimension but cut don't consider the training and guard cell boundaries
       
        % Set the noise back to zero for every new iteration (on the CUT)
        noise_level = zeros(1,1);
        
        % Inner loop to iterate over every single cell surrounding
        % the CUT in order to sum up the noise 
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                % If a guard cell has been reched, skip the noise level
                % calculation. Otherwise, transform the 2D FFT matrix value
                % into linear space and sum it up
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        
        % Normalize the noise value, transform it back to logarithmic space and add the offset
        numTrainingCells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1); 
        threshold = pow2db(noise_level/numTrainingCells);
        threshold = threshold + offset;
        
        % Compare CUT against the calculated threshold
        CUT = RDM(i,j);
        if (CUT < threshold)
            CFAR_Map(i,j) = 0; 
        else
           CFAR_Map(i,j) = 1; 
        end
    end
end

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis, CFAR_Map);
%colorbar;

