function IQ=process_IQ(Data,Properties)
% process_IQ process the RF data in order to obtain IQ signals. To do
% that, an IQ demodulation is apply to the RF data in three steps:
% Down-mixing: The real valued RF-signal is multiplied (“mixed”) with a
% complex sinusoidal signal:x_IQ(t) = x_RF(t)×exp(-i2pi×Fe×t)
%                                   = x_RF(t)×(cos(wt)-i×sin(wt))
% where t is the time along the beam. The down-mixing can be thought of as 
% mixing the RF-signal with two sinusoid signals with 90° phase difference
% Low-pass filtering: After down mixing, the complex signal is low-pass 
% filtered to remove the negative frequency spectrum and noise outside the 
% desired bandwidth
% Decimation: We can reduce the sampling frequency, the smallest integer 
% fraction has to be larger than twice the filter cut-off frequency.
% INPUTS: 
% Data:              variable in which the RF is stored
% Properties:        properties in which you can find the PRF,sampling freq  
%                    central freq...etc
% 
% OUTPUT:
% IQ:                IQ data that has the in-phase and quadrature signals
% 
% Author: Juvenal Ormachea Quispe - March 2014
%         j.ormachea@pucp.pe
% Laboratorio de imagenes medicas-PUCP-PERU


% data information
% PRF= Properties.dr;            % pulse repetition frequency
Fe = Properties.txf;           % emition frequency MHz
Fs = Properties.sf;            % sampling frequency MHz
ensemble = Properties.extra;   % number of RF acquisition 
LineLength = Properties.h;     % RF Line length
numOfLine = Properties.w;      % number of scan lines
factor = 5;
% factor = 1;
%%%%%%%%%%%%%%%%%%% Down-mixing %%%%%%%%%
% Sine and Cosine Table calculation
% Ts = 1/Fs;      % sampling period
c  = 1540;      % speed of sound m/s
D  = Properties.h*c/(Fs*2);    % Depth of window [meter]
T  = D/(c/2);                   % Time of each window [second]
t=linspace(0,T,LineLength);        % Time indixes for each window
SinDem = sin(2*pi*(Fe)*t); % sine table
CosDem = cos(2*pi*(Fe)*t); % cosine table
I     = zeros(LineLength, ensemble,numOfLine,Properties.nframes);
Q     = zeros(LineLength, ensemble,numOfLine,Properties.nframes);
filtI = zeros(LineLength, ensemble,numOfLine,Properties.nframes);
filtQ = zeros(LineLength, ensemble,numOfLine,Properties.nframes);
decI  = zeros((LineLength/factor), ensemble,numOfLine,Properties.nframes);
decQ  = zeros((LineLength/factor), ensemble,numOfLine,Properties.nframes);
Data=reshape( Data, [Properties.h,Properties.extra, Properties.w,  Properties.nframes]);
for frame=1:Properties.nframes
    for ii=1:ensemble
        I(:,ii,:,frame)=squeeze(Data(:,ii,:,frame)).*(CosDem'*ones(1, numOfLine));
        Q(:,ii,:,frame)=squeeze(Data(:,ii,:,frame)).*(-SinDem'*ones(1, numOfLine));
    end
end
%%%%%%%%%%%%%%%%%%% Low-pass filtering %%%%%%%%%
filterOrder = 2;            % order of  filter 
Wn = (1.2*Fe/2)/(Fs/2);     % 20% mas alla de 3.3MHz      
B = fir1(filterOrder,Wn);

for frame=1:Properties.nframes
    for ii=1:ensemble        
        filtI(:,ii,:,frame)=sqrt(2)*filtfilt(B,1,I(:,ii,:,frame));
        filtQ(:,ii,:,frame)=sqrt(2)*filtfilt(B,1,Q(:,ii,:,frame));
    end
end

% for frame=1:Properties.nframes
%     for ee=1:ensemble
%         for jj=1:numOfLine
%             decI(:,ee,jj,frame)=decimate(filtI(:,ee,jj,frame),5);
%             decQ(:,ee,jj,frame)=decimate(filtQ(:,ee,jj,frame),5);
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%% Decimation %%%%%%%%%
for frame=1:Properties.nframes
    for ee=1:ensemble
        A1=squeeze(filtI(:,ee,:,frame));
        A2=squeeze(filtQ(:,ee,:,frame));
        a1=reshape(A1,[],1);
        a2=reshape(A2,[],1);
        aa1=decimate(a1,factor);
        aa2=decimate(a2,factor);
        decI(:,ee,:,frame)=reshape(aa1,[(size(A1,1)/factor),size(A1,2)]);
        decQ(:,ee,:,frame)=reshape(aa2,[(size(A1,1)/factor),size(A1,2)]);
    end
end
IQ=decI+1i*decQ;
end