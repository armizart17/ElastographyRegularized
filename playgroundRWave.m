%% Loading filtered sono signal
clear, clc
%baseDir = 'C:\Users\sebas\Documents\MATLAB\Elastography';
baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets';
sonoPath = [baseDir,'\heterogeneo_sono'];
swsRange = [2,8];  
move = 'left';

% Pre-processing
load([sonoPath, '\Image9\sono.mat']);
Properties.pitch = 3.0800e-04;
[sono_filt,~,~] = process_sono_data(sono,Properties,move,swsRange);
Prop = Properties;

%% Selecting line
[Nz,Nx,Nt] = size(sono_filt);
depth = floor(Nz*0.4);
frame = floor(Nt/2);
x = (1:Nx)*Prop.pitch * 1000;
sonoLine = squeeze(sono_filt(depth,:,frame));
plot(x,sonoLine)
xlabel('Lateral position [mm]')

%% Peak finding
[b] = peakfinder(sonoLine/max(sonoLine(:)),0);
figure,subplot(2,1,1)
plot(x,sonoLine)
hold on
plot(x(b),sonoLine(b),'o')

if(b(1)==1),b=b(2:end);end
if(b(end)==size(sono_filt,2)),b=b(1:end-1);end
subplot(2,1,2)
plot(x,sonoLine)
hold on
plot(x(b),sonoLine(b),'o')

%% Estimation of SWS
% Se genera un vector con la sws de cada espacio entre picos
lamdas = diff(b);
sws_l = lamdas*2* Prop.pitch * Prop.VibFreq;

SWS = zeros(size(sonoLine));
for i=1:size(sws_l,2) %esto se hace la cantidad de secciones que tengas. Ej: 13 veces
    SWS(b(i)+1:b(i+1))=sws_l(i);
end

%plot(x,SWS)

%% Estimation of weighted coefficients
A = zeros(1,size(sono_filt,2));
n = 1;
clear B;
for i=1:size(sws_l,2) %esto se hace la cantidad de secciones que tengas. Ej: 13 veces
    %%
    %A = zeros(1,size(sono_filt,2)); n = 1; i = 3;

    % Analyzing each wavelength
    y=hilbert(sonoLine(b(i)+1:b(i+1)));     % Analytic signal
    angul=unwrap(atan2(imag(y),real(y)));   % Phase recovery
    [~,p_b]=max(diff(diff(angul)));         % Inflexion point
    p_b = p_b + 1;

    % Plotting
    iCycle = b(i)+1:b(i+1);
    figure, subplot(211), 
    plot(iCycle,[real(y)',imag(y)'])
    xlabel('Samples')
    ylabel('Complex sono signal')
    axis tight
    subplot(212),
    plot(iCycle,angul)
    hold on, xline(iCycle(p_b),'k--'), hold off
    xlabel('Samples')
    ylabel('Phase of sono signal')
    axis tight
    legend({'','Phase inflexion point'}, 'Location','northwest')

    % Weight computation
    weight1 = angul(p_b)/max(angul(:))/p_b;
    weight2 = (max(angul(:))-angul(p_b))/max(angul(:))/(lamdas(i)-p_b);
    A(n,b(i)+1:b(i)+p_b) = weight1;
    A(n,b(i)+p_b+1:b(i+1)) = weight2;

    % SWS values for each wavelength
    B(n)=sws_l(i); n=n+1;
end

%% A and B matrices
close all
figure, 
subplot(1,2,1), stem(B)
xlabel('# of wavelength')
ylabel('SWS')
title('b vector')
subplot(1,2,2), imagesc(A)
title('Weight matrix A')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar


%% Extrapolation of Avg SWS values
% Sample points (SWS values from inside a wavelength)
temp = SWS(b(1)+1:b(end));
x_temp = b(1)+1:b(end);       % Index vector
x1_temp = horzcat(1:(b(1)),(b(end)+1):size(sono_filt,2)); % Query points
boundary = interp1(x_temp,temp,x1_temp,'pchip','extrap');
figure,
plot(x_temp,temp)
hold on
plot(x1_temp,boundary,'b.')

%% Weights of extrapolated values
% The weights must sum 1, so they are 1/(length of the extended part)
A(n,1:b(1))=1/(b(1));   % Weight left side
B(n) = boundary(1);     % Just takes the further boundary
n = n+1;
A(n,(b(end)+1):end) = 1/(size(A,2)-(b(end)+1)+1); % a los ultimos que no tienen valor se les pone ese valor
B(n)=boundary(end);     % no entiendo bien
n=n+1;

% Plotting matrices
figure, 
subplot(1,2,1), stem(B)
xlabel('# of wavelength')
ylabel('SWS')
subplot(1,2,2), imagesc(A)
title('Weights')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar

%% Computing matrix in all frames
clear B;
% Find peaks
A = zeros(1,size(sono_filt,2));
SWS = zeros(size(sono_filt,2),1);
n=1;
for frame=1:size(sono_filt,3)
    % Se extraen los indices de los picos de cada linea (1x128)
    sonoLine=squeeze(sono_filt(depth,:,frame));
    [b]=peakfinder(sonoLine/max(sonoLine(:)),0);

    % Obviamos los picos si estan al inicio o al final
    if size(b,2)<2,continue;end
    if(b(1)==1),b=b(2:end);end
    if(b(end)==size(sono_filt,2)),b=b(1:end-1);end
    if size(b,2)<2,continue;end

    % Se genera un vector con la sws de cada espacio entre picos
    lamdas = diff(b);
    sws_l =lamdas*2* Prop.pitch * Prop.VibFreq;

    % Estimation of weighted coefficients
    for i=1:size(sws_l,2) %esto se hace la cantidad de secciones que tengas. Ej: 13 veces
        SWS(b(i)+1:b(i+1))=sws_l(i);
        y=hilbert(sonoLine(b(i)+1:b(i+1))); % real + img del vector LATERAL
        angul=unwrap(atan2(imag(y),real(y))); % fase recovering, le saca el angulito a y
        [~,p_b]=max(diff(diff(angul))); %posición del máximo valor de diff(diff... -> 2da derivada-> change index of phase velocity
        A(n,b(i)+1:b(i)+p_b)=angul(p_b)/max(angul(:))/p_b; %del inicio hasta la posición p_b lo llenamos de un valor
        A(n,b(i)+p_b+1:b(i+1))=(max(angul(:))-angul(p_b))/... %lo que sigue de p_b hasta antes del sgt pico lo llenamos de otro valor
            max(angul(:))/(lamdas(i)-p_b);
        B(n)=sws_l(i); n=n+1; %va llenando B con los valorsde sws_l
    end

    % Extrapolation of Avg SWS values
    temp=SWS(b(1)+1:b(end)); % mete a esa matriz temporal los valores de sws desde el 7:125 (los no vacíos)
    x_temp=b(1)+1:b(end); %crea un vector desde el numero 7 al 125 (los no vacíos)
    x1_temp=horzcat(1:(b(1)),(b(end)+1):size(sono_filt,2)); %crea un vector temporal concatenando los numeros vacíos -> 1,2,..,6,126,127,128
    boundary=interp1(x_temp,temp,x1_temp,'pchip','extrap'); %coge el primer valor y el ultimo, longitud igual que x1 temp (valores vacios)
    A(n,1:b(1))=1/(b(1));%a los primeros que no tienen valor, le ponen 1/primer posicion de pico
    B(n)=boundary(1); % no entiendo bien
    n=n+1;
    A(n,(b(end)+1):end)=1/(size(A,2)-(b(end)+1)+1); % a los ultimos que no tienen valor se les pone ese valor
    B(n)=boundary(end); % no entiendo bien
    n=n+1;
end

% Plotting matrices
figure, 
subplot(2,1,1), stem(B)
xlabel('# of wavelength')
ylabel('SWS')
xlim([1,50])
subplot(2,1,2), imagesc(A)
set(gca,'cLim',[0.04 0.08])
title('Weight matrix A')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar
ylim([1,50])

%% Adding valleys
% Here we repeat the process but add the measurements from valley to valley
for frame=1:size(sono_filt,3)
    sonoLine=squeeze(-1*sono_filt(depth,:,frame));
    [b]=peakfinder(sonoLine/max(sonoLine(:)),0);
    if size(b,2)<2,continue;end
    if(b(1)==1),b=b(2:end);end
    if(b(end)==size(sono_filt,2)),b=b(1:end-1);end
    if size(b,2)<2,continue;end
    lamdas=diff(b); sws_l=lamdas*2*Prop.pitch*Prop.VibFreq;
    % Estimation of weighted coefficients
    for i=1:size(sws_l,2)
        SWS(b(i)+1:b(i+1))=sws_l(i);
        y=hilbert(sonoLine(b(i)+1:b(i+1)));
        angul=unwrap(atan2(imag(y),real(y)));
        [~,p_b]=max(diff(diff(angul)));
        A(n,b(i)+1:b(i)+p_b)=angul(p_b)/max(angul(:))/p_b;
        A(n,b(i)+p_b+1:b(i+1))=(max(angul(:))-angul(p_b))/...
            max(angul(:))/(lamdas(i)-p_b);
        B(n)=sws_l(i); n=n+1;
    end
    % Extrapolation of Avg SWS values
    temp=SWS(b(1)+1:b(end)); %creo un vector temporal con valores llenos de sws, 1x119
    x_temp=b(1)+1:b(end); %creo un vector temporal con los  indices de los valores llenos de SWS, ej: 7,8,...,125. 1x119
    x1_temp=horzcat(1:(b(1)),(b(end)+1):size(sono_filt,2)); % creo un vector temporal, unión de los indices con SWS vacios, 1,2,...,126,127,128,129. 1x9
    boundary=interp1(x_temp,temp,x1_temp,'pchip','extrap'); % extrapolación de los valores fuera de rango
    A(n,1:b(1))=1/(b(1));
    B(n)=boundary(1);
    n=n+1;
    A(n,(b(end)+1):end)=1/(size(A,2)-(b(end)+1)+1);
    B(n)=boundary(end); n=n+1;
end

%% Iterative algorithm
% k value of 0.1 gives a sharper image but misses some borders
RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;

WA = itreg(A,B,RW_Param);
figure, plot(WA)
ylim([3.6,5.4])

%% Computing whole image
tic
SWS = R_WAVE(sono_filt,Prop,RW_Param);
toc

% Displaying image
SWS_im_range = [2 5.5];
x = 1000*Prop.Width_S;
z = 1000*Prop.Depth_S;
figure('Position', [200 200 320 320]);
imagesc(x,z,SWS,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS R-WAVE, Freq = ' num2str(Prop.VibFreq) ' Hz'])
ax = gca; ax.FontSize = 12;

%% ------------------------------- TESTING -------------------------------
%% Trying new weights
clear B;
A = zeros(1,size(sono_filt,2));
n=1;
for frame=1:size(sono_filt,3)
    % Finding indices for each peak in each line (1x128)
    sonoLine = squeeze(sono_filt(depth,:,frame));
    [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
    
    % Ignoring if peaks are in the start, the end, or there are none
    if size(iPeaks,2)<2, continue; end
    if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end          
    if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
    if size(iPeaks,2)<2, continue; end

    % SWS for each wavelength
    lamdaSamples = diff(iPeaks);
    swsLamda = lamdaSamples*2* Prop.pitch * Prop.VibFreq;

    % Extrapolation at the start (nearest neighbour)
    A(n,1:iPeaks(1)) = 1/(iPeaks(1));
    B(n) = swsLamda(1);
    n = n+1;

    % Estimation of weighted coefficients
    for i=1:size(swsLamda,2) % For each wavelength
        A(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
        B(n) = swsLamda(i); 
        n = n+1;
    end
    
    % Extrapolation at the end
    A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1); 
    B(n) = swsLamda(end); 
    n = n+1;
end

if RW_Param.dual
    % Find valleys too
    for frame=1:size(sono_filt,3)
        sonoLine=squeeze(-1*sono_filt(depth,:,frame));

        % Repeating the same procedure
        [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
        if size(iPeaks,2)<2, continue; end
        if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end
        if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
        if size(iPeaks,2)<2, continue; end
        lamdaSamples = diff(iPeaks);
        swsLamda = lamdaSamples*2* Prop.pitch * Prop.VibFreq;
        A(n,1:iPeaks(1)) = 1/(iPeaks(1));
        B(n) = swsLamda(1);
        n = n+1;
        for i=1:size(swsLamda,2) % For each wavelength
            A(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
            B(n) = swsLamda(i);
            n = n+1;
        end
        A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1);
        B(n) = swsLamda(end);
        n = n+1;
    end
end

% Plotting matrices
figure, 
subplot(2,1,1), stem(B)
xlabel('# of wavelength')
ylabel('SWS')
xlim([1,50])
subplot(2,1,2), imagesc(A)
set(gca,'cLim',[0.04 0.08])
title('Weight matrix A')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar
ylim([1,50])

%% Plotting side profile
RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;

WA = itreg(A,B,RW_Param);
figure, plot(WA)
ylim([3.6,5.4])
%% New function (simple weights and extrapolation)
RW_Param.alpha=2;
tic
SWS2 = rWave2(sono_filt,Prop,RW_Param);
toc

SWS_im_range = [2 5.5];
x = 1000*Prop.Width_S;
z = 1000*Prop.Depth_S;
fig = figure('Position', [200 200 320 320]);
imagesc(x,z,SWS2,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS R-WAVE, Freq = ' num2str(Prop.VibFreq) ' Hz'])
ax = gca; ax.FontSize = 12;

%% Version 3, calculating frequency instead of SWS
clear B;
A = zeros(1,size(sono_filt,2));
n=1;
for frame=1:size(sono_filt,3)
    % Finding indices for each peak in each line (1x128)
    sonoLine = squeeze(sono_filt(depth,:,frame));
    [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
    
    % Ignoring if peaks are in the start, the end, or there are none
    if size(iPeaks,2)<2, continue; end
    if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end          
    if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
    if size(iPeaks,2)<2, continue; end

    % SWS for each wavelength
    lamdaSamples = diff(iPeaks);
    lamdaLengths = lamdaSamples * Prop.pitch;

    % Estimation of weighted coefficients
    for i = 1:size(lamdaLengths,2) % For each wavelength
        A(n,iPeaks(i)+1:iPeaks(i+1)) = 1/lamdaSamples(i);
        B(n) = 1/lamdaLengths(i); 
        n = n+1;
    end
    
    % Extrapolation
    A(n,1:iPeaks(1)) = 1/(iPeaks(1));
    B(n) = 1/lamdaLengths(1);
    n = n+1;
    A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1); 
    B(n) = 1/lamdaLengths(end); 
    n = n+1;
end

if RW_Param.dual
    % Find valleys too
    for frame=1:size(sono_filt,3)
        sonoLine=squeeze(-1*sono_filt(depth,:,frame));

        % Repeating the same procedure
        [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
        if size(iPeaks,2)<2, continue; end
        if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end
        if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
        if size(iPeaks,2)<2, continue; end
        lamdaSamples = diff(iPeaks);
        lamdaLengths = lamdaSamples * Prop.pitch;
        for i=1:size(lamdaLengths,2) % For each wavelength
            spFreq(iPeaks(i)+1:iPeaks(i+1)) = 1/lamdaLengths(i); 
            A(n,iPeaks(i)+1:iPeaks(i+1)) = 1/lamdaSamples(i);
            B(n) = 1/lamdaLengths(i); 
            n = n+1;
        end
        A(n,1:iPeaks(1)) = 1/(iPeaks(1));
        B(n) = 1/lamdaLengths(1);
        n = n+1;
        A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1); 
        B(n) = 1/lamdaLengths(end); 
        n = n+1;
    end
end

% Plotting matrices
figure, 
subplot(2,1,1), stem(B)
xlabel('# of wavelength')
ylabel('SWS')
xlim([1,50])
subplot(2,1,2), imagesc(A)
set(gca,'cLim',[0.04 0.08])
title('Weights')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar
ylim([1,50])

%%
% k value of 0.1 gives a sharper image but misses some borders
RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha = 12;

% Plotting side profile
WA =  2* Prop.VibFreq ./ itreg(A,B,RW_Param);
figure, plot(WA)
ylim([3.6,5.4])

%% Trying new function
RW_Param.alpha = 12;
tic
SWS = rWave3(sono_filt,Prop,RW_Param);
toc

SWS_im_range = [2 5.5];
x = 1000*Prop.Width_S;
z = 1000*Prop.Depth_S;
fig = figure('Position', [200 200 320 320]);
imagesc(x,z,SWS,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS R-WAVE, Freq = ' num2str(Prop.VibFreq) ' Hz'])
ax = gca; ax.FontSize = 12;

%% Functions
function [x,discr,cost]=itreg(A,B,RW_Param)

% Choice of regularization function
switch (RW_Param.operator)
    % Penalizing greater SWS values
    case 'I'
        L=speye(size(A,2));
    % Penalizing greater gradients
    case 'G'
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]);
        L(end,:)=0;
end
B=B';
x0=zeros(size(A,2),1);
err=1;
while err > RW_Param.tolerance
    Lx = L*x0; %Tikhonov matrix
    W = diag( RW_Param.k/2*( abs(Lx.^2+RW_Param.beta).^(RW_Param.k/2 - 1) ) );
    %fprintf("A:%dx%d, B:%d, W:%d\n ",size(A,1),size(A,2), length(B), length(W))

    x = ((A'*A+RW_Param.alpha^2 *L'*W*L)\A') *B;
    err=norm(x-x0)^2/norm(x)^2; %%error relativo
    x0=x;
end
end