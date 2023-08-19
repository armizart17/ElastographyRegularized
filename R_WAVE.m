% Regularized Wavelength Velocity Estimator (R-WAVE) using iterative
% Tikhonov regularization for the inversion problem of each depth
% Inputs:
%
%   sono_filt:      normalized sonoelasticity data, 3D matrix
%
%   Prop.
%       Freq:      vibration freq for CrW generation           [Hz] 
%       Pitch:     spacing between transducer elements         [m]
%
%   RW_Param.
%       dual:       TRUE = valleys included                    [boolean]
%       k:          norm order                                 [ ]
%       beta:       offset to ensure convex solution           [ ]
%       tolerance:  stopping criteria                          [ ]
%       operator:   I = Identity, G = Gradient                 [ ]
%       alpha:      regularization coefficient                 [ ]
%   
% Outputs:
%
%   WA:             Shear Wave speed image                      % [m/s]
%
% Author:           Eduardo A. Gonzalez

function WA=R_WAVE(sono_filt,Prop,RW_Param)
WA = zeros(size(sono_filt,1),size(sono_filt,2));

for depth=1:size(sono_filt,1)
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

    if RW_Param.dual
    %     % Find valleys too
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
    end
    %B=B;
    WA(depth,:)=itreg(A,B,RW_Param);
    %depth
end

end
%%
function [x,discr,cost]=itreg(A,B,RW_Param)
switch (RW_Param.operator)
    case 'I'
        L=speye(size(A,2));
    case 'G'
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]);
        L(end,:)=0;
end
B=B';
x0=zeros(size(A,2),1);
err=1;
while err > RW_Param.tolerance
    W = zeros(size(A,2));
    Lx = L*x0; %Tikhonov matrix
    for i=1:numel(x0) %el segundo sumando de xsombrero j
        W(i,i)=RW_Param.k/2*(abs(Lx(i))^2+RW_Param.beta)^((RW_Param.k/2)-1);
    end
    x = ((A'*A+RW_Param.alpha^2 *L'*W*L)\A') *B;
    err=norm(x-x0)^2/norm(x)^2; %%error relativo
    x0=x;
end

end