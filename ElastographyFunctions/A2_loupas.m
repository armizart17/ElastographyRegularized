%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particle velocity estimates using Loupas' algorithm
%
% IQ: The complex IQ data.
% ROI_ax, ROI_lat, ROI_fr: The coordinates of the ROI cuboid.
% =========================================================================
function [ u, disp_axis, ROI ] = A2_loupas( IQ, M, N, dW )

%%%  Acquisition parameters
% c0 = evalin('base','Resource.Parameters.speedOfSound');
% fc = evalin('base','Trans.frequency')*1e6;
% dz = evalin('base','dinf.dz');
c0 = 1540;
fc = 7240000;
dz = 3e-04;


% =========================================================================
%%%  Compute mean axial displacement (with adjustable window shift dW)

%  Define lateral ROI
ROI_lat =  1 : 1 : size(IQ,2);

%  Define temporal ROI
ROI_fr  =  1 : 1 : size(IQ,3)-N;

%  Define axial ROI
if mod(M,2)  % M is odd
    M1 = floor(M/2);
    M2 = floor(M/2);
    ax_s = floor(M/2)+1;           % axial index start
    ax_e = size(IQ,1)-floor(M/2);  % axial index end
else         % M is even
    M1 = M/2;
    M2 = M/2-1;
    ax_s = M/2+1;
    ax_e = size(IQ,1)-M/2+1;
end
ROI_ax  = ax_s : dW : ax_e;


%  ROI
ROI.ROI_ax  = ROI_ax;
ROI.ROI_lat = ROI_lat;
ROI.ROI_fr  = ROI_fr;


%  Define displacement axis in meters
disp_axis = ROI_ax * dz;


%  Define disp matrix
u = zeros( length(ROI_ax), length(ROI_lat), length(ROI_fr) );


%  Compute disp
A = c0/(4*pi*fc);  % (calc once to reduce multiplications)
for m = ROI_ax(1) : ROI_ax(end)     % axial (range gate) index
    
    for n = ROI_fr(1) : ROI_fr(end) % temporal index
        
        %  Extract I and Q components for the N and M ranges
        I = real( IQ( m-M1:m+M2, ROI_lat, n:n+N-1 ) );
        Q = imag( IQ( m-M1:m+M2, ROI_lat, n:n+N-1 ) );
        
        %  Compute expressions within the sigmas
        uu = Q(1:M,   :, 1:N-1).*I(1:M, :, 2:N) - I(1:M  , :, 1:N-1).*Q(1:M, :, 2:N);
        ud = I(1:M,   :, 1:N-1).*I(1:M, :, 2:N) + Q(1:M  , :, 1:N-1).*Q(1:M, :, 2:N);
        du = Q(1:M-1, :, 1:N  ).*I(2:M, :, 1:N) - I(1:M-1, :, 1:N  ).*Q(2:M, :, 1:N);
        dd = I(1:M-1, :, 1:N  ).*I(2:M, :, 1:N) + Q(1:M-1, :, 1:N  ).*Q(2:M, :, 1:N);
        
        %  Compute the sigmas
        uu = sum(sum(uu,1),3);
        ud = sum(sum(ud,1),3);
        du = sum(sum(du,1),3);
        dd = sum(sum(dd,1),3);
        
        %  Pinton (2006) (mean disp)
        num = atan2(uu,ud);
        den = 1 + atan2(du,dd) / (2*pi);
		u( m-ax_s+1, :, n ) = A * num./den ;
    end
end



