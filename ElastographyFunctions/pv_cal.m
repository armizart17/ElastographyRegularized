function [u,dinf] = pv_cal(IQ,dinf,nframes) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particle velocity estimates using Loupas' algorithm for a chose region of
% interest (ROI)
% 
% Inputs:  IQData - IQ volume data obtained with Verasonics scan
%          dinf   - struct data containing spatial resolution
%          PData  - struct data containing ROI size
%           cmpd  - struct data that has # angles for beamforming
% Outputs:      u - 3D matrix of particle velocity estimates
%          dinf   - struct data with additional columns and rows info from
%                    particle velocity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Obtain IQ data
% IQ = squeeze(cell2mat(IQData(2)));	%  size(IQ)
% RF = squeeze(cell2mat(RcvData(2)));


% % =========================================================================
% %%%  Crop IQData to define ROI for Cs map
% % =======================================
% %  Define axial ROI
% dep1 = dinf.z(1);
% dep2 = dinf.z(end);
% row1 = round(dep1/dinf.dz);
% row2 = round(dep2/dinf.dz);
% dinf.Bscan_dep  = [dep1 dep2];
% dinf.Bscan_rows = [row1 row2];
% 
% %  Define lateral ROI
% wid1 = 0e-2;
% wid2 = size(IQ,2)*dinf.dx;
% col1 = round(wid1/dinf.dx)+1;
% col2 = round(wid2/dinf.dx);
% dinf.Bscan_wid  = [wid1 wid2];
% dinf.Bscan_cols = [col1 col2];
	

%  Crop IQData
% if ndims(IQ)==4
% 	IQ = IQ( row1:row2, col1:col2, :, : );
% else
% 	IQ = IQ( row1:row2, col1:col2, : );
% end

% =========================================================================
%%%  compound_angle_IQ
% =========================================================================
if nframes > 2
	
	%  If more than one IQ data set
	if ndims(IQ)==4
		tmp = [];
		for frm = 1 : size(IQ,4)
			tmp(:,:,:,frm) = A1_compound_angle_IQ ( squeeze(IQ(:,:,:,frm)), nframes );
		end
		IQ = tmp;
		clear tmp
		
	%  If only one IQ data set
	else
		IQ = A1_compound_angle_IQ ( IQ, nframes );
	end
end


% =========================================================================
%%%  Displacement estimation
% ==========================
disp('Calculating particle velocity signals...')
M  = round(1.0e-3/(dinf.dz));  % range gate length (ax) 0.5mm
N  = 2;   % ensemble length (fr) N=2 when i used it first
dW = 2;   % window shift (axial) dW = 1 when i used it first

%  If more than one IQ data set
if ndims(IQ)==4
	U = [];
	for frm = 1 : size(IQ,4)
		disp(['IQ frame ' num2str(frm) ' out of ' num2str(size(IQ,4))])
		[ U(:,:,:,frm), ~, ROI ] = A2_loupas( IQ(:,:,:,frm), M, N, dW );
	end
	u = sum(U,4)/size(U,4);		% size(u)
	clear U
	
%  If only one IQ data set
else
	[ u, ~, ROI ] = A2_loupas( IQ, M, N, dW );
    %[ u, disp_axis, ROI ] = A2_kasai( IQ, N );
end

% u = u*dinf.PRFe;
% dinf.u_rows = [dinf.Bscan_rows(1) + (ROI.ROI_ax(1)-1) ...
% 	           dinf.Bscan_rows(1) + (ROI.ROI_ax(1)-1) + (numel(ROI.ROI_ax)-1)];           
% dinf.u_cols = [dinf.Bscan_cols(1) + (ROI.ROI_lat(1)-1) ...
% 	           dinf.Bscan_cols(1) + (ROI.ROI_lat(1)-1) + (numel(ROI.ROI_lat)-1)];
           
end