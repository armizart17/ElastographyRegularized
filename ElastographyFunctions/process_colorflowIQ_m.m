function [sono] = process_colorflowIQ_m(IQ,Properties,wf)
%%%%%%%%%%%%%%%%%% Process IQ data using autocorrelation-based %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% spectral variance estimation %%%%%%%%%%%%%%%%%%%%%%%%%%% 
     PRF = Properties.dr;
     N = Properties.extra; % number of samples in packet
     N_wf = N-1;           % number of samples in wallfiltered signal packet
     sono=zeros(size(IQ,1),size(IQ,3),size(IQ,4));
     for frame = 1:Properties.nframes
       s2var=reshape(IQ(:,:,:,frame),[size(IQ,1) size(IQ,2) size(IQ,3)]);
       s_wfvar=wf(1)*s2var(:,2:N,:)+wf(2)*s2var(:,1:N-1,:);
       R0var=single(2*N_wf)\squeeze(sum(abs(s_wfvar(:,1:N_wf-1,:)).^2+abs(s_wfvar(:,2:N_wf,:)).^2,2));
       RTvar=single(N_wf)\squeeze(sum(s_wfvar(:,1:N_wf-1,:).*conj(s_wfvar(:,2:N_wf,:)),2));
       R0var=R0var+eps;
       sono(:,:,frame)=(2*PRF^2)*(1-(abs(RTvar)./R0var));
%          figure;imagesc(sono(:,:,frame),[0 max(sono(:))]); set(gcf,'Colormap',sonomap)
     end
%     close all;
end