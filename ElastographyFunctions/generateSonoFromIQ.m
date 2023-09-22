function [sono] = generateSonoFromIQ(IQ,PRF)
%%%%%%%%%%%%%%%%%% Process IQ data using autocorrelation-based %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% spectral variance estimation %%%%%%%%%%%%%%%%%%%%%%%%%%% 

[~,~,Nens,Nframes] = size(IQ);
N_wf = Nens-1;           % number of samples in wallfiltered signal packet
wf = [1 0];
sono=zeros(size(IQ,1),size(IQ,2),size(IQ,4));
for frame = 1:Nframes
    s2var = squeeze(IQ(:,:,:,frame));
    s_wfvar=wf(1)*s2var(:,:,2:Nens)+wf(2)*s2var(:,:,1:Nens-1);
    R0var=single(2*N_wf)\squeeze(sum(abs(s_wfvar(:,:,1:N_wf-1)).^2+abs(s_wfvar(:,:,N_wf)).^2,3));
    RTvar=single(N_wf)\squeeze(sum(s_wfvar(:,:,1:N_wf-1).*conj(s_wfvar(:,:,N_wf)),3));
    R0var=R0var+eps;
    sono(:,:,frame)=(2*PRF^2)*(1-(abs(RTvar)./R0var));
end

end