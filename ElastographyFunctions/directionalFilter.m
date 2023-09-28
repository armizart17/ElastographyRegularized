function [sono_filt_mov,filt_mov] = directionalFilter(sono_filt,FR,pitch,f,df,sd,mov)
index_p=0;
[f1,f2] = freqspace([size(sono_filt,2) size(sono_filt,3)],'meshgrid');

[~,a]=min((abs(f1(1,:))-df/(FR)*2).^2);
[~,b]=min((abs(f1(1,:))-df/(FR+1)*2).^2);

if (b-a)>0 && a+index_p<=b
    index=a+index_p;
elseif index_p<0
    index=a+index_p;
else
    index=a;
end

temp_mask=(abs(f1)==abs(f1(1,index)));
        
spat_mask=(abs(f2)<=(f/sd(1)*pitch*2 + 1/size(sono_filt,2)*2))...
            &(abs(f2)>=(f/sd(2)*pitch*2 - 1/size(sono_filt,2)*2))...
        &(abs(f2)>0);
mask=and(temp_mask,spat_mask);
L=bwlabel(mask);
if strcmp(mov,'right')
    L(L==1)=0;
    L(L==4)=0;
    %L(L==3)=0; % JEJE
elseif strcmp(mov,'left')
    L(L==2)=0;
    L(L==3)=0;
    %L(L==4)=0; % JEJE
else
    disp('Incorrect direction');
end
mask=L>0;
win = fspecial('gaussian',[size(sono_filt,2) round(size(sono_filt,3)-1)],0.5);
filt_mov=conv2(double(mask),win,'same'); 
filt_mov = filt_mov./max(filt_mov(:)); 
sono_filt_mov=zeros(size(sono_filt));

filt_mov_shift = ifftshift(filt_mov);
for i=1:size(sono_filt,1)
    sono_filt_mov(i,:,:)=real(ifft2(fft2(squeeze(sono_filt(i,:,:))).*filt_mov_shift));
    %sono_filt_mov(i,:,:)=ifft2(ifftshift(fftshift(fft2(squeeze(sono_filt(i,:,:)))).*filt_mov));
end

end