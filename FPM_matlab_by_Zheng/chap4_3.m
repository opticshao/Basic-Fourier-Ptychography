%% input high-resolution color image
imsize=256;
imtemp=double(imread('baboon.jpg'));
objectIntensity_r = imresize(imtemp(:,:,1),[imsize imsize]);
objectIntensity_g = imresize(imtemp(:,:,2),[imsize imsize]);
objectIntensity_b = imresize(imtemp(:,:,3),[imsize imsize]);
%% setup the incoherent imaging system 
waveLength_r = 0.63e-6;
waveLength_g = 0.53e-6;
waveLength_b = 0.47e-6;
psize = 0.5e-6; % sampling pixel size 
NA = 0.15;
k0_r = 2*pi/waveLength_r;
k0_g = 2*pi/waveLength_g;
k0_b = 2*pi/waveLength_b;
cutoffFrequency_r = NA * k0_r;
cutoffFrequency_g = NA * k0_g;
cutoffFrequency_b = NA * k0_b;
[m n] = size(objectIntensity_r);
kx = -pi/psize:2*pi/(psize*(n-1)):pi/psize;
ky = -pi/psize:2*pi/(psize*(n-1)):pi/psize;
[kxm kym] = meshgrid(kx,ky);
CTF_r = double(((kxm.^2+kym.^2)<cutoffFrequency_r^2)); 
CTF_g = double(((kxm.^2+kym.^2)<cutoffFrequency_g^2)); 
CTF_b = double(((kxm.^2+kym.^2)<cutoffFrequency_b^2)); 
cpsf_r = fftshift(ifft2(ifftshift(CTF_r))); % coherent PSF
cpsf_g = fftshift(ifft2(ifftshift(CTF_g))); 
cpsf_b = fftshift(ifft2(ifftshift(CTF_b))); 
ipsf_r = (abs(cpsf_r)).^2; % incoherent PSF
ipsf_g = (abs(cpsf_g)).^2;
ipsf_b = (abs(cpsf_b)).^2;
OTF_r = abs(fftshift(fft2(ifftshift(ipsf_r)))); % incoherent transfer function 
OTF_r = OTF_r./max(max(OTF_r));
OTF_g = abs(fftshift(fft2(ifftshift(ipsf_g)))); % incoherent transfer function 
OTF_g = OTF_g./max(max(OTF_g));
OTF_b = abs(fftshift(fft2(ifftshift(ipsf_b)))); % incoherent transfer function 
OTF_b = OTF_b./max(max(OTF_b));
%% define illumination pattern
pattern_r = imnoise(ones(m,n),'speckle',0.5);
pattern_g = imnoise(ones(m,n),'speckle',0.5);
pattern_b = imnoise(ones(m,n),'speckle',0.5);
patternNum = 169; % number of speckle patterns
% define the scanning path of the speckle pattern
arraySize = 15; 
xLocation = zeros(1,arraySize^2);
yLocation = zeros(1,arraySize^2);
scanStep = 2; 
for i=1:arraySize 
    xLocation(1,1+arraySize*(i-1):15+arraySize*(i-1)) = (-(arraySize-1)/2:1:(arraySize-1)/2)*scanStep;
    yLocation(1,1+arraySize*(i-1):15+arraySize*(i-1)) = ((arraySize-1)/2-(i-1))*scanStep;
end
snakeSeq = gseq(arraySize);
patternSeq_r = zeros(m,n,patternNum);
patternSeq_g = zeros(m,n,patternNum);
patternSeq_b = zeros(m,n,patternNum);
% shift one speckle pattern to generate a sequence of scanning speckle patterns
for i=1:patternNum
    patternSeq_r(:,:,i) = circshift(pattern_r,[xLocation(snakeSeq(i)) yLocation(snakeSeq(i))]);         
    patternSeq_g(:,:,i) = circshift(pattern_g,[xLocation(snakeSeq(i)) yLocation(snakeSeq(i))]);   
    patternSeq_b(:,:,i) = circshift(pattern_b,[xLocation(snakeSeq(i)) yLocation(snakeSeq(i))]);   
end 
%% forward imaging model
imSeqLowRes_r = zeros(m,n,patternNum);
imSeqLowRes_g = zeros(m,n,patternNum);
imSeqLowRes_b = zeros(m,n,patternNum);
imSeqLowRes_t = zeros(m,n,patternNum);
for i=1:patternNum
    lowResFT_r = OTF_r.*fftshift(fft2(objectIntensity_r.*patternSeq_r(:,:,i)));
    imSeqLowRes_r(:,:,i) = abs(ifft2(ifftshift(lowResFT_r)));
    lowResFT_g = OTF_g.*fftshift(fft2(objectIntensity_g.*patternSeq_g(:,:,i)));
    imSeqLowRes_g(:,:,i) = abs(ifft2(ifftshift(lowResFT_g)));
    lowResFT_b = OTF_b.*fftshift(fft2(objectIntensity_b.*patternSeq_b(:,:,i)));
    imSeqLowRes_b(:,:,i) = abs(ifft2(ifftshift(lowResFT_b)));
    imSeqLowRes_t(:,:,i) = imSeqLowRes_r(:,:,i) + imSeqLowRes_g(:,:,i) + imSeqLowRes_b(:,:,i);
end 
%% pattern-illuminated FP recovery 
objectRecover_r = sum(imSeqLowRes_t,3)/(3*patternNum); % initial guess
objectRecover_g = sum(imSeqLowRes_t,3)/(3*patternNum); 
objectRecover_b = sum(imSeqLowRes_t,3)/(3*patternNum); 
for loopnum=1:5
    for i=1:patternNum        
    object_pattern_r=objectRecover_r.*patternSeq_r(:,:,i);
    object_pattern_g=objectRecover_g.*patternSeq_g(:,:,i);
    object_pattern_b=objectRecover_b.*patternSeq_b(:,:,i);    
    object_pattern2_r=object_pattern_r;
    object_pattern2_g=object_pattern_g;
    object_pattern2_b=object_pattern_b;    
    object_patternFT_r=fftshift(fft2(object_pattern_r));
    object_patternFT_g=fftshift(fft2(object_pattern_g));
    object_patternFT_b=fftshift(fft2(object_pattern_b));    
    lowResFT1_r=OTF_r.*object_patternFT_r;
    lowResFT1_g=OTF_r.*object_patternFT_g;
    lowResFT1_b=OTF_r.*object_patternFT_b;    
    im_lowRes_r=ifft2(ifftshift(lowResFT1_r));
    im_lowRes_g=ifft2(ifftshift(lowResFT1_g));
    im_lowRes_b=ifft2(ifftshift(lowResFT1_b));        
    intensity = im_lowRes_r + im_lowRes_g + im_lowRes_b;     
    im_lowRes_r=imSeqLowRes_t(:,:,i).*im_lowRes_r./intensity.*exp(1i.*angle(im_lowRes_r));
    im_lowRes_g=imSeqLowRes_t(:,:,i).*im_lowRes_g./intensity.*exp(1i.*angle(im_lowRes_g));
    im_lowRes_b=imSeqLowRes_t(:,:,i).*im_lowRes_b./intensity.*exp(1i.*angle(im_lowRes_b));     
    lowResFT2_r=fftshift(fft2(im_lowRes_r));
    lowResFT2_g=fftshift(fft2(im_lowRes_g));
    lowResFT2_b=fftshift(fft2(im_lowRes_b));    
    object_patternFT_r=object_patternFT_r + (OTF_r)./(max(max((abs(OTF_r)).^2))).*(lowResFT2_r-lowResFT1_r);
    object_patternFT_g=object_patternFT_g + (OTF_g)./(max(max((abs(OTF_g)).^2))).*(lowResFT2_g-lowResFT1_g);
    object_patternFT_b=object_patternFT_b + (OTF_b)./(max(max((abs(OTF_b)).^2))).*(lowResFT2_b-lowResFT1_b);
    object_pattern_r=ifft2(ifftshift(object_patternFT_r));
    object_pattern_g=ifft2(ifftshift(object_patternFT_g));
    object_pattern_b=ifft2(ifftshift(object_patternFT_b));     
    objectRecover_r=objectRecover_r + patternSeq_r(:,:,i).*(object_pattern_r - object_pattern2_r)./(max(max(patternSeq_r(:,:,i)))).^2;
    objectRecover_g=objectRecover_g + patternSeq_g(:,:,i).*(object_pattern_g - object_pattern2_g)./(max(max(patternSeq_g(:,:,i)))).^2;
    objectRecover_b=objectRecover_b + patternSeq_b(:,:,i).*(object_pattern_b - object_pattern2_b)./(max(max(patternSeq_b(:,:,i)))).^2;
    end    
end
figure;imshow(abs(objectRecover_r),[]); title('Multiplexed-SI recovery (red)');
figure;imshow(abs(objectRecover_g),[]); title('Multiplexed-SI recovery (green)');
figure;imshow(abs(objectRecover_b),[]); title('Multiplexed-SI recovery (blue)');