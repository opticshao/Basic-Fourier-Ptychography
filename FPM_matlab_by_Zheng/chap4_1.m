%% input high-resolution image
objectIntensity = double(imread('cameraman.tif'));
figure;imshow(abs(objectIntensity),[]);
%% setup the incoherent imaging system 
waveLength = 0.63e-6;
psize = 0.5e-6; % sampling pixel size 
NA = 0.15;
k0 = 2*pi/waveLength;
cutoffFrequency = NA*k0;
[m n] = size(objectIntensity);
kx = -pi/psize:2*pi/(psize*(n-1)):pi/psize;
ky = -pi/psize:2*pi/(psize*(n-1)):pi/psize;
[kxm kym] = meshgrid(kx,ky);
CTF = double(((kxm.^2+kym.^2)<cutoffFrequency^2)); % the coherent transfer function
cpsf = fftshift(ifft2(ifftshift(CTF))); % coherent PSF
ipsf = (abs(cpsf)).^2; % incoherent PSF
OTF = abs(fftshift(fft2(ifftshift(ipsf)))); % incoherent transfer function 
OTF = OTF./max(max(OTF));
figure;imshow(abs(OTF),[]);title('Incoherent transfer function in the Fourier domain');
%% define illumination pattern
pattern = imnoise(ones(m,n),'speckle',0.5);
patternNum = 49; % number of speckle patterns
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
patternSeq = zeros(m,n,patternNum);
% shift one speckle pattern to generate a sequence of scanning speckle patterns
for i=1:patternNum
    patternSeq(:,:,i) = circshift(pattern,[xLocation(snakeSeq(i)) yLocation(snakeSeq(i))]);         
end 
%% forward imaging model
imSeqLowRes = zeros(m,n,patternNum);
for i=1:patternNum
    lowResFT = OTF.*fftshift(fft2(objectIntensity.*patternSeq(:,:,i)));
    imSeqLowRes(:,:,i) = abs(ifft2(ifftshift(lowResFT)));     
end 
%% pattern-illuminated FP recovery 
objectRecover = sum(imSeqLowRes,3)/patternNum; % initial guess
figure;imshow(objectRecover,[]);title('diffraction-limited image');
for loopnum=1:10
    for i=1:patternNum        
    object_pattern=objectRecover.*patternSeq(:,:,i);
    object_pattern2=object_pattern;
    object_patternFT=fftshift(fft2(object_pattern));
    lowResFT1=OTF.*object_patternFT;
    im_lowRes=ifft2(ifftshift(lowResFT1));
    im_lowRes=imSeqLowRes(:,:,i).*exp(1i.*angle(im_lowRes));
    lowResFT2=fftshift(fft2(im_lowRes));
    object_patternFT=object_patternFT+conj(OTF)./(max(max((abs(OTF)).^2))).*(lowResFT2-lowResFT1);
    object_pattern=ifft2(ifftshift(object_patternFT));
    objectRecover=objectRecover+patternSeq(:,:,i).*(object_pattern-object_pattern2)./(max(max(patternSeq(:,:,i)))).^2;
    end    
end
figure;imshow(abs(objectRecover),[]); title('Pattern-illuminated FP recovery');
objectRecoverFT=fftshift(fft2(objectRecover));
figure;imshow(log(abs(objectRecoverFT)),[]);