%% input high-resolution image
objectIntensity = double(imread('cameraman.tif'));
figure;imshow(abs(objectIntensity),[]);
%% setup the incoherent imaging system 
waveLength = 0.63e-6;
psize = 0.3e-6; % sampling pixel size 
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
%% define the illumination patterns
NA_ill=0.14; % illumination NA for the sin patterns
sita=(0:2*pi/3:2*pi-2*pi/3);% three illumination patterns with three orientations 
n=size(objectIntensity,1);
for i=1:3
      kx_relative(i) = NA_ill.*sin(sita(i));
      ky_relative(i) = NA_ill.*cos(sita(i)); % create kx, ky wavevectors 
      kx = k0*kx_relative(i);
      ky = k0*ky_relative(i);
      x = -n/2*psize:psize:(n/2-1)*psize;
      y = -n/2*psize:psize:(n/2-1)*psize;
      [X, Y]=meshgrid(x,y);
      patternSeq(:,:,i)=(1.1+1.*sin(kx.*X+ky.*Y));% change the background and modulation here
end 
patternSeq(:,:,4) = 2.2-patternSeq(:,:,1);% pattern 4 is the complementary pattern of pattern 1
%% forward imaging model
imSeqLowRes = zeros(m,n,4);
for i=1:4
    lowResFT = OTF.*fftshift(fft2(objectIntensity.*patternSeq(:,:,i)));
    imSeqLowRes(:,:,i) = abs(ifft2(ifftshift(lowResFT)));     
end 
imshow(imSeqLowRes(:,:,1),[])
%% pattern-illuminated FP recovery
objectRecover = imSeqLowRes(:,:,1)+ imSeqLowRes(:,:,4); % initial guess
figure;imshow(objectRecover,[]);title('diffraction-limited image');
objectRecoverFT = fftshift(fft2(objectRecover));
for loop=1:50
    for i=1:4
    object_pattern = objectRecover.*patternSeq(:,:,i);
    object_pattern2 = object_pattern;
    object_patternFT = fftshift(fft2(object_pattern));
    lowResFT = OTF.*object_patternFT;
    im_lowRes = ifft2(ifftshift(lowResFT));
    object_pattern = object_pattern + deconvwnr(imSeqLowRes(:,:,i)-im_lowRes,ipsf,0.01);    
    objectRecover = objectRecover + (patternSeq(:,:,i)).*(object_pattern-object_pattern2)./((max(max(patternSeq(:,:,i)))).^2);% update the high resolution image in spatial 
    end
end
figure;imshow(abs(objectRecover),[]);title('Recovered high resolution image');
