%% simulate the forward imaging process of Fourier ptychography
%% simulate the high resoluiton complex object 
objectAmplitude = double(imread('cameraman.tif'));
phase = double(imread('westconcordorthophoto.png'));    %用一张图表示相位，说明相位可恢复
phase = pi*imresize(phase,[256 256])./max(max(phase));
object = objectAmplitude.*exp(1i.* phase);
imshow(abs(object),[]);title('Input complex object');

%% setup the parameters for the coherent imaging system
waveLength = 0.63e-6;
k0 = 2*pi/waveLength;
spsize = 2.75e-6; % sampling pixel size of the CCD
psize = spsize / 4; % final pixel size of the reconstruction
NA = 0.08;

%% create the wave vectors for the LED illumiantion 
arraysize = 15;
xlocation = zeros(1,arraysize^2);
ylocation = zeros(1,arraysize^2);
LEDgap = 4; % 4mm between adjacent LEDs
LEDheight = 90; % distance bewteen the LED matrix and the sample

for i=1:arraysize % from top left to bottom right
    xlocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
    ylocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = ((arraysize-1)/2-(i-1))*LEDgap;
end;
kx_relative = -sin(atan(xlocation/LEDheight));  
ky_relative = -sin(atan(ylocation/LEDheight)); 

%% generate the low-pass filtered images==低通滤波器是光瞳函数的振幅
[m,n] = size(object); % image size of high resolution object
m1 = m/(spsize/psize);n1 = n/(spsize/psize); % image size of the final output 
imSeqLowRes = zeros(m1, n1, arraysize^2); % the final low resolution image sequence
kx = k0 * kx_relative;
ky = k0 * ky_relative;
dkx = 2*pi/(psize*n);
dky = 2*pi/(psize*m);
cutoffFrequency = NA * k0;
kmax = pi/spsize;
[kxm kym] = meshgrid(-kmax:kmax/((n1-1)/2):kmax,-kmax:kmax/((n1-1)/2):kmax);
CTF = ((kxm.^2+kym.^2)<cutoffFrequency^2); % pupil function circ(kmax); no aberration
% z = 10e-6; kzm = sqrt(k0^2-kxm.^2-kym.^2);
% pupil = exp(1i.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm)));
% aberratedCTF = pupil.*CTF;
%%
close all;clc;
objectFT = fftshift(fft2(object));
for tt =1:arraysize^2
    kxc = round((n+1)/2+kx(1,tt)/dkx);
    kyc = round((m+1)/2+ky(1,tt)/dky);
    kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
    kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
    imSeqLowFT = (m1/m)^2 * objectFT(kyl:kyh,kxl:kxh).*CTF;
    imSeqLowRes(:,:,tt) = abs(ifft2(ifftshift(imSeqLowFT)));
%     A = histeq(imSeqLowRes(:,:,tt));
    imwrite( imSeqLowRes(:,:,tt),['B:\ZJU\2019年12月汇报_徐豪\basic\数据来源\','image' num2str(tt) '.jpg']);
end;
figure;imshow(imSeqLowRes(:,:,1),[]);title('test')

%% recover high resolution image
seq = gseq(arraysize);
objectRecover = ones(m,n);
objectRecoverFT = fftshift(fft2(objectRecover));
loop = 5;
for tt=1:loop
    for i3=1:arraysize^2
        i2=seq(i3);
        kxc = round((n+1)/2+kx(1,i2)/dkx);
        kyc = round((m+1)/2+ky(1,i2)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        lowResFT = (m1/m)^2 * objectRecoverFT(kyl:kyh,kxl:kxh).*CTF;
        im_lowRes = ifft2(ifftshift(lowResFT));
        im_lowRes = (m/m1)^2 * imSeqLowRes(:,:,i2).*exp(1i.*angle(im_lowRes)); %replace the amplitude
        lowResFT=fftshift(fft2(im_lowRes)).*CTF;
        objectRecoverFT(kyl:kyh,kxl:kxh)=(1-CTF).*objectRecoverFT(kyl:kyh,kxl:kxh) + lowResFT;  %update the spectrum 
%         alpha=0:pi/50:2*pi;
%         R = cutoffFrequency;
%         x=R*cos(alpha)+kxc; 
%         y=R*sin(alpha)+kyc; 
%         imshow(log(objectRecoverFT),[]);
%         hold on  %保持
%         plot(x,y,'-','LineWidth',10,'color','r') 
%         hold on  %保持
%         axis equal
     end;
end;
objectRecover=ifft2(ifftshift(objectRecoverFT));
imshow(abs(objectRecover),[]);
figure;imshow(angle(objectRecover),[]);
figure;imshow(log(objectRecoverFT),[]);





