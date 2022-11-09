function [RECON_image, R_HOLO, I_HOLO] = OSPR_RECONS(IMG,N,depth)
%N=OSPR方法使用的图片张数
%depth=1 for z1,depth=2 for z2

IMG=im2double(IMG);
IMG=IMG/max(IMG(:));    %Normalize input image to range [0, 1]
lamda=520e-9;   %Wavelength of light =520nm
psize=6.4e-6;   %Pixel size=6.4um
z1=0.16;        %Distance between image and hologram plane
z2=0.18;        %The second depth
H_DIM=1024;     %Size of hologram
I=zeros(H_DIM,H_DIM);   %A canvas of similar size as the hologram

[n m]=size(IMG(:,:));   %Size of source image
%Copy source image in the center of a larger canvas that is identical
%in size as the hologram
sx=(H_DIM-m)/2;ex=(H_DIM+m)/2-1;
sy=(H_DIM-n)/2;ey=(H_DIM+n)/2-1;

%进行OSPR方法的Hologram处理:对多个重建后的图片取平均
RECON_image=zeros(H_DIM,H_DIM);
if depth==1
    z=z1;
elseif depth==2
    z=z2;
end

for i=1:N
    IMG_RNA1=RNA(IMG);
    I(sy:ey,sx:ex)=IMG_RNA1;
    H_Up1=Holo_gen(I,z1,psize,lamda);
    H_Down1=Holo_gen(I,z2,psize,lamda);
    H_Double1=zeros(H_DIM,H_DIM);
    H_Double1(1:H_DIM/2,:)=H_Up1(1:H_DIM/2,:);
    H_Double1((H_DIM/2+1):H_DIM,:)=H_Down1((H_DIM/2+1):H_DIM,:);
    RECON_image1=abs(Holo_recon(H_Double1,z,psize,lamda));%深度z取决于输入depth
    RECON_image=RECON_image+RECON_image1;  
end
RECON_image=RECON_image/N;%F步骤的重建图

%SHOW real/imag components and recons image
R_HOLO=real(H_Double1);
I_HOLO=imag(H_Double1);%处理OSPR 时最后一个Hologram图片的实部和虚部
R_HOLO=mat2gray(R_HOLO(sy:ey,sx:ex));
I_HOLO=mat2gray(I_HOLO(sy:ey,sx:ex));%仅展示中间512*512的部分

RECON_image=mat2gray(RECON_image(sy:ey,sx:ex));%这是需要展示的重建图
end

%Generating a digital Fresnel hologram
function H=Holo_gen(I,z,psize,lamda)
[N M]=size(I(:,:));
wn=2*pi/lamda;      %Wave number

%Generate Fresnel Zone Plate
z2=z*z; %square of the depth value z
[kx ky]=meshgrid(-M/2:M/2-1,-N/2:N/2-1);
kx=kx*psize;    %Physical horizontal position
ky=ky*psize;    %Physical vertical position
fzp=exp(j*wn*sqrt(kx.^2+ky.^2+z2));

%Compute hologram of source image
A=fft2(I);  %Fourier transform of the hologram 
B=fft2(fzp);  %Fourier transform of the free space impulse response
H=fftshift(ifft2(A.*B));

ANG=angle(H);
H=cos(ANG)+j*sin(ANG);%Extract phase component(which means magnitude=1 for all the points)

end

%Numerical reconstruction of hologram 
function IR=Holo_recon(H,z,psize,lamda)
[N M]=size(H(:,:)); %Size of hologram

wn=2*pi/lamda;      %Wave number

%Generate Fresnel Zone Plate
z2=z*z; %square of the depth value z
[kx ky]=meshgrid(-M/2:M/2-1,-N/2:N/2-1);
kx=kx*psize;    %Physical horizontal position
ky=ky*psize;    %Physical vertical position
r=sqrt(kx.^2+ky.^2+z2);
fzp_c=exp(-j*wn*r)./r; %Conjugate of FZP
%fzp_c=exp(-j*wn*r);

A=fft2(H);  %Fourier transform of the hologram 
B=fft2(fzp_c);  %Fourier transform of the free space impulse response
%Compute reconstructed image 
IR=fftshift(ifft2(A.*B));
end

function I_down=grid_cross(Img,t)
[X,Y]=size(Img(:,:));
S=zeros(X,Y);
for i=1:X
    for j=1:Y
        a=mod(i,t);
        b=mod(j,t);
        if((a==0)||(b==0)||(a==b)||(a+b==t))
            S(i,j)=1;
        end
    end
end
I_down=Img.*S;
end


function IN=RNA(img)
[X,Y]=size(img(:,:));
R=rand(X,Y);
N=exp(j*2*pi*R);
IN=img.*N;
end