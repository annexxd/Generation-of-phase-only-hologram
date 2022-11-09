function [RECON_image, R_HOLO, I_HOLO] = TwoDepth_RECONS(IMG,select,depth,t)
%select=1 for RNA, select=2 for down_sampling
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

%输入三个参数时，t为默认值6；输入四个参数时，t为用户定义的值
if nargin==3
    IMG_down=grid_cross(IMG,6);
elseif nargin==4
    IMG_down=grid_cross(IMG,t);
end

IMG_RNA=RNA(IMG);

if select==0
    I(sy:ey,sx:ex)=IMG;%选择原图
elseif select==1
    I(sy:ey,sx:ex)=IMG_RNA;%选择RNA处理
elseif select==2
    I(sy:ey,sx:ex)=IMG_down; %选择下采样处理
end

%Generate a digital complex-valued Fresnel hologram
H_Up=Holo_gen(I,z1,psize,lamda);%这里是E步得two-depth Hologram的过程
H_Down=Holo_gen(I,z2,psize,lamda);
H_Double=zeros(H_DIM,H_DIM);
H_Double(1:H_DIM/2,:)=H_Up(1:H_DIM/2,:);
H_Double((H_DIM/2+1):H_DIM,:)=H_Down((H_DIM/2+1):H_DIM,:);%分成上下两个部分的H，得two_depth H

R_HOLO=real(H_Double);%处理two-depth Hologram图片时的实部和虚部
I_HOLO=imag(H_Double);
R_HOLO=mat2gray(R_HOLO(sy:ey,sx:ex));%仅展示中间512*512的部分
I_HOLO=mat2gray(I_HOLO(sy:ey,sx:ex));


%Reconstruct hologram to image, extract magnitude component
if depth==1
    RECON_image=abs(Holo_recon(H_Double,z1,psize,lamda));%E步骤展示z1深度的重建图
elseif depth==2
    RECON_image=abs(Holo_recon(H_Double,z2,psize,lamda));%E步骤展示z2深度的重建图
end

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