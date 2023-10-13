
clear;
clc;
%%
nm = 1e-9;
lamda = 627*nm;
%高度要注意调节，越高则倏逝波越小
d = 40*nm;
k0 = 2*pi/lamda;
n = 1;
NA = n*sin(pi/2);

kx = linspace(-n*k0,n*k0,2001);
ky = linspace(-n*k0,n*k0,2001);
[KX,KY] = meshgrid(kx,ky);
K0 = ones(2001)*k0;
KZ = sqrt(K0.^2-KX.^2-KY.^2);
KZ2 = sqrt(n^2*K0.^2-KX.^2-KY.^2);
C = exp(1j*KZ*d).*sqrt(n^2*K0.^2-KX.^2-KY.^2)./KZ;
TP = 2*n*KZ./(KZ2 + n^2*KZ);
TS = 2*KZ./(KZ + KZ2);
%p = [-0.9797-0.2003j,0,0.0352-0.2234j];

for ampz = (0:0.025:0.5)*2*2
    %这是颗粒散射偶极矩
    p = [1,1j,ampz*1j]; %ampz = (-0.5:0.1:0.5)*2
    %p = [1,1j,0.2*exp(1j*ampz)]; %ampz = (-0.5:0.1:0.5)*2*pi
    %p = [0.76,0,1*1j];
    
    %上半区域
    pEfp = (p(1)*C.*KX.*KZ)./(sqrt(KX.^2+KY.^2).*K0) + (p(2)*C.*KY.*KZ)./(sqrt(KX.^2+KY.^2).*K0) - p(3)*C.*sqrt(KX.^2+KY.^2)./K0;
    pEfs = (-1*p(1)*C.*KY./sqrt(KX.^2+KY.^2)) + p(2)*C.*KX./sqrt(KX.^2+KY.^2);
    %下半区域
    nEfp = (-1*p(1)*C.*KX.*KZ)./(sqrt(KX.^2+KY.^2).*K0) + (-1*p(2)*C.*KY.*KZ)./(sqrt(KX.^2+KY.^2).*K0) - p(3)*C.*sqrt(KX.^2+KY.^2)./K0;
    nEfs = (-1*p(1)*C.*KY./sqrt(KX.^2+KY.^2)) + p(2)*C.*KX./sqrt(KX.^2+KY.^2);

    
    pI = abs(pEfs).^2 + abs(pEfp).^2;
    nI = abs(nEfs).^2 + abs(nEfp).^2;
    
    sinphi = KY./sqrt(KX.^2+KY.^2);
    cosphi = KX./sqrt(KX.^2+KY.^2);
    %上半区域
    pEx = -1*pEfs.*sinphi + pEfp.*cosphi;
    pEy =    pEfs.*cosphi + pEfp.*sinphi;
    pElp = (pEx + 1j*pEy)/sqrt(2);
    pErp = (pEx - 1j*pEy)/sqrt(2);
    pIrp =abs(pErp).^2;
    pIlp =abs(pElp).^2;
%     %下半区域%%%%%%%%%%%%%%%%%%%%有问题
%     nEx = -1*nEfs.*cosphi - nEfp.*sinphi;
%     nEy =    nEfs.*sinphi - nEfp.*cosphi;
    %下半区域%%%%%%%%%%%%%%%%%%%%有问题
    nEx = -1*nEfs.*sinphi - nEfp.*cosphi;
    nEy =    nEfs.*cosphi - nEfp.*sinphi;
    nElp = (nEx - 1j*nEy)/sqrt(2);
    nErp = (nEx + 1j*nEy)/sqrt(2);
    nIrp =abs(nErp).^2;
    nIlp =abs(nElp).^2;

    %画图
    pf = figure(1);
    pf.Position(1:2) = [100 200];
    pf.Position(3:4) = [920 800];
   
    sgtitle(num2str(ampz))
    subplot(2,2,1)
    imagesc(kx/k0,ky/k0,pI);title('z>0')
    colormap("jet")
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(2,2,2)
    imagesc(kx/k0,ky/k0,abs(pErp).^2);title('Irp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(2,2,3)
    imagesc(kx/k0,ky/k0,abs(pElp).^2);title('Ilp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(2,2,4)
    imagesc(kx/k0,ky/k0,angle(pElp));title('phase(left)')
%     imagesc(kx/k0,ky/k0,angle(pErp));title('phase(right)')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%下半空间的
    nf = figure(2);
    nf.Position(1:2) = [1100 200];
    nf.Position(3:4) = [920 800];

    subplot(2,2,1)
    imagesc(kx/k0,ky/k0,nI);title('z<0')
    colormap("jet")
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(2,2,2)
    %imagesc(kx/k0,ky/k0,log(abs(nErp).^2));title('Irp')
    imagesc(kx/k0,ky/k0,abs(nErp).^2);title('Irp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off
    
    subplot(2,2,3)
    imagesc(kx/k0,ky/k0,abs(nElp).^2);title('Ilp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(2,2,4)
%     imagesc(kx/k0,ky/k0,angle(nElp));title('phase(left)')
    imagesc(kx/k0,ky/k0,angle(nErp));title('phase(right)')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    pause(3)
end

%%
imagesc(kx/k0,ky/k0,nIrp+nIlp-nI);title('Irp')
colorbar
axis xy
imagesc(kx/k0,ky/k0,imag(KZ));title('Irp')
colorbar
axis xy

%     %Erp = 1j*Ex*sin(theta) + Ey*cos(theta);
%     %Erp = 1j*Efp*sin(theta) + Efs*cos(theta);
%     %Erp = Ex*cos(theta) + Ey*sin(theta);
