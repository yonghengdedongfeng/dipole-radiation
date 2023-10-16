%博士论文，格林函数法计算角谱
% Jorg S.Eiamann
% C:\Users\Administrator\Desktop\程序示例及修改\斯托克斯参量\iscat\二维
%还有一个随着dipole ，位置的改变，中心涡旋也应该改变位置
clear;
clc;

%%
nm = 1e-9;
lamda = 533*nm;
%the focusing length of lens
f1 = 5e6*nm;
f2 = 5e6*nm;
%the position of dipole
xj = 40*nm;
yj = 0*nm;
zj = 0*nm;
ptj = [xj,yj,zj];
k0 = 2*pi/lamda;
%真空介电常数
eps = 8.854e-12;
n1 = 1;
n2 = 1;
%%
%对于自由空间
k1 = n1*k0;
k2 = n2*k0;
NA1 = 1;
%最大横向波矢分量和波矢比值
kxy1 = NA1/n1;
kx = linspace(-kxy1*k1,kxy1*k1,2001);
ky = linspace(-kxy1*k1,kxy1*k1,2001);
[KX,KY] = meshgrid(kx,ky);
K1 = ones(2001)*k1;
K2 = ones(2001)*k2;
D1 = K1.^2.*exp(1j*K1*f1)/(4*pi*eps*f1);

KZ1 = sqrt(K1.^2-KX.^2-KY.^2);
KZ2 = sqrt(K2.^2-KX.^2-KY.^2);
D2 = K1.^2.*exp(1j*K2*f1).*KZ2./KZ1/(4*pi*eps*f2);

TP = 2*n1*n2*KZ1./(n1^2*KZ2 + n2^2*KZ1);
TS = 2*KZ1./(KZ1 + KZ2);

%%
phase = exp(1j*(-KX*xj-KY*yj+KZ1*zj));
p = [1,1j,0.5j];
%下半部分
nEfp = phase.*((-1*p(1)*D1.*KX.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) + (-1*p(2)*D1.*KY.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) - p(3)*D1.*sqrt(KX.^2+KY.^2)./K1);
nEfs = phase.*((-1*p(1)*D1.*KY./sqrt(KX.^2+KY.^2)) + p(2)*D1.*KX./sqrt(KX.^2+KY.^2));

%将波矢周围赋零
% for i=1:2001
%     for j=1:2001
%         if KX(i,j)^2+KY(i,j)^2 > 1.01*k1^2
%             nEfp(i,j) = 0;
%             nEfs(i,j) = 0;
%         end
%     end
% end

sinphi = KY./sqrt(KX.^2+KY.^2);
cosphi = KX./sqrt(KX.^2+KY.^2);
nEx = -1*nEfs.*sinphi - nEfp.*cosphi;
nEy =    nEfs.*cosphi - nEfp.*sinphi;
nElp = (nEx - 1j*nEy)/sqrt(2);
nErp = (nEx + 1j*nEy)/sqrt(2);
nI = abs(nEfs).^2 + abs(nEfp).^2;
nIrp =abs(nErp).^2;
nIlp =abs(nElp).^2;



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
%imagesc(kx/k0,ky/k0,angle(nElp));title('phase(left)')
imagesc(kx/k0,ky/k0,angle(nErp));title('phase(right)')
colorbar
axis xy
hold on
k01 = ones(1,100)*k0;
x=linspace(-k0,k0,100);
plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
hold off

%%
%Z>0,上半部分，可以添加基底

p = [1,1j,0.5j];

zint = 40*nm;
pphase =exp(1j*(-1*KX*xj-KY*yj-KZ1*(zj-zint)-KZ2*zint));

pEfp = pphase.*(((p(1)*D2.*KX.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) + (p(2)*D2.*KY.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) - p(3)*D2.*sqrt(KX.^2+KY.^2)./K1));
pEfs = pphase.*(((-1*p(1)*D2.*KY./sqrt(KX.^2+KY.^2)) + p(2)*D2.*KX./sqrt(KX.^2+KY.^2)));
% pEfp = phase.*((p(1)*D1.*KX.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) + (p(2)*D1.*KY.*KZ1)./(sqrt(KX.^2+KY.^2).*K1) - p(3)*D1.*sqrt(KX.^2+KY.^2)./K1);
% pEfs = phase.*((-1*p(1)*D1.*KY./sqrt(KX.^2+KY.^2)) + p(2)*D1.*KX./sqrt(KX.^2+KY.^2));

sinphi = KY./sqrt(KX.^2+KY.^2);
cosphi = KX./sqrt(KX.^2+KY.^2);
pEx = -1*pEfs.*sinphi + pEfp.*cosphi;
pEy =    pEfs.*cosphi + pEfp.*sinphi;
pElp = (pEx + 1j*pEy)/sqrt(2);
pErp = (pEx - 1j*pEy)/sqrt(2);
pI = abs(pEfs).^2 + abs(pEfp).^2;
pIrp =abs(pErp).^2;
pIlp =abs(pElp).^2;


nf = figure(3);
nf.Position(1:2) = [1100 200];
nf.Position(3:4) = [920 800];

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
%imagesc(kx/k0,ky/k0,angle(pElp));title('phase(left)')
imagesc(kx/k0,ky/k0,angle(pErp));title('phase(right)')
colorbar
axis xy
hold on
k01 = ones(1,100)*k0;
x=linspace(-k0,k0,100);
plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
hold off











%%
phase = exp(1j*(-KX*40*nm-KY*40*nm+KZ1*zj));
imagesc(kx/k0,ky/k0,angle(phase));
colorbar

phase = exp(1j*(KZ1*40*nm));
imagesc(kx/k0,ky/k0,angle(phase));
colorbar

phase = exp(1j*(-KX*40*nm));
imagesc(kx/k0,ky/k0,angle(phase));
colorbar

phase = exp(1j*(-KY*40*nm));
imagesc(kx/k0,ky/k0,angle(phase));
colorbar

imagesc(kx/k0,ky/k0,angle(pphase));
colorbar

imagesc(kx/k0,ky/k0,real(D2));
colorbar

image(real(KX))
colorbar


% nexttile
% sphere(16)
% shading interp
% title('Interpolated Shading')