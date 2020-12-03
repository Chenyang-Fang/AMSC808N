function MakeEmojiData()
close all
N = 40;
face = zeros(N);
t = linspace(0,1,N);
[x,y] = meshgrid(t,t);
R = 0.45; % face radius
reye = 0.1; % eye radius
deye = R/2;
aleye = 5*pi/6;
scleye = [cos(aleye);sin(aleye)];
areye = pi/6;
screye = [cos(areye);sin(areye)];
c = [0.5;0.5]; % face center;
% mouth
mx1 = c(1) - R/2;
mx2 = c(1) + R/2;
my0 = c(2) - R/2;
mysin =@(a,x1,x2,y0,x)y0+a*sin(pi*(x-x1)/(x2-x1));
hm = R/5;
% eye
cleft = c + deye*scleye;
cright = c + deye*screye;
ae = pi/3;
xl1 = cleft(1) - reye*sin(ae);
xl2 = cleft(1) + reye*sin(ae);
xr1 = cright(1) - reye*sin(ae);
xr2 = cright(1) + reye*sin(ae);
y0 = cleft(2) - reye*cos(ae);
he = reye*cos(ae);
% draw face
nt = 8;
nu = 4;
t = linspace(-1,1,nt);
u = linspace(0,1,nu);
I = (1 : N)';
bb = 0.1;
figure;
for j = 1 : nu
    bu = 0;
    fu = 1;
    for i = 1 : nt
        face = bu*ones(N);
        ind = find((x-c(1)).^2 + (y-c(2)).^2 < R^2);
        face(ind) = fu;
        si = @(x)mysin(-t(i)*he,xl1,xl2,y0,x);
        ind = find((x-cleft(1)).^2 + (y-cleft(2)).^2 <= reye^2 & si(x) <= y);
        face(ind) = bu;
        si = @(x)mysin(-t(i)*he,xr1,xr2,y0,x);
        ind = find((x-cright(1)).^2 + (y-cright(2)).^2 <= reye^2 & si(x) <= y);
        face(ind) = bu;
        si = @(x)mysin(t(i)*hm,mx1,mx2,my0,x);
        if t(i) >= 0
            ind = find(si(x) >= y & y >= my0 & x >= mx1 & x <= mx2);
        else
            ind = find(si(x) <= y & y <= my0 & x >= mx1 & x <= mx2);
        end    
        face(ind) = bu;
        % smooth
        subplot(nu,nt,(j-1)*nt + i);
        for k = 1 : 2 + 4*j
            face = bb*(face(circshift(I,[-1,0]),:)+face(circshift(I,[1,0]),:)+...
                face(:,circshift(I,[-1,0]))+face(:,circshift(I,[1,0])))+(1-4*bb)*face;
        end
        imagesc(face);
        caxis([0,1]);
        axis xy
        axis off
        colormap gray
        drawnow
    end
end
    
% make data
nuold = nu;
nt = 32;
nu = 32;
t = linspace(-1,1,nt);
u = linspace(0,0.4,nu);
I = (1 : N)';
n = nt*nu;
count = 0;
X = zeros(n,N^2);
del = bb*(nuold/nu);
col = parula(1000);
c = zeros(n,3);
for j = 1 : nu
    bu = 0;
    fu = 1;
    for i = 1 : nt
        count = count + 1;
        face = bu*ones(N);
        ind = find((x-c(1)).^2 + (y-c(2)).^2 < R^2);
        face(ind) = fu;
        si = @(x)mysin(-t(i)*he,xl1,xl2,y0,x);
        ind = find((x-cleft(1)).^2 + (y-cleft(2)).^2 <= reye^2 & si(x) <= y);
        face(ind) = bu;
        si = @(x)mysin(-t(i)*he,xr1,xr2,y0,x);
        ind = find((x-cright(1)).^2 + (y-cright(2)).^2 <= reye^2 & si(x) <= y);
        face(ind) = bu;
        si1 = @(x)mysin(t(i)*hm,mx1,mx2,my0,x);
        si2 = @(x)mysin((1-t(i)^2)*hm,mx1,mx2,my0,x);
        if t(i) >= 0
            ind = find(si1(x) >= y & y >= si2(x) & x >= mx1 & x <= mx2);
        else
            ind = find(si1(x) <= y & y <= si2(x) & x >= mx1 & x <= mx2);
        end    
        face(ind) = bu;
        % smooth
        for k = 1 : 4*(j+1)
            face = del*(face(circshift(I,[-1,0]),:)+face(circshift(I,[1,0]),:)+...
                face(:,circshift(I,[-1,0]))+face(:,circshift(I,[1,0])))+face*(1-4*del);
        end
        X(count,:) = face(:)';
        c(count,:) = col(ceil(1000*sqrt((i^2+j^2)/(nt^2+nu^2))),:);
    end
end
data3 = X;
save('FaceData.mat','data3','c');
end
    




