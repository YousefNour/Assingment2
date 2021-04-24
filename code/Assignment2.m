%Author: Yousef Nour
%Elec 4607-Modelling of integrated Device
%Assignmetn 2 

set(0,'DefaultFigureWindowStyle', 'docked')
set(0,'defaultaxesfontsize', 20)
set(0, 'defaultaxesfontname', 'Times New Roman')
set(0,'DefaultLineLineWidth', 2);

close all
clear all 
clc

%% 1 Use the Finite Difference Method (I strongly suggest you use the matrix form
...of the problem GV = F) to solve for the electrostatic potential in the rectangular
...region L * W. 

% a) Solve the simple case where V = V0 at x = 0 and V = 0 at x = L. Note that in
...this case the top/bottom boundary conditions are not fixed. assume dV/dy = 0 for the BC

% grid size 
numX = 100;
numY = 100;
L = 3;
W = 2;

% Initialize a nonzero G matrix using Sparse function
sp = numX*numY;
G = sparse(sp,sp);  

%Initialize a zeros F matrix 
F = zeros(sp,1);
V0 = 1;

% g-matrix scanning
for ii = 1 : numY %rows 
    for jj = 1 : numX % columns 
        nMap = (ii-1)*numX + jj; % mapping of crrect values in G-matrix   
        
        numG1 = (ii-2)*numX + jj; numG2 = (jj-1) + (ii-1)*numX;
        numGy = (ii*numX) + jj; numGx = (jj+1) + (ii-1)*numX;
        
        if ii == 1
            G(nMap,nMap) = 1;
        elseif ii == numY
            G(nMap,nMap) = 1;
        elseif jj == 1
            G(nMap,nMap) = 1;
        elseif jj == numX 
            G(nMap,nMap) = 1;   
        else
            G(nMap,nMap) = -4;
            G(nMap,numG1) = 1; G(nMap,numG2) = 1;
            G(nMap,numGx) = 1; G(nMap,numGy) = 1; 
        end
        
        if ii < numX && ii > (numY-(numY-1))
            G(ii,ii + numY) = -1;
            G(ii + (numX - 1)*numY,ii + (numX - 2)*numY) = -1;   
        end
    end
end

for kk = 1: numX
    F(1+(kk-1)*numX,1) = V0;
end

V = G\F;

for ii = 1:numY
    for jj = 1:numX
        nMap = jj + (ii - 1)*numX;
        voltMap(ii,jj) = V(nMap);
    end
end

Dx = linspace(0,W,numY);
Dy = linspace(0,L,numX);

figure(1)
surf(Dx,Dy,voltMap)
title('Voltage Plot, 101046991')
xlabel('X position')
ylabel('Y position')
view(0,270)

%1B) Solve the case where V = V0 at x = 0, x = L and V = 0 at y = 0, y = W.
...Compare the solution of a bunch of mesh sizes to the analytical series solution
clear 
clc

% grid size 
numX = 100;
numY = 100;
L = 3;
W = 2;

% Initialize a nonzero G matrix using Sparse function
sp = numX*numY;
G = sparse(sp);

%Initialize a zeros F matrix 
F = zeros(sp,1);
V0 = 1;

% g-matrix scanning
for ii = 1 : numX %rows 
    for jj = 1 : numY % columns 
        nMap = (ii*numY-numY) + jj; % mapping of crrect values in G-matrix   
        
        numG1 = (ii*numY-2*numY) + jj; numG2 = (jj-1) + (ii*numY-numY);
        numGy = (ii*numY-numY) + (jj+1); numGx = jj + (ii*numY);
        
        if ii == 1
            G(nMap,nMap) = 1;
        elseif ii == numY
            G(nMap,nMap) = 1;
        elseif jj == 1
            G(nMap,nMap) = 1;
        elseif jj == numX 
            G(nMap,nMap) = 1;   
        else
            G(nMap,nMap) = -4;
            G(nMap,numG1) = 1; G(nMap,numG2) = 1;
            G(nMap,numGx) = 1; G(nMap,numGy) = 1; 
        end
    end
end

for kk = 2:(numY-1)
    F(1+(kk*numY-numY),1) = V0;
    F(numX+(kk*numY-numY),1) = V0;
end

V = G\F;

for ii = 1:numX
    for jj = 1:numY
        nMap = jj + (ii*numY-numY);
        voltMap(ii,jj) = V(nMap);
    end
end

Dx = linspace(0,W,numY);
Dy = linspace(0,L,numX);

figure(2)
surf(Dx,Dy,voltMap)
title('Voltage Plot, 101046991')
xlabel('X position')
ylabel('Y position')
zlabel('Voltage')

%Compare the solution of a bunch of mesh sizes to the analytical series
AnalyticalSer = zeros(numX,numY);

aa = linspace(-L/2,L/2,numX);
bb = linspace(0,W,numY);
[V,y] = meshgrid(aa, bb);  

for iter = 1:100
    %n = 1,3,5..., where n must be odd.
    n = 2*iter - 1;
    
   % Analytical Solution:
    eqution_a = W;
    eqution_b = L/2;
    const = (4*V0/pi);
    num = cosh((n*pi).*V./eqution_a);
    den = cosh((n*pi).*eqution_b./eqution_a);
    multip = sin((n*pi).*y./eqution_a);
    AnalyticalSer = AnalyticalSer + const .* (n^-1) .* (num./ den) .* multip;
     
     pause(0.01)
end

%3D plot of the Analytical Solution 
figure(3)
surf(Dx,Dy,AnalyticalSer)
title({'Anaytical Series Solution,101046991'})
xlabel('X')
ylabel('Y')
zlabel('Z')

%% 2 Use the Finite Diference Method to solve for the current flow in the rectangular region
...LXW shown in the given Figure.

% close all
% clear all
% clc
% 
% numX = 100;
% numY = 100*(3/2);
% 
% % Initialize a nonzero G matrix using Sparse function
% sp = numX*numY;
% G = sparse(sp);  
% 
% %Initialize a zeros F matrix 
% F = zeros(1,sp);
% V0 = 1;
% % sigma matrix and sigma intialization
% SigMatrx = zeros (numX, numY)
% sigma1 = 1;
% sigma2 = 10e-2
% 
% BoxDim = [(numX*2/5) (numX*3/5) (numY*2/5) (numY*3/5)];
% 
% for ii = 1: numX 
%     for jj = 1: numY
%         if ii > BoxDim(1)
%             SigMatrx(ii, jj) = sigma2;
%         elseif ii < BoxDim(2)
%             SigMatrx(ii, jj) = sigma2;
%         elseif jj < BoxDim(3)
%             SigMatrx(ii, jj) = sigma2;    
%         elseif jj > BoxDim(4)
%             SigMatrx(ii, jj) = sigma2;    
%         else
%             SigMatrx(ii, jj) = sigma1;
%         end
%     end
% end

