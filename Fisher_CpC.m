%% Test against exact solution of Fisher Eq. in 1D, CpC
clear; 
tic
Mh=15; %How many different timesteps we want 
D=1; x0=0; xf=4; tf=2; TIME=tf; %Diffusion coefficient, lenght of the space and time interval, x parameter
N=401;  alf=2.5;   % initial number of cells, coefficient of nonlin term
p0=0.5; dp=0.2; Np=3;
MaxD=zeros(Mh,20);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1);%%Errors and axis for the timestep-sizes
h=TIME/4; Dx=(xf-x0)/(N-1); %%initial timestep and Space-step
xaxis=zeros(N,1); for i=1:N xaxis(i)=x0+(i-1)*Dx; end  %x axis
M1=zeros(N,N); ee1=zeros(N,1);
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
E=D/(Dx)^2; %%matrix elements
for i=1:N 
 u0(i)=1/(1+exp(sqrt(alf/6)*xaxis(i)))^2; 
 uex(i)=1/(1+exp(sqrt(alf/6)*xaxis(i)-5/6*alf*TIME))^2;  %%Initial and exact final results
end
M1(1,1)=-E;M1(1,2)=E; M1(N,N)=-E;M1(N,N-1)=E;
for i=2:N-1     M1(i,i-1)=E; M1(i,i)=-2*E; M1(i,i+1)=E; end
Neb1=zeros(N,1);Neb2=zeros(N,1); 
for i=1:N
    if (i>1) Neb1(i)=M1(i,i-1); end  %bal szomszéd
    if (i<N) Neb2(i)= M1(i,i+1); end %jobb szomszéd
   end
%%%%%%%%%%%%
EG = eig(M1);Emax=0; E0=-1; Emin=-2; %%Eigenvalues
for i=1:N if EG(i)<Emax Emax=EG(i); end
    if EG(i)>E0 
       E0=EG(i); end 
end
for i=1:N if EG(i)>Emin
      if EG(i)<E0 
           Emin=EG(i); end 
    end
end
hMAX= -2/Emax  %% max timestep for Expl Euler
StiffRatio=Emax/Emin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
T = ceil(TIME/h-0.1)+1;
taxis=zeros(T,1); %%physical time
for t=1:T
      taxis(t)=(t-1)*h;
end
if(h>0.00001 && h<0.0001) rfun(ih)=2000000000*h*h*h; end
UC=zeros(N,1);UCx=zeros(N,1);UTEMP1=zeros(N,1); UTEMP2=zeros(N,1);  %%Final U and temporary results
a0=zeros(N,1);
Axhstep(ih)=h; %%Timestep axis
r=h*E;r2=2*h*E; ee=exp(-r2); eehalf=exp(-r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UC(:)=u0(:); %%%% 1 Stage CN
UC(1)=1/(1+exp(sqrt(alf/6)*xaxis(1)-5/6*alf*taxis(1)))^2;
UC(N)=1/(1+exp(sqrt(alf/6)*xaxis(N)-5/6*alf*taxis(1)))^2;
tic
for t=1:T-1
    for i=2:N-1  %%%Neighbours are constans
        a0(i)=(UC(i-1)+UC(i+1))/2;
        UTEMP1(i)=UC(i)*ee+ a0(i)*(1-ee)+0*alf*UC(i)*(1-UC(i))*h;
    end
    UTEMP1(1)=1/(1+exp(sqrt(alf/6)*xaxis(1)-5/6*alf*taxis(t+1)))^2;
    UTEMP1(N)=1/(1+exp(sqrt(alf/6)*xaxis(N)-5/6*alf*taxis(t+1)))^2;
 UC(:)=UTEMP1(:);
 for i=2:N-1 
     UC(i)=(1+alf*h)*UC(i)/(1+alf*h*UC(i));
    end
end
toc
%%Constant errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UC(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,1)=MaxKul; 
fprintf('Const h: %f \nMaxK: %f \n',h,MaxKul)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ip=0:Np-1  %% Big loop for the parameter
p=p0+ip*dp;eep=exp(-p*r2);
UCx=u0(:); %%%% CpC method
tic
for t=1:T-1
    for i=2:N-1  %%%first stage, neighbours are constans
        a0(i)=(UCx(i-1)+UCx(i+1))/2;
        UTEMP1(i)=UCx(i)*eep+ a0(i)*(1-eep);
    end
    UTEMP1(1)=1/(1+exp(sqrt(alf/6)*xaxis(1)-5/6*alf*taxis(t+1)))^2;
    UTEMP1(N)=1/(1+exp(sqrt(alf/6)*xaxis(N)-5/6*alf*taxis(t+1)))^2;
 UTEMP1(:)=(1-1/2/p)*UCx(:)+1/2/p*UTEMP1(:);
  for i=2:N-1 %%Second stage
     a0(i)=(UTEMP1(i-1)+UTEMP1(i+1))/2;
     UCx(i)=UCx(i)*ee+ a0(i)*(1-ee); 
     UCx(i)=(1+alf*h)*UCx(i)/(1+alf*h*UCx(i));
  end
    UCx(1)=1/(1+exp(sqrt(alf/6)*xaxis(1)-5/6*alf*taxis(t+1)))^2;
    UCx(N)=1/(1+exp(sqrt(alf/6)*xaxis(N)-5/6*alf*taxis(t+1)))^2; 
end
toc
%% A2 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UCx(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,ip+2)=MaxKul;
plot(xaxis,u0,xaxis,uex,xaxis,UCx,'--o','LineWidth',2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=h/2;

end

    orange        = [1 0.5 0];
    dark_red      = [0.75, 0.0780, 0.1840]; 
    dark_yellow   = [0.9290, 0.6940, 0.1250]; 
    dark_green    = [0, 0.5, 0];
    dark_purple   = [0.55, 0.0840, 0.5560];
    dark_blue     = [0, 0.4470, 0.7410];
 
	%% figure 1:
  figure('Name', 'Errors as a fuction of time step h');   

  plot(Axhstep(:),MaxD(:,1), '-v', 'Color', orange, 'LineWidth', 1.3); %% cn 
	hold on;
    	plot(Axhstep(:),MaxD(:,2), '--o',  'Color', dark_purple, 'LineWidth', 1.2); %% cncn p=1/3
	hold on;
	plot(Axhstep(:),MaxD(:,3), ':*',  'Color', dark_blue, 'LineWidth', 2.6); %% cncn p=1/2
	hold on;
	plot(Axhstep(:),MaxD(:,4), '-x',  'Color', dark_yellow, 'LineWidth', 2.4); %% cncn p=2/3
	hold on;
 	plot(Axhstep(:),MaxD(:,5), '-.s',  'Color', dark_green, 'LineWidth', 1); %% cncn p=3/4
	hold on;
	hold off;
	
XLabel = xlabel('Time step size h');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 16); 
set(YLabel, 'FontWeight', 'bold');

Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'CN','CpC p=0.5','CpC p=0.7','CpC p=0.9'},'Location','southeast');
set(Legend, 'FontSize', 14); 
set(Legend, 'FontWeight', 'bold');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;

