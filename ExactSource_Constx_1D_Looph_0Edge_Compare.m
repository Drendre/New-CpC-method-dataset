%% Test against exact solution in 1D
%% loop for the timesteps 
clear; fprintf('\n')
Mh=14; %How many different timestep we want 
D=1; xf=1; TIME=0.02; p=0.5; %Diffusion coefficient, lenght of the space and time interval, x parameter
N=51; c1=5; c2=17; lamb=16; k0=-8; % initial number of cells, coefficients of sin modes, wavelength
MaxD=zeros(Mh,6);Axhstep=zeros(Mh,1); %%Errors and axis for the timestep-sizes
h=TIME/2; Dx=xf/(N-1); %%initial timestep and Space-step
axis=zeros(N,1); for i=1:N axis(i)=(i-1)*Dx; end  %x axis
M1=zeros(N,N); 
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
uee=zeros(N,1);utr=zeros(N,1);umi=zeros(N,1);ura=zeros(N,1);udf=zeros(N,1);u2=zeros(N,1);uk1=zeros(N,1); du=zeros(N,1);
E=D/(Dx)^2; %%matrix elements
k=zeros(N,1); %%sources, now zero
for i=1:N 
    k(i)=k0*pi^2*sin(pi*axis(i));
u0(i)=c1*sin(pi*axis(i))+c2*sin(lamb*pi*axis(i)); 
     uex(i)=(c1*exp(-pi^2*TIME)+k0*(1-exp(-pi^2*TIME)))*sin(pi*axis(i)) +c2*sin(lamb*pi*axis(i))*exp(-lamb^2*pi^2*TIME); end  %%Initial and exact final results
M1(1,1)=-E;M1(1,2)=E; M1(N,N)=-E;M1(N,N-1)=E;
for i=2:N-1     M1(i,i-1)=E; M1(i,i)=-2*E; M1(i,i+1)=E; end
Neb1=zeros(N,1);Neb2=zeros(N,1); Sum=zeros(N,1); %Dufort Frankel
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
 h1=p*h; T = uint32(TIME/h)+1;
 Axhstep(ih)=h; %%Timestep axis
 UC=zeros(N,1); Ux=zeros(N,1); UTEMP=zeros(N,1); UTEMP2=zeros(N,1); a0=zeros(N,1); %%Final U and temporary results
 b=zeros(N,1);ee=zeros(N,1); 
for i=1:N  b(i)=-M1(i,i);
ee(i)=exp(-b(i)*h); ee1(i)=exp(-b(i)*h1);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UC(:)=u0(:); Ux(:)=u0(:);  utr(:)=u0(:);umi(:)=u0(:);ura(:)=u0(:);udf(:)=u0(:); uee(:)=u0(:);%%initial values
%  GITC=it*(T-1); 
%% Loop for time
tic
for t=1:T-1
    for i=2:N-1  %%%Neighbours are constans
        a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*UC(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*UC(i+1); end     %%rigth neighbour
       UTEMP(i)=UC(i)*ee(i)+ a0(i)/b(i)*(1-ee(i));
    end
    UTEMP(1)=0;UTEMP(N)=0;
 UC(:)=UTEMP(:); 
end
toc
%% 2nd order Constant-x starts here
tic
for t=1:T-1
  for i=2:N-1  %%%Neighbours are constans, first stage
         a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*Ux(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*Ux(i+1); end     %%rigth neighbour  
       UTEMP(i)=Ux(i)*ee1(i)+ a0(i)/b(i)*(1-ee1(i));
  end
     UTEMP(1)=0;UTEMP(N)=0; %%Boundary conditions
 % UTEMP(:)=(1-1/2/p)*U(:,1)+1/2/p*UTEMP(:);
 for i=2:N-1  %%%Neighbours are constans, second stage
     a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*UTEMP(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*UTEMP(i+1); end     %%rigth neighbour    
       UTEMP2(i)=Ux(i)*ee(i)+ a0(i)/b(i)*(1-ee(i));
  end
  Ux(:)=UTEMP2(:);   Ux(1)=0;Ux(N)=0;
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
%% Const-x errors
MaxKul=0; 
for i=1:N
    Kul(i)=Ux(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,2)=MaxKul; 
fprintf('Const-x h: %f \nMaxK: %f \n',h,MaxKul)
%%%%%%%%%%%%%%%%%%%%%%%%Explicit Euler FTCS
if h<hMAX
tic
for t=1:T-1
    for i=1:N
        du(i)=M1(i,i)*uee(i);
     if(i>1) du(i)=du(i)+M1(i,i-1)*uee(i-1); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*uee(i+1); end     %%rigth neighbour 
       uk1(i)=uee(i)+h*(du(i)+k(i));
    end
uee(:)=uk1(:);uee(1)=0;uee(N)=0;
end
toc
%% EE errors
MaxKul=0; 
for i=1:N
    Kul(i)=uee(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,3)=MaxKul;
fprintf('\n EE h: %f \nMaxK: %f \n',h,MaxKul)
end
if h<1.1*hMAX
tic %%%Trapez
for t=1:T-1
    for i=1:N
        du(i)=M1(i,i)*utr(i);
     if(i>1) du(i)=du(i)+M1(i,i-1)*utr(i-1); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*utr(i+1); end     %%rigth neighbour 
        uk1(i)=utr(i)+h*(du(i)+k(i));
    end
    uk1(1)=0;uk1(N)=0;
    for i=1:N
        du(i)=M1(i,i)*(utr(i)+uk1(i));
     if(i>1) du(i)=du(i)+M1(i,i-1)*(utr(i-1)+uk1(i-1)); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*(utr(i+1)+uk1(i+1)); end     %%rigth neighbour        
        UTEMP(i)=utr(i)+h*(du(i)/2+k(i));
    end
    utr(:)=UTEMP(:);
utr(1)=0;utr(N)=0;
end
toc
%% Trapez errors
MaxKul=0; 
for i=1:N
    Kul(i)=utr(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,4)=MaxKul;
fprintf('\n Trap h: %f \nMaxK: %f \n',h,MaxKul)
tic %%Midpoint
for t=1:T-1
    for i=1:N
        du(i)=M1(i,i)*umi(i);
     if(i>1) du(i)=du(i)+M1(i,i-1)*umi(i-1); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*umi(i+1); end     %%rigth neighbour 
        uk1(i)=umi(i)+h*(du(i)+k(i))/2;
    end
    uk1(1)=0;uk1(N)=0;
    for i=1:N
        du(i)=M1(i,i)*uk1(i);
     if(i>1) du(i)=du(i)+M1(i,i-1)*uk1(i-1); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*uk1(i+1); end     %%rigth neighbour        
        umi(i)=umi(i)+h*(du(i)+k(i));
    end
umi(1)=0;umi(N)=0;
end
toc
%% Midpoint errors
MaxKul=0; 
for i=1:N
    Kul(i)=umi(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,5)=MaxKul;
fprintf('\n Midpoint h: %f \nMaxK: %f \n',h,MaxKul)
tic %%%Ralston
for t=1:T-1
    for i=1:N
        du(i)=M1(i,i)*ura(i);
     if(i>1) du(i)=du(i)+M1(i,i-1)*ura(i-1); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*ura(i+1); end     %%rigth neighbour 
        uk1(i)=ura(i)+2/3*h*(du(i)+k(i));
    end
    uk1(1)=0;uk1(N)=0;
    for i=1:N
        du(i)=M1(i,i)*(ura(i)+3*uk1(i));
     if(i>1) du(i)=du(i)+M1(i,i-1)*(ura(i-1)+3*uk1(i-1)); end  %%left neighbour
     if(i<N) du(i)=du(i)+M1(i,i+1)*(ura(i+1)+3*uk1(i+1)); end     %%rigth neighbour        
        UTEMP(i)=ura(i)+h*(du(i)/4+k(i));
    end
ura(:)=UTEMP(:); ura(1)=0;ura(N)=0;
end
toc
%% Ralston errors
MaxKul=0; 
for i=1:N
    Kul(i)=ura(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,6)=MaxKul;
end
%% DufortFrankel
u2(:)=u0(:); uk1(:)=0; Sum(:)=h*(Neb1(:)+Neb2(:));
for i=1:N  udf(i)=sin(pi*axis(i))*(c1*exp(-pi^2*h)+k0*(1-exp(-pi^2*h)))+c2*sin(lamb*pi*axis(i))*exp(-lamb^2*pi^2*h); end %%Previous values
tic
for t=2:T-1
     for i=2:N-1  %%Dufort Frankel
      a0=k(i); 
        a0=a0+Neb1(i)*udf(i-1);  %bal szomszéd
        a0=a0+ Neb2(i)*udf(i+1); %jobb szomszéd
         uk1(i)=(2*h*a0 +(1-Sum(i))*u2(i)) / (1+Sum(i));
     end
    uk1(1)=0;uk1(N)=0; 
u2(:)=udf(:);
udf(:)=uk1(:);
end
toc
%% DufortFrank errors
MaxKul=0; 
for i=1:N
    Kul(i)=udf(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,7)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  remh=rem(ih,3);
%      if (remh==2) h=2*h/5;
%      else h=h/2; end
h=h/2;
end
% plot(axis,u0,axis,uex,axis,UC,axis,udf,'--o','LineWidth',2)
plot(Axhstep(:),MaxD(:,1),Axhstep(:),MaxD(:,2),Axhstep(:),MaxD(:,3),Axhstep(:),MaxD(:,4),Axhstep(:),MaxD(:,5),Axhstep(:),MaxD(:,6),Axhstep(:),MaxD(:,7),'--o','LineWidth',2)
set(gca,'xScale','log')
set(gca,'yScale','log')
grid on
