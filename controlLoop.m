

clear functions
clear all
%close all




load gasifierModel3.mat


N = 400;
u = zeros(1,N);
% u(11:end)=NaN % mostra só amostras até k no plot do sinal de controlo 
y = zeros(1,N);



sp=800
reference=sp*ones(1,N);
reference(floor(end/3):end)=750;
reference(floor(end/1.4):end)=850;


%Parametros
pH=5%prevision horizon  
cH=1 %control horizon
p=100% 

n=0.1; %resposta a perturbações (Zé)
kp=0.1;          
 
%Inicializações
R=zeros(pH,1) %reference Horizon
Y=zeros(pH,1) %prediction horizon
E=zeros(pH,1) %error horizon
U_inc=zeros(cH,1) %control increment horizon
U=zeros(cH,1) %control increment horizon



%Perturbações

%inVar=0.35/1000
%inDisturb=sqrt(inVar)*randn(1,N);
ld=0.3/10;
inDisturb=zeros(1,N);
inDisturb(floor(50):end)=ld;
% inDisturb(floor(100):end)=ld/2;
% inDisturb(floor(150):end)=ld/3;
outVar=4; %erro de 2 graus 
outDisturb=outVar*(rand(1,N)-0.5);
% 

%Saturação
uHigh=0.35*ones(1,cH);
uLow=0.1*ones(1,cH);


h =  1;


%Controlador só começa em k=10
y=zeros(1,N-10);
y(11:end)=NaN;
y(1:10)=out(1:10);
u(1:10)=in(1:10);
predictions=NaN(1,N-10);
controlPredictions=NaN(1,N-10);
X=con2seq(u);
T=con2seq(y);

% Control signal limits

%Error filter
 fs=1
 Tao=10
fc=1/Tao;
Tmin=3;
fu=1/Tmin;
Whp=(1/Tao)/fs;
num=[(1-Whp) -exp(-1)];
den=[1 -exp(-1)]


% wc=fc/fs%*2*pi
%  [num,den]=fir1(40,wc,'low')
% [num,den]=fir1(40,wc,'high')
% [num,den]=fir1(41,[wc 0.9])


%%% Net init
% plant=ctrlModelOpen;
ctrlModelOpen=plant;
ctrlModelClosed=closeloop(ctrlModelOpen);
ctrlModelClosed=removedelay(ctrlModelClosed,1);


d=4; %numero delays da rede

[x1,xio,aio] = preparets(plant,X(1:9),{},T(1:9));
[y1,xfo,afo] = plant(x1,xio,aio);
plantO=plant;
[plant,xic,aic] = closeloop(plant,xfo,afo);
plant=removedelay(plant,1);

global prevPrediction;
global outHorizon;
global inputHorizon;
 prevPrediction=out(9);

dm=0;
count=0;
dmProf=zeros(1,N);
DmProf=zeros(1,N);
wd=zeros(1,N)
b=zeros(1,N)
e=zeros(1,N);
ee=zeros(1,N);
PM = 100*eye(4);
Theta= 0.5*ones(1, 4);
eFilt=zeros(1,N);
Dm=zeros(1,pH);
disp('Running...      ')
load distNet;


for k = 10 : N-10

       % fprintf(1,'\b\b\b\b\b\b  %.2d %%',floor(k/N*100)) %progress status
        
    
        %%%Plant           
        [y2,xfc,aic] = plant(X(k),X(k-(d-1):k-1),aic);
        y(k)=cell2mat(y2);% + outDisturb(k);
        %y=filter(num,den,y);
        
        
        
 
        
        %%%
        %Controller
        R=reference(k:k+pH-1)';
        Uprev=U(1);   
        
        %Disturbance compensation
        e(k)=y(k)-prevPrediction;
%         eFiltDelayed=filter(num,den,e);
       % eFilt=[eFiltDelayed(11:end) zeros(1,10)];


      
        fprintf('plant=%.2f  pred=%.2f',y(k),prevPrediction)
        
        %%Compensador da dissertação do Zé
%       eFilt(k)=num(1)*e(k)+num(2)*e(k-1)+den(2)*eFilt(k-1);
%       b(k)=b(k-1)+n*eFilt(k)+kp*(eFilt(k)-eFilt(k-1));
%       wd(k)=wd(k-1)+n*(eFilt(k)^2);
%       dm=wd(k)*eFilt(k)+b(k);
       
       
       %DMC
%       dm=e(k);
        
        
        %Nada
        dm=0;
        dmProf(k)=dm;
        Dm=Dm*0+dm; %perturbação constante no horizonte

        
%Modelo de perturbações com rede neuronal
%           distNetC=closeloop(distNet);
%           [x1,xid,aid] = preparets(distNet,{},{},con2seq(dmProf(k-4:k)));
%             x1=cell2mat(x1);
%             xid=cell2mat(xid);
%           for j=1:pH
%           %[y3,aid] = distModel(aid)    
%           [x1,xid] = distModel(x1,xid);
%           Dm(j)=x1;
%           end





    %RLS 4th Order AR
%           Phi=dmProf(k-4:k-1);
%           
%           npar=4;
%           ForgFactor=0.99;
%           Dm_ext=[dmProf(k-3:k) Dm];
%           [Theta, PM, PredError] = rls_am(Theta, Phi, PM, dm, npar, ForgFactor)
%           for j=5:5+pH-1
%           Dm_ext(j)=Theta(1)*Dm_ext(j-4)+Theta(2)*Dm_ext(j-3)+Theta(3)*Dm_ext(j-2)+Theta(4)*Dm_ext(j-1);
%           end
%           Dm=Dm_ext(5:end);
%           DmProf(k+1)=Dm(1);
% 
%           Dm=Dm*0;




        
        %Nota: A maneira correta é criar uma handle:
        % h = @(x0) costOnly(x0,,'plant','R','k','X','cH',...);
        % Mas o save é mais confortavel para depurar.(É só fazer load na função de custo)
        save('vars.mat','aic','plant','R','k','X','cH','pH','p','Uprev','ctrlModelOpen','T','d','n','kp','dm','u','Dm','ctrlModelClosed','dmProf','y')

        %Minimização
%        x0=U;
          x0=[U(2:end) U(end)]';
%          x0=(uHigh'-uLow')/2;
       
        options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
        U = lsqnonlin(@costOnly,x0,uLow',uHigh',options)

         X{k+1}=U(1);%+inDisturb(k+1);
         u(k+1)=U(1);
         
        %PLOTS      
        predictions=[y(1:k) outHorizon' NaN(1,N-10-(k+pH))];
        controlPredictions=[u(1:k) inputHorizon' NaN(1,N-10-(k+cH))];
        disturbancePredictions=[e(1:k) Dm NaN(1,N-10-(k+pH))];
         figure(10)
        subplot(2,1,1)
        plot([predictions(1:N-10)' y(1:N-10)' reference((1:N-10))'])
        %plot([y' out((1:N-10))'])
        legend('Prediction horizon','Plant Output','Reference')
%          axis([1 N-10 720 900])
        title([' TRR MemDepth=',num2str(d),' pH=',num2str(pH),' cH=',num2str(cH),' p=',num2str(p),' kp=',num2str(kp),' n=',num2str(n)])
        subplot(2,1,2)
        plantInput=cell2mat(X(1:N-10));
        controlSignal=u(1:N-10);
        plot([controlPredictions(1:N-10)' plantInput(1:N-10)' controlSignal(1:N-10)' ])
        legend('Control Horizon','Plant input','Computed Control Signal')
 
        figure(11)
        subplot(2,1,1)
        plot([dmProf(1:N-10)'  b(1:N-10)' wd(1:N-10)' ])
        legend('dm','b','wd')
        title(['Disturbance model fc=1/',num2str(Tao),' kp=',num2str(kp),' n=',num2str(n)])
        subplot(2,1,2)
        plot([disturbancePredictions(1:N-10)' DmProf(1:N-10)' e(1:N-10)' eFilt(1:N-10)'])
        legend('Disturbance predictions','Step ahead prediction','Disturbance','Filtered Disturbance')
        
        
        drawnow 
      
end


