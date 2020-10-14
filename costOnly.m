
function j= costOnly(U)
                    global prevPrediction;
                    global outHorizon;
                    global inputHorizon;

                    load vars.mat
                    
                    U_long=[U ;zeros(pH-cH,1)];
                    for n=cH:pH
                        U_long(n)=U(end);
                    end    

             
%                   Actualizar estados da rede com amostras observadas.(Matlab)
                     [x1_tmp,xictrl,aictrl] = preparets(ctrlModelClosed,con2seq(u(k-d:k)),{},con2seq(y(k-d:k)));
                     Y=cell2mat(ctrlModelClosed(con2seq(U_long'),X(k-(d-2):k),aictrl));

                     
                     
%                   Calcular saídas uma de cada vez, aplicando dm entre predições. Necessária retirar dm para calcular o erro 
%                      tt=T(k-4:k);
%                     uu=u(k-2:k);
%                     for z=1:pH
%                         [x1,xictrl,aictrl] = preparets(ctrlModelClosed,con2seq(u(k-4:k)),{},tt);
%                          Y(z) = ctrlModel(U_long(z),uu,cell2mat(aictrl(2,:)))+dm;
%                          tt=[tt(2:end) Y(z)];
%                          uu=[uu(2:end) U_long(z)];
%                      end
                     
                
                    Y=Y';

                    prevPrediction=Y(1);%-dm;
                    Y=Y+Dm';
                    
                    outHorizon=Y; %Plots
                    inputHorizon=U;
                    
                    E=R-Y;
                    U_inc=[U(1)-Uprev; diff(U)];
                    

                    j=zeros(1,cH);

                    for i=1:cH
                       j(i)=sqrt(E(i)^2+p*U_inc(i)^2);
                   
                    end               
                    if(pH>cH)
                        j(cH)=sqrt(E(cH:end)'*E(cH:end)+p*U_inc(cH:end)'*U_inc(cH:end));
                    end    
%                     
                    
                    
end