function [K_total,P_total,thetamin_total,thetamax_total,q] = setcover(alpha,A_total,B_total,thetamin,thetamax)

% This function covers a polytopic set of systems
% with a finite number of subsets that are individually 
% alpha-quadratically stabilizable (alpha-QS) and polytopic 

j=1;
q=0;

good=[];
K_total=[];
P_total=[];
thetamin_total=[];
thetamax_total=[];

while q~=j
    
    q=j
    
    for i=1:q
        
        if double(ismember(i,good))==0
        
        thetamid=(thetamax(:,i)+thetamin(:,i))/2;
                
        [check,P,K]=qs_vertices(alpha,A_total,B_total,thetamin(:,i),thetamax(:,i)); % Check for alpha-QS
        
        if check==1
        
        P_total=[P_total P];
        K_total=[K_total K];
        thetamin_total=[thetamin_total thetamin];
        thetamax_total=[thetamax_total thetamax];
        
        good=[good;i] % This is set \cal{I} in the paper
        
        else
            
            [~,jmax]=max(thetamax(:,i)-thetamin(:,i));
            
            thetamax(:,j+1)=thetamax(:,i);
            thetamin(:,j+1)=thetamin(:,i);
            
            thetamax(jmax,i)=thetamid(jmax);
            thetamin(jmax,j+1)=thetamid(jmax);
            
            j=j+1;
            
        end
        
        end
        
    end
    
end


end

