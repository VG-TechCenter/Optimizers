
function [Silverback_Score,Silverback,convergence_curve]=GTO(pop_size,max_iter,lower_bound,upper_bound,variables_no,fobj)


% initialize Silverback
Silverback=[];
Silverback_Score=inf;

%Initialize the first random population of Gorilla
X=initialization(pop_size,variables_no,upper_bound,lower_bound);


convergence_curve=zeros(max_iter,1);

for i=1:pop_size   
    Pop_Fit(i)=fobj(X(i,:));%#ok
    if Pop_Fit(i)<Silverback_Score 
            Silverback_Score=Pop_Fit(i); 
            Silverback=X(i,:);
    end
end


GX=X(:,:);
lb=ones(1,variables_no).*lower_bound; 
ub=ones(1,variables_no).*upper_bound; 

%%  Controlling parameter

p=0.03;
Beta=3;
w=0.8;

%%Main loop
for It=1:max_iter 
    
    a=(cos(2*rand)+1)*(1-It/max_iter);
    C=a*(2*rand-1); 

%% Exploration:

    for i=1:pop_size
        if rand<p    
            GX(i,:) =(ub-lb)*rand+lb;
        else  
            if rand>=0.5
                Z = unifrnd(-a,a,1,variables_no);
                H=Z.*X(i,:);   
                GX(i,:)=(rand-a)*X(randi([1,pop_size]),:)+C.*H; 
            else   
                GX(i,:)=X(i,:)-C.*(C*(X(i,:)- GX(randi([1,pop_size]),:))+rand*(X(i,:)-GX(randi([1,pop_size]),:))); %ok ok 

            end
        end
    end       
       
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation 
    for i=1:pop_size
         New_Fit= fobj(GX(i,:));          
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
    
%% Exploitation:  
    for i=1:pop_size
       if a>=w  
            g=2^C;
            delta= (abs(mean(GX)).^g).^(1/g);
            GX(i,:)=C*delta.*(X(i,:)-Silverback)+X(i,:); 
       else
           
           if rand>=0.5
              h=randn(1,variables_no);
           else
              h=randn(1,1);
           end
           r1=rand; 
           GX(i,:)= Silverback-(Silverback*(2*r1-1)-X(i,:)*(2*r1-1)).*(Beta*h); 
           
       end
    end
   
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation    
    for i=1:pop_size
         New_Fit= fobj(GX(i,:));
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
             
convergence_curve(It)=Silverback_Score;
% fprintf("In Iteration %d, best estimation of the global optimum is %4.4f \n ", It,Silverback_Score );
         
end 
end
%%
function [ X ]=initialization(N,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(N,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end
%%
function [ X ] = boundaryCheck(X, lb, ub)

    for i=1:size(X,1)
            FU=X(i,:)>ub;
            FL=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
    end
end