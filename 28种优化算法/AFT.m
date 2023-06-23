
%% Ali baba and the Forty Thieves (AFT) algorithm source code version 1.0
%
%  Developed in MATLAB R2018a
%
%  Programmed by Malik Braik
%  Al-Balqa Applied University (BAU) %
%  e-Mail: m_fjo@yahoo.com
%          mbraik@bau.edu.au

% -------------------------------------------------
% This demo implements a standard version of AFT algorithm for minimization problems 
% of a standard test function on MATLAB (R2018).
% -------------------------------------------------	
% Note:
% Due to the stochastic nature of meta-heuristcs, 
% different runs may slightly produce different results.

function [fitness,gbest,ccurve]=AFT(noThieves,itemax,lb,ub,dim,fobj)

%% Convergence curve
ccurve=zeros(1,itemax);

% f1 =  figure (1);
% set(gcf,'color','w');
% hold on
% xlabel('Iteration','interpreter','latex','FontName','Times','fontsize',10)
% ylabel('fitness value','interpreter','latex','FontName','Times','fontsize',10); 
% grid;

% Fit dimension 
if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

%% Start AFT  
% % Generation of initial solutions (position of thieves)

xth=zeros(noThieves, dim);

for i=1:noThieves 
    for j=1:dim
        xth(i,j)=  lb(j)-rand()*(lb(j)-ub(j)) ; % Position of the thieves in the space
    end
end

%% Evaluate the fitness of the initial population (thieves)

fit=zeros(noThieves, 1);

for i=1:noThieves
     fit(i,1)=fobj(xth(i,:));
end

%% Initalization of the parameters of AFT

fitness=fit; % Initial the fitness of the random positions of the thieves
 
%Calculate the initial fitness of the thieves

[sorted_thieves_fitness,sorted_indexes]=sort(fit);

for index=1:noThieves
    Sorted_thieves(index,:)=xth(sorted_indexes(index),:);
end

gbest=Sorted_thieves(1,:); % initialization of the global position of the thieves
fit0=sorted_thieves_fitness(1);

best = xth;          % initialization of the best position of the thieves (based on Marjaneh's astute plans)
xab=xth;             % initialization of the position of Ali Baba


%% Start running AFT
for ite=1:itemax

    % Two essential parameters for the AFT algorithm
    Pp=0.1*log(2.75*(ite/itemax)^0.1); % Perception potential 
    Td = 2* exp(-2*(ite/itemax)^2);   % Tracking distance

    % Generation of random candidate followers (thieves (chasing))
    a=ceil((noThieves-1).*rand(noThieves,1))'; %        
%%  % Generation of a new position for the thieves 
for i=1:noThieves
     if (rand>=0.5)% In  this case, the thieves know where to search (TRUE); 
            if rand>Pp 
              xth(i,:)=gbest +(Td*(best(i,:)-xab(i,:))*rand+Td*(xab(i,:)-best(a(i),:))*rand)*sign(rand-0.50);% Generation of a new position for the thieves (case 1) 
            else
                
                for j=1:dim
                xth(i,j)= Td*((ub(j)-lb(j))*rand+lb(j)); % Generation of a new position for the thieves (case 3)
                end
                
            end
            % Generation of a new position for the thieves (case 2)
     else  %In  this case, the thieves do NOT know where to search (based on Marjaneh's astute plans), but the thieves can GUESS to look at several opposite directions (Tricks)
       
            for j=1:dim
                xth(i,j)=gbest(j) -(Td*(best(i,j)-xab(i,j))*rand+Td*(xab(i,j)-best(a(i),j))*rand)*sign(rand-0.50);          
                
            end
     end
   
end
    
%% Update the global, best, position of the thieves and Ali Baba as well
  
 for i=1:noThieves 
        
  % Evaluation of the fitness

       fit(i,1)=fobj(xth(i,:));

       % Handling the boundary conditions
                
       if and (~(xth(i,:)-lb<= 0),~(xth(i,:) - ub>= 0))            
           xab(i,:)=xth(i,:); % Update the position of Ali Baba based on the movements of the thieves and Marjaneh plans
            
            if fit(i)<fitness(i)
                 best(i,:) = xth(i,:); % Update Marjaneh astute plans based on the best position of the thieves
                 fitness(i)=fit(i);    % Update the fitness
            end
        
      % Finding out the best positions of the thieves   
      % % Updating gbest solutions based on thieves' position and Marjaneh trickes

            if fitness(i)<fit0
               fit0=fitness(i);
               gbest=best(i,:); % Update the global best positions of the thieves 
            end
        end
 end
%% 
% Show the results
% disp(['Iteration# ', num2str(ite) , '  Fitness= ' , num2str(fit0)]);
  
ccurve(ite)=fit0; % Best found value until iteration ite
 
%  if ite>2
%      set(0, 'CurrentFigure', f1)
% 
%         line([ite-1 ite], [ccurve(ite-1) ccurve(ite)],'Color','b'); 
%         title({'Convergence characteristic curve'},'interpreter','latex','FontName','Times','fontsize',12);
%         xlabel('Iteration');
%         ylabel('Best score obtained so far');
%         drawnow 
%  end 
end

bestThieves=find(fitness==min(fitness)); 
gbestSol=best(bestThieves(1),:);  % Solutin of the problem

fitness =fobj(gbestSol); 

end