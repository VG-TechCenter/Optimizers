%_________________________________________________________________________%
%  Grasshopper Optimization Algorithm (GOA) source codes demo V1.0        %
%                                                                         %
%  Developed in MATLAB R2016a                                             %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%  Main paper: S. Saremi, S. Mirjalili, A. Lewis                          %
%              Grasshopper Optimisation Algorithm: Theory and Application %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004   %
%                                                                         %
%_________________________________________________________________________%

% The Grasshopper Optimization Algorithm
function [TargetFitness,TargetPosition,Convergence_curve,Trajectories,fitness_history, position_history]=GOA(N, Max_iter, lb,ub, dim, fobj)

% disp('GOA is now estimating the global optimum for your problem....')

flag=0;
if size(ub,1)==1
    ub=ones(dim,1).*ub;
    lb=ones(dim,1).*lb;
end

if (rem(dim,2)~=0) % this algorithm should be run with a even number of variables. This line is to handle odd number of variables
    dim = dim+1;
    ub(dim)=100;
    lb(dim)=100;
%     ub = [ub,100];
%     lb = [lb,-100];
    flag=1;
end

%Initialize the population of grasshoppers
GrassHopperPositions=initialization(N,dim,ub,lb);
GrassHopperFitness = zeros(1,N);

fitness_history=zeros(N,Max_iter);
position_history=zeros(N,Max_iter,dim);
Convergence_curve=zeros(1,Max_iter);
Trajectories=zeros(N,Max_iter);

cMax=1;
cMin=0.00004;
%Calculate the fitness of initial grasshoppers

for i=1:size(GrassHopperPositions,1)
    if flag == 1
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
    else
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
    end
    fitness_history(i,1)=GrassHopperFitness(1,i);
    position_history(i,1,:)=GrassHopperPositions(i,:);
    Trajectories(:,1)=GrassHopperPositions(:,1);
end

[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);

% Find the best grasshopper (target) in the first population 
for newindex=1:N
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end

TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);

% Main loop
l=2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
while l<Max_iter+1
    
    c=cMax-l*((cMax-cMin)/Max_iter); % Eq. (2.8) in the paper
    
    for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
        for k=1:2:dim
            S_i=zeros(2,1);
            for j=1:N
                if i~=j
                    Dist=distance(temp(k:k+1,j), temp(k:k+1,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(k:k+1,j)-temp(k:k+1,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((ub(k:k+1) - lb(k:k+1))*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total(k:k+1, :) = S_i;
            
        end
        
        X_new = c * S_i_total'+ (TargetPosition); % Eq. (2.7) in the paper      
        GrassHopperPositions_temp(i,:)=X_new'; 
    end
    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp;
    
    for i=1:size(GrassHopperPositions,1)
        % Relocate grasshoppers that go outside the search space 
        Tp=GrassHopperPositions(i,:)>ub';Tm=GrassHopperPositions(i,:)<lb';GrassHopperPositions(i,:)=(GrassHopperPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        
        % Calculating the objective values for all grasshoppers
        if flag == 1
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
        else
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
        end
        fitness_history(i,l)=GrassHopperFitness(1,i);
        position_history(i,l,:)=GrassHopperPositions(i,:);
        
        Trajectories(:,l)=GrassHopperPositions(:,1);
        
        % Update the target
        if GrassHopperFitness(1,i)<TargetFitness
            TargetPosition=GrassHopperPositions(i,:);
            TargetFitness=GrassHopperFitness(1,i);
        end
    end
        
    Convergence_curve(l)=TargetFitness;
%     disp(['In iteration #', num2str(l), ' , target''s objective = ', num2str(TargetFitness)])
    
    l = l + 1;
end


if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end

end
%%
function d = distance(a,b)
d=sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
end
%%
function [X]=initialization(N,dim,up,down)

if size(up,1)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
end
%%
function o=S_func(r)
f=0.5;
l=1.5;
o=f*exp(-r/l)-exp(-r);  % Eq. (2.3) in the paper
end
