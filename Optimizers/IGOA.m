function [TargetFitness,TargetPosition,Convergence_curve,Trajectories,fitness_history, position_history]=IGOA(N, Max_iter, lb,ub, dim, fobj)

% tic
% disp('GOA is now estimating the global optimum for your problem....')

flag=0;
if size(ub,2)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
else
    ub=ub';
    lb=lb';
end

if (rem(dim,2)~=0) % this algorithm should be run with a even number of variables. This line is to handle odd number of variables
    dim = dim+1;
%     ub = [ub'; 100];
%     lb = [lb'; -100];
    ub(dim)=100;
    lb(dim)=-100;
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
    
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
       % for k=1:2:dim  
            S_i=zeros(dim,1);
            for j=1:N
                if i~=j
                    Dist=distance(temp(:,j), temp(:,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(:,j)-temp(:,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((ub - lb).*c/2).*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total = S_i;
            
      %  end
        
        X_new = c * S_i_total'+ (TargetPosition); % Eq. (2.7) in the paper      
        GrassHopperPositions_temp(i,:)=X_new'; 
      end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% time=toc
end
%%
function d = distance(a,b)
% DISTANCE - computes Euclidean distance matrix
%
% E = distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Oct 29 16:35:48 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
% Thanx    : Nikos Vlassis

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

if (nargin ~= 2)
   error('Not enough input arguments');
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
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