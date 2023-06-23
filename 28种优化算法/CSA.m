%% Chameleon Swarm Algorithm (CSA) source codes version 1.0
%
%  Developed in MATLAB R2018a
%
%  Author and programmer: Malik Braik
%
%         e-Mail: m_fjo@yahoo.com
%                 mbraik@bau.edu.au
%

%   Main paper:
%   Malik Sh. Braik,
%   Chameleon Swarm Algorithm: A Bio-inspired Optimizer for Solving Engineering Design Problems
%   Expert Systems with Applications
%   DOI: https://doi.org/10.1016/j.eswa.2021.114685
%____________________________________________________________________________________

function [fmin0,gPosition,cg_curve]=CSA(searchAgents,iteMax,lb,ub,dim,fobj)

%%%%* 1
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
 
% if size(ub,1)==1
%     ub=ones(dim,1)*ub;
%     lb=ones(dim,1)*lb;
% end

%% Convergence curve
cg_curve=zeros(1,iteMax);
% 
% f1 =  figure (1);
% set(gcf,'color','w');
% hold on
% xlabel('Iteration','interpreter','latex','FontName','Times','fontsize',10)
% ylabel('fitness value','interpreter','latex','FontName','Times','fontsize',10); 
% grid;

%% Initial population

chameleonPositions=initialization(searchAgents,dim,ub,lb);% Generation of initial solutions 

%% Evaluate the fitness of the initial population

fit=zeros(searchAgents,1);


for i=1:searchAgents
     fit(i,1)=fobj(chameleonPositions(i,:));
end

%% Initalize the parameters of CSA
fitness=fit; % Initial fitness of the random positions
 
[fmin0,index]=min(fit);

chameleonBestPosition = chameleonPositions; % Best position initialization
gPosition = chameleonPositions(index,:); % initial global position

v=0.1*chameleonBestPosition;% initial velocity
 
v0=0.0*v;

%% Start CSA 
% Main parameters of CSA
rho=1.0;
p1=2.0;  
p2=2.0;  
c1=2.0; 
c2=1.80;  
gamma=2.0; 
alpha = 4.0;  
beta=3.0; 
 

 %% Start CSA
for t=1:iteMax
a = 2590*(1-exp(-log(t))); 
omega=(1-(t/iteMax))^(rho*sqrt(t/iteMax)) ; 
p1 = 2* exp(-2*(t/iteMax)^2);  % 
p2 = 2/(1+exp((-t+iteMax/2)/100)) ;
        
mu= gamma*exp(-(alpha*t/iteMax)^beta) ;

ch=ceil(searchAgents*rand(1,searchAgents));
%% Update the position of CSA (Exploration)
for i=1:searchAgents  
             if rand>=0.1
                  chameleonPositions(i,:)= chameleonPositions(i,:)+ p1*(chameleonBestPosition(ch(i),:)-chameleonPositions(i,:))*rand()+... 
                     + p2*(gPosition -chameleonPositions(i,:))*rand();
             else 
                 for j=1:dim
                   chameleonPositions(i,j)=   gPosition(j)+mu*((ub(j)-lb(j))*rand+lb(j))*sign(rand-0.50) ;
                 end 
              end   
end       
 %% Rotation of the chameleons - Update the position of CSA (Exploitation)

%%% Rotation 180 degrees in both direction or 180 in each direction
%  

% [chameleonPositions] = rotation(chameleonPositions, searchAgents, dim);
 
 %%  % Chameleon velocity updates and find a food source
     for i=1:searchAgents
               
        v(i,:)= omega*v(i,:)+ p1*(chameleonBestPosition(i,:)-chameleonPositions(i,:))*rand +.... 
               + p2*(gPosition-chameleonPositions(i,:))*rand;        

         chameleonPositions(i,:)=chameleonPositions(i,:)+(v(i,:).^2 - v0(i,:).^2)/(2*a);
     end
    
  v0=v;
  
 %% handling boundary violations
 for i=1:searchAgents
     if chameleonPositions(i,:)<lb
        chameleonPositions(i,:)=lb;
     elseif chameleonPositions(i,:)>ub
            chameleonPositions(i,:)=ub;
     end
 end
 
 %% Relocation of chameleon positions (Randomization) 
for i=1:searchAgents
    
    ub_=sign(chameleonPositions(i,:)-ub)>0;   
    lb_=sign(chameleonPositions(i,:)-lb)<0;
       
    chameleonPositions(i,:)=(chameleonPositions(i,:).*(~xor(lb_,ub_)))+ub.*ub_+lb.*lb_;  %%%%%*2
 
  fit(i,1)=fobj (chameleonPositions(i,:)) ;
      
      if fit(i)<fitness(i)
                 chameleonBestPosition(i,:) = chameleonPositions(i,:); % Update the best positions  
                 fitness(i)=fit(i); % Update the fitness
      end
 end


%% Evaluate the new positions

[fmin,index]=min(fitness); % finding out the best positions  


% Updating gPosition and best fitness
if fmin < fmin0
    gPosition = chameleonBestPosition(index,:); % Update the global best positions
    fmin0 = fmin;
end

%% Print the results
%   outmsg = ['Iteration# ', num2str(t) , '  Fitness= ' , num2str(fmin0)];
%   disp(outmsg);

%% Visualize the results

   cg_curve(t)=fmin0; % Best found value until iteration t

%     if t>2
%      set(0, 'CurrentFigure', f1)
% 
%         line([t-1 t], [cg_curve(t-1) cg_curve(t)],'Color','b'); 
%         title({'Convergence characteristic curve'},'interpreter','latex','FontName','Times','fontsize',12);
%         xlabel('Iteration');
%         ylabel('Best score obtained so far');
%         drawnow 
%     end 
end
 
ngPosition=find(fitness== min(fitness)); 
g_best=chameleonBestPosition(ngPosition(1),:);  % Solutin of the problem
fmin0 =fobj(g_best);

end
%%
function pos=initialization(searchAgents,dim,u,l)

% This function initialize the first population of search agents
Boundary_no= size(u,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both u and l
if Boundary_no==1
    u_new=ones(1,dim)*u;
    l_new=ones(1,dim)*l;
else
     u_new=u;
     l_new=l;   
end

% If each variable has a different l and u
    for i=1:dim
        u_i=u_new(i);
        l_i=l_new(i);
        pos(:,i)=rand(searchAgents,1).*(u_i-l_i)+l_i;
    end
end
%%
function answer = get_orthonormal(m,n)
% Produces an m x n set of orthonormal vectors, 
% (thats n vectors, each of length m)
% 
% Inputs should be two scalars, m and n, where n is smaller than 
% or equal to m.
%
% Example: >> get_orthonormal(5,4)
%
% ans =
%     0.1503   -0.0884   -0.0530    0.8839
%    -0.4370   -0.7322   -0.1961   -0.2207
%    -0.3539    0.3098    0.7467   -0.0890
%     0.7890   -0.1023    0.0798   -0.3701
%    -0.1968    0.5913   -0.6283   -0.1585



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK USER INPUT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ( (nargin==2) && (m>n) && (isnumeric(m)*isnumeric(n)) )
    
elseif ( nargin==1 && isnumeric(m) && length(m)==1 )
    
    n=m;
    
else
   error('Incorrect Inputs. Please read help text in m-file.')
end

% to get n orthogonal vectors (each of size m), we will first get a larger mxm 
% set of orthogonal vectors, and then just trim the set so it is 
% of size mxn,
%
% to get an mxm set of orthogonal vectors, 
% we can exploit the fact that the eigenvectors
% from distinct eigenspaces (corresponding to different eigenvalues) of a
% symmetric matrix are orthogonal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the orthonormal vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
while (count==0)

    % generate an mxm matrix A, then make a symmetric mxm matrix B
    A=rand(m);
    B=A'*A ;

    % now take the eigenvectors of B, 
    % eigenvectors from different eigenspaces (corresponding to different
    % eigenvalues) will be orthogonal

    % there is a chance that there will be repeat eigenvalues, which would give
    % rise to non-orthogonal eigenvectors (though they will still be independent)
    % we will check for this below
    % if this does happen, we will just start the loop over again 
    % and hope that the next randomly created symmetric matrix will not
    % have repeat eigenvalues,
    % (another approach would be to put the non-orthogonal vectors
    % through the gram-schmidt process and get orthogonal vectors that way)

    % since matlab returns unit length eigenvectors, they will also be
    % orthonormal

    [P,D] = eig(B) ;

    % can double check the orthonormality, by taking the difference between 
    % P'P and I, if it is non-zero, then there is an error (caused by repeat
    % eigenvalues), repeat the loop over again

    if ((P'*P - eye(m))>eps) 
        % error, vectors not orthonormal, repeat the random matrix draw again
        count=0;
    else
        % we want the first n of these orthonormal columns
        answer=P(:,1:n) ;
        count=1;
    end


end
end

%%
function [chameleonPositions]=rotation(chameleonPosition, searchAgents, dim)
for i=1:searchAgents      
          if (dim>2) 
              xmax=1;xmin=-1;
              th=round(xmin+rand(1,1)*(xmax-xmin));
              vec=get_orthonormal(dim,2);
              vecA=vec(:,1);
              vecB=vec(:,2);
              theta=(th*rand()*180)*(pi/180) ;
              Rot = RotMatrix(theta,vecA, vecB) ;
             if (theta~=0)
                V=[chameleonPosition(i,:) ]; 
                V_centre=mean(V,1); %Centre, of line
                Vc=V-ones(size(V,1),1)*V_centre; %Centering coordinates

                Vrc=[Rot*Vc']'; %Rotating centred coordinates
%                 Vruc=[Rot*V']'; %Rotating un-centred coordinates
                Vr=Vrc+ones(size(V,1),1)*V_centre; %Shifting back to original location
                 chameleonPosition(i,:)=((Vr)/1); 
 
             end
         else
              xmax=1;xmin=-1;
              th=round(xmin+rand(1,1)*(xmax-xmin));
              theta=th*rand()*180*(pi/180);
              Rot = RotMatrix(theta);
              
               if (theta~=0)
                V=[chameleonPosition(i,:) ];  
                V_centre=mean(V,1); %Centre, of line
                Vc=V-ones(size(V,1),1)*V_centre; %Centering coordinates

                Vrc=[Rot*Vc']'; %Rotating centred coordinates
                Vr=Vrc+ones(size(V,1),1)*V_centre; %Shifting back to original location
                chameleonPosition(i,:)=((Vr)/1);
               end
          end
end  
   chameleonPositions=chameleonPosition;
end
%%
function R = RotMatrix(alpha, u, v)
% RotMatrix - N-dimensional Rotation matrix
% R = RotMatrix(alpha, u, v)
% INPUT:
%   alpha: Angle of rotation in radians, counter-clockwise direction.
%   u, v:  Ignored for the 2D case.
%          For the 3D case, u is the vector to rotate around.
%          For the N-D case, there is no unique axis of rotation anymore, so 2
%          orthonormal vectors u and v are used to define the (N-1) dimensional
%          hyperplane to rotate in.
%          u and v are normalized automatically and in the N-D case it is cared
%          for u and v being orthogonal.
% OUTPUT:
%   R:     Rotation matrix.
%          If the u (and/or v) is zero, or u and v are collinear, The rotation
%          matrix contains NaNs.
%
% REFERENCES:
% analyticphysics.com/Higher%20Dimensions/Rotations%20in%20Higher%20Dimensions.htm
% en.wikipedia.org/wiki/Rotation_matrix
% application.wiley-vch.de/books/sample/3527406204_c01.pdf
%
% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
if numel(alpha) ~= 1
   error('JSimon:RotMatrrix:BadInput1', ...
      'Angle of rotation must be a scalar.');
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
s = sin(alpha);
c = cos(alpha);

% Different algorithms for 2, 3 and N dimensions:
switch nargin
   case 1
      % 2D rotation matrix:
      R = [c, -s;  s, c];
      
   case 2
      if numel(u) ~= 3
         error('JSimon:RotMatrrix:BadAxis2D', ...
            '3D: Rotation axis must have 3 elements.');
      end
      
      % Normalized vector:
      u = u(:);
      u = u ./ sqrt(u.' * u);
            
      % 3D rotation matrix:
      x  = u(1);
      y  = u(2);
      z  = u(3);
      mc = 1 - c;
      R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];
         
      % Alternative 1 (about 60 times slower):
      % R = expm([0, -z,  y; ...
      %           z,  0, -x; ...
      %          -y,  x,  0]  * alpha);
      
      % Alternative 2:
      % R = [ 0, -z, y; ...
      %       z, 0, -x; ...
      %      -y, x,  0] * s + (eye(3) - u * u.') * c + u * u.';
              
   case 3
      n = numel(u);
      if n ~= numel(v)
         error('JSimon:RotMatrrix:BadAxes3D', ...
            'ND: Axes to define plane of rotation must have the same size.');
      end
      
      % Normalized vectors:
      u = u(:);
      u = u ./ sqrt(u.' * u);
      
      % Care for v being orthogonal to u:
      v = v(:);
      v = v - (u.' * v) * u;
      v = v ./ sqrt(v.' * v);
      
      % Rodrigues' rotation formula:
      R = eye(n) + ...
         (v * u.' - u * v.') * s + ...
         (u * u.' + v * v.') * (c - 1);
      
   otherwise
      error('JSimon:RotMatrrix:BadNInput', ...
            '1 to 3 inputs required.');
end

end

