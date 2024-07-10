function [bestnest,fmin,bb]=KMNCEC2014(Fun,Dim,n,N_iter,Lb,Ub);
% [Dim,Fun,Ub,Lb]  = Select_Functions(function_name);
% (Lb,Ub,Dim,Fun,Iterations,PopSize);
global initial_flag
% n=20;
a1=0.95;
a2=0.05;
Q =0;
%pa=a2+((a1-a2).*exp(-(t/(N_iter/10))));
% nd=30;
Lb=Lb.*ones(1,Dim);
Ub=Ub.*ones(1,Dim);
% disp('Version 2.0 CS');
for i=1:n,
    nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end
% iter=1;
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,Fun);
% N_iter=500;
for t=1:N_iter

    %%% -----------Adapt I-----------RPO---------------%
    % I =rand();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %  r1=rand(); % r1 is a random number in [0,1]
           %  r2=rand(); % r2 is a random number in [0,1]
           % I=2-t*((2)/N_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
            % %%%%%% weight-1 %%%%%%Natural Exponent Inertia Weight %%%%%%%%%%
                       alpha_min= rand;%0.25;
                       alpha_max= rand; %0.8;
                       if alpha_min==alpha_max
                            alpha_max=0.9;
                          end
                  I= alpha_min+(alpha_max-alpha_min)*exp((-t)/(N_iter/10));

            %%%%%%%%%%%%  oscillating inertia weight%%%%%%%%%%%%%%%
                          %  alpha_min= 0.9;
                          %  alpha_max= 0.3;
                          % T=(2*10000)/(3+(2*0.3));
                          % I=((alpha_max+alpha_min)/2)+((alpha_min-alpha_max)/2)*cos((2*pi*t)/T);


    %%% -----------Adapt R---------------RPO------------%
% % 1%     
                 R = rand();
% % %  2   %%% **** Chaotic Inertia weight *******%%%%%
           %       alpha_min=0.5;%0.25;
           %       alpha_max=0.9; k=rand;
           %      k=4*k*(1-k);
           % R=((alpha_max-alpha_min)*((N_iter-t)/N_iter))+(alpha_max)*k;
% % %   3  % weight-1 %%%%%%Natural Exponent Inertia Weight %%%%%%%%%%
         %        alpha_min=rand;%0.25;
         %        alpha_max=rand;%0.8;
         %         if alpha_min==alpha_max
         %            alpha_max=0.9;
         %          end
         % R=alpha_min+(alpha_max-alpha_min)*exp((-t)/(N_iter/10));




    
%   1)  %%% -----------Adapt Lambda --------------NMRA-----------%
%     alpha_min=rand;%0.25;
%     alpha_max=rand;%0.8;
%     k=rand;
%     lambda=alpha_max+(((alpha_max-alpha_min)*k)/N_iter);
% %    2)        %%%%%%%%%%%%%%%%Simulated Annealing-orignal %%%%%%%%%%%%
         alpha_min=0.5;%0.25;
         alpha_max=0.9; k=rand; p=0.95;
        lambda=alpha_min+((alpha_max-alpha_min)*p^(k-1));

% 3) Sigmoidal decreasing iw (sig)   
        %       ha=rand;%0.25;
        %        ue=0.51;
        % lambda=(0.5-0.9)/(1+exp(((-ue)*(t-ha*51))))+0.9;



    %%%%%--------------Adapt pa ----------------CS-GWO--------%%%%%%
%%% Value of PA
% % 1) weight-1 %%%%%%Natural Exponent Inertia Weight %%%%%%%%%% Exponential decreasing iw (exp)
                 alpha_min=rand;%0.25;
                 alpha_max=rand;%0.8;
                 if alpha_min==alpha_max
                     alpha_max=0.9;
                  end
           pa=alpha_min+(alpha_max-alpha_min)*exp((-t)/(N_iter/10));
% % % 2) %%%%%%%%%%%%  oscillating inertia weight%%%%%%%%%%%%%%% (Oscillating decreasing)
%                             T=(2*500)/(3+(2*0.3));
%                             pa=((0.3+0.9)/2)+((0.9-0.3)/2)*cos((2*pi*t)/T); 
% % 3) Logarithmic decreasing iw (log)
                            % pa=0.9+(0.25-0.9)*log(1+(10*t/N_iter));
                                

    if t>N_iter/2
        new_nest=get_cuckoos(nest,bestnest,Lb,Ub,t); %CS-GWO
        [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Fun);
        new_nest=empty_nests(best,nest,Lb,Ub,pa, R, I) ; % RPO
        [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Fun);
        if fnew<fmin,
            fmin=fnew;
            bestnest=best;
        end
    else
        new_nest=get_different_cuckoos(nest,bestnest,Lb,Ub, Q);
        [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Fun);
        new_nest=new_empty_nests(nest,Lb,Ub,pa, lambda); % NMRA
        [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Fun);
        if fnew<fmin,
            fmin=fnew;
            bestnest=best;
        end
    end
    bb(t)=fmin;
    t=t+1;
    initial_flag=1;
end
%% --------------- All subfunctions are list below ------------------
%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Fun)
% Evaluating all new solutions
for j=1:size(nest,1),
    %fnew=Fun(newnest(j,:));
    fnew=Fun(new_nest(j,:));
   % fnew=feval(Fun,new_nest(j,:)',func_no)';
    if fnew<=fitness(j),
        fitness(j)=fnew;
        nest(j,:)=new_nest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness);
best=nest(K,:);

function nest=get_cuckoos(nest,best,Lb,Ub,t)
n=size(nest,1);
beta=1;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
%maxiter=10000;
for j=1:n/2,
          %%%%%%%%%%%% 1) oscillating inertia weight%%%%%%%%%%%%%%%
                           T=(2*10000)/(3+(2*0.3));
                           a=((0.3+0.9)/2)+((0.9-0.3)/2)*cos((2*pi*t)/T);

            % % 2)  weight-1 %%%%%%Natural Exponent Inertia Weight %%%%%%%%%%
            %             alpha_min= rand;%0.25;
            %             alpha_max= rand; %0.8;
            %                 if alpha_min==alpha_max
            %                 alpha_max=0.9;
            %                 end
            %        a= alpha_min+(alpha_max-alpha_min)*exp((-t)/(10000/10));

                % % 3)  % % % % % a decreases linearly fron 2 to 0
                %     a=2-j*((2)/10000); % a decreases linearly fron 2 to 0
    s=nest(j,:);
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    stepsize=0.01*step.*(s-best);
    s=s+stepsize.*randn(size(s));
end
%---------------------------------------------------------------
% ----- Grey Wolf Optimization ----- %%%
%---------------------------------------------------------------
for j=((n/2)+1):n
    r1=rand(); % r1 is a random number in [0,1]
    r2=rand(); % r2 is a random number in [0,1]
    A1=2*a*r1-a; % Equation (3.3)
    C1=2*r2;% Equation (3.4)
    D_alpha=abs(C1*best-nest(j,:)); % Equation (3.5)-part 1
    X1=best-A1*D_alpha;% Equation (3.6)-part 1
    r1=rand();
    r2=rand();
    A2=2*a*r1-a; % Equation (3.3)
    C2=2*r2; % Equation (3.4)
    D_beta=abs(C2*best-nest(j,:)); % Equation (3.5)-part 2
    X2=best-A2*D_beta; % Equation (3.6)-part 2
    r1=rand();
    r2=rand();
    A3=2*a*r1-a; % Equation (3.3)
    C3=2*r2; % Equation (3.4)
    D_delta=abs(C3*best-nest(j,:)); % Equation (3.5)-part 3
    X3=best-A3*D_delta; % Equation (3.5)-part 3
    nest(j,:)=(X1+X2+X3)/3;% Equation (3.7)        %     pause
    s=nest(j,:);
end
for j=1:n
    % Apply simple bounds/limits
    s=nest(j,:);
    nest(j,:)=simplebounds(s,Lb,Ub);
end

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(best,nest,Lb,Ub,pa, R, I)
n=size(nest,1);
%--------------------------------------------------------------
% ----- Red Panda Optimization ----- %%%
%---------------------------------------------------------------
for i=1:n
    K=rand(size(nest))>pa;
    stepsize= R *(nest(randperm(n),:))- I * (nest(randperm(n),:));
    new_nest=nest+stepsize.*K;
end
for j=1:size(new_nest,1)
    s=new_nest(j,:);
    new_nest(j,:)=simplebounds(s,Lb,Ub);
end
%--------------------------------------------------------------
% ----- Meerkat Optimization Algorithm ----- %%%
%---------------------------------------------------------------
function new_nest=new_empty_nests(nest,Lb,Ub,pa,lambda) %NMRA
n=size(nest,1);
K=rand(size(nest))>pa;
stepsize= nest(randperm(n),:) + (2.* rand -1 * (nest(randperm(n),:))+ rand * lambda);
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
    new_nest(j,:)=simplebounds(s,Lb,Ub);
end

%--------------------------------------------------------------
% ----- Keplers Optimization Algorithm ----- %%%
%---------------------------------------------------------------
function nest=get_different_cuckoos(nest,best,Lb,Ub, Q)
%% Computer Q that represents the Euclidian distance between the best-so-far solution and the ith solution
n=size(nest,1);
%% A randomly-assigned binary vector
rd = rand();
r = rand();
M=(rand.*(1-r)+r); %% Eq.(16)
U1=rd<r; %% Eq.(21)
sum=0;
for i=1:n
    sum=sum+(nest(i,:)-max(nest(i,:)));
    Q=sqrt((min(nest(i,:)) - nest(i,:)).^2); %% Eq.(7) of KOA
    MS=rand*(best-max(nest(i,:)))/(sum); %% Eq.(8)
    m=(nest(i,:)-max(nest(i,:)))/(sum); %% Eq.(9)
end



%K=rand(size(nest))>pa;
%beta=3/2;
%sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
for j=1:n,
    %% Step 3: Calculating an object velocity
    % A flag to opposite or leave the search direction of the current planet
    if rand<0.5 %% Eq.(18)
        f=1;
    else
        f=-1;
    end

    Rnorm=(Q-min(Q))/(max(Q)-min(Q)); %% The normalized R (Eq.(24))
    MSnorm=(MS-min(MS))/(max(MS)-min(MS)); %% The normalized MS
    Mnorm=(m-min(m))/(max(m)-min(m)); %% The normalized m
    MSnorm = min(MSnorm);
    Mnorm = min(Mnorm);
    Rnorm= min(Rnorm);
    Fg=rand().*M.*((MSnorm.*Mnorm)/(Rnorm.*Rnorm+eps))+(rand); %% Eq.(6)
    L=M*(MS+m).*abs((2/(Rnorm+eps))-(1/(rand()+eps)))^(0.5); %% Eq.(15)
    U=rd>rand(); %% A binary vector
    if Rnorm<0.5 %% Eq.(13)
        l=L*M*U; %% Eq.(14)
        Mv=(rand*(1-rd)+rd); %% Eq.(20)
        l1=L.*Mv.*(1-U);%% Eq.(19)
        V=l.*(2*rand*nest(j,:)-nest(randperm(n),:))+l1.*(nest(randperm(n),:)-nest(randperm(n),:))+(1-Rnorm)*f*U1.*rand().*(Ub-Lb); %% Eq.(13a)
    else
        U2=rand>rand; %% Eq. (22)
        V=rand.*L.*(nest(randperm(n),:)-nest(j,:))+(1-Rnorm)*f*U2*rand().*(rand*Ub-Lb);  %% Eq.(13b)
    end %% End IF

    %% Step 4: Escaping from the local optimum
    % Update the flag f to opposite or leave the search direction of the current planet
    if rand<0.5 %% Eq.(18)
        f=1;
    else
        f=-1;
    end
    %% Step 5
    nest(j,:)=((nest(j,:)+V(j,:).*f)+(Fg+abs(randn))*U.*(best-nest(j,:))); %% Eq.(25)



    % s=nest(j,:);
    % u=randn(size(s))*sigma;
    % v=randn(size(s));
    % step=u./abs(v).^(1/beta);
    % stepsize=0.01*step.*(s-best);
    % s=s+stepsize.*randn(size(s));
    % nest(j,:)=simplebounds(s,Lb,Ub);
end
% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;

%% You can replace the following by your own functions
% A d-dimensional objective function
% function z=rosenbrock(x)
% D= 30;
% sum1 = 0;
%         for i = 1:1:D
%             sum1 = sum1 + i*x(i)^4;
%         end
%         z = sum1 + rand();