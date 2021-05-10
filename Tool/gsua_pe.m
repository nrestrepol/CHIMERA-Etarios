function [Table,res,IP] = gsua_pe(Table,xdata,ydata,varargin)
% Parameter estimation function
% [T,res] = gsua_pe(T,xdata,ydata,N)
% Parameters:
% T     <-- summary table from gsua_dataprep function
% xdata <-- array of points where the model will be evaluated
% ydata <-- array with expected model output
% Outputs:
% T   <-- summary table with parameter estimation results
% res <-- cost functions for each estimation
% Additional paired features:
% 'N',N     <-- number of parameter estimations
% 'Multistart',k <-- activate multistart feature for lsqcurvefit optimizer,
% perform k parameter estimations in each cycle
% 'solver',{'lsqc','lsqn','ga','particle','psearch','surrogate','annealing'}
% <-- allows to choose among several matlab optimizers. Default: lsqc
% 'opt',optimoptions(optimizer,...) <-- allows to configure the respective
% matlab optimizer, optimizer must match with the real matlab optimizer
% name.
% 'ipoint',point <-- allows to give to the optimizer the initial points for
% estimations. point must have the same number of columns as factors to
% estimate and the same number of rows as N.
% 'Show',{'off','on'}. If activated, the function plot optimizer results
% Default: off.
% [T,res] = gsua_pe(T,xdata,ydata,'Show','on','solver','particle','N',3)



p=inputParser;

dfN=1;
dfMS=1;
dfIP=[];
dfOS='lsq';
dfOpt=[];

validOS={'lsqc','lsqn','ga','particle','psearch','surrogate','annealing','fmincon'};
checkOS = @(x) any(validatestring(x,validOS));
checkIP = @(x) (isnumeric(x) || isempty(x));
defaultShow='off';
validShow={'off' 'on' 'on2'};
checkShow = @(x) any(validatestring(x,validShow));
dfmarg=0.1;
dfSave=true;

addRequired(p,'Table',@istable);
addRequired(p,'xdata',@isnumeric);
addRequired(p,'ydata',@isnumeric);
addParameter(p,'N',dfN,@isnumeric);
addParameter(p,'Multistart',dfMS,@isnumeric);
addParameter(p,'ipoint',dfIP,checkIP);
addParameter(p,'solver',dfOS,checkOS);
addParameter(p,'opt',dfOpt);
addParameter(p,'Show',defaultShow,checkShow);
addParameter(p,'A',dfOpt);
addParameter(p,'B',dfOpt);
addParameter(p,'Aeq',dfOpt);
addParameter(p,'Beq',dfOpt);
addParameter(p,'nonlcon',dfOpt);
addParameter(p,'margin',dfmarg,@isnumeric);
addParameter(p,'save',dfSave,@islogical);

parse(p,Table,xdata,ydata,varargin{:})
T=p.Results.Table;
xdata=p.Results.xdata;
ydata=p.Results.ydata;
N=p.Results.N;
MS=p.Results.Multistart;
IP=p.Results.ipoint;
OS=p.Results.solver;
Opt=p.Results.opt;
show=p.Results.Show;
A=p.Results.A;
B=p.Results.B;
Aeq=p.Results.Aeq;
Beq=p.Results.Beq;
nonlcon=p.Results.nonlcon;
margin=p.Results.margin+1;
saver=p.Results.save;

lb=T.Range(:,1);
ub=T.Range(:,2);
Np=size(T,1);

if MS>1
    OS='Multistart';
    disp('only lsqc supports Multistart ... Running lsq optimization')
end
if (isempty(IP) || size(IP,1)~=N) && any(strcmp(OS,{'lsqc','Multistart','lsqn','fmincon'}))
    disp('Generating a valid matrix for estimations')
    [IP,~]=gsua_dmatrix(T,N);
end

x=zeros(N,Np);
res=zeros(N,1);
% try
%     TP=T.Properties.CustomProperties;   
% catch
%     TP=load('ATable.mat');
%     TP=TP.Table2;
% end
% if any(strcmp(OS,{'lsqc','lsqn'})) && N>1
%     
%     parfor i=1:N
%         try
%             switch OS
%                 case 'lsqc'
%                      [x(i,:),res(i)] = lsqcurvefit(@(pars,xaxis) gsua_deval(pars,TP,xaxis),IP(i,:),xdata,ydata,lb,ub,Opt);
%                 case 'lsqn'
%                     obj=@(pars) gsua_deval(pars,TP,xdata)-ydata;  
%                     [x(i,:),res(i)] = lsqnonlin(obj,IP(i,:),lb,ub,Opt);
%             end
%         catch ME
%             warning(ME.message)  
%         
%         end
%     end    
% else
prog=N>1;
if prog
    progressbar
end
    for i=1:N
        disp(['Estimation ' num2str(i)])
        done=true;
        attempt=1;
    while done
    try
    switch OS
        case 'Multistart'
            problem = createOptimProblem(OS,'x0',IP(i,:),'objective',@(pars,xaxis) gsua_deval(pars,T,xaxis),...
            'lb',lb,'ub',ub,'xdata',xdata,'ydata',ydata);
            ms = MultiStart(UseParallel',1);
            [x(i,:),res(i)] = run(ms,problem,multi);           
        case 'lsqc'
            [x(i,:),res(i)] = lsqcurvefit(@(pars,xaxis) gsua_deval(pars,T,xaxis),IP(i,:),xdata,ydata,lb,ub,Opt);
        case 'lsqn' 
             obj=@(pars) (gsua_deval(pars,T,xdata)-ydata)/std(ydata);   
             [x(i,:),res(i)] = lsqnonlin(obj,IP(i,:),lb,ub,Opt);
        case {'ga','particle','psearch','surrogate','annealing','fmincon'}
             inputs=size(ydata,1);
             len=length(xdata);
             regulator=sum((ydata-ydata*margin).^2,2)/len;
             obj=@(pars) gsua_costf(inputs,regulator,len,ydata,gsua_deval(pars,T,xdata));              
            switch OS
                case 'ga'
                    [x(i,:),res(i)] = ga(obj,Np,A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
                case 'particle'
                    [x(i,:),res(i)] = particleswarm(obj,Np,lb,ub,Opt);    
                case 'psearch'
                    [x(i,:),res(i)] = patternsearch(obj,Np,A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
                case 'surrogate'
                    [x(i,:),res(i)] = surrogateopt(obj,lb,ub,Opt);
                case 'annealing'
                    [x(i,:),res(i)] = simulannealbnd(obj,lb,ub,Opt);
                case 'fmincon'
                    [x(i,:),res(i)] = fmincon(obj,IP(i,:),A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
            end
    end
    done=false;
    if saver
        save('Estimations.mat','x','res');
    end
    catch ME
        switch ME.identifier
            case 'optimlib:lsqcurvefit:YdataSizeMismatchFunVal'
            if attempt==1
                ydata=ydata';
                attempt=2;
            else
                throw(ME)
            end
            case 'MATLAB:warning:invalidMessageType'
                throw(ME)
            otherwise
                warning(ME.message)
                IP(i,:)=gsua_dmatrix(T,1); 
        end
    end
    end
    if prog
        progressbar(i/N)
    end
    end
% end


[res,index] = sort(res);
%res=res(res ~= 0);
res=res';
x=x(index(1:size(res,2)),:);
eval(strcat('Table.Est',OS,"=x';"));
if strcmp(show,'on')||strcmp(show,'on2')
    figure(1);
    clf
    [~]=gsua_eval(x',Table,xdata,ydata);
    if strcmp(show,'on2')
        gsua_ia(Table,x')
    end
end
end