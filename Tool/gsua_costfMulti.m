function cos=gsua_costfMulti(ydata,yfunction,margin,parallel)
if nargin<3
    margin=0.1;
    parallel=false;
end
if parallel
    parforArg=Inf;
else
    parforArg=0;
end
margin=1+margin;

[reps,len,inputs]=size(ydata);
if size(yfunction,3)>1
    yfunction=squeeze(yfunction)';
end
regulator=sum((yfunction-yfunction*margin).^2,2)/len;
cos=zeros(1,reps);

parfor (j=1:reps,parforArg) 
%for j=1:reps
    cost=zeros(1,inputs);
    for i=1:inputs
        icorre=(2-corr2(ydata(j,:,i),yfunction(i,:)));
        if isnan(icorre)
            icorre=3;
        end
        cost(i)=(icorre*(sum((ydata(j,:,i)-yfunction(i,:)).^2,2)/len)/regulator(i))^2;
    end
    cos(j)=sum(cost)/inputs;
%     if isnan(cos(j))
%         disp('fail')
%     end
end
end