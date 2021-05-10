function T=gsua_load(T)

prop=T.Properties.CustomProperties;
func=prop.creator.func;
Range=prop.creator.Range;
try
    names=prop.creator.names;
catch
    names=[];
end
nominal=prop.creator.nominal;
vectorized=prop.creator.vectorized;
domain=prop.Domain;
rMethod=prop.rMethod;
output=prop.output;
vars=prop.Vars;
opt=prop.copt;



T=gsua_dataprep(func,Range,'names',names,'domain',domain,'rMethod',rMethod,...
    'nominal',nominal,'output',output,'out_names',vars,'vectorized',vectorized,...
    'opt',opt);


end