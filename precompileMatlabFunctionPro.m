function precompileMatlabFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileMatlabFunctionPro generates a m-file with 'matlabFunction()' 
%   in the bin directory from a SymbolicMatrix for a runtime execution.
%
%   precompileMatlabFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated file NameMF.m
%       SymbolicMatrix: symbolic matrix
%       Params: cell array of vectors of symbolic elements.
%               This elements are contained in the SymbolicMatrix.
%               This symbolic elements will be replaced with the elements in SymbolicList.
%           example:
%           {[q1,q2,q3],[dq1,dq2,dq3]} 
%           will be replaced with
%           q(1),q(2),q(3),dq(1),dq(2)....
%
%       SymbolList: cell array of string that contains the name for the
%                   symbolic variables
%           example:
%           {'q','dq'}
%
%	example:
%
% 	syms t1 t2 t3 d1 d2;
% 	M=[sin(t1)*sin(t2+t3)-d1*d2,d1+d2;t1*t2*t3,sin(t1+t2+t3)-d1*d2];
% 	mkdir('bin');
% 	precompileMatlabFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a file MatrixMF.m
%   
%   addpath('bin');
%   MatrixMF([1,2,3],[2,3]);
%   ans =
%
%       -6.8069    5.0000
%        6.0000   -6.2794  
%
%   PhD Gastone Pietro Rosati Papini
%   Scuola Superiore Sant'Anna Pisa 
%   Percro Laboratory
%   $Revision: 1.0 $  $Date: 2014/03/25 15:23:00 $

    nlist = length(Params);
    if(nlist~=length(SymbolList))
         exception = MException('MATLAB:InconsistentDataType','SymbolList and Params number are different');
         throw(exception);
    end
    cd('bin');
    fprintf(['Build ',Name,' with MatlabFunction to M-Function']);
    [s1,s2]=size(SymbolicMatrix);
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('.');
            matlabFunction(SymbolicMatrix(indr,indc), 'file', [Name,'_',num2str(indr),'_',num2str(indc),'_MF'],'vars',Params);
        end
    end
    fun_name = strcat(Name,'MF.m');
    fid = fopen(fun_name, 'wt');
    args = '';
    for i = 1 : nlist
        args = strcat(args,SymbolList{i});
        if (i ~= nlist)
            args = strcat(args,',');
        end
    end
    str = strcat('function M = ',Name,'MF(',args,')\n');
    fprintf(fid, str);
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('-');
            var_mat = strcat('\n\tM(',num2str(indr),',',num2str(indc),')=',Name,'_',num2str(indr),'_',num2str(indc),'_MF(',args,');');
            fprintf(fid, var_mat);
        end
    end
    fprintf(fid,'\n\nend');
    fclose(fid);
    fprintf(['\nBuilt: ',Name,'MF.m\n']);
    cd('..');
end