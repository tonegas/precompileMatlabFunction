function precompileMFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileMFunctionPro generates a m-file in the bin directory from a 
%   SymbolicMatrix for a runtime execution.
%
%   precompileMFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated file NameM.m
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
% 	precompileMFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a file MatrixM.m
%   
%   addpath('bin');
%   MatrixM([1,2,3],[2,3]);
%   ans =
%
%       -6.8069    5.0000
%        6.0000   -6.2794  
%
%   PhD Gastone Pietro Rosati Papini
%   Scuola Superiore Sant'Anna Pisa 
%   Percro Laboratory
%   $Revision: 1.0 $  $Date: 2014/03/25 15:23:00 $

    nlist = length(Params);     %Numero di parametri
 
    if(nlist~=length(SymbolList))
         exception = MException('MATLAB:InconsistentDataType','SymbolList and Params number are different');
         throw(exception);
    end

    %Converto la matrice simbolica in una stringa
    fprintf(['Build ',Name,' to M-Function']);
    M = char(vpa(SymbolicMatrix,10));
    M = strrep(M,'matrix([','\t\t');
    M = strrep(M,']])',']');
    
    %Nome del file su cui scrivere la funzione precompilata
    fun_name = strcat('bin/',Name,'M.m');
    fid = fopen(fun_name, 'wt');
    
    %Argomenti della funzione precompilata
    args = '';
    for i = 1 : nlist
        args = strcat(args,SymbolList{i});
        if (i ~= nlist)
            args = strcat(args,',');
        end
    end
    
    str = strcat('function M = ',Name,'M(',args,')\n\n');
    
    for indlist=1:nlist
        lunghezza = length(Params{indlist});     %Numero elementi
        
        str=strcat(str,['assert(length(',SymbolList{indlist},')==',num2str(lunghezza),',''','Elements in ',SymbolList{indlist},' must be ',num2str(lunghezza),''')\n']);

        %Converto i parametri in formato par(x)
        for x = 1:lunghezza
            fprintf('.');
            M = strrep(M,char(Params{indlist}(x)),strcat(SymbolList{indlist},'(',char(sym(x)),')'));
        end
    end
    M = strrep(M,'], [',';\n\t\t ');
    M = strrep(M,'],[',';\n\t\t ');
    
    %Costruisco la funzione precompilata
    fprintf('-');
    str = strcat(str,'\nM = ');
    str = strcat(str,M);
    str = strcat(str,';\n\nend');

    fprintf(fid, str);
    fclose(fid);
    fprintf(['Built: ',Name,'M.m\n']);
    
end