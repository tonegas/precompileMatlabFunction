function precompileCFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileCFunctionPro generates a C++ source MEX-Files in the bin directory from a 
%   SymbolicMatrix than compile it for a runtime execution.
%
%   precompileCFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated file NameC.cpp
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
% 	precompileCFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a file MatrixC.cpp and MatrixC.mexw64
%   
%   addpath('bin');
%   M=MatrixC([1,2,3],[2,3]);
%   M =
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
    
    fprintf(['Build ',Name,' to Mex-Function']);
    fun_name = strcat('bin/',Name,'C.cpp');
    fid = fopen(fun_name, 'wt');
    
    [s1,s2]=size(SymbolicMatrix);
    lunghezza{nlist}=[];
    costrue{nlist}=[];
    sintrue{nlist}=[];
    max=0;
    for indlist=0:nlist-1
        lunghezza{indlist+1} = length(Params{indlist+1});     %Numero elementi
        if lunghezza{indlist+1}>max
            max=lunghezza{indlist+1};
        end
    end
    
    fprintf(fid,'#include "mex.h"\n');
    fprintf(fid,'#include <math.h>\n\n\n');
       
    fprintf(fid,'void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){\n\n');
    
    fprintf(fid,['\tif(nrhs!=',num2str(nlist),'){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","',num2str(nlist),' inputs required.");}\n']);
    fprintf(fid,'\tif(nlhs!=1){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 outputs required.");}\n\n');
    
    if max > 1
        powsintrue{nlist,max-1}=[];
        powcostrue{nlist,max-1}=[];
    end
    
    for indlist=0:nlist-1
        lunghezza{indlist+1} = length(Params{indlist+1});     %Numero elementi
        costrue{indlist+1} = zeros(1,lunghezza{indlist+1});
        sintrue{indlist+1} = zeros(1,lunghezza{indlist+1});
        fprintf(fid,'\tif(mxGetM(prhs[%d])!=1){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","The dimensions of %s must be 1x%d");}\n',indlist,SymbolList{indlist+1},lunghezza{indlist+1});
        fprintf(fid,'\tif(mxGetN(prhs[%d])!=%d){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","The dimensions of %s must be 1x%d");}\n',indlist,lunghezza{indlist+1},SymbolList{indlist+1},lunghezza{indlist+1});
        for pot = 2:lunghezza{indlist+1}
            powcostrue{indlist+1,pot-1} = zeros(1,lunghezza{indlist+1});
            powsintrue{indlist+1,pot-1} = zeros(1,lunghezza{indlist+1});
        end
    end
    
    fprintf(fid,['\n\tplhs[0] = mxCreateDoubleMatrix(',num2str(s1),',',num2str(s2),',mxREAL);\n']);
    fprintf(fid,'\tdouble\t*M = mxGetPr(plhs[0]);\n\n');
    
    for indlist=0:nlist-1
        fprintf(fid,['\tdouble\t*',SymbolList{indlist+1},' = mxGetPr(prhs[',num2str(indlist),']);\n']);
    end
    fprintf(fid,'\n');

    indvar=1;
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('.');
            M = ccode(vpa(SymbolicMatrix(indr,indc),10));
            M = strrep(M,'~','');
            M = strrep(M,'t0 = ','');
            while (isempty(strfind(M, strcat('MapleGenVar',num2str(indvar))))==0)
                fprintf(fid,strcat('\tdouble\tMapleGenVar',num2str(indvar),';\n'));
                indvar=indvar+1;
            end
            for indlist=1:nlist
                for x = 1:lunghezza{indlist}
                    M = strrep(M,char(Params{indlist}(x)),strcat(SymbolList{indlist},'[',num2str(x-1),']'));
                    if(costrue{indlist}(x)==0)
                        if(isempty(strfind(M, strcat('cos(',SymbolList{indlist},'[',num2str(x-1),'])')))==0)
                            costrue{indlist}(x)=1;
                            fprintf(fid,strcat('\tdouble\tc',num2str(SymbolList{indlist}),num2str(x),...
                                '=cos(',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                        end
                    end

                    if(sintrue{indlist}(x)==0)
                        if(isempty(strfind(M, strcat('sin(',SymbolList{indlist},'[',num2str(x-1),'])')))==0)
                            sintrue{indlist}(x)=1;
                            fprintf(fid,strcat('\tdouble\ts',num2str(SymbolList{indlist}),num2str(x),...
                                '=sin(',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                            M = strrep(M,strcat('sin(',SymbolList{indlist},'[',num2str(x-1),'])'),...
                                strcat('s',num2str(SymbolList{indlist}),num2str(x)));
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(powcostrue{indlist,pot-1}(x)==0)
                            if(isempty(strfind(M, strcat('pow(cos(',SymbolList{indlist},'[',num2str(x-1),']),',num2str(pot),'.0)')))==0)
                                powcostrue{indlist,pot-1}(x)=1;
                                fprintf(fid,strcat('\tdouble\tpc',num2str(SymbolList{indlist}),...
                                    num2str(x),'_',num2str(pot),'=pow(c',SymbolList{indlist},num2str(x),',',num2str(pot),'.0);\n'));
                                M = strrep(M,strcat('pow(c',SymbolList{indlist},num2str(x),',',num2str(pot),'.0)'),...
                                    strcat('pc',SymbolList{indlist},num2str(x),'_',num2str(pot)));
                            end
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(powsintrue{indlist,pot-1}(x)==0)
                            if(isempty(strfind(M, strcat('pow(sin(',SymbolList{indlist},'[',num2str(x-1),']),',num2str(pot),'.0)')))==0)
                                powsintrue{indlist,pot-1}(x)=1;
                                fprintf(fid,strcat('\tdouble\tps',num2str(SymbolList{indlist}),...
                                    num2str(x),'_',num2str(pot),'=pow(s',SymbolList{indlist},num2str(x),',',num2str(pot),'.0);\n'));
                                M = strrep(M,strcat('pow(s',SymbolList{indlist},num2str(x),',',num2str(pot),'.0)'),...
                                    strcat('ps',SymbolList{indlist},num2str(x),'_',num2str(pot)));
                            end
                        end
                    end
                    if(costrue{indlist}(x)==1)
                        M = strrep(M,strcat('cos(',SymbolList{indlist},'[',num2str(x-1),'])'),...
                            strcat('c',num2str(SymbolList{indlist}),num2str(x)));
                    end
                    if(sintrue{indlist}(x)==1)
                        M = strrep(M,strcat('sin(',SymbolList{indlist},'[',num2str(x-1),'])'),...
                            strcat('s',num2str(SymbolList{indlist}),num2str(x)));
                    end
                    for pot = 2:lunghezza{indlist}
                        if(powcostrue{indlist,pot-1}(x)==1)
                            M = strrep(M,strcat('pow(c',SymbolList{indlist},num2str(x),',',num2str(pot),'.0)'),...
                                strcat('pc',SymbolList{indlist},num2str(x),'_',num2str(pot)));
                        end
                    end
                    for pot = 2:lunghezza{indlist}
                        if(powsintrue{indlist,pot-1}(x)==1)
                            M = strrep(M,strcat('pow(s',SymbolList{indlist},num2str(x),',',num2str(pot),'.0)'),...
                                strcat('ps',SymbolList{indlist},num2str(x),'_',num2str(pot)));
                        end
                    end
                end
            end
            M = strcat('\tM[',num2str((indr-1)+(indc-1)*s1),'] =',M,'\n');
            fprintf(fid, M);
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'}\n\n\n');
    fclose(fid);
    
    cd('bin');
    Vai = [Name,'C.cpp'];
    fprintf('\n*\n');
    eval(['mex ',Vai]);
    cd('..');
    fprintf(['Built: ',Name,'C.mexw64\n']);

end
