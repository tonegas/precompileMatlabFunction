function precompileDistCFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileDistCFunctionPro generates NxM cpp-files and one C++ source MEX-Files 
%   in the bin directory from a SymbolicMatrix(NxM) then compiles all files in 
%   one it for a runtime execution.
%
%   precompileDistCFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated files Name_1_1_C.cpp,
%             Name_1_2_C.cpp,... (each for each elements) and NameC.cpp, extern_NameC.h
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
% 	precompileDistCFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a files Matrix_1_1_C.cpp, Matrix_1_2_C.cpp, MatrixC.cpp, extern_MatrixC.h and MatrixC.mexw64
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
         exception = MException('MATLAB:InconsistentDataType','ListSimboli and Params number are different');
         throw(exception);
    end
    
    fprintf(['Build Elements of ',Name,' to Mex-Function']);
    fun_name = strcat('bin/',Name,'C.cpp');
    fid = fopen(fun_name, 'wt');
    
    [s1,s2]=size(SymbolicMatrix);
    
    lunghezza{nlist}=[];
    costrue{nlist}=[];
    sintrue{nlist}=[];
    call_fun{s1,s2}=[];
    max=0;
    for indlist=0:nlist-1
        lunghezza{indlist+1} = length(Params{indlist+1});     %Numero elementi
        if lunghezza{indlist+1}>max
            max=lunghezza{indlist+1};
        end
    end
    
    extern_name = strcat('bin/extern_',Name,'C.h');
    fid_extern = fopen(extern_name, 'wt');
    
    fprintf(fid,'#include "mex.h"\n');
    fprintf(fid,strcat('#include "extern_',Name,'C.h"\n'));
    fprintf(fid,'#include <math.h>\n\n\n');
       
    fprintf(fid,'void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){\n\n');
    
    fprintf(fid,['\tif(nrhs!=',num2str(nlist),'){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","',num2str(nlist),' inputs required.");}\n']);
    fprintf(fid,'\tif(nlhs!=1){mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 outputs required.");}\n\n');   
    
    powsintrue{nlist,max-1}=[];
    powcostrue{nlist,max-1}=[];
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

    fprintf(fid,['\tplhs[0] = mxCreateDoubleMatrix(',num2str(s1),',',num2str(s2),',mxREAL);\n']);
    fprintf(fid,'\tdouble\t*M = mxGetPr(plhs[0]);\n\n');
    
    for indlist=0:nlist-1
        fprintf(fid,['\tdouble\t*',SymbolList{indlist+1},' = mxGetPr(prhs[',num2str(indlist),']);\n']);
    end
    fprintf(fid,'\n');

    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            cd('bin');
            mat_name = strcat(Name,'_',num2str(indr),'_',num2str(indc),'_C.cpp');
            fid_mat = fopen(mat_name, 'wt');
            fprintf('.');
            M = ccode(vpa(SymbolicMatrix(indr,indc),10));
            M = strrep(M,'~','');
            M = strrep(M,'t0 = ','');
            
            fprintf(fid_mat,'#include <math.h>\n\n');
            name_fun = strcat('double\t',Name,'_',num2str(indr),'_',num2str(indc),'_(');
            call_fun{indr,indc} = strcat(Name,'_',num2str(indr),'_',num2str(indc),'_(');
            
            name_fun=[name_fun,'double &',char(Params{1}(1))];
            call_fun{indr,indc}=strcat(call_fun{indr,indc},strcat(SymbolList{1},'[0]'));
            for x = 2:lunghezza{1}
            	name_fun=[name_fun,',double &',char(Params{1}(x))];
                call_fun{indr,indc}=strcat(call_fun{indr,indc},',',strcat(SymbolList{1},'[',num2str(x-1),']'));
            end
            
            for indlist=2:nlist
                for x = 1:lunghezza{indlist}
                    name_fun=[name_fun,',double &',char(Params{indlist}(x))];
                    call_fun{indr,indc}=strcat(call_fun{indr,indc},',',strcat(SymbolList{indlist},'[',num2str(x-1),']'));
                end
            end
            
            for indlist=1:nlist
                for x = 1:lunghezza{indlist}
                    if(isempty(strfind(M, strcat('cos(',char(Params{indlist}(x)),')')))==0)
                        name_fun=strcat(name_fun,',double &c',char(Params{indlist}(x)));
                        call_fun{indr,indc}=strcat(call_fun{indr,indc},',c',char(Params{indlist}(x)));
                        M = strrep(M,strcat('cos(',char(Params{indlist}(x)),')'),...
                            strcat('c',char(Params{indlist}(x))));
                        if(costrue{indlist}(x)==0)
                            fprintf(fid,strcat('\tdouble\tc',char(Params{indlist}(x)),...
                                '=cos(',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                            costrue{indlist}(x)=1;
                        end
                    end
                    
                    if(isempty(strfind(M, strcat('sin(',char(Params{indlist}(x)),')')))==0)
                        name_fun=strcat(name_fun,',double &s',char(Params{indlist}(x)));
                        call_fun{indr,indc}=strcat(call_fun{indr,indc},',s',char(Params{indlist}(x)));
                        M = strrep(M,strcat('sin(',char(Params{indlist}(x)),')'),...
                            strcat('s',char(Params{indlist}(x))));
                        if(sintrue{indlist}(x)==0)
                            fprintf(fid,strcat('\tdouble\ts',char(Params{indlist}(x)),...
                                '=sin(',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                            sintrue{indlist}(x)=1;
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(isempty(strfind(M, strcat('pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0)')))==0)
                            name_fun=strcat(name_fun,',double &pc',char(Params{indlist}(x)),'_',num2str(pot));
                            call_fun{indr,indc}=strcat(call_fun{indr,indc},',pc',char(Params{indlist}(x)),'_',num2str(pot));
                            M = strrep(M,strcat('pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0)'),...
                                strcat('pc',char(Params{indlist}(x)),'_',num2str(pot)));
                            if(powcostrue{indlist,pot-1}(x)==0)
                                fprintf(fid,strcat('\tdouble\tpc',char(Params{indlist}(x)),...
                                    '_',num2str(pot),'=pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0);\n'));
                                powcostrue{indlist,pot-1}(x)=1;
                            end
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(isempty(strfind(M, strcat('pow(s',char(Params{indlist}(x)),',',num2str(pot),'.0)')))==0)
                            name_fun=strcat(name_fun,',double &ps',char(Params{indlist}(x)),'_',num2str(pot));
                            call_fun{indr,indc}=strcat(call_fun{indr,indc},',ps',char(Params{indlist}(x)),'_',num2str(pot));
                            M = strrep(M,strcat('pow(s',char(Params{indlist}(x)),',',num2str(pot),'.0)'),...
                                strcat('ps',char(Params{indlist}(x)),'_',num2str(pot)));
                            if(powsintrue{indlist,pot-1}(x)==0)
                                fprintf(fid,strcat('\tdouble\tps',char(Params{indlist}(x)),...
                                    '_',num2str(pot),'=pow(s',char(Params{indlist}(x)),',',num2str(pot),'.0);\n'));
                                powsintrue{indlist,pot-1}(x)=1;
                            end
                        end
                    end
                end
            end
            
            name_fun=strcat(name_fun,')');
            call_fun{indr,indc}=strcat(call_fun{indr,indc},');\n');
            
            fprintf(fid_extern,['extern ',name_fun,';\n']);
            fprintf(fid_mat,[name_fun,'{\n']);
            
            indvar=1;
            while (isempty(strfind(M, strcat('MapleGenVar',num2str(indvar))))==0)
                fprintf(fid_mat,strcat('\tdouble\tMapleGenVar',num2str(indvar),';\n'));
                indvar=indvar+1;
            end
            M = strcat('\treturn ',M,'\n}\n');
            fprintf(fid_mat, M);
            fclose(fid_mat);
            
            fprintf('*');
            eval(strcat('mex -c ./',mat_name));
            cd('..');
        end
    end
    
    fprintf(fid,'\n');
    file = [Name,'C.cpp '];
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('-');
            var_mat = strcat('\tM[',num2str((indr-1)+(indc-1)*s1),']=',call_fun{indr,indc});
            file = [file,Name,'_',num2str(indr),'_',num2str(indc),'_C.obj '];
            fprintf(fid, var_mat);
        end
    end
    fprintf(fid,'}\n\n\n');
    fclose(fid);
    fclose(fid_extern);
    %fprintf(file);
    
    cd('bin');
    fprintf('\n*\n');
    eval(['mex ',file]);
    cd('..');
    fprintf(['Built: ',Name,'C.mexw64\n']);

end