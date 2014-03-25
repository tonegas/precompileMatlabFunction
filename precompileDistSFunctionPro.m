function precompileDistSFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileDistSFunctionPro generates NxM+1 S-function cpp-files in the bin directory from a 
%   SymbolicMatrix(NxM) then compiles all files in one it for a runtime execution with S-function.
%
%   precompileDistSFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated files Name_1_1_S.cpp,
%             Name_1_2_S.cpp,... (each for each elements) and NameS.cpp, extern_NameS.h
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
% 	precompileDistSFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a files Matrix_1_1_S.cpp, Matrix_1_2_S.cpp, MatrixS.cpp, extern_Matrix.h and MatrixS.mexw64
%   use it in simulink with S-function.
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
    
    fprintf(['Build Elements of ',Name,' to S-Function']);
    fun_name = strcat('bin/',Name,'S.cpp');
    fid = fopen(fun_name, 'wt');
    
    [s1,s2]=size(SymbolicMatrix);
    
    nlist = length(Params);
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
    
    extern_name = strcat('bin/extern_',Name,'S.h');
    fid_extern = fopen(extern_name, 'wt');
    
    def_name = strcat('#define S_FUNCTION_NAME\t',Name,'S\n');
    fprintf(fid,def_name);
    fprintf(fid,'#define S_FUNCTION_LEVEL 2\n');
    fprintf(fid,'#include "simstruc.h"\n');
    fprintf(fid,strcat('#include "extern_',Name,'S.h"\n'));
    fprintf(fid,'#include <math.h>\n\n\n');
       
    fprintf(fid,'static void mdlInitializeSizes(SimStruct *S)\n{\n');
    fprintf(fid,'\tssSetNumSFcnParams(S, 0);\n');
    fprintf(fid,'\tif (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))\n\t{\n');
    fprintf(fid,'\t\treturn; /* Parameter mismatch will be reported by Simulink */\n\t}\n\n');

    fprintf(fid,strcat('\tif (!ssSetNumInputPorts(S,',num2str(nlist),')) return;\n'));
    
    powsintrue{nlist,max-1}=[];
    powcostrue{nlist,max-1}=[];
    for indlist=0:nlist-1
        lunghezza{indlist+1} = length(Params{indlist+1});     %Numero elementi
        fprintf(fid,strcat('\tssSetInputPortWidth(S,',num2str(indlist),',',num2str(lunghezza{indlist+1}),');\n'));
        fprintf(fid,strcat('\tssSetInputPortDirectFeedThrough(S,',num2str(indlist),',1);\n\n'));
        costrue{indlist+1} = zeros(1,lunghezza{indlist+1});
        sintrue{indlist+1} = zeros(1,lunghezza{indlist+1});
        for pot = 2:lunghezza{indlist+1}
            powcostrue{indlist+1,pot-1} = zeros(1,lunghezza{indlist+1});
            powsintrue{indlist+1,pot-1} = zeros(1,lunghezza{indlist+1});
        end
    end
 
    fprintf(fid,'\tif (!ssSetNumOutputPorts(S,1)) return;\n');
    fprintf(fid,strcat('\tssSetOutputPortMatrixDimensions(S,0,',num2str(s1),',',num2str(s2),');\n\n'));

    fprintf(fid,'\tssSetNumSampleTimes(S, 1);\n\n');

    fprintf(fid,'\t/* Take care when specifying exception free code - see sfuntmpl.doc */\n');
    fprintf(fid,'\tssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);\n}\n\n\n');

    fprintf(fid,'static void mdlInitializeSampleTimes(SimStruct *S)\n{\n');
    fprintf(fid,'\tssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);\n');
    fprintf(fid,'\tssSetOffsetTime(S, 0, 0.0);\n}\n\n\n');
    
    fprintf(fid,'static void mdlOutputs(SimStruct *S, int_T tid)\n{\n');
    fprintf(fid,'\treal_T\t*M\t=\tssGetOutputPortRealSignal(S,0);\n\n');
    
    for indlist=0:nlist-1
        fprintf(fid,strcat('\tInputRealPtrsType\t',SymbolList{indlist+1},'\t=\tssGetInputPortRealSignalPtrs(S,',num2str(indlist),');\n'));
    end
    fprintf(fid,'\n');

    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('.');
            cd('bin');
            mat_name = strcat(Name,'_',num2str(indr),'_',num2str(indc),'_S.cpp');
            fid_mat = fopen(mat_name, 'wt');
            M = ccode(vpa(SymbolicMatrix(indr,indc),10));
            M = strrep(M,'~','');
            M = strrep(M,'t0 = ','');
            
            fprintf(fid_mat,'#include "simstruc.h"\n');
            fprintf(fid_mat,'#include <math.h>\n\n');
            name_fun = strcat('real_T\t',Name,'_',num2str(indr),'_',num2str(indc),'_(');
            call_fun{indr,indc} = strcat(Name,'_',num2str(indr),'_',num2str(indc),'_(');
            
            name_fun=[name_fun,'real_T\t&',char(Params{1}(1))];
            call_fun{indr,indc}=strcat(call_fun{indr,indc},strcat('(real_T &)(*',SymbolList{1},'[0])'));
            for x = 2:lunghezza{1}
            	name_fun=[name_fun,',real_T\t&',char(Params{1}(x))];
                call_fun{indr,indc}=strcat(call_fun{indr,indc},',',strcat('(real_T &)(*',SymbolList{1},'[',num2str(x-1),'])'));
            end
            
            for indlist=2:nlist
                for x = 1:lunghezza{indlist}
                    name_fun=[name_fun,',real_T\t&',char(Params{indlist}(x))];
                    call_fun{indr,indc}=strcat(call_fun{indr,indc},',',strcat('(real_T &)(*',SymbolList{indlist},'[',num2str(x-1),'])'));
                end
            end
            
            for indlist=1:nlist
                for x = 1:lunghezza{indlist}
                    if(isempty(strfind(M, strcat('cos(',char(Params{indlist}(x)),')')))==0)
                        name_fun=strcat(name_fun,',real_T\t&c',char(Params{indlist}(x)));
                        call_fun{indr,indc}=strcat(call_fun{indr,indc},',c',char(Params{indlist}(x)));
                        M = strrep(M,strcat('cos(',char(Params{indlist}(x)),')'),...
                            strcat('c',char(Params{indlist}(x))));
                        if(costrue{indlist}(x)==0)
                            fprintf(fid,strcat('\treal_T\tc',char(Params{indlist}(x)),...
                                '=cos(*',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                            costrue{indlist}(x)=1;
                        end
                    end
                    
                    if(isempty(strfind(M, strcat('sin(',char(Params{indlist}(x)),')')))==0)
                        name_fun=strcat(name_fun,',real_T\t&s',char(Params{indlist}(x)));
                        call_fun{indr,indc}=strcat(call_fun{indr,indc},',s',char(Params{indlist}(x)));
                        M = strrep(M,strcat('sin(',char(Params{indlist}(x)),')'),...
                            strcat('s',char(Params{indlist}(x))));
                        if(sintrue{indlist}(x)==0)
                            fprintf(fid,strcat('\treal_T\ts',char(Params{indlist}(x)),...
                                '=sin(*',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                            sintrue{indlist}(x)=1;
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(isempty(strfind(M, strcat('pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0)')))==0)
                            name_fun=strcat(name_fun,',real_T\t&pc',char(Params{indlist}(x)),'_',num2str(pot));
                            call_fun{indr,indc}=strcat(call_fun{indr,indc},',pc',char(Params{indlist}(x)),'_',num2str(pot));
                            M = strrep(M,strcat('pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0)'),...
                                strcat('pc',char(Params{indlist}(x)),'_',num2str(pot)));
                            if(powcostrue{indlist,pot-1}(x)==0)
                                fprintf(fid,strcat('\treal_T\tpc',char(Params{indlist}(x)),...
                                    '_',num2str(pot),'=pow(c',char(Params{indlist}(x)),',',num2str(pot),'.0);\n'));
                                powcostrue{indlist,pot-1}(x)=1;
                            end
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(isempty(strfind(M, strcat('pow(s',char(Params{indlist}(x)),',',num2str(pot),'.0)')))==0)
                            name_fun=strcat(name_fun,',real_T\t&ps',char(Params{indlist}(x)),'_',num2str(pot));
                            call_fun{indr,indc}=strcat(call_fun{indr,indc},',ps',char(Params{indlist}(x)),'_',num2str(pot));
                            M = strrep(M,strcat('pow(s',char(Params{indlist}(x)),',',num2str(pot),'.0)'),...
                                strcat('ps',char(Params{indlist}(x)),'_',num2str(pot)));
                            if(powsintrue{indlist,pot-1}(x)==0)
                                fprintf(fid,strcat('\treal_T\tps',char(Params{indlist}(x)),...
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
                fprintf(fid_mat,strcat('\treal_T\tMapleGenVar',num2str(indvar),';\n'));
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
    
    file = [Name,'S.cpp '];
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('-');
            var_mat = strcat('\tM[',num2str((indr-1)+(indc-1)*s1),']=',call_fun{indr,indc});
            file = [file,Name,'_',num2str(indr),'_',num2str(indc),'_S.obj '];
            fprintf(fid, var_mat);
        end
    end
    fprintf(fid,'}\n\n\n');
    
    fprintf(fid,'static void mdlTerminate(SimStruct *S){}\n\n');

    fprintf(fid,'#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */\n');
    fprintf(fid,'#include "simulink.c"      /* MEX-file interface mechanism */\n');
    fprintf(fid,'#else\n');
    fprintf(fid,'#include "cg_sfun.h"       /* Code generation registration function */\n');
    fprintf(fid,'#endif\n');

    fclose(fid);
    fclose(fid_extern);
    %fprintf(file);
    
    cd('bin');
    fprintf('\n*\n');
    eval(['mex ',file]);
    cd('..');
    fprintf(['Built: ',Name,'S.mexw64\n']);

end