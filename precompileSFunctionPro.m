function precompileSFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%precompileSFunctionPro generates a S-function cpp-file in the bin directory from a 
%   SymbolicMatrix than compile it for a runtime execution with S-function.
%
%   precompileSFunctionPro(Name,SymbolicMatrix,Params,SymbolList)
%       Name: name of the generated file NameS.cpp
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
% 	precompileSFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});
%
%   result is a file MatrixS.cpp and MatrixS.mexw64
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

    fprintf(['Build ',Name,' to S-Function']);
    fun_name = strcat('bin/',Name,'S.cpp');
    fid = fopen(fun_name, 'wt');
    
    [s1,s2]=size(SymbolicMatrix);
    
    def_name = strcat('#define S_FUNCTION_NAME\t',Name,'S\n');
    fprintf(fid,def_name);
    fprintf(fid,'#define S_FUNCTION_LEVEL 2\n');
    fprintf(fid,'#include "simstruc.h"\n');
    fprintf(fid,'#include <math.h>\n\n\n');
       
    fprintf(fid,'static void mdlInitializeSizes(SimStruct *S)\n{\n');
    fprintf(fid,'\tssSetNumSFcnParams(S, 0);\n');
    fprintf(fid,'\tif (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))\n\t{\n');
    fprintf(fid,'\t\treturn; /* Parameter mismatch will be reported by Simulink */\n\t}\n\n');

    fprintf(fid,strcat('\tif (!ssSetNumInputPorts(S,',num2str(nlist),')) return;\n'));
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

    indvar=1;
    for indr=1:s1
        fprintf('\n');
        for indc=1:s2
            fprintf('.');
            M = ccode(vpa(SymbolicMatrix(indr,indc),10));
            M = strrep(M,'~','');
            M = strrep(M,'t0 = ','');
            while (isempty(strfind(M, strcat('MapleGenVar',num2str(indvar))))==0)
                fprintf(fid,strcat('\treal_T\tMapleGenVar',num2str(indvar),';\n'));
                indvar=indvar+1;
            end
            for indlist=1:nlist
                for x = 1:lunghezza{indlist}
                    M = strrep(M,char(Params{indlist}(x)),strcat('(*',SymbolList{indlist},'[',num2str(x-1),'])'));
                    if(costrue{indlist}(x)==0)
                        if(isempty(strfind(M, strcat('cos((*',SymbolList{indlist},'[',num2str(x-1),']))')))==0)
                            costrue{indlist}(x)=1;
                            fprintf(fid,strcat('\treal_T\tc',num2str(SymbolList{indlist}),num2str(x),...
                                '=cos(*',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                        end
                    end

                    if(sintrue{indlist}(x)==0)
                        if(isempty(strfind(M, strcat('sin((*',SymbolList{indlist},'[',num2str(x-1),']))')))==0)
                            sintrue{indlist}(x)=1;
                            fprintf(fid,strcat('\treal_T\ts',num2str(SymbolList{indlist}),num2str(x),...
                                '=sin(*',SymbolList{indlist},'[',num2str(x-1),']);\n'));
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(powcostrue{indlist,pot-1}(x)==0)
                            if(isempty(strfind(M, strcat('pow(cos((*',SymbolList{indlist},'[',num2str(x-1),'])),',num2str(pot),'.0)')))==0)
                                powcostrue{indlist,pot-1}(x)=1;
                                fprintf(fid,strcat('\treal_T\tpc',num2str(SymbolList{indlist}),...
                                    num2str(x),'_',num2str(pot),'=pow(c',SymbolList{indlist},num2str(x),',',num2str(pot),'.0);\n'));
                            end
                        end
                    end

                    for pot = 2:lunghezza{indlist}
                        if(powsintrue{indlist,pot-1}(x)==0)
                            if(isempty(strfind(M, strcat('pow(sin((*',SymbolList{indlist},'[',num2str(x-1),'])),',num2str(pot),'.0)')))==0)
                                powsintrue{indlist,pot-1}(x)=1;
                                fprintf(fid,strcat('\treal_T\tps',num2str(SymbolList{indlist}),...
                                    num2str(x),'_',num2str(pot),'=pow(s',SymbolList{indlist},num2str(x),',',num2str(pot),'.0);\n'));
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
    
    fprintf(fid,'static void mdlTerminate(SimStruct *S){}\n\n');

    fprintf(fid,'#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */\n');
    fprintf(fid,'#include "simulink.c"      /* MEX-file interface mechanism */\n');
    fprintf(fid,'#else\n');
    fprintf(fid,'#include "cg_sfun.h"       /* Code generation registration function */\n');
    fprintf(fid,'#endif\n');

    fclose(fid);
    
    cd('bin');
    Vai = [Name,'S.cpp'];
    fprintf('\n*\n');
    eval(['mex ',Vai]);
    cd('..');
    fprintf(['Built: ',Name,'S.mexw64\n']);

end