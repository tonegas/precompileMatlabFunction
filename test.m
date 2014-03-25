%% Build
mkdir('bin');
load('B');
disp('----------------------------Start Build---------------------------');
syms t1 t2 t3 t4 t5 t6
precompileMFunctionPro('B',B,{[t1,t2,t3,t4,t5,t6]},{'q'});
precompileMatlabFunctionPro('B',B,{[t1 t2 t3 t4 t5 t6]},{'q'});
precompileCFunctionPro('B',B,{[t1,t2,t3,t4,t5,t6]},{'q'});
precompileDistCFunctionPro('BD',B,{[t1,t2,t3,t4,t5,t6]},{'q'});
precompileSFunctionPro('B',B,{[t1,t2,t3,t4,t5,t6]},{'q'});
precompileDistSFunctionPro('BD',B,{[t1,t2,t3,t4,t5,t6]},{'q'});
disp('----------------------------Stop Build----------------------------');

%% Test
disp('--------------------------Start Test Sub--------------------------');
%test with subs
load('B');
tic
M=eval(subs(B,{'t1','t2','t3','t4','t5','t6'},{1,2,3,4,5,6}));
disp('For one clicle of subs the');
toc
clear B;
disp('--------------------------Stop Test Sub---------------------------');

%% Test Function
disp('----------------------------Start Test----------------------------');
clear B;
addpath('bin');
%test M-function
tic
N=1000;
for i=1:N
    M=BM([i,i*2,i*3,i*4,i*5,i*6]);
end
disp('For 1000 clicle M-function the');
toc

%test M-function
tic
N=1000;
for i=1:N
    M=BMF([i,i*2,i*3,i*4,i*5,i*6]);
end
disp('For 1000 clicle MatlabFunction the');
toc

%test M-function
tic
N=1000;
for i=1:N
    M=BC([i,i*2,i*3,i*4,i*5,i*6]);
end
disp('For 1000 clicle Mex-Function the');
toc

%test M-function
tic
N=1000;
for i=1:N
    M=BDC([i,i*2,i*3,i*4,i*5,i*6]);
end
disp('For 1000 clicle Mex-Function distinct element the');
toc
disp('----------------------------Stop Test-----------------------------');