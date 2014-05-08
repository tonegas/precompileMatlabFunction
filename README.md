precompileMatlabFunction
========================
Precompile FunctionPro package
is a set of function for generate a fast execution (mfile, s-function, mex-function) file from a symbolic Matlab matrix.

The package contains:
precompileMFunctionPro to generate the m-file function from a symbolic matrix like 'matlabFunction()'.
	
	precompileMFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameM.m


precompileMatlabFunctionPro to generate the m-file function using a 'matlabFunction()'.
	
	precompileMatlabFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameMF.m

precompileCFunctionPro to generate C++ source MEX-Files from a symbolic matrix, 100 times fast in execution respect NameMF.m.
	
	precompileCFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameC.mexw64

precompileDistCFunctionPro will generate C++ source MEX-Files from a symbolic matrix, 100 times fast in execution respect NameMF.m, the difference from precompileCFunctionPro is that this function create one file for each element of the matrix and then will compile all file together.
	
	precompileDistCFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameC.mexw64

precompileSFunctionPro to generate S-Function from a symbolic matrix.
	
	precompileSFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameS.mexw64

precompileDistSFunctionPro will generate S-Function from a symbolic matrix the difference from precompileSFunctionPro is that this function create one file for each element of the matrix and then will compile all file together.
	
	precompileDistSFunctionPro(Name,SymbolicMatrix,Params,SymbolList) -> NameS.mexw64


example:

	syms t1 t2 t3 d1 d2;
	M=[sin(t1)*sin(t2+t3)-d1*d2,d1+d2;t1*t2*t3,sin(t1+t2+t3)-d1*d2];
	mkdir('bin');
	precompileCFunctionPro('Matrix',M,{[t1,t2,t3],[d1,d2]},{'t','d'});

	% result is a files Matrix_1_1_C.cpp, Matrix_1_2_C.cpp, MatrixC.cpp, extern_MatrixC.h and MatrixC.mexw64
   
	addpath('bin');
	M=MatrixC([1,2,3],[2,3])
	
	% M =
	%    -6.8069    5.0000
	%     6.0000   -6.2794  

usage 

	test.m 
	
for a total benchmark.

PhD Gastone Pietro Rosati Papini, Percro Laboratory Scuola Superiore Sant'Anna Pisa 

