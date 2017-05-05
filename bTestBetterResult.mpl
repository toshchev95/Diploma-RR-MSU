$include "C:\\Kursovay\\maple\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Generate.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

# function getBetterResultMatrix
getBetterResultMatrix := proc(listResultMatrix)
  local i, resultMatrix, tempOpMatrix;

  if nops(listResultMatrix) > 1 then
  	resultMatrix := cmpOreMatrix(listResultMatrix[1], listResultMatrix[2]);

  	for i from 3 to nops(listResultMatrix) do
  		tempOpMatrix := resultMatrix;
  	  resultMatrix := compareOreMatrix(tempOpMatrix, listResultMatrix[i]);
  	end do;
  else
  	resultMatrix := listResultMatrix[1];
  end if;

  return resultMatrix;
end proc:

# function compareOreMatrix
compareOreMatrix := proc(opMatrixA, opMatrixB)
	local i,j, bCompare, listNumbersDiffRows, A, B;

	A := matrixOreWithoutGCD(matrixOreWithoutDenom(A));
	B := matrixOreWithoutGCD(matrixOreWithoutDenom(B));
	#print(A,B);
	bCompare := true;
	listNumbersDiffRows := getListNumbersDiffRows(A, B);


	if bCompare = true then
		return opMatrixA;
	else
		return opMatrixB;
	end if;
end proc:

# function getListNumbersDiffRows
getListNumbersDiffRows := proc(A, B)
	local i,j,listNumbersDiffRows, size;

	listNumbersDiffRows := list();
	size := op(1, A)[1];

	for i to size do
		for j to size do
	#MultiplyPolyOnListOre
		end do;
	end do;

	return listNumbersDiffRows;	
end proc: