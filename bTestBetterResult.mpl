$include "C:\\Kursovay\\maple\\Git\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Git\\Generate.mpl"
$include "C:\\Kursovay\\maple\\Git\\ImprovingRR.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

# function getListResultMatrixLessIteration
getListResultMatrixIteration := proc(list_UID_Result, numberIteration) # [resultMatrix, numberIteration]
  local resultMatrix, i, j, size, listResultMatrixLessIteration;
  listResultMatrixLessIteration := list();

  for i to nops(list_UID_Result) do
    if numberIteration = list_UID_Result[2] then
      listResultMatrixLessIteration := [op(listResultMatrixLessIteration), list_UID_Result[1]];
    end if;
  end do;

  return listResultMatrixLessIteration;
end proc:

# function getBetterResultMatrix
getBetterResultMatrix := proc(list_UID_Result, iteration) # [resultMatrix, numberIteration]
  local i, resultMatrix, tempOpMatrix, listResult, listResultMatrix;
  global UID_betterNumber;

  listResultMatrix := getListResultMatrixLessIteration(list_UID_Result, iteration);

  if nops(listResultMatrix) > 1 then
  	UID_betterNumber := 1;
  	listResult := compareOreMatrix(listResultMatrix[1], listResultMatrix[2]);
  	resultMatrix := listResult[1];
  	if listResult[2] = false then
  		UID_betterNumber := 2;
  	end if;

  	for i from 3 to nops(listResultMatrix) do
  		tempOpMatrix := resultMatrix;
  	  
  	  listResult := compareOreMatrix(tempOpMatrix, listResultMatrix[i]);
  	  resultMatrix := listResult[1];

  	  if listResult[2] = false then
  			UID_betterNumber := i;
  		end if;
  	end do;
  else
  	resultMatrix := listResultMatrix[1];
  end if;

  #print(resultMatrix, UID_betterNumber);
  return resultMatrix;
end proc:

# function compareOreMatrix
compareOreMatrix := proc(opMatrixA, opMatrixB)
	local i,j, bCompare, bResult, listNumbersDiffRows, A, B,
		size, countNumbersDiffRows, bMatrixRows, listNumbersDiffRowsA, listNumbersDiffRowsB;
  global UID_equalMatrix;

  #print(whattype(opMatrixA), opMatrixA);
  A := matrixOreWithoutGCD(matrixOreWithoutDenom(opMatrixA));
  B := matrixOreWithoutGCD(matrixOreWithoutDenom(opMatrixB));
	
	print("Compare Matrix");
	print(A);
	print(B);
	bCompare := true;
	size := op(1, opMatrixA)[1];

	listNumbersDiffRows := getListNumbersDiffRows(A, B);
	print(listNumbersDiffRows);
	listNumbersDiffRowsA := listNumbersDiffRows[1];
	listNumbersDiffRowsB := listNumbersDiffRows[2];
	countNumbersDiffRows := nops(listNumbersDiffRowsA);

	if countNumbersDiffRows <> 0 then
		bMatrixRows := Matrix(countNumbersDiffRows);

		for i to countNumbersDiffRows do 
			for j to countNumbersDiffRows do
				#if evalb(convert(A[i],list) = convert(B[j],list)) or evalb(convert(A[i],list) = MultiplyPolyOnListOre(-1, convert(B[j],list))) then
				#	bMatrixRows[i, j] := 0;
				#else
				#end if;
					# compute
					bResult := compareRowsOpMatrix2( convert(A[listNumbersDiffRowsA[i]], list), convert(B[listNumbersDiffRowsB[j]],list));
					if bResult = true then
						bMatrixRows[i, j] := 1;
					else
						bMatrixRows[i, j] := -1;
					end if;
			end do;
		end do;

		if op(1, bMatrixRows)[1] <> 1 then
			bCompare := processDataCollection(bMatrixRows);
		elif bMatrixRows[1,1] = 1 then
			bCompare := true;
		elif bMatrixRows[1,1] = -1 then
			bCompare := false;
		end if;

		print("bMatrixRows=",bMatrixRows, bCompare);
	end if;


	if bCompare = true then
		return [A, true];
	else
		return [B, false];
	end if;
end proc:

# function compareRowsOpMatrix2 (rowA, rowB)
compareRowsOpMatrix2 := proc(rowA, rowB)
  local i,j, m_resRowInfoA, m_resRowInfoB, bResult,
    countRowA, countRowB, sumOrdersPolynomialsRowA, sumOrdersPolynomialsRowB,sumNumbersAdditionsRowA,sumNumbersAdditionsRowB,
    highDiffA,highDiffB, sumDiffOrdersA, sumDiffOrdersB, listRowA, listRowB, numberDiffRowA,numberDiffRowB,
    bMatrix;
  global bRepeatCompareRows, size;

  size := nops(convert(rowA, list));

  m_resRowInfoA := OrePoly(op(getNonNullList(getSumOrePolyInRows(rowA))));
  m_resRowInfoB := OrePoly(op(getNonNullList(getSumOrePolyInRows(rowB))));

  # Сравним строки по критериям сравнения 2х элементов 1 - 5
  # 1
  highDiffA := getHighDifferRow2(rowA);
  highDiffB := getHighDifferRow2(rowB);
  # 2
  listRowA := convert(map(proc (ore) options operator, arrow; [seq(x, `in`(x, ore))] end proc, rowA), list);
  listRowB := convert(map(proc (ore) options operator, arrow; [seq(x, `in`(x, ore))] end proc, rowB), list);
  sumDiffOrdersA := sum('map(proc (ore) options operator, arrow; getSumDiffOrders(ore) end proc, listRowA)[k]', k = 1 .. size);
  sumDiffOrdersB := sum('map(proc (ore) options operator, arrow; getSumDiffOrders(ore) end proc, listRowB)[k]', k = 1 .. size);
  # 3
  numberDiffRowA := getNumberDiffRow(listRowA);
  numberDiffRowB := getNumberDiffRow(listRowB);
  # 4
  sumNumbersAdditionsRowA := getSumNumbersAdditionsRow(listRowA);
  sumNumbersAdditionsRowB := getSumNumbersAdditionsRow(listRowB);
  # 5
  sumOrdersPolynomialsRowA := getSumOrdersPolynomialsRow(listRowA);
  sumOrdersPolynomialsRowB := getSumOrdersPolynomialsRow(listRowB);
  if highDiffA < highDiffB then
    return true;
  elif highDiffA > highDiffB then
    return false;
  elif sumDiffOrdersA < sumDiffOrdersB then
    return true;
  elif sumDiffOrdersA > sumDiffOrdersB then
    return false;
  elif numberDiffRowA < numberDiffRowB then  
    return true;
  elif numberDiffRowA > numberDiffRowB then
    return false;
  elif sumNumbersAdditionsRowA < sumNumbersAdditionsRowB then
    return true;
  elif sumNumbersAdditionsRowA > sumNumbersAdditionsRowB then
    return false;
  elif sumOrdersPolynomialsRowA < sumOrdersPolynomialsRowB then
    return true;
  elif sumOrdersPolynomialsRowA > sumOrdersPolynomialsRowB then
    return false;
  else
    print("In func compareRowsOpMatrx: execute special comparing algorithm");
  end if;  

  # a) сумма всех порядков дифференцирования (<)
  # b) кол-во термов порядков диф-ния (>)
  # c) степени у множителей коэффициентов полиномов Оре (<)
  # d) кол-во термов f(x) (<)
  # e) числа у множителей (<)
  
  bRepeatCompareRows := false;  
  bResult := compareOrePoly(m_resRowInfoA,m_resRowInfoB, rowA, rowB);
  if bRepeatCompareRows = true then
    #print("!");

    bMatrix := Matrix(size);
    for i to size do # rowA
      for j to size do # rowB
        bRepeatCompareRows := false;
        bResult := compareOrePoly(rowA[i], rowB[j], rowA, rowB);
        if bRepeatCompareRows = false then
          if bResult = true then
            bMatrix[i,j] := 1;
          else
            bMatrix[i,j] := -1;
          end if;
        end if;
      end do;

    end do;

    #print(rowA,rowB);
    bResult := processDataCollection(bMatrix);
  end if;

  return bResult;
end proc:


# function getListNumbersDiffRows
getListNumbersDiffRows := proc(A, B)
	local i,j,listNumbersDiffRowsA,listNumbersDiffRowsB, size;

	size := op(1, A)[1];
	listNumbersDiffRowsA := [seq(i, i=1..size)];
	listNumbersDiffRowsB := [seq(i, i=1..size)];

	for i to size do
		for j to size do
			if evalb(convert(A[i],list) = convert(B[j],list)) or evalb(convert(A[i],list) = MultiplyPolyOnListOre(-1, convert(B[j],list))) then
				listNumbersDiffRowsA[i] := 0;
				listNumbersDiffRowsB[j] := 0;
				break;
			end if;
		end do;
	end do;

	listNumbersDiffRowsA := sort(ListTools:-MakeUnique(listNumbersDiffRowsA));
	if nops(listNumbersDiffRowsA) = 1 and listNumbersDiffRowsA[1] = 0 then
		listNumbersDiffRowsA := list();
	elif listNumbersDiffRowsA[1] = 0 then
		listNumbersDiffRowsA := listNumbersDiffRowsA[2.. nops(listNumbersDiffRowsA)];
	end if;

	listNumbersDiffRowsB := sort(ListTools:-MakeUnique(listNumbersDiffRowsB));
	if nops(listNumbersDiffRowsB) = 1 and listNumbersDiffRowsB[1] = 0 then
		listNumbersDiffRowsB := list();
	elif listNumbersDiffRowsB[1] = 0 then
		listNumbersDiffRowsB := listNumbersDiffRowsB[2.. nops(listNumbersDiffRowsB)];
	end if;

	return [listNumbersDiffRowsA, listNumbersDiffRowsB];	
end proc: