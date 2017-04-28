$include "C:\\Kursovay\\maple\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Generate.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

# function printRankNullSpaceOpMatrix
printRankNullSpaceOpMatrix := proc(opMatrix::Matrix) 
  local m_front; m_front := getFrontMatrix(opMatrix); 
  print("frontMatrix=", m_front); 
  print("Rank = ", LinearAlgebra:-Rank(m_front)); 
  print("NullSpace=", LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(m_front)));
end proc:

# function printList_opIterMatrix
printList_opIterMatrix := proc(listMatrix)
  local i;
  print(nops(listMatrix));
  for i to nops(listMatrix) do
    print(i,listMatrix[i]);
  end do;
end proc:

# function modifyRR
modifyRR := proc(opMatrix::Matrix)
  local A,B,step,i, saved,nullSpace_indexRowOrderDiff, firstListNumberHighDiffForUniMat;
  local vector_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff, uni, m_newRow;
  local height := op(1, opMatrix)[1], 
        width := op(1, opMatrix)[2];
  global front, m_matrix, fullNullSpace,m_deg_rows, m_infoOnMatrix, m_matrixInfo,
    listIterMatrix, size := height;

  if height <> width then
    error "Func modifyRR: wrong scale opMatrix <- height =/= width";
  end if;

  # init
  step:=1; listIterMatrix := [opMatrix];

  m_matrix := matrixOreWithoutDenom(opMatrix);
  front := getFrontMatrix(m_matrix);

  # Statistics:-Count
  #print(A, LinearAlgebra:-Rank(A),LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(A)));

  while (LinearAlgebra:-Rank(front) < LinearAlgebra:-RowDimension(front)) do 

    # Data collection
    fullNullSpace := nullspaceWithoutDenom(front);
    m_matrixInfo := Matrix(size);
    m_deg_rows := vector(size, 0);
    for i to size do
      m_deg_rows[i] := getHighDifferRow(m_matrix[i]);
    end do;

    firstListNumberHighDiffForUniMat := getListNumberHighDiffForUniMat(fullNullSpace[1], m_deg_rows);
    if nops(fullNullSpace) = 1 and nops(firstListNumberHighDiffForUniMat) = 1 then
      m_nullSpace := fullNullSpace[1];
      m_indexRowOrderDiff := firstListNumberHighDiffForUniMat[1];
    else
      # Info:Data collection
  
      m_infoOnMatrix := list();
      getInfoOnOpMatrix();
  
      # issues modify
      nullSpace_indexRowOrderDiff := estimations();
      m_nullSpace := nullSpace_indexRowOrderDiff[1];
      m_indexRowOrderDiff := nullSpace_indexRowOrderDiff[2];
      m_newRow := nullSpace_indexRowOrderDiff[3];
    end if;

    uni := getUnimodulMatrix(m_matrix, height, m_nullSpace, m_deg_rows, m_indexRowOrderDiff);
    #saved := getReverseLUMatrix(m_matrix, uni);
    #saved := matrixOreWithoutDenom(saved);

    saved := m_matrix;
    saved[m_indexRowOrderDiff] := m_newRow;

    print(step, m_nullSpace, "index=",m_indexRowOrderDiff,uni,step+1);
    listIterMatrix := [op(listIterMatrix), saved];

    # check again
    front := getFrontMatrix(saved);
    m_matrix := saved;
  end do;  

  printList_opIterMatrix(listIterMatrix);

  return m_matrix;
end proc:



# Функция, собирающая все данные операторной матрицы, 
# строк унимод и полученных на следующем шаге матриц
# ! доработать
# function getInfoOnOpMatrix
getInfoOnOpMatrix := proc()
  local i,j,k, listNumberRowsForUniMatrix, numbersRowsNullSpace, temp, maxOrd, highDiffMatrix,
    listNumbersRowsColsInfo, listOfOrderDiff, listOfPoly, listEmpty, listOrderDiffRowUniMat;
  global m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, size, # global for modifyRR()
    m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix, nullSpace;

  highDiffMatrix := max(convert(m_deg_rows, list));
  maxOrd := 0;

  # Список номеров строк для неоднозначности при получении унимод матриц
  listNumberRowsForUniMatrix := list(); # [[1,2],[3,4]]
  # Список номеров строк, столбцов операторной матрицы, по к. есть неоднозначности uniMatrix
  listNumbersRowsColsInfo := list(); # [1,2,3,4]
  for i to nops(fullNullSpace) do
    nullSpace := fullNullSpace[i];

    numbersRowsNullSpace := getListNumberHighDiffForUniMat(nullSpace, m_deg_rows);
    # ?
    listNumbersRowsColsInfo := [op(listNumbersRowsColsInfo),op(numbersRowsNullSpace)];
    listNumberRowsForUniMatrix := [op(listNumberRowsForUniMatrix),numbersRowsNullSpace];
  end do;
  listNumbersRowsColsInfo := ListTools:-MakeUnique(listNumbersRowsColsInfo);
  print("listNumberRowsForUniMatrix=",listNumberRowsForUniMatrix);
  print("listNumbersRowsColsInfo=",listNumbersRowsColsInfo);

  # Вычислим матрицу, где каждый элемент содержит инфу об OrePoly в виде 2-х списков:
  # listOfOrderDiff, listOfPoly
  # + 1 за нулевой порядок дифференцирования
  listEmpty := convert(vector(highDiffMatrix + 1,0),list); # [0,0,...,0]

  for i to size do
    for j to size do # Списки можно перевернуть (????)
      # compute
      listOfPoly := [seq(x, `in`(x, m_matrix[i,j]))];

      # compute
      listOfOrderDiff := list();
      for k to nops(listOfPoly) do 
        if listOfPoly[k] <> 0 then 
          temp := 1;
        else 
          temp := 0;
        end if; 
        listOfOrderDiff := [op(listOfOrderDiff), temp]; 
      end do;
      #listOfOrderDiff := ListTools:-Reverse(listOfOrderDiff);

      m_matrixInfo[i,j]:= [sumListDifferentLength(listOfOrderDiff, listEmpty), 
                    sumListDifferentLength(listOfPoly,listEmpty)];
    end do;

  end do;


  # А также вычислим общую информацию по строкам и столбцам opMatrix, индексы к. нам известны
  # Вычислим инфу по строкам и столбцам opMatrix
  m_RowsInfo := convert(vector(size,0),list);
  m_ColsInfo := convert(vector(size,0),list);

  # compute Rows
  for i to size do
    m_RowsInfo[i] := getInfoMatrixRow(m_matrixInfo[i]);
  end do;

  # compute Cols
  for i to size do
    m_ColsInfo[i] := getInfoMatrixCol( LinearAlgebra:-Column(m_matrixInfo,i) );
  end do;

  # Вычислим строки унимодулярных операторов по спец номерам rowHighOrderDiff
  # using listNumberRowsForUniMatrix, m_deg_rows, fullNullSpace
  m_listRowsInfoUniMatrix := convert(vector(size,0),list);
  for i to size do
    listOrderDiffRowUniMat := convert(vector(size,m_deg_rows[i]),list) - m_deg_rows;
    m_listRowsInfoUniMatrix[i]:= [listOrderDiffRowUniMat, max(listOrderDiffRowUniMat)];       # , info_2, info_3, ...]
    # if member(i,listNumbersRowsColsInfo) then
    # end if;
  end do;

  # m_infoOnMatrix - структура, содержащая доп инфу
  m_infoOnMatrix := [m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix,listNumberRowsForUniMatrix];

end proc:

# function getSumOrePolyInRows
getSumOrePolyInRows := proc (orePoly) 
  local i, listPolynomials, maxHighDiffOrder, listEmpty; 
  listPolynomials := convert(map(proc (ore) options operator, arrow; [seq(x, `in`(x, ore))] end proc, orePoly), list); 
  maxHighDiffOrder := max(getHighDifferRow(orePoly)); 
  listEmpty := convert(vector(maxHighDiffOrder+1, 0), list); 
  for i to nops(listPolynomials) do 
    listEmpty := sumListDifferentLength(listEmpty, listPolynomials[i]);
  end do; 

  return listEmpty;
end proc:

# function getListNumberHighDiffForUniMat
getListNumberHighDiffForUniMat := proc(nullSpace, m_deg_rows)
  local j, numbersRowsNullSpace, maxOrd;
  
  maxOrd := 0;
  numbersRowsNullSpace := list();
  for j to nops(nullSpace) do
    if nullSpace[j] <> 0 then
      if m_deg_rows[j] > maxOrd then
        maxOrd := m_deg_rows[j];
        numbersRowsNullSpace := [j];
      elif m_deg_rows[j] = maxOrd then
        numbersRowsNullSpace := [op(numbersRowsNullSpace),j];
      end if;
    end if;

  end do;

  return numbersRowsNullSpace;
end proc:

# function sumListDifferentLength
sumListDifferentLength := proc(list1, list2) 
  local i, j, listBig, listSmall, sumList; 
  if nops(list1) < nops(list2) then 
    listSmall := list1; 
    listBig := list2;
  else 
    listSmall := list2; 
    listBig := list1;
  end if; 
  sumList := [op(listSmall), op(convert(vector(nops(listBig)-nops(listSmall), 0), list))]; 
  return sumList+listBig; 
end proc:

# function getInfoMatrixRow
getInfoMatrixRow := proc(Row)
  local i, listInfo;
  listInfo := Row[1];

  for i from 2 to nops(Row) do
    listInfo := listInfo + Row[i];
  end do;

  return listInfo;
end proc:

# function getInfoMatrixCol
getInfoMatrixCol := proc(Col)
  local i, listInfo;
  listInfo := Col[1];

  for i from 2 to nops(Col) do
    listInfo := listInfo + Col[i];
  end do;

  return listInfo;
end proc:


##################################
# Optional estimations are here! #
##################################

# function estimations
estimations := proc()
  local i,j, 
    m_nullSpace, m_indexRowOrderDiff, listNumberRowsForUniMatrix, newRow, listNumberRowsForNullSpace, numberRowForNullSpace,
    nullSpace, listPar_s, tempPar_s, betterParameters, m_newRow;
  global size,m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, front;

  listNumberRowsForUniMatrix := m_infoOnMatrix[4];

  # Сформируем список параметров спец образом, по которым существуют разветвления.
  # в последовательном порядке расположения nullSpace в fullNullSpace и номеров строк HighOrderDiff
  listPar_s := list();

  for i to nops(fullNullSpace) do
    listNumberRowsForNullSpace := listNumberRowsForUniMatrix[i];
    nullSpace := fullNullSpace[i];
    for j to nops(listNumberRowsForNullSpace) do
      numberRowForNullSpace := listNumberRowsForNullSpace[j];
      newRow := getReplaceRowMatrixRR(m_matrix, size, nullSpace, m_deg_rows, numberRowForNullSpace);

      # fill out -> tempPar_s[4]
      tempPar_s := [nullSpace, numberRowForNullSpace, newRow, getHighDifferRow(newRow)];

      listPar_s := [op(listPar_s), tempPar_s];
    end do;
  end do;

  betterParameters := cmpParameters(listPar_s[1], listPar_s[2]);
  for i from 3 to nops(listPar_s) do
    betterParameters := cmpParameters(betterParameters, listPar_s[i]);
  end do;

  # Result
  m_nullSpace := betterParameters[1];
  m_indexRowOrderDiff := betterParameters[2];
  m_newRow := betterParameters[3];
  return [m_nullSpace, m_indexRowOrderDiff, m_newRow];
end proc:

# function cmpParameters ( listParameters := [nullSpace, numberRowForNullSpace, newRow, getHighDifferRow(newRow)] )
cmpParameters := proc(listPar_A, listPar_B)
  local m_listRowsInfoUniMatrix,listOrderDiffRowUniMat,
    listPar_temp, bCompareRows, polyNullSpace_A, polyNullSpace_B;
  global size,m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, front,
    indexA, indexB;
  
  m_listRowsInfoUniMatrix := m_infoOnMatrix[3]; # [listOrderDiffRowUniMat, max(listOrderDiffRowUniMat)];
  listOrderDiffRowUniMat := m_listRowsInfoUniMatrix[1];

    listPar_temp := [0,0,0,0]; # (!):(delete)
    bCompareRows := false;
    indexA := listPar_A[2];
    indexB := listPar_B[2];

  # Будем применять эвристики по убыванию
  # if then elif then elif then else end if;
  # 1) Выберем min( delta(i) + HighOrderDiff(opMatrix in non-zero row's nullSpace) )
  # Выберем минимальную сумму разности порядков дифференцирования по выбранной строке и 
  # максимальный порядок дифф. в строках, кот. соответствуют не нулевые значения nullSpace
  # При выборе учесть, что в новой строке HighOrderDiff м.б. = 0.
  # (!) Note: достаточно вычислить новую строку
  if evalb(listPar_A[3] <> listPar_B[3]) then

    polyNullSpace_A := listPar_A[1][indexA];
    polyNullSpace_B := listPar_B[1][indexB];

    # maxHighDiffOrder
    if listPar_A[4] < listPar_B[4] then
      listPar_temp := listPar_A;
    elif listPar_A[4] > listPar_B[4] then
      listPar_temp := listPar_B;
    
    # nops polyNullSpace
    elif nops(polyNullSpace_A) < nops(polyNullSpace_B) then # # listPar_A[4] == listPar_B[4]
      listPar_temp := listPar_A;
    elif nops(polyNullSpace_A) > nops(polyNullSpace_B) then
      listPar_temp := listPar_B;

    # # degree polyNullSpace
    # elif degree(polyNullSpace_A) < degree(polyNullSpace_B) then
    #   listPar_temp := listPar_A;
    # elif degree(polyNullSpace_A) > degree(polyNullSpace_B) then
    #   listPar_temp := listPar_B;

    elif indexA = indexB then
      if compareRowsOpMatrx(listPar_A[3], listPar_B[3]) = true then
        listPar_temp := listPar_A;
      else
        listPar_temp := listPar_B;
      end if;
    else 
      bCompareRows := true;
    end if;

  # 2) Если на шаге 1) все значения равны, то смотри на строки в оп. матрицах,
  # отличающиеся между собой параметрами исследуемых строк: (по убыванию)
  else # listPar_A[3] == listPar_B[3]
    bCompareRows := true;
  end if;

  if bCompareRows = true then
    if compareRowsOpMatrx(m_matrix[indexB], m_matrix[indexA]) = true then
      listPar_temp := listPar_A;
    else
      listPar_temp := listPar_B;
    end if;
  
  end if;

  return listPar_temp;
end proc:

# function compareRowsOpMatrx (rowA, rowB)
compareRowsOpMatrx := proc(rowA, rowB)
  local i,j, m_resRowInfoA, m_resRowInfoB, bResult,
    countRowA, countRowB, bListA, bListB;
  global size,m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, front,
    indexA, indexB, bRepeatCompareRows,
    bMatrix;

  print(indexA, indexB);
  if evalb(rowA = m_matrix[indexB]) and evalb(rowB = m_matrix[indexA]) then
    # m_infoOnMatrix := [m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix,listNumberRowsForUniMatrix];
    # m_RowsInfo := [sumListDifferentLength(listOfOrderDiff), sumListDifferentLength(listOfPoly)];
    # listOfOrderDiff - бесполезен
    print("if then");
    m_resRowInfoA := OrePoly(op(getNonNullList(m_infoOnMatrix[1][2][indexB])));
    m_resRowInfoB := OrePoly(op(getNonNullList(m_infoOnMatrix[1][2][indexA])));
  else
    print("else", rowA);
    m_resRowInfoA := OrePoly(op(getNonNullList(getSumOrePolyInRows(rowA))));
    m_resRowInfoB := OrePoly(op(getNonNullList(getSumOrePolyInRows(rowB))));
  end if;

  # a) сумма всех порядков дифференцирования (<)
  # b) кол-во термов порядков диф-ния (>)
  # c) степени у множителей коэффициентов полиномов Оре (<)
  # d) кол-во термов f(x) (<)
  # e) числа у множителей (<)
  
  bRepeatCompareRows := false;  
  bResult := compareOrePoly(m_resRowInfoA,m_resRowInfoB, rowA, rowB);
  if bRepeatCompareRows = true then
    print("!");
    countRowA := 0; 
    countRowB := 0;
    bListA := [seq(true,k=1..size)];
    bListB := [seq(false,k=1..size)];
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
          #bListA[j] := bListA[j] and bResult;
          #bListB[j] := bListB[j] or bResult;
        end if;
      end do;

    end do;

    #for i to size do
    #  if bListA[i] = false then
    #    countRowA := countRowA + 1;
    #  end if;
    #  if bListB[i] = true then
    #    countRowB := countRowB + 1;
    #  end if;
    #end do;
    #bResult := evalb(countRowA > countRowB);

    bResult := processDataCollection();
  end if;

  return bResult;
end proc:

# function processDataCollection
processDataCollection := proc()
  local i,j,list_NumberA, list_NumberB;
  global size, bMatrix; # {-1,0,1}
  countRowA := 0;
  countRowB := 0;
  list_NumberA := computeOptimalVector(bMatrix); # {-1,0,1}
  bMatrix := - bMatrix;
  list_NumberB := computeOptimalVector(bMatrix); # {-1,0,1}

  if list_NumberA[1] < list_NumberB[1] then
    return true;
  elif list_NumberA[1] > list_NumberB[1] then
    return false;
  elif list_NumberA[3] > list_NumberB[3] then
    return true;
  elif list_NumberA[3] < list_NumberB[3] then
    return false;
  else
    print("!!: processDataCollection");
  end if;
end proc:

# function computeVectors
computeVectors := proc (m::Matrix) 
  local i, j, listLists, minorLists, count, size; 
  listLists := list(); 
  size := LinearAlgebra:-RowDimension(m); 
  if size = 1 then 
    return [m[1, 1]];
  end if; 
  for i to size do 
    count := m[i][1]; 
    minorLists := computeVectors(LinearAlgebra:-Minor(m, i, 1, output = ['matrix'])); 
    for j to size-1 do 
      listLists := [op(listLists), [count, op(minorLists[j])]]; 
    end do;
  end do; 
  return listLists;
end proc:

# function computeOptimalVector
computeOptimalVector := proc (m::Matrix) 
  local i, j, listNumberLists, listValues, listTemp, listEmpty, resultVector, count, size, listSort, temp; 
  size := LinearAlgebra:-RowDimension(m); 
  listNumberLists := computeVectors(m); 
  print(listNumberLists); 
  listEmpty := convert(vector(size, 0), list); 
  listValues := list(); 
  for i to nops(listNumberLists) do 
    listTemp := listEmpty; 
    for j to size do 
      if listNumberLists[i][j] = -1 then 
        listTemp[1] := listTemp[1]+1;
      elif listNumberLists[i][j] = 0 then 
        listTemp[2] := listTemp[2]+1;
      elif listNumberLists[i][j] = 1 then 
        listTemp[3] := listTemp[3]+1;
      end if;
    end do; 
    listValues := [op(listValues), listTemp];
  end do; print(listValues); 
  listSort := list(); 
  for i to nops(listValues) do 
    temp := insertValuesInList(listValues[i], listSort); 
    listSort := temp;
  end do; 
  print(listSort); 
  resultVector := listSort[1]; 
  for i to nops(listValues) do 
    if resultVector = listValues[i] then 
      count := i; 
      break; 
    end if;
  end do; 
  return listNumberLists[count];
end proc:

# function insertValuesInList
insertValuesInList := proc (values, listSort) 
  local i, j, newListSort; 
  if nops(listSort) = 0 then 
    return [values];
  end if; 
  for i to nops(listSort) do 
    if listSort[i] <> values then 
      if compareValues(values, listSort[i]) then 
        if i = 1 then 
          return [values, op(listSort)];
        else 
          return [op(listSort[1 .. i-1]), values, op(listSort[i+1, nops(listSort)])];
        end if;
      elif i = nops(listSort) then 
        return [op(listSort), values];
      end if;
    else 
      return listSort;
    end if;
  end do;
end proc: 

# function compareValues
compareValues := proc (valueA, valueB) 
  if valueA[1] < valueB[1] then 
    return true;
  elif valueB[1] < valueA[1] then 
    return false; 
  elif valueB[3] < valueA[3] then 
    return true;
  elif valueA[3] < valueB[3] then 
    return false;
  end if;
end proc:

# function compareOrePoly
compareOrePoly := proc (oreA, oreB, rowA, rowB) 
  local i, j, listA, listB, sizeA, sizeB, sumPolyA, sumPolyB, listIsDiffA, listIsDiffB, 
        numberDiffA, numberDiffB, sumPlusCoeffsA, sumPlusCoeffsB; 
  global bRepeatCompareRows;
  listA := [seq(x, `in`(x, oreA))]; 
  listB := [seq(x, `in`(x, oreB))]; 

  print(oreA,oreB);
  sizeA := nops(listA); 
  sizeB := nops(listB); 
  sumPolyA := sum('listA[k]', k = 1 .. sizeA); 
  sumPolyB := sum('listB[k]', k = 1 .. sizeB); 

  listIsDiffA := map(proc (x) options operator, arrow; if x = 0 then return 0 else return 1 end if end proc, listA); 
  listIsDiffB := map(proc (x) options operator, arrow; if x = 0 then return 0 else return 1 end if end proc, listB); 
  numberDiffA := sum('listIsDiffA[k]', k = 1 .. sizeA); 
  numberDiffB := sum('listIsDiffB[k]', k = 1 .. sizeB); 

  sumPlusCoeffsA := sum('seq(abs(c), `in`(c, coeffs(sumPolyA,x)))[k]', k = 1 .. nops(sumPolyA)); 
  sumPlusCoeffsB := sum('seq(abs(c), `in`(c, coeffs(sumPolyB,x)))[k]', k = 1 .. nops(sumPolyB)); 

  if nops(oreA) < nops(oreB) then 
    return true 
  elif nops(oreB) < nops(oreA) then 
    return false 
  elif numberDiffA < numberDiffB then 
    return true 
  elif numberDiffB < numberDiffA then 
    return false 
  elif degree(sumPolyA) < degree(sumPolyB) then 
    return true 
  elif degree(sumPolyB) < degree(sumPolyA) then 
    return false 
  elif nops(sumPolyA) < nops(sumPolyB) then 
    return true 
  elif nops(sumPolyB) < nops(sumPolyA) then 
    return false 
  elif sumPlusCoeffsA < sumPlusCoeffsB then 
    bRepeatCompareRows := true;
    return true;
  elif sumPlusCoeffsB < sumPlusCoeffsA then 
    bRepeatCompareRows := true;
    return false; 
  else 
    # using sumListDifferentLength(), coeffFull()
    print("In func compareOrePoly: else"); 
    bRepeatCompareRows := true;
    return true;
  end if 
end proc:


# function coeffFull
coeffFull := proc (poly) 
  local i, listCoeff; 
  listCoeff := list(); 
  for i from 0 to degree(poly) do 
    listCoeff := [op(listCoeff), coeff(poly, x, i)];
  end do; 

  return listCoeff;
end proc:



getNonNullList := proc (listA) 
  local i, tempList, countZero; 
  countZero := 0; 
  tempList := ListTools:-Reverse(listA); 
  for i to nops(listA) do 
    if tempList[i] = 0 then 
      countZero := countZero+1; 
    else 
      break;
    end if;
  end do; 

  return listA[1 .. nops(listA)-countZero]; 
end proc:

##########################################
###### Init functions RR algorithm #######
##########################################

# function MultiplyOrePolyOnListOre
MultiplyOrePolyOnListOre := proc (orePoly, listOrePoly) 
  return [seq(OreTools:-Multiply(orePoly, polyOre, R), `in`(polyOre, listOrePoly))];
end proc: 

# function MultiplyPolyOnListOre
MultiplyPolyOnListOre := proc (poly, listOrePoly) 
  return [seq(OreTools:-Multiply(OrePoly(poly), polyOre, R), `in`(polyOre, listOrePoly))];
end proc: 

# function AddListsOrePoly
AddListsOrePoly := proc (listA, listB) 
  local i, listSum := list(); 
  for i to nops(listA) do 
    listSum := [op(listSum), OreTools:-Add(listA[i], listB[i])]; 
  end do; 

  return listSum;
end proc:

# function getReplaceMatrixRR
getReplaceMatrixRR := proc (opMatrix, height, nullSpace, m_deg_rows, m_indexRowOrderDiff) 
  local m_matrix; 

  m_matrix := opMatrix; 
  m_matrix[m_indexRowOrderDiff] := getReplaceRowMatrixRR(opMatrix, height, 
                      nullSpace, m_deg_rows, m_indexRowOrderDiff);
  return m_matrix;
end proc:

# function getReplaceRowMatrixRR
getReplaceRowMatrixRR := proc(opMatrix, height, nullSpace, m_deg_rows, m_indexRowOrderDiff) 
  local i, j, m_matrix, size, poly, rowResultSum, row, oreDiffOrder, polyOre; 
  m_matrix := opMatrix; 
  size := height; 
  rowResultSum := [seq(OrePoly(x), `in`(x, convert(vector(size, 0), list)))]; 
  
  for i to size do 
    if nullSpace[i] <> 0 then 
      row := convert(opMatrix[i], list); 
      oreDiffOrder := m_deg_rows[m_indexRowOrderDiff]-m_deg_rows[i]; 
      polyOre := LinearOperators[DEToOrePoly](diff(y(x), [`$`(x, oreDiffOrder)]), y(x)); 
      row := MultiplyOrePolyOnListOre(polyOre, row); 
      poly := nullSpace[i]; 
      row := MultiplyPolyOnListOre(poly, row); 
      rowResultSum := AddListsOrePoly(rowResultSum, row); 
    end if; 
  end do; 
  
  return eqLCMinMatrix(convert(rowResultSum,list)); 
end proc:

# function getUnitMatrix
getUnitMatrix := proc (size) 
  local i, j, matrix; 
  matrix := Matrix(size); 
  for i to size do 
    for j to size do 
      if i <> j then 
        matrix[i, j] := OrePoly(0);
      else 
        matrix[i, j] := OrePoly(1);
      end if; 
    end do; 
  end do; 
  
  return matrix;
end proc: