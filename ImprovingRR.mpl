$include "C:\\Kursovay\\maple\\bTest.mpl"
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

# function modifyRR
modifyRR := proc(opMatrix::Matrix)
  local A,B,i, saved,nullSpace_indexRowOrderDiff, firstListNumberHighDiffForUniMat;
  local vector_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff, uni;
  local height := op(1, opMatrix)[1], 
        width := op(1, opMatrix)[2];
  global front, m_matrix, fullNullSpace,m_deg_rows, m_infoOnMatrix, m_matrixInfo,
    size := height;

  if height <> width then
    error "Func modifyRR: wrong scale opMatrix <- height =/= width";
  end if;

  # init
  front := getFrontMatrix(opMatrix);
  m_matrix := opMatrix;

  # Statistics:-Count
  #print(A, LinearAlgebra:-Rank(A),LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(A)));

  while (LinearAlgebra:-Rank(front) < LinearAlgebra:-RowDimension(front)) do 

    # Data collection
    fullNullSpace := LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(front));
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
      m_nullSpace := vector_indexRowOrderDiff[1];
      m_indexRowOrderDiff := vector_indexRowOrderDiff[2];
    end if;

    uni := getUnimodulMatrix(m_matrix, height, m_nullSpace, m_deg_rows, m_indexRowOrderDiff);
    saved := getReverseLUMatrix(m_matrix, uni);

    # check again
    front := getFrontMatrix(saved);
    m_matrix := saved;
  end do;  

  return m_matrix;
end proc:



# Функция, собирающая все данные операторной матрицы, 
# строк унимод и полученных на следующем шаге матриц
# ! доработать
# function getInfoOnOpMatrix
getInfoOnOpMatrix := proc()
  local i,j, listNumberRowsForUniMatrix, numbersRowsNullSpace, temp, maxOrd, highDiffMatrix,
    listNumbersRowsColsInfo, listOfOrderDiff, listOfPoly, listEmpty, listOrderDiffRowUniMat;
  global m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, size, # global for modifyRR()
    m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix, nullSpace;

  highDiffMatrix := max(m_deg_rows);
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
      for i to nops(listOfPoly) do 
        if listOfPoly[i] <> 0 then 
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

# function getListNumberHighDiffForUniMat
getListNumberHighDiffForUniMat := proc(nullSpace, m_deg_rows)
  local j, numbersRowsNullSpace, maxOrd;
  
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
    nullSpace, listPar_s, tempPar_s, betterParameters;
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
      newRow := getReplaceRowMatrixRR(opMatrix, size, nullSpace, m_deg_rows, numberRowForNullSpace);

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
  return [m_nullSpace, m_indexRowOrderDiff];
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
  if listPar_A[3] <> listPar_A[3] then

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
  local i,j, m_resRowInfoA, m_resRowInfoB;
  global size,m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, front,
    indexA, indexB;

  # m_infoOnMatrix := [m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix,listNumberRowsForUniMatrix];
  m_resRowInfoA := getNonNullList(m_infoOnMatrix[1][indexB]);
  m_resRowInfoB := getNonNullList(m_infoOnMatrix[1][indexA]);

  # a) сумма всех порядков дифференцирования (<)
  # b) кол-во термов порядков диф-ния (>)
  # c) степени у множителей коэффициентов полиномов Оре (<)
  # d) кол-во термов f(x) (<)
  # e) числа у множителей (<)
  
  return compareOrePoly(m_resRowInfoA,m_resRowInfoB);
end proc:

# function compareOrePoly
compareOrePoly := proc (oreA, oreB) 
  local i, j, listA, listB, sizeA, sizeB, sumPolyA, sumPolyB, listIsDiffA, listIsDiffB, 
        numberDiffA, numberDiffB, sumPlusCoeffsA, sumPlusCoeffsB; 
  listA := [seq(x, `in`(x, oreA))]; 
  listB := [seq(x, `in`(x, oreB))]; 
  sizeA := nops(listA); 
  sizeB := nops(listB); 
  sumPolyA := sum('listA[k]', k = 1 .. sizeA); 
  sumPolyB := sum('listB[k]', k = 1 .. sizeB); 

  listIsDiffA := map(proc (x) options operator, arrow; if x = 0 then return 0 else return 1 end if end proc, listA); 
  listIsDiffB := map(proc (x) options operator, arrow; if x = 0 then return 0 else return 1 end if end proc, listB); 
  numberDiffA := sum('listIsDiffA[k]', k = 1 .. sizeA); 
  numberDiffB := sum('listIsDiffB[k]', k = 1 .. sizeB); 

  sumPlusCoeffsA := sum('seq(abs(c), `in`(c, coeffs(sumPolyA)))[k]', k = 1 .. nops(sumPolyA)); 
  sumPlusCoeffsB := sum('seq(abs(c), `in`(c, coeffs(sumPolyB)))[k]', k = 1 .. nops(sumPolyB)); 

  if nops(oreA) < nops(oreB) then 
    return true 
  elif nops(oreB) < nops(oreA) then 
    return false 
  elif numberDiffA < nuberDiffB then 
    return true 
  elif nuberDiffB < numberDiffA then 
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
    return true 
  elif sumPlusCoeffsB < sumPlusCoeffsA then 
    return false 
  else 
    # using sumListDifferentLength(), coeffFull()
    print("In func compareOrePoly: else"); 
    return true 
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
                      nullSpace, m_deg_rows, m_indexRowOrderDiff) ;
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
  
  return convert(rowResultSum, Vector); 
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