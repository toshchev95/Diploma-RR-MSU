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
        width := op(1, opMatrix)[2],
        size := height;
  global front, m_matrix, fullNullSpace,m_deg_rows, m_infoOnMatrix, m_matrixInfo;

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
    listNumbersRowsColsInfo, size, listOfOrderDiff, listOfPoly, listEmpty, listOrderDiffRowUniMat;
  global m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, m_matrixInfo, # global for modifyRR()
    m_RowsInfo, m_ColsInfo,m_listRowsInfoUniMatrix, nullSpace;

  size := LinearAlgebra:-RowDimension(m_matrix);
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

  # Вычислим матрицу, где элемент содержит инфу об OrePoly в виде 2-х списков:
  # listOfOrderDiff, listOfPoly
  listEmpty := convert(vector(highDiffMatrix,0),list); # [0,0,...,0]

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




#Optional estimations are here!


# function estimations
estimations := proc()
  local size,i,j, estResult, nullSpace_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff;
  global front, m_matrix, fullNullSpace,m_deg_rows;

  size := Matrix(op(1, opMatrix));

  # обработка всех имеющихся параметров
  # process();
  # выбор вектора ЛЗ и индекс максимального дифференциального оператора
  # m_nullSpace := selectNullspace();
  # m_indexRowOrderDiff := selectIndexDiffOrderRows();

  nullSpace_indexRowOrderDiff := [m_nullSpace, m_indexRowOrderDiff];
  return nullSpace_indexRowOrderDiff;
end proc:

# function selectNullspace
selectNullspace := proc()
  local nullSpace, i;
  global front, m_matrix, fullNullSpace;


  return nullSpace;
end proc:


# function selectIndexDiffOrderRows
selectIndexDiffOrderRows := proc()
  local m_index, i;
  global front, m_matrix, fullNullSpace;


  return m_index;
end proc:

##########################################
###### Init functions RR algorithm ######
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
  
  m_matrix[m_indexRowOrderDiff] := convert(rowResultSum, Vector); 
  return m_matrix;
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