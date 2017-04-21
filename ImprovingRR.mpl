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
  local A,B,i, saved,nullSpace_indexRowOrderDiff;
  local vector_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff, uni;
  local height := op(1, opMatrix)[1], 
        width := op(1, opMatrix)[2],
        size := height;
  global front, m_matrix, fullNullSpace,m_deg_rows, m_infoOnMatrix;

  if height <> width then
    error "Func modifyRR: wrong scale opMatrix <- height =/= width";
  end if;

  # init
  front := getFrontMatrix(opMatrix);
  m_matrix := opMatrix;

  #estimation := estimations(opMatrix, 0);
  # Statistics:-Count
  #print(A, LinearAlgebra:-Rank(A),LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(A)));

  while (LinearAlgebra:-Rank(front) < LinearAlgebra:-RowDimension(front)) do
   #  and (estimation) then
    m_deg_rows := vector(size, 0);
    fullNullSpace := LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(front));

    for i to size do
      m_deg_rows[i] := getHighDifferRow(m_matrix[i]);
    end do;

    # Info
    m_infoOnMatrix := list();
    getInfoOnOpMatrix();

    # issues modify
    nullSpace_indexRowOrderDiff := estimations();
    m_nullSpace := vector_indexRowOrderDiff[1];
    m_indexRowOrderDiff := vector_indexRowOrderDiff[2];

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

# function getInfoOnOpMatrix
getInfoOnOpMatrix := proc()
  local i,j, listNumberRowsForUniMatrix, numbersRowsNullSpace, 
    listNumbersRowsColsInfo;
  global m_matrix, m_infoOnMatrix, m_deg_rows,fullNullSpace, nullSpace;

  # Список номеров строк для неоднозначности для получения унимод матриц
  listNumberRowsForUniMatrix := list();
  # Список номеров строк, столбцов операторной матрицы, по к. есть неоднозначности uniMatrix
  listNumbersRowsColsInfo := list();
  for i to nops(fullNullSpace) do
    nullSpace := fullNullSpace[i];
    numbersRowsNullSpace := list();
    for j to nops(nullSpace) do
      if nullSpace[j] <> 0 then
        numbersRowsNullSpace := [op(numbersRowsNullSpace),j];
        listNumbersRowsColsInfo := [op(listNumbersRowsColsInfo),j];
      end if;
    end do;
    listNumberRowsForUniMatrix := [op(listNumberRowsForUniMatrix),numbersRowsNullSpace];
  end do;
  listNumbersRowsColsInfo := ListTools:-MakeUnique(listNumbersRowsColsInfo);

  # Вычислим матрицу, где элемент 
end proc:

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


#Optional estimations are here!


