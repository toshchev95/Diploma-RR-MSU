
# function modifyRR
modifyRR := proc(opMatrix::Matrix)
  local A,B,i, saved,m_deg_rows;
  local vector_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff, uni;
  local height := op(1, opMatrix)[1], 
        width := op(1, opMatrix)[2],
        size := height;
  global front, m_matrix, fullNullSpace;

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

    # issues modify
    nullSpace_indexRowOrderDiff := estimations();
    m_nullSpace := vector_indexRowOrderDiff[1];
    m_indexRowOrderDiff := vector_indexRowOrderDiff[2];

    for i to size do
      m_deg_rows[i] := getHighDifferRow(m_matrix[i]);
    end do;

    uni := getUnimodulMatrix(m_matrix, height, m_nullSpace, m_deg_rows, m_indexRowOrderDiff);
    saved := getReverseLUMatrix(m_matrix, uni);

    # check again
    front := getFrontMatrix(saved);
    m_matrix := saved;
  end do;  

  return m_matrix;
end proc:


# function estimations
estimations := proc();
  local size,i,j, estResult, nullSpace_indexRowOrderDiff, m_nullSpace, m_indexRowOrderDiff;
  global front, m_matrix, fullNullSpace;

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

