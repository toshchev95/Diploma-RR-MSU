$include "C:\\Kursovay\\maple\\bTest.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

L := Matrix([[OrePoly(x, 0, 0, 1),   OrePoly(0, 0, 2),    OrePoly(x^2+x)], 
             [OrePoly(0, 1),         OrePoly(0, 0, x),    OrePoly(2*x^2+1)], 
             [OrePoly(0, 0, 1),      OrePoly(0, x),       OrePoly(1)]]):

M  := Matrix([[OrePoly(0, 1, 3, 7),    OrePoly(0, 0, 11),      OrePoly(2*x^2+1),   OrePoly(1)], 
              [OrePoly(x, 0, x),       OrePoly(0, 0, 2, 6),    OrePoly(0, x^2+x),  OrePoly(1)],
              [OrePoly(0, 3),          OrePoly(0, x, 8, 3),    OrePoly(1),         OrePoly(1)],
              [OrePoly(0, 2),          OrePoly(0, x, 8, 3),    OrePoly(1),         OrePoly(1)] ]):

P := Matrix([[OrePoly((1/2)*x^2), OrePoly(1, -(1/2)*x)], [OrePoly(-3, -x), OrePoly(0, 0, 1)]]):

Bug1 := Matrix(4, 4, 
          {(1, 1) = 0, (1, 2) = 1, (1, 3) = 0, (1, 4) = 0, 
          (2, 1) = 0, (2, 2) = 15, (2, 3) = 0, (2, 4) = 0, 
          (3, 1) = 5, (3, 2) = 0, (3, 3) = 10, (3, 4) = 6, 
          (4, 1) = 0, (4, 2) = 0, (4, 3) = 0, (4, 4) = 0}):

#	Неоднозначности алгоритма:
#	
#	
#	1) Выбор вектора коэффициентов линейной зависимости строк вырожденной фронтальной матрицы
#	help function s: frontGenerate, hasNull, nullspaceWithoutDenom
#	
#	
#	2) Выбор элемента вектора коэффициентов линейной зависимости строк дифференциального оператора
#	
#	
#	


# function frontGenerate(size, countVectors, highBound, bNullVectors).
# size - Row x Column
# countVectors - count vectors from nullspace
# highBound - high bound of values from matrix
# bNullVectors - if true then include to work hasNull, otherwise hasNull doesn't work

frontGenerate := proc(size, countVectors, highBound::integer, bNullVectors)
  local i, j, front, delta, stepenPoly,polyRand;
  global maxStepenPoly;

  if countVectors > size - 1 or countVectors < 0 then 
    error "there isnot count of vectors for matrix!" 
  end if;

  # init
  maxStepenPoly := 2;
  delta := 1.0 / (countVectors + 1);  
  front := RandomMatrix(size, density = delta, generator = 0 .. 1);

  for i to size do
    for j to size do

      if front[i,j] <> 0 then
      
        stepenPoly := Generate(integer(range=1 .. maxStepenPoly));
        front[i,j] := Generate(polynom(integer(range= 0 - highBound .. highBound),x,degree=stepenPoly));
      end if;
    od;
  od;

  #while LinearAlgebra:-Rank(front) > LinearAlgebra:-RowDimension(front) - countVectors 
  while countVectors <> size - LinearAlgebra:-Rank(front)
        or hasNull(front, bNullVectors) do
    front := RandomMatrix(size, density = delta, generator = 0 .. 1);
    
    for i to size do
      for j to size do

        if front[i,j] <> 0 then
          stepenPoly := Generate(integer(range=1 .. maxStepenPoly));
          front[i,j] := Generate(polynom(integer(range= 0 - highBound .. highBound),x,degree=stepenPoly));
        end if;
      od;
    od;
  end do;

  return front;
end proc:


# hasNull - функция,проверяющая есть ли ненулевые строки в фронтальной матрице?
# function hasNull(matrix, worked)
hasNull := proc(matrix, worked)
  local i, row;
  if not(worked) then return false;
  end if;
  # 1 - Column
  # 2 - row
  # 3 - nothing
  row := ArrayTools[AddAlongDimension](matrix, 2);

  for i to nops(row) do
    if row[i] = 0 then
      return true;
    end if;
  od;
  return false;
end proc:

# function nullspaceWithoutDenom without denominator (знаменатилей)
nullspaceWithoutDenom := proc(front)
  local i, setOfVectors,vector, newVec, sSet;
  sSet := {};
  setOfVectors := LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(front));
  # Избавляемся от знаменатилей в векторах nullspace
  for i to nops(setOfVectors) do
    vector := setOfVectors[i];
    newVec := map(proc (x, y) options operator, arrow; x*y end proc, vector, lcm(seq(x, `in`(x, map(denom, vector)))));
    newVec := simplify(newVec);
    sSet:={op(sSet), newVec};
  od;
  return sSet;
end proc:

# function maxValuesRowsFront(frontMatrix)
maxValuesRowsFront := proc(front)
  local i, listIndex, listValue, list, vMax, indMax;
  listValue := list();
  listIndex := list();
  for i to LinearAlgebra:-RowDimension(front) do
    list := front[i];
    vMax := max(list);
    member(max(list), list, 'indMax');
    listValue := [op(listValue), vMax];
    listIndex := [op(listIndex), indMax];
  od;
  return [listIndex, listValue];
end proc:


# Цикл по каждому вектору из globalNullspace. К завершению цикла выбираем ненулевые коэффициенты вектора,
# номера строк соответствующие этим коэффициентам в векторе удаляем из globalDelta. В это время
# Заполним элементы операторной матрицы тем порядком дифференцирования, который мы захотим выбрать.
# При обработке других векторов коэффициентов ЛЗ фронтальной матрицы порядок дифференцирования
# выбираем меньшего размера, если коэффициентов в globalDelta осталось больше одного.

# function matrixOperatorGenerate
matrixOperatorGenerate := proc(size, countVectors, bUnimodular, highDiff, highBound)
  local front, i, globalList, listIndex, listValue,
        globalNullspace, operatorMatrix, globalDelta, highDifferChangingRow,
        countNonValue, sumDontUseRows, bUnimodular_1;
  local j, k1, y1, vec, high, genStepenMainPoly, genCoeff;
  global vectorUnEnableDiff;
  # \/_init_\/
  operatorMatrix := Matrix(size); # create empty matrix

    front := frontGenerate(size,countVectors, highBound,true);#x1;
    globalDelta := Vector[row](size, 1);
    vectorUnEnableDiff := list();
    
    ######
    #front := Bug1; size := 4;
    #highDiff := 4; highBound := 15;
    ######
    # \/_processing_nullspace_\/
#    if countVectors <> 0 then
#      globalNullspace := [Vector[column](size, 1)];
#      bUnimodular := false;
#    else
#      globalNullspace := nullspaceWithoutDenom(front);
#    end if;
    if countVectors <> 0 then
      globalNullspace := nullspaceWithoutDenom(front);
      for i to nops(globalNullspace) do
        high := getDiffEquation(highDiff);
        vec:=globalNullspace[i];
        countNonValue := 0;
        listIndex := list();
    
        for j to nops(vec) do
          if vec[j] <> 0 then 
            countNonValue:=countNonValue + 1;
            listIndex:=[op(listIndex), j];
          end if;
        od;
  
        # process rows operator Matrix
        sumDontUseRows := sum('globalDelta[k]', k = 1 .. size);
        #print(vec);      #print(listIndex, countNonValue, sumDontUseRows);
  
        if countNonValue > 1 and sumDontUseRows > 1 then
          # conditional about наличие строк операторной матрицы
          # high надо изменить на рандомный
          vectorUnEnableDiff := [op(vectorUnEnableDiff), high];
          if bUnimodular = true then
            bUnimodular_1 := false;
          else
            bUnimodular_1 := true;
          end if;
  
          for j to nops(listIndex) do
            if bUnimodular_1 = true or bUnimodular = true then
              if bUnimodular = false then
                bUnimodular_1 := false;
              end if;
            y1 := listIndex[j];
            print("j= ",j);
            for k1 to size do 
              operatorMatrix[y1,k1] := LinearOperators[DEToOrePoly](front[y1,k1]*high, y(x)); # ??? = 0
              genStepenMainPoly := Generate(integer(range= 0 .. highDiff-1));
              genCoeff := Generate(integer(range= -highBound .. highBound));
              operatorMatrix[y1,k1] := OreTools:-Add(generateOrePoly(genStepenMainPoly, highBound, true, genCoeff), operatorMatrix[y1,k1]);
            od;
            globalDelta[y1] := 0;
  
            end if;
          od;
        end if;
  
      end do;
    end if;
  
    # fill out 
    print(front);
    for i to size do
      if globalDelta[i] <> 0 then
        print(i);
        operatorMatrix[i] := generateRowOrePoly(front[i], highBound, highDiff);
        #highDifferChangingRow := getHighDifferRow(operatorMatrix[i]);
        #vectorUnEnableDiff := [op(vectorUnEnableDiff), highDifferChangingRow];
      end if;
    od;

  return operatorMatrix;
end proc:

# function generateRowOrePoly(vecFront, highBound, highDiff)
# vectorEnableDiff - vector of enable diff elements (non zero)
# highBound - bound of polynomial

generateRowOrePoly := proc(vecFront, highBound, highDiff)
  local i, j, vectorRow, elem, size, stepenDiff, boundCoef;
  global vectorUnEnableDiff, bExistNullDiff;

  stepenDiff := HasDiff(highDiff);
  size := nops(vecFront);
  vectorRow := Vector[row](size);
  
  for i to size do
    
    elem := vecFront[i];
    if (elem <> 0) then
      vectorRow[i] := generateOrePoly(stepenDiff, highBound, true, vecFront[i]);
    else
      vectorRow[i] := generateOrePoly(stepenDiff, highBound, false, 0);
    end if;

  od;
  return vectorRow;
end proc:

# function generateOrePoly(highDiff, highBound, isHighDiffer, optionalElem)
# isHighDiffer - flag: if true then: generating OrePoly stepenDiff = highDiff
#                         otherwise: generating OrePoly stepenDiff < highDiff !

generateOrePoly := proc(highDiff, highBound, isHighDiffer, optionalElem)
  local i, j, stepenDiff, rndom, high, polyOre, highDiff0;
  polyOre := OrePoly(0);

  if not(isHighDiffer) then
    highDiff0 := Generate(integer(range= 0 .. highDiff - 1));
  else
    highDiff0 := highDiff;
  end if;

  for i from 0 to highDiff0 do
 
    stepenDiff := getDiffEquation(i);
    if i = highDiff0 and isHighDiffer then
      rndom := optionalElem;
    else
      rndom := Generate(integer(range= -highBound .. highBound));
    end if;
    polyOre := OreTools:-Add(polyOre, LinearOperators[DEToOrePoly](rndom*stepenDiff, y(x)));
  od;
  
  return polyOre;
end proc:

# function HasDiff(highDiff)
# globalFillDiff - вектор из элементов наивысших порядков дифференцирования
# highDiff - порядок дифференцирования

HasDiff := proc(highDiff)
  local i,j, genValue, flag, size;
  global vectorUnEnableDiff;

  genValue := -1;
  size := nops(vectorUnEnableDiff);
  flag := true;
  
  if (highDiff = size - 1) then
    error "List of values is fill!" 
  end if;
  
  while evalb(flag) do
    genValue := Generate(integer(range= 1 .. highDiff));
    #for i to size do
    if not(member(genValue, vectorUnEnableDiff)) then
      flag := false;
      break;
    end if;
    #od;
  od;

  if genValue = -1 then
    error "func HasDiff: no exist right order of differential";
  end if;

  vectorUnEnableDiff := [op(vectorUnEnableDiff), genValue];
  return genValue;
end proc:


# function getDiffEquation
getDiffEquation := proc(highDiff)
  return diff(y(x), [`$`(x, highDiff)]);
end proc:

# function getHighDifferRow(list)
getHighDifferRow := proc(list)
  local i, ord;
  ord := 0;
    for i from 1 to nops(list) do
      if nops(OreTools:-Normalize(list[i]))>ord then
        ord := nops(OreTools:-Normalize(list[i]));
      end if; 
    end do; 
  return ord - 1;
end proc:



# @parametrs:
# 0. Параметр по размеру матрицы size: size * size
# 1. Параметр по итерациям: 
# 2. Параметры по выбору не однозначности (int, int), 
#    2ой параметр - количество векторов линейной зависимости фронтальной матрицы
#  3ий параметр - количество строк с равным максимальным порядком дифференцирования в операторной матрице
#  Максимальный порядок дифференцирования варьируется в пределах от 1 .. 6
# @return: 
#  Операторная матрица (дифференциальная)

# function GenerateMatrixRR
GenerateMatrixRR := proc(size, iter, highDiff::integer)      #counter_vectors::integer) option overload; #counter_strings
  local i, m_Gen, m_prevGen, m_stepenPoly, m_highDiff, m_UniMatrix, two_vectors,
    m_indexMaxDiffOrderInList,m_getMatrix, bFirstIterWhile;
  global m_vector_LZ, m_list_MaxDiffOrder;

  # Сначала генерируем матрицу, которая невырождена (/= 0, ранк == размера)
  m_highDiff := highDiff;

  m_Gen := matrixOperatorGenerate(size, 0, false, m_highDiff, 10);
  m_prevGen := m_Gen;
  print(m_Gen, getFrontMatrix(m_Gen), nullspaceWithoutDenom(getFrontMatrix(m_Gen)));

  for i to iter do
    m_getMatrix := m_Gen;
    bFirstIterWhile := true;

    while LinearAlgebra:-Rank(getFrontMatrix(m_getMatrix)) = size
      or LinearAlgebra:-Equal(m_getMatrix,m_prevGen) or bFirstIterWhile do

      bFirstIterWhile := false;
      print("while");
      
      # Выполняем алгоритм Row-Reduction в обратном порядке с учётом обратной унимодулярной матрицы.
      # Сгенерируем обратную унимодулярную матрицу (хз как) и умножим на невырожденную матрицу
      # Для этого мы создадим два списка размера size, 
      # 1. один список - вектор коэффициентов ЛЗ строк фронтальной матрицы
      # 2. второй - вектор, элементами которого являются максимальные порядки дифференцирования строк линейного оператора
      m_stepenPoly := Generate(integer(range = 0 .. 1));
      m_vector_LZ := Generate(list(polynom(integer(range = -10 .. 10), x, degree = m_stepenPoly), size));
      m_list_MaxDiffOrder := Generate(list(integer(range= 0..m_highDiff), size));
      
      # Построим вектор без однозначностей
      m_indexMaxDiffOrderInList := simultaneouslyProcessLists();

      # Создадим обратную операторную унимодулярную матрицу 
      m_UniMatrix := getReverseUnimodularMatrix(m_indexMaxDiffOrderInList);
  
      # Compute L' = m_UniMatrix*m_Gen
      # Аналог function getReverseLUMatrix(L,U) = U*L
      m_getMatrix := getReverseLUMatrix(m_Gen,m_UniMatrix);
    end do;
    
    m_prevGen := m_Gen;
    m_Gen := m_getMatrix;
    print(i,m_Gen,LinearAlgebra:-Equal(m_Gen,m_prevGen));

    # Пока не понятно, как обработать параметры
    # <!> Возникает сложность, связанная с получением неоднозначности, у которой количество векторов ЛЗ будет больше одного
  od;

  return m_Gen;
  #return m_getMatrix;
end proc:

# function simultaneouslyProcessLists()
simultaneouslyProcessLists := proc()
  local i,j,listNumber,maxDiffOrder,size, countMaxDiffOrder, list_MaxDiffOrder, vector_LZ,
    indexMaxDiffOrderInList,listEmptyNumber;
  global m_vector_LZ, m_list_MaxDiffOrder;

  # := formal parametrs
#  list_MaxDiffOrder := m_list_MaxDiffOrder;
#  vector_LZ := m_vector_LZ;

  # init
  size := nops(m_vector_LZ);
  listNumber := list();
  listEmptyNumber := list();
  maxDiffOrder := 0;
  countMaxDiffOrder := 0;
  indexMaxDiffOrderInList := 0;

  for i to size do
    if m_vector_LZ[i] <> 0 then

      # по логике алгоритма нам не критично выбирать максимальный порядок из всех 
      # максимальных порядков ненулевых строк операторной матрицы
      if maxDiffOrder < m_list_MaxDiffOrder[i] then 
        maxDiffOrder := m_list_MaxDiffOrder[i];
        countMaxDiffOrder := 1;
        indexMaxDiffOrderInList := i;
      else
        countMaxDiffOrder := countMaxDiffOrder + 1;
      end if;
      
      listNumber := [op(listNumber), i];
    else
      listEmptyNumber := [op(listEmptyNumber), i];
    end if;
  od;

  if countMaxDiffOrder <= 1 then
#  while countMaxDiffOrder <= 1 do 
    for i to nops(listNumber) do
      j := listNumber[i];
      if m_list_MaxDiffOrder[j] <> maxDiffOrder then

        #m_vector_LZ[j] := Generate(integer(range=-10..10));

        m_list_MaxDiffOrder[j] := maxDiffOrder;
        #break; # Можно самим установить количетсво одинаковых максимальных порядков строк
      end if;

    od;

    if countMaxDiffOrder <= 1 then
      for i to nops(listEmptyNumber) do

        j := listEmptyNumber[i];
        m_vector_LZ[j] := Generate(integer(range=-10..10));
        m_list_MaxDiffOrder[j] := maxDiffOrder;
        countMaxDiffOrder := countMaxDiffOrder + 1;
        
        if countMaxDiffOrder > 1 then
          break;
        end if;

      end do;
    end if;

#  end do;
  end if;

  return indexMaxDiffOrderInList;
  #return [vector_LZ, list_MaxDiffOrder, indexMaxDiffOrderInList];
end proc:

# function getReverseUnimodularMatrix()
getReverseUnimodularMatrix := proc(m_indexMaxDiffOrderInList::integer)
  local i,j, m_reverseUniMatrix,size;
  global m_vector_LZ, m_list_MaxDiffOrder;

  # init
  size := nops(m_vector_LZ);
  m_reverseUniMatrix := Matrix(LinearAlgebra:-ScalarMatrix(1, size));

  for i to size do

    for j to size do
      if m_indexMaxDiffOrderInList <> i then
        
        if i <> j then
          m_reverseUniMatrix[i,j] := OrePoly(0);
        else
          m_reverseUniMatrix[i,j] := OrePoly(1);
        end if;
      else
        # formul
        #m_reverseUniMatrix[m_indexMaxDiffOrderInList,j] := OrePoly
    m_reverseUniMatrix[m_indexMaxDiffOrderInList,j] := coefForReverseConversation(
                                      convert(m_vector_LZ,Vector), 
                                      j, 
                                      m_indexMaxDiffOrderInList, 
                                      convert(m_list_MaxDiffOrder,vector));
      end if;
    end do;
  end do;

  return m_reverseUniMatrix;
end proc:

