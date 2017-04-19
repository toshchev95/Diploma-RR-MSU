with(OreTools):
with(LinearAlgebra):
R := SetOreRing(x, 'differential'):

L := Matrix([[OrePoly(x, 0, 0, 1),   OrePoly(0, 0, 2),    OrePoly(x^2+x)], 
             [OrePoly(0, 1),         OrePoly(0, 0, x),    OrePoly(2*x^2+1)], 
             [OrePoly(0, 0, 1),      OrePoly(0, x),       OrePoly(1)]]):

L1 := Matrix([[OrePoly(0, x, 0, 0, 1),   OrePoly(0, 0, 0, 2),    OrePoly(0, x^2+x)], 
              [OrePoly(x, 0, 0, 1),   OrePoly(0, 0, 2),    OrePoly(x^2-x)],
              [OrePoly(0, 0, 1),      OrePoly(0, x),       OrePoly(1)]]):

M  := Matrix([[OrePoly(0, 1, 3, 7),    OrePoly(0, 0, 11),      OrePoly(2*x^2+1),   OrePoly(1)], 
              [OrePoly(x, 0, x),       OrePoly(0, 0, 2, 6),    OrePoly(0, x^2+x),  OrePoly(1)],
              [OrePoly(0, 3),          OrePoly(0, x, 8, 3),    OrePoly(1),         OrePoly(1)],
              [OrePoly(0, 2),          OrePoly(0, x, 8, 3),    OrePoly(1),         OrePoly(1)] ]):

M1 := Matrix([[OrePoly(0, 1, 3, 7),    OrePoly(0, 0, 11),      OrePoly(2*x^2+1),   OrePoly(1)], 
              [OrePoly(x, 0, x),       OrePoly(0, x^2+x),    OrePoly(0, 2, 6),  OrePoly(1)],
              [OrePoly(0, 3),          OrePoly(0, x, 8, 3),    OrePoly(1),         OrePoly(1)],
              [OrePoly(0, x, 8, 3),          OrePoly(x),    OrePoly(1),         OrePoly(1)] ]):

P := Matrix([[OrePoly((1/2)*x^2), OrePoly(1, -(1/2)*x)], [OrePoly(-3, -x), OrePoly(0, 0, 1)]]):

getMatrix := proc()
  return L1;
end proc:

getNewMatrix := proc()
  return M1;
end proc:


#######################
# cal frontal matrix  #
#######################
# function getFrontMatrix
getFrontMatrix := proc(A::Matrix)
  local r, i, front, width, height, ord;
  height, width := op(1, A);  # size matrix of A 
  if height <> width then 
    error "not square matrix!" 
  end if;
  front := Matrix(height, width);
  for r from 1 to height do
    ord := 0;
    for i from 1 to width do
      if nops(OreTools:-Normalize(A[r,i]))>ord then
        ord := nops(OreTools:-Normalize(A[r,i]));
      end if; 
    end do; 
    for i from 1 to width do
      if nops(OreTools:-Normalize(A[r,i]))=ord then
        front[r,i] := op(ord, A[r,i])
      end if;
    end do;
  end do;
  return front;
end proc:

###############################################
### get list Problem of lists Unimod matrix ###
###############################################
# function getListOfProblemUnimodulMatrix
getListOfProblemUnimodulMatrix := proc(K::Matrix, numberOpMatrix::integer) 
  local vlist := list(), 
        height := op(1, K)[1], 
        width := op(1, K)[2];
  local front, nullSpace, fullNullSpace;
  local count, problemList, listUniMatrixs, numberNullSpace := 1;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs;
   
  if height <> height then 
    error "not square matrix!" 
  else
    front := getFrontMatrix(K);
  end if;
  #print("front: ", front);
  fullNullSpace := LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(front));
  #print(Statistics:-Count([op(e)]));
  count := 1;
  problemList := list();
  try

    while member(fullNullSpace[count], fullNullSpace) do 
      #инструкции 
      nullSpace := fullNullSpace[count];
      nullSpace := map(proc (x, y) options operator, arrow; x*y end proc, nullSpace, lcm(seq(x, `in`(x, map(denom, nullSpace)))));
      nullSpace := simplify(nullSpace);

      # UID
      UID_vector := [op(UID_vector), [numberOpMatrix, nullSpace]];

      listUniMatrixs := getListUnimodulMatrix(K, nullSpace, numberOpMatrix);
      problemList := [op(problemList), [nullSpace, listUniMatrixs]];
      #problemList := [op(problemList), op(listUniMatrixs)];

      count := count + 1;

      #print("=============================================================");
    end do;
    print("QA END1");
  catch:
  end try;
  
  if problemList = list() then
    error "problemList is empty";
  end if;

  return problemList;
end proc:



##############################
### get list Unimod matrix ###
##############################
# function getListUnimodulMatrix
getListUnimodulMatrix := proc(K::Matrix, nullSpace::Vector, numberOpMatrix::integer) 
  local x2, indexNullSpace, ord2, j,i,j2, posX, ord1, ordTemp;

  # init
  local vlist := list(), 
        height := op(1, K)[1], 
        width := op(1, K)[2];
  local front, 
        repeatRow := -1, 
        repeatCount := 0;
  local uni, count, problemList, listUniMatrixs, sizeList;
  local porydok, numberList := list();
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs;
  
  with(OreTools):
  x2 := 0;
  ord2 := 0;

  sizeList := Statistics:-Count(nullSpace);
  porydok := vector(sizeList, 0);
  
  for indexNullSpace from 1 to sizeList do

    if nullSpace[indexNullSpace] <> 0 then
      
      posX := indexNullSpace;
      ord1 := 0;
      for j from 1 to height do
        ordTemp := nops(OreTools:-Normalize(K[posX][j]));
        if ordTemp > ord1 then
          ord1 := nops(OreTools:-Normalize(K[posX][j]));
        end if;
      end do;
      porydok[indexNullSpace] := ord1 - 1;
      if ord1 > ord2 then
        ord2 := ord1;
        repeatCount := 0;
        numberList := [indexNullSpace];
      elif ord1 = ord2 then
        repeatCount := repeatCount + 1;
        numberList := [op(numberList), indexNullSpace];
      end if;
    else
      porydok[indexNullSpace] := -1;
    end if;
  end do;

  ord2 := ord2 - 1; #!!!
  repeatCount := repeatCount + 1;
  
  #####################
  #####_CONDITION_#####
  #####################
  for i from 1 to repeatCount do
    uni := getUnimodulMatrix(K, height, nullSpace, porydok, numberList[i]);
    vlist := [op(vlist), uni];

    # UID
    UID_uniMatrix := [op(UID_uniMatrix), [numberOpMatrix, nullSpace, uni]];
  od;
 
  return vlist;
end proc:

# function getUnimodulMatrix
getUnimodulMatrix := proc(K::Matrix, height, nullSpace, porydokDiffRows, ord2)
  local uni, x2, j2;
  
  if ord2 <> 0 then
    
    uni := Matrix(LinearAlgebra:-ScalarMatrix(1, height));
    
    for x2 from 1 to height do
        uni[ord2, x2] := coef(nullSpace, x2, ord2, porydokDiffRows);
    od;
    #print("uni[ord2] = ",uni[ord2], "porydokDiffRows = ", porydokDiffRows);
    for x2 from 1 to height do
      for j2 from 1 to height do
        if uni[x2, j2] = 0 then
          uni[x2, j2] := OrePoly(0);
        elif uni[x2, j2] = 1 then
          uni[x2, j2] := OrePoly(1);
        end if;
      od;
    od;

    return uni;
  else
    error "ERROR !"
  end if;
end proc:

# function coef
coef := proc(C::Vector, x2::integer, ord2::integer, porydok::vector)
  local i,j;
  local x1, n, ode;
  #print("nullSpace = ", C, "row x2 = ", x2,    "ord2 = ",ord2 ,"porydok = ",porydok);
  x1 := C[][x2]; 
  #x1 := C[][x2] / C[][ord2];
  
  if x1=0 then
    return OrePoly(0);
  else
    n := porydok[ord2] - porydok[x2];
    ode := diff(y(x), [`$`(x, n)]);
    return LinearOperators[DEToOrePoly](x1*ode, y(x)); 
  end if;
end proc:

# function getReverseLUMatrix
getReverseLUMatrix := proc(K::Matrix, uni::Matrix)

  local i, j, r;
  local temp; 
  local H := Matrix(op(1, uni)[1], op(1, K)[2]);

  for i to op(1, uni)[1] do 
    for j to op(1, K)[2] do 
      H[i,j] := OrePoly(0);
      for r to op(1, K)[2] do 
        temp := OreTools:-Multiply(uni[i, r], K[r, j], R);
        H[i, j] := OreTools:-Add(temp, H[i,j]);
      end do;
    end do;
  end do;
  return H;
end proc:



############################
##########_REPLACE_#########
############################

#listIK := []; print(RR(P, listIK)):

# function RR
RR := proc(opMatrix::Matrix, listG::list, numberOpMatrix::integer)
  local A,B, saved,saved_poly, nextRR, count;
  local uni, uniList, nextNumber, nextNumberList, estimation;
  local listChain, listGlobal, listGain;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs;


  # init
  A := getFrontMatrix(opMatrix);
  listChain := [];
  #listGain := [];
  listGlobal := [];

  #estimation := estimations(opMatrix, 0);
  count := Statistics:-Count(listG);
  if count > 0 then
    listGain := [op(listG), opMatrix];
  else
    listGain := [opMatrix];
  end if;

  #print(A, LinearAlgebra:-Rank(A),LinearAlgebra[NullSpace](LinearAlgebra:-Transpose(A)));

  if (LinearAlgebra:-Rank(A) < LinearAlgebra:-RowDimension(A)) then
   #  and (estimation) then

    print(0);
    uniList:= getListOfProblemUnimodulMatrix(opMatrix, numberOpMatrix);
    #uni:= getListOfProblemUnimodulMatrix(opMatrix, numberOpMatrix);
    #print(uniList);

    for nextNumberList from 1 to nops(uniList) do

      uni := uniList[nextNumberList][2]; # список унимодулярных матриц зависящих от nullSpace = uniList[nextNumberList][1]

      for nextNumber from 1 to Statistics:-Count(uni) do
      
        saved := getReverseLUMatrix(opMatrix, uni[nextNumber]);
  
        # UID
        UID_opMatrix := [op(UID_opMatrix), saved];
        List_UIDs := [op(List_UIDs), [numberOpMatrix, uniList[nextNumberList][1], uni[nextNumber], nops(UID_opMatrix)]];
  
        nextRR := RR(saved, listGain, nops(UID_opMatrix));
        listChain := nextRR;
        listGlobal := [op(listGlobal), listChain];
      end do;
    end do;

    return op(listGlobal);
    
  else
    return listGain;
  end if;  
end proc:

# function outputRR()
outputRR := proc(opMatrix::Matrix)
  local i,j, listResultOpMatrix;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs;

  # init # list ? vector
  listResultOpMatrix := list();
  UID_opMatrix := list();
  UID_vector := list();
  UID_uniMatrix := list();
  List_UIDs := list();

  UID_opMatrix := [op(UID_opMatrix), opMatrix];
  listResultOpMatrix := RR(opMatrix, [], 1 );

  #convertRR();

  print(UID_opMatrix);

  print(List_UIDs);
  
  #saveAs("C:\\output.txt");
  #return listResultOpMatrix;
end proc:


# function convertRR
convertRR:= proc()
  local i;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs;

  # UID_opMatrix
  for i to nops(UID_opMatrix) do
    UID_opMatrix[i] := convertMatrixOrePolyToPoly(UID_opMatrix[i]);
  end do;

  # UID_uniMatrix
  for i to nops(UID_uniMatrix) do
    UID_uniMatrix[i][3] := convertMatrixOrePolyToPoly(UID_uniMatrix[i][3]);
  end do;

  # List_UIDs
  for i to nops(List_UIDs) do
    List_UIDs[i][3] := convertMatrixOrePolyToPoly(List_UIDs[i][3]);
  end do;

end proc:

# function convertMatrixOrePolyToPoly
convertMatrixOrePolyToPoly := proc(matrix::Matrix) # delta #∂
  return map(proc (x) options operator, arrow; OreTools[Converters]:-FromOrePolyToPoly(x, delta) end proc, matrix);
end proc:

# filename := "C:\\Kursovay\\maple\\Git\\output.txt"
# function saveAs(path) 
SaveAs:=proc(filename)
  local i,j, fd;
  try
    fd := fopen(filename,WRITE,TEXT);
    writedata(fd,UID_opMatrix);
    #writedata(!);
    #ExportMatrix(!)
  catch:
    error "func SaveAs: wrong saved";
  end try;

  fclose(fd);
end proc:

############################
##########_REPLACE_#########
############################

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
