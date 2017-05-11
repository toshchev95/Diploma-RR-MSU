$include "C:\\Kursovay\\maple\\Git\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Git\\Generate.mpl"
$include "C:\\Kursovay\\maple\\Git\\ImprovingRR.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

# function getAbramovStyle
getAbramovStyle := proc (opMatrix::Matrix) 
  local i, j, k, size, listMatrixDiff, temp, m_deg_rows, maxDiff, emptyList, explicitListsMat, polyList; 
  size := op(1, opMatrix)[1]; 
  m_deg_rows := list(); 
  for i to size do 
    temp := getHighDifferRow2(opMatrix[i]); 
    m_deg_rows := [op(m_deg_rows), temp] 
  end do; 
  maxDiff := max(m_deg_rows);
  print(maxDiff); 
  emptyList := convert(vector(maxDiff, 0), list); 
  listMatrixDiff := convert(vector(maxDiff + 1, 0), list); 
  explicitListsMat := convertOreMatrixToExplicitMatrix(opMatrix, maxDiff); 
  for k from 0 to maxDiff do 
    temp := LinearAlgebra[Copy](Matrix(size)); 
    for i to size do 
      for j to size do 
        polyList := explicitListsMat[i, j]; 
        temp[i, j] := polyList[k+1];
      end do; 
    end do; 
  listMatrixDiff[k+1] := temp; 
  end do; 
  return Matrix(ListTools:-Reverse(listMatrixDiff));
end proc:

# function convertOreMatrixToExplicitMatrix
convertOreMatrixToExplicitMatrix := proc (opMatrix::Matrix, maxDiff) 
  local i, j, k, size, emptyList, explicitMat, polyList; 
  size := op(1, opMatrix)[1]; 
  explicitMat := Matrix(size); 
  emptyList := convert(vector(maxDiff + 1, 0), list); 
  for i to size do 
    for j to size do 
      polyList := [seq(x, `in`(x, opMatrix[i, j]))]; 
      explicitMat[i, j] := sumListDifferentLength(emptyList, polyList); 
    end do; 
  end do; 
  return explicitMat; 
end proc:

# function getOreStyle
getOreStyle := proc (explMat::Matrix, r::integer) # r - DiffOrder
  local i, j, k, size, temp, oreMat, polyElem, length, urez; 
  size := op(1, explMat)[2]/(r+1); 
  if frac(size) <> 0 then 
    print(size, frac(size)); 
    error "In func getOreStyle: parameter r is wrong"; 
  elif size <> op(1, explMat)[1] then 
    print(size, op(1, explMat)[1]); 
    error "In func getOreStyle: row and column dimensions are wrong"; 
  end if; 
  oreMat := Matrix(size); 
  for i to size do 
    for j to size do 
      length := 0; 
      temp := vector(r+1, 0); 
      for k from 0 to r do 
        temp[k+1] := explMat[i][(r-k)*size+j]; 
      end do; 
      length := getNops(temp); 
      if length > 0 then
        urez := convert(temp, list)[1 .. length]; 
      else
        urez := [0];
      end if;
      oreMat[i, j] := OrePoly(op(convert(urez, list))); 
    end do; 
  end do; 
  return oreMat; 
end proc: 

# function getNops
getNops := proc (struct0) 
  local i, size, list0, temp; temp := 0; 
  list0 := convert(struct0, list); 
  size := Statistics:-Count(list0); 
  for i to size do 
    if list0[size-i+1] <> 0 then 
      break; 
    end if; 
    temp := temp+1; 
  end do; 
  return size-temp; 
end proc: