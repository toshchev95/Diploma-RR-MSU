$include "C:\\Kursovay\\maple\\Git\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Git\\Generate.mpl"
$include "C:\\Kursovay\\maple\\Git\\ImprovingRR.mpl"
$include "C:\\Kursovay\\maple\\Git\\bTestBetterResult.mpl"
$include "C:\\Kursovay\\maple\\Git\\AbramovOreStyle.mpl"
$include "C:\\Kursovay\\maple\\Git\\egrr.mpl"
with(EGRR):
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

read("C:\\Kursovay\\maple\\Git\\bTest.mpl");
read("C:\\Kursovay\\maple\\Git\\Generate.mpl");
read("C:\\Kursovay\\maple\\Git\\ImprovingRR.mpl");
read("C:\\Kursovay\\maple\\Git\\bTestBetterResult.mpl");
read("C:\\Kursovay\\maple\\Git\\AbramovOreStyle.mpl");
read("C:\\Kursovay\\maple\\Git\\egrr.mpl");
# read("C:\\Kursovay\\maple\\Git\\Testing.mpl");

# Bugs
# Error, (in getOreStyle) In func getOreStyle: row and column dimensions are wrong
# 

# function test
test := proc(m, r, iter, bOptionalRandMat)
	local i,j,k, m_Gen, m_Opt, m_matrix, m_start, m_start_explicit, m_better,
		rr_matrix, triangle_rr_matrix, rr_matrix_explicit, triangle_rr_matrix_explicit, rr_explicit, triangle_rr_explicit,
		bUsing, countRR, countTriangleRR, countMyRR, countTests, countUsingAlgoritm, countLessIteration, countIteration,
		highDiff;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs, UID_Results, UID_iteration,UID_bIteration, 
  	UID_betterNumber, UID_equalMatrix, UID_using;

  countRR := 0; 
  countTriangleRR := 0; 
  countMyRR := 0;
  countUsingAlgoritm := 0;
  countLessIteration := 0;
  countTests := 20;
  countIteration := 0;

  for i to countTests do
  	if bOptionalRandMat = true then
  		m_Opt := matrixOperatorGenerate(m, Generate(integer(range = 1 .. trunc(m/2)+1 )), true, r, 30):
  		m_Gen := GenerateMatrixRR(m, iter, r, bOptionalRandMat, m_Opt):
  	else
  		m_Gen := matrixOperatorGenerate(m, Generate(integer(range = 1 .. trunc(m/2)+1 )), true, r, 30):
  	end if;

	 	m_start := copy(matrixOreWithoutGCD(matrixOreWithoutDenom(m_Gen))):
	 	highDiff := getOrderDiffMatrix(m_start);
#	 	print("m_start=",m_start,highDiff);
	
	 	# myself
	 	UID_bIteration := false;
	 	outputRR(m_start):
	 	if nops(UID_Results) = 1 then
	 		i := i - 1;
	 		print("Repeat Test");
	 		next;
	 	end if;
	 	if UID_bIteration = true then
	 		countLessIteration := countLessIteration + 1;
	 	end if;
    countIteration := countIteration + nops(UID_Results);
	
  	m_better := getBetterResultMatrix(UID_Results, UID_iteration):
#  	print("m_better=",m_better);
  	UID_using := false;
  	m_matrix := modifyRR(m_start):
#    print("m_matrix=",m_matrix);
  	if UID_using = true then
  		countUsingAlgoritm := countUsingAlgoritm + 1;
  	end if;

  	# RR, TriangleRR
  	#print("Abramov");
  	#m_start_explicit := getAbramovStyle(m_start);
		#print("m_start_explicit");
  	
  	#rr_matrix_explicit := RR(m_start_explicit, highDiff + 1, x);
  	#print("rr_matrix_explicit");
  	#triangle_rr_matrix_explicit := TriangleRR(m_start_explicit, highDiff + 1, x);
  	#print("triangle_rr_matrix_explicit");
		#rr_explicit := getOreStyle(rr_matrix_explicit[1], highDiff);
		#triangle_rr_explicit := getOreStyle(triangle_rr_matrix_explicit[1], highDiff);
  	#rr_matrix := matrixOreWithoutGCD(matrixOreWithoutDenom(rr_explicit));
  	#print("rr_matrix=",rr_matrix);
  	#triangle_rr_matrix := matrixOreWithoutGCD(matrixOreWithoutDenom(triangle_rr_explicit));
  	#print("triangle_rr_matrix=",triangle_rr_matrix);

  	# check
  	if LinearAlgebra[Equal](m_matrix, m_better) = false or UID_bIteration = true then
  		###printUID_opMatrix();
  		###print(List_UIDs);
  		###print("End Results");
  		###printList(UID_Results, true, 1);
  		###print("numberIteration=",UID_iteration);
  		###print("Better matrix is ", m_better, UID_betterNumber);
  	else
  		countMyRR := countMyRR + 1;
		end if;
	
		#UID_equalMatrix := false;
		#compareOreMatrix(m_better, rr_matrix);
		#if UID_equalMatrix = true then
		#	countRR := countRR + 1;
		#end if;
	
		#UID_equalMatrix := false;
		#compareOreMatrix(m_better, triangle_rr_matrix);
		#if UID_equalMatrix = true then
		#	countTriangleRR := countTriangleRR + 1;
		#end if;
	end do;

	print("percent MyRR= ", evalf(countMyRR/countTests));
	#print("percent AbrRR= ", evalf(countRR/countTests));
	#print("percent AbrTrRR= ", evalf(countTriangleRR/countTests));
  print("percent different number Iteration= ", evalf(countLessIteration/countTests));
	print("percent UsingAlgoritm= ", evalf(countUsingAlgoritm/countTests));
  print("count end result matrixes=", evalf(countIteration/countTests));
  
end proc:

# function getOrderDiffMatrix
getOrderDiffMatrix := proc(opMatrix::Matrix)
	local i,size, order, temp, m_deg_rows;
  size := op(1, opMatrix)[1]; 
  m_deg_rows := list(); 
  for i to size do 
      temp := getHighDifferRow2(opMatrix[i]); 
      m_deg_rows := [op(m_deg_rows), temp] 
  end do; 
  #print("m_deg_rows",m_deg_rows);
  return max(m_deg_rows);
end proc:


# function testFindEvent
testFindEvent := proc(m, r, iter, bOptionalRandMat)
  local i,j,k, m_Gen, m_Opt, m_matrix, m_start, m_start_explicit, m_better,
    rr_matrix, triangle_rr_matrix, rr_matrix_explicit, triangle_rr_matrix_explicit, rr_explicit, triangle_rr_explicit,
    bUsing, countRR, countTriangleRR, countMyRR, countTests, countUsingAlgoritm, countLessIteration, countIteration,
    highDiff, temp;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs, UID_Results, UID_iteration,UID_bIteration, 
    UID_betterNumber, UID_equalMatrix, UID_using;

  countTests := 1;

  for i to countTests do
    if bOptionalRandMat = true then
      m_Opt := matrixOperatorGenerate(m, 2, true, r, 10):
      m_Gen := GenerateMatrixRR(m, iter, r, bOptionalRandMat, m_Opt):
    else
      m_Gen := matrixOperatorGenerate(m, 2, true, r, 10):
    end if;

    m_start := copy(matrixOreWithoutGCD(matrixOreWithoutDenom(m_Gen))):
    highDiff := getOrderDiffMatrix(m_start);
#   
  
    # myself
    UID_bIteration := false;
    temp := outputRR(m_start):
    if nops(UID_Results) = 1 then
      i := i - 1;
      print("Repeat Test");
      next;
    end if;

    if UID_bIteration = true then
      next;
    end if;

    print("m_start=",m_start,UID_iteration, nops(UID_Results), nops(UID_vector));
    graphVisualisation(createListEdgesFromNumberOpMatrix());


    if true or UID_iteration = 2 and evalb(nops(UID_Results) >= 3) and nops(UID_vector) >= 4 
      #and ( LinearAlgebra[Equal](convert(UID_vector[1][2], vector), convert(UID_vector[2][2], vector))
      #or LinearAlgebra[Equal](convert(UID_vector[2][2], vector), convert(UID_vector[3][2], vector))
      #or LinearAlgebra[Equal](convert(UID_vector[1][2], vector), convert(UID_vector[3][2], vector)) ) 
      then
      

      printUID_opMatrix();
      #print(List_UIDs);
      print("End Results");
      printList(UID_Results, true, 1);
      print("numberIteration=",UID_iteration);
      #print(UID_vector);
    
      break;
    else
      next;
    end if;


    
  end do;

  
end proc:

# function outputTest
outputTest := proc()
      printUID_opMatrix();
      print(List_UIDs);
      print("End Results");
      printList(UID_Results, true, 1);
      print("numberIteration=",UID_iteration);
    
      graphVisualisation(createListEdgesFromNumberOpMatrix());
end proc: