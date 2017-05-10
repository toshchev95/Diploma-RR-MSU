$include "C:\\Kursovay\\maple\\Git\\bTest.mpl"
$include "C:\\Kursovay\\maple\\Git\\Generate.mpl"
$include "C:\\Kursovay\\maple\\Git\\ImprovingRR.mpl"
$include "C:\\Kursovay\\maple\\Git\\bTestBetterResult.mpl"
$include "C:\\Kursovay\\maple\\Git\\AbramovOreStyle.mpl"
$include "C:\\Kursovay\\maple\\CA\\egrr.mpl"
with(OreTools):
with(LinearAlgebra):
with(RandomTools):
with(ArrayTools):
with(ListTools):
R := SetOreRing(x, 'differential'):

# function test
test := proc(m, r, iter, bOptionalRandMat)
	local i,j,k, size, m_Gen, m_Opt, m_matrix, m_start, m_start_explicit, m_better,
		rr_matrix, triangle_rr_matrix, rr_matrix_explicit, triangle_rr_matrix_explicit,
		bUsing, countRR, countTriangleRR, countMyRR, countTests, countUsingAlgoritm;
  global UID_opMatrix, UID_vector, UID_uniMatrix, List_UIDs, UID_Results, UID_iteration,UID_bIteration, 
  	UID_betterNumber, UID_equalMatrix, UID_using;

  countRR := 0; 
  countTriangleRR := 0; 
  countMyRR := 0;
  countUsingAlgoritm := 0
  countTests := 10;

  for i to countTests do
  	if bOptionalRandMat = true then
  		m_Opt := matrixOperatorGenerate(size, Generate(integer(range = 1 .. trunc(m/2)+1 )), true, r, 30);
  	end if;
	
	  m_Gen := GenerateMatrixRR(m, iter, r, bOptionalRandMat, m_Opt);
	
	 	m_start := matrixOreWithoutGCD(matrixOreWithoutDenom(m_Gen));
	 	print("m_start=",m_start);
	
	 	# myself
	 	outputRR(m_start);
	
  	m_better := getBetterResultMatrix(UID_Results, UID_iteration);
  	UID_using := false;
  	m_matrix := modifyRR(m_start);
  	if UID_using = true then
  		countUsingAlgoritm := countUsingAlgoritm + 1;
  	end if;
	
  	# RR, TriangleRR
  	m_start_explicit := getAbramovStyle(m_start);
  	rr_matrix_explicit := RR(m_start_explicit, r + 1, x);
  	triangle_rr_matrix_explicit := TriangleRR(m_start_explicit, r + 1, x);
	
  	rr_matrix := matrixOreWithoutGCD(matrixOreWithoutDenom(getOreStyle(rr_matrix_explicit[1], r)));
  	triangle_rr_matrix := matrixOreWithoutGCD(matrixOreWithoutDenom(getOreStyle(triangle_rr_matrix_explicit[1], r)));

  	# check
  	if LinearAlgebra[Equal](m_matrix, m_better) = false or UID_bIteration = true then
  		printUID_opMatrix();
  		print(List_UIDs);
  		print("End Results");
  		printList(UID_Results, true, 1);
  		print("numberIteration=",UID_iteration);
  		print("Better matrix is ", m_better, UID_betterNumber);
  	else
  		countMyRR := countMyRR + 1;
		end if;
	
		UID_equalMatrix := false;
		compareOreMatrix(m_better, rr_matrix);
		if UID_equalMatrix = true then
			countRR := countRR + 1;
		end if;
	
		UID_equalMatrix := false;
		compareOreMatrix(m_better, triangle_rr_matrix);
		if UID_equalMatrix = true then
			countTriangleRR := countTriangleRR + 1;
		end if;
	end do;

	print("percent MyRR= ", countMyRR/countTests);
	print("percent AbrRR= ", countRR/countTests);
	print("percent AbrTrRR= ", countTriangleRR/countTests);
	print("percent UsingAlgoritm= ", countUsingAlgoritm/countTests);
end proc: