using FinancialMonteCarlo
path1=joinpath(dirname(pathof(FinancialMonteCarlo)),"..","test")
test_listTmp=readdir(path1);
BlackList=["REQUIRE","runtests.jl","runner.jl","test_black_cuda.jl","cuda"]
if(!DIFFEQ_MONTECARLO_ACTIVE_FLAG)
	BlackList=vcat("test_diffeq_log.jl","test_diffeq.jl","test_diff_eq_montecarlo_array.jl",BlackList)
end
test_list=[test_element for test_element in test_listTmp if !Bool(sum(test_element.==BlackList))]
println("Running tests:\n")
for (current_test,i) in zip(test_list,1:length(test_list))
	println("------------------------------------------------------------")
    println("  * $(current_test) *")
    include(joinpath(path1,current_test))
	println("------------------------------------------------------------")
	if (i<length(test_list))
		println("")
	end
end