path1=joinpath(Pkg.dir("FinancialMonteCarlo"),"benchmark")
test_listTmp=readdir(path1);
BlackList=["benchmarks.jl","runner.jl","cuda","af"]
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