path1=joinpath(Pkg.dir("FinancialMonteCarlo"),"benchmark")
test_listTmp=readdir(path1);
BlackList=["REQUIRE","benchmarks.jl","runner.jl","cuda","af"]
test_list=[test_element for test_element in test_listTmp if !Bool(sum(test_element.==BlackList))]
println("Running tests:\n")
i=1;
for current_test in test_list
	println("------------------------------------------------------------")
    println("  * $(current_test) *")
    include(current_test)
	println("------------------------------------------------------------")
	if (i<length(test_list))
		println("")
	end
	i+=1;
end