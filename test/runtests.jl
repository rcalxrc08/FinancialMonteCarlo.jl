using FinancialMonteCarlo
path1=joinpath(dirname(pathof(FinancialMonteCarlo)),"..","test")
test_listTmp=readdir(path1);

BlackList=["REQUIRE",
"runtests.jl",
"Project.toml",
"Manifest.toml",
"runner.jl",
"test_black_cuda.jl",
"cuda",
"af",
"wip",
"test_black_rev_diff.jl",
"test_diff_eq_montecarlo_array.jl",
"test_black_pre.jl"];

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