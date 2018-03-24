path1=joinpath(Pkg.dir("MonteCarlo"),"test")
test_listTmp=readdir(path1);
test_list=[test_element for test_element in test_listTmp if test_element!="REQUIRE"&&test_element!="runtest.jl"]
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