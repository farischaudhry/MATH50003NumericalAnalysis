#####
# labs
####

import Literate

Literate.notebook("src/labs/lab1s.jl", "labs/")
write("src/labs/lab1.jl", replace(replace(read("src/labs/lab1s.jl", String), r"## SOLUTION(.*?)## END"s => ""), r"@test" => "@test_broken"))
