#####
# labs
####

import Literate

Literate.notebook("src/labs/lab1s.jl", "labs/"; execute=false)
write("src/labs/lab1.jl", replace(replace(read("src/labs/lab1s.jl", String), r"## SOLUTION(.*?)## END"s => "")))
Literate.notebook("src/labs/lab1.jl", "labs/"; execute=false)
