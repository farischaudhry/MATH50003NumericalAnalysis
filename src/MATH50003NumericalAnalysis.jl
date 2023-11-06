#####
# notes
#####

using Weave


function replacetheorem(path, thrm, Thrm)
    str = read(path, String)
    str = replace(str, Regex("\\\\textbf{$Thrm \\(([^}]*)\\)}(.*?)\\\\textbf{Proof}", "s") => SubstitutionString("\\\\begin{$thrm}[\\1]\\2\\\\end{$thrm}\\n\\\\textbf{Proof}"))
    write(path, str)
end
function compilenotes(filename)
    weave("src/notes/$filename.jmd"; out_path="notes/", doctype="md2tex", template="src/notes/template.tpl")
    path = "notes/$filename.tex"
    replacetheorem(path, "theorem", "Theorem")
    replacetheorem(path, "lemma", "Lemma")
    # work around double newline before equation
    write(path, replace(read(path, String), "\n\n\\[" => "\n\\["))
end

compilenotes("I.1.RectangularRule")

#####
# labs
#####

import Literate

write("src/labs/lab1.jl", replace(replace(read("src/labs/lab1s.jl", String), r"## SOLUTION(.*?)## END"s => "")))
Literate.notebook("src/labs/lab1.jl", "labs/"; execute=false)
Literate.notebook("src/labs/lab1s.jl", "labs/"; execute=false)
