using MATH50003NumericalAnalysis


####
# Notes
####

compilenotes("I.1.RectangularRule")
compilenotes("I.2.DividedDifferences")
compilenotes("I.3.DualNumbers")
compilenotes("I.4.NewtonMethod")

compilenotes("II.1.Integers")
compilenotes("II.2.Reals")
compilenotes("II.3.Arithmetic")
compilenotes("II.4.Intervals")

compilenotes("III.1.StructuredMatrices")
compilenotes("III.2.DifferentialEquations")
compilenotes("III.3.Cholesky")
compilenotes("III.4.Regression")
compilenotes("III.5.OrthogonalMatrices")
compilenotes("III.6.QR")


compilenotes("A.Asymptotics")
compilenotes("A.Permutations")


####
# Sheets
####

for k = 1:7
    compilesheet(k)
end

for k = 1:5
    compilesheetsolution(k)
end

####
# Labs
####
for k = 1:7
    compilelab(k)
end


compilelabsolution(1)
compilelabsolution(2)
compilelabsolution(3)
compilelabsolution(4)
compilelabsolution(5)



#####
# extras
#####

using Weave

nkwds = (out_path="notes/", jupyter_path="$(homedir())/.julia/conda/3/x86_64/bin/jupyter", nbconvert_options="--allow-errors")
notebook("src/notes/A.Julia.jmd"; nkwds...)