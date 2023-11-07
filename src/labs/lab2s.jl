# # MATH50003 (2023–24)
# # Lab 2: I.3 Dual Numbers and I.4 Newton's Method

# In this lab we explore an alternative approach to computing derivatives:
# using _dual numbers_. This is a special mathematical object akin to complex numbers
# that allows us to compute derivatives to very high accuracy in an automated fashion,
# i.e. an example of [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)
# that is extremely important in Machine Learning and other computational applications.
# To realise dual numbers on a computer we need to introduce the notation of a "type"
# and create a customised type to represent dual numbers, which is what we discuss first.
# As an application of computing derivatives we consider root finding via [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).


# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 1. Definition of dual numbers and functions applied dual numbers.
# 2. Approximating second derivatives using second-order divided differences. 
# 3. Newton's method for root finding.
#
# Coding knowledge:
#
# 1. The notion of a type and how to make your own type.
# 2. Defining functions for specific types.
# 3. Overloading functions like `+`, `*`, and `exp` for a custom type.


# We load the `Test` and `Plots` packages to be used below.
# We also define the function `nanabs` which is useful for logarithmically scaled
# plots. For brevity we use a shorthand `cond ? expr1 : expr2` which just means
# ```julia
# if cond
#     expr1
# else
#     expr2
# end
# ```

using Plots, Test
nanabs(x) = x == 0 ? NaN : abs(x)


# ## Types in Julia


# In Julia everything has a "type". The function `typeof` can be used to determine the type of,
# for example, a number.
# By default when we write an Int (e.g. `-123`) it is of type `Int`:

typeof(5)

# On a 64-bit machine this will print `Int64`, where the `64` indicates it is using precisely 64 bits
# to represent the number (a topic we will come back to next lab). If we write something with
# a decimal point it represents a "real" number, but the actual type is `Float64`:

typeof(5.3)

# This is called a "floating point" number, and again the `64` indicates it is using precisely
# 64-bits to represent this number. (We will see this is why computations like divided differences 
# have error: because we are limiting the number of "digits" to represent numbers we need to
# round our computations.) Note that some operations involving `Int`s return `Float64`s:

1/5 # 1 and 5 are Int but output is a Float64

# It is possible to have functions behave differently depending on the input value. 
# To do so we can add a restriction denoted `::Int` or `::Float64` to the function "signature".





# Types allow for combining multiple numbers (or other types) to represent a more complicated
# object. That is, while a computer can only apply functions on $p$-bits at a time,
# these functions can be combined to perform more complicated operations on types
# that require more than $p$-bits. A simple example of this is a complex number, 
# which stores two real numbers $x$ and $y$ (either `Int` or `Float64` or indeed other real number types not yet discussed)
# to represent the complex number $x + {\rm i} y$. In Julia ${\rm i}$ is denoted `im` and
# hence we can create a complex number like $1+2{\rm i}$ as follows:

z = 1 + 2im

# This complex number has two "fields": the real and imaginary part. Accessing the fields is done
# using a `.`, here we display the real and imaginary parts as a "tuple":

z.re, z.im

# When we ask  its type we see it is a `Complex{Int}`:

typeof(z)

# That is, it is of a type `Complex` and the `{Int}` indicates that each of the fields is an `Int`.
# Note we can add, subtract, multiply, or apply functions like `exp` to complex numbers:

exp(2z^2 + 3im)

# -----
# **Problem 1(a)** Use `typeof` to determine the type of `1.2 + 2.3im`.

## TODO: What is the type of `1.2 + 2.3im`?
## SOLUTION
typeof(1.2 + 2.3im)
## `ComplexF64` is short hand for `Complex{Float64}`
## END

# **Problem 1(b)** Write a function `runge(z)` that computes $1/(25z^2+1)$.

## TODO: Write a function `runge(z)`
## SOLUTION
runge(z) = 1/(25z^2+1)
## END
@test runge(0.1) ≈ 0.8
# ------ 

# **Problem 2(a)** Consider the Taylor series approximation to the exponential:
# $$ 
# \exp z ≈ ∑_{k=0}^n {z^k \over k!}
# $$
# Complete the function `exp_t(z, n)` that computes this and returns a
# `Complex{Float64}` if the input is complex and a `Float64`. 
# Do not use the inbuilt `factorial` function.
# Hint: Floating point numbers cope much better with the factorial than Ints
# do. It might help to think inductively: for $s_k = z^k/k!$ we have 
# $$
#   s_{k+1}  = {z \over k+1} s_k.
# $$

function exp_t(z, n)
    ## TODO: Compute the first (n+1)-terms of the Taylor series of exp evaluated at z
    ret = one(z) # the function one a "1" of the same type as z
    s = one(z) # Use s to represent s_k at each step, starting with s_0
    ## SOLUTION
    for k = 1:n
        s = s/k * z
        ret = ret + s
    end
    ret
    ## END
end

@test exp_t(1, 10) isa Float64 # isa is used to test the type of a result
@test exp_t(im, 10) isa ComplexF64 # isa is used to test the type of a result

@test exp_t(1, 100) ≈ exp(1)

# **Problem 2(b)** Plot the error for `n = 1:1000` of `exp_t(z, n)` for `z = 1, im, -5`, and `-100`,
# scaling the y-axis logarithmically.
# Does the method appear to converge for all values of $z$? 

## TODO: plot the error for the Taylor series approximation.

## SOLUTION
plot(1:1000, [nanabs(exp_t(1, n) - exp(1)) for n = 1:1000]; yscale=:log10, label="z = 1")
plot!(1:1000, [nanabs(exp_t(im, n) - exp(im)) for n = 1:1000]; yscale=:log10, label="z = im")
plot!(1:1000, [nanabs(exp_t(-10, n) - exp(-10)) for n = 1:1000]; yscale=:log10, label="z = -10")
plot!(1:1000, [nanabs(exp_t(-100, n) - exp(-100)) for n = 1:1000]; yscale=:log10, label="z = -100")

## It appears to converge to a fixed constant. But this constant is growing exponentially with $z$ giving
## very inaccurate results for `z = -100`.
## END

# ------


# One of the powerful parts of Julia is its very easy to make our own types. Lets begin with a simple
# implementation of a rational function $p/q$ where $p$ and $q$ are Ints.  Thus we want to create a new
# type called `Rat` with two fields `p` and `q` to represent the numerator and denominator, respectively.
# (For simplicity  we won't worry about restricting $p$ and $q$ to be `Int`.)
# We can construct such a type using the `struct` keyword:

struct Rat
    p
    q
end

# A new instance of `Rat` is created via e.g. `Rat(1, 2)` represents 1/2
# where the first argument specifies `p` and the second argument `q`.
# The fields are accessed by `.`:

x = Rat(1, 2) # Rat(1, 2) creates an instance with fields equal to the input
@test x.p == 1
@test x.q == 2

# Unfortunately we can't actually do anything with this type, yet:

x + x

# The error is telling us to overload the `+` function when the inputs are both `Rat`.
# Do to do this we need to "import" the `+` function and then we can overload it like any
# other function:

import Base: + # allows us to overload +

+(x::Rat, y::Rat) = Rat(x.p * y.q + y.p * x.q, x.q * y.q)

Rat(1,2) + Rat(3,4) # 1/2 + 3/4 == 10/8 (== 5/4)

# We can support mixing `Rat` and `Int` by adding additional functionality:

Rat(p::Int) = Rat(p,1) # an Int is converted to p/1
+(x::Rat, y::Int) = x + Rat(y) # To add a Rat to an Int we convert the Int into a Rat and use the previously defined +

Rat(1,2) + 1  # 1 + 1/2 == 3/2

# -----

# **Problem 3** Support `*`, `-`, `/`, and `==` for `Rat`.

## We import `-`, `*`, `/` so we can "overload" these operations specifically for `Rat`.
import Base: +, -, *, /, ==

## The ::Rat means the following version of `==` is only called if both arguments
## are Rat
function ==(x::Rat, y::Rat)
    ## TODO: implement equality, making sure to check the case where
    ## the numerator/denominator are possibly reducible
    ## Hint: `gcd` and `div` may be useful. Use `?` to find out what they do

    ## SOLUTION
    xg = gcd(x.p, x.q)
    yg = gcd(y.p, y.q)
    div(x.p, xg) == div(y.p, yg) && div(x.q, xg) == div(y.q, yg)
    ## END
end

## We can also support equality when `x isa Rat` and `y isa Int`
function ==(x::Rat, y::Int)
    ## TODO: implement
    ## SOLUTION
    x == Rat(y, 1)
    ## END
end

## TODO: implement ==(x::Int, y::Rat)
## SOLUTION
function ==(x::Int, y::Rat)
    ## TODO: implement
    ## SOLUTION
    Rat(x,1) == y
    ## END
end

## END

@test Rat(1, 2) == Rat(2, 4)
@test Rat(1, 2) ≠ Rat(1, 3)
@test Rat(2,2) == 1
@test 1 == Rat(2,2)

## TODO: implement +, -, *, and /, 
## SOLUTION

+(x::Rat, y::Rat) = Rat(x.p * y.q + y.p * x.q, x.q * y.q)
-(x::Rat, y::Rat) = Rat(x.p * y.q - y.p * x.q, x.q * y.q)
*(x::Rat, y::Rat) = Rat(x.p * y.p, x.q * y.q)
/(x::Rat, y::Rat) = x * Rat(y.q, y.p)

## END

@test Rat(1, 2) + Rat(1, 3) == Rat(5, 6)
@test Rat(1, 3) - Rat(1, 2) == Rat(-1, 6)
@test Rat(2, 3) * Rat(3, 4) == Rat(1, 2)
@test Rat(2, 3) / Rat(3, 4) == Rat(8, 9)


# ## I.3 Dual Numbers
# 
# We now consider implementing a type `Dual` to represent the dual number `a + b*ϵ`,
# in a way similar to `Complex` or `Rat`. For simplicity we don't restrict the types of `a` and `b`
# but for us they will usually be `Float64`. We create this type very similar to `Rat` above:

struct Dual
    a
    b
end

# We can easily support addition of dual numbers as in `Rat` using the formula
# $$
# (a+bε) + (c+dε) = (a+c) + (b+d)ε
# $$

function +(x::Dual, y::Dual)
    a,b = x.a, x.b # x == a+bε. This gets out a and b
    c,d = y.a, y.b # y == c+dε. This gets out c and d
    Dual(a+c, b+d)
end

Dual(1,2) + Dual(3,4) # just adds each argument

# For multiplication weuse the formula
# $$
# (a+bε)*(c+dε) = ac +(bc+ad)ε
# $$

import Base: * # we want to also overload *

function *(x::Dual, y::Dual) 
    a,b = x.a, x.b # x == a+bε. This gets out a and b
    c,d = y.a, y.b # y == c+dε. This gets out c and d
    Dual(a*c, b*c + a*d)
end

# We can already differentiate simple polynomials:

f = x -> x*x*x + x
f(Dual(2,1)) # (2^3 + 2) + (3*2^2+1)*ε

# A polynomial like `x^3 + 1` is not yet supported.
# To support this we need to add addition of `Dual` with `Int` or `Float64`.
# Note that both of these are "subtypes" of `Real` and so restricting on `Real`
# will support both at the same time.
# We can overload the appropriate functions as follows:

import Base: ^

Dual(a::Real) = Dual(a, 0) # converts a real number to a dual number with no ε

+(x::Real, y::Dual) = Dual(x) + y
+(x::Dual, y::Real) = x + Dual(y)

# a simple recursive function to support x^2, x^3, etc.
function ^(x::Dual, n::Int)
    if n < 0
        error("Not implemented") # don't support negative n, yet
    end
    if n == 1
        x # Just return the input
    else
        ret = x
        for k = 1:n-1
            ret = ret*x
        end
        ret # returns the last argument
    end
end

f = x -> x^3 + 1
f(Dual(2,1))  # 2^3+1 + 3*2^2*ε



# Algebraic operationds for duals
-(x::Dual) = Dual(-x.a, -x.b)
-(x::Dual, y::Dual) = Dual(x.a - y.a, x.b - y.b)


exp(x::Dual) = Dual(exp(x.a), exp(x.a) * x.b)


# We can also try it on the two polynomials as above:

f = x -> 1 + x + x^2
g = x -> 1 + x/3 + x^2
f(ϵ).b, g(ϵ).b

# The first example exactly computes the derivative, and the
# second example is exact up to the last bit rounding!
# It also works for higher order polynomials:

f = x -> 1 + 1.3x + 2.1x^2 + 3.1x^3
f(Dual(0.5,1))

# It is indeed "accurate to (roughly) 16-digits", the best we can hope for 
# using floating point.

# We can use this in "algorithms" as well as simple polynomials.
# Consider the polynomial $1 + … + x^n$:

function s(n, x)
    ret = 1 + x # first two terms
    for k = 2:n
        ret += x^k
    end
    ret
end
s(10, 0.1 + ϵ).b

# This matches exactly the "true" (up to rounding) derivative:

sum((1:10) .* 0.1 .^(0:9))


# Finally, we can try the more complicated example:

f = x -> exp(x^2 + exp(x))
f(1 + ϵ)


# What makes dual numbers so effective is that, unlike divided differences, they are not
# prone to disasterous growth due to round-off errors. 


# ------

# **Problem 1(a)** Add support for `cos`, `sin`, and `/` to the type `Dual`
# by replacing the `# TODO`s in the below code.

function cos(x::Dual)
    ## TODO: implement cos for Duals
    ## SOLUTION
    Dual(cos(x.a), -sin(x.a) * x.b)
    ## END
end

function sin(x::Dual)
    ## TODO: implement sin for Duals
    ## SOLUTION
    Dual(sin(x.a), cos(x.a) * x.b)
    ## END
end

function /(x::Dual, y::Dual)
    ## TODO: implement division for Duals
    ## SOLUTION
    if iszero(y.a)
        error("Division for dual numbers is ill-defined when denonimator real part is zero.")
    end
    return Dual(x.a / y.a, (y.a * x.b - x.a * y.b) / y.a^2)
    ## END
end

x = 0.1
@test cos(sin(x+ϵ)/(x+ϵ)).b ≈ -((cos(x)/x - sin(x)/x^2)sin(sin(x)/x))


# **Problem 1(b)** Use dual numbers to compute the derivatives to
# $$
# \exp(\exp x \cos x + \sin x), ∏_{k=1}^{1000} \left({x \over k}-1\right), \hbox{ and } f^{\rm s}_{1000}(x).
# $$
# Compare with divided differences to give evidence that your implementation is correct.

## TODO: Use dual numbers to compute the derivatives of the 3 functions above.
## SOLUTION

## Define the functions
f = x -> exp(exp(x)cos(x) + sin(x))
g = x -> prod([x] ./ (1:1000) .- 1)
function cont(n, x)
    ret = 2*one(x)
    for k = 1:n-1
        ret = 2 + (x-1)/ret
    end
    1 + (x-1)/ret
end

## With the previous problems solved, this is as simple as running

fdual = f(0.1+ϵ)
fdual.b
#
gdual = g(0.1+ϵ)
gdual.b
#
contdual = cont(1000,0.1+ϵ)
contdual.b
## END



# ------
# ## I.4 Newton's method