# # MATH50003 (2023–24)
# # Lab 1: Introduction to Mathematical Computing

# Numerical analysis primarily studies the mathematical construction and analysis of algorithms
# for solving continuous problems like computing integrals or solving differential equations.
# It is fundamental to the area to also understand how to implement such algorithms
# in software. In year 1 you learned basic programming concepts such as loops, conditions,
# functions, etc. and in this first lab we will employ these concepts in the implementation
# of some basic algorithms you have already seen. In particular, we will look at implementing
# the rectangular and triangular rules for approximating integrals. 
#
# We will use the Julia programming language which is in some ways similar to Python.
# Julia is a _compiled language_ whereas Python is interpreted. It is also more adapted to the
# implementation of algorithms in numerical analysis and scientific computing. 
# Being a compiled language means it will help us later on in the module understand exactly how
# the computer functions when performing numerical calculations.
#
# We have included exercises interspersed with the material which are highly recommended for
# preparation for the computer-based exam later this term. Note each exercise comes with a
# "unit-test". 

## 1. Rectangular rules

# One possible definition for an integral is the limit of a Riemann sum, for example:
# $$
#   ∫_0^1 f(x) {\rm d}x = \lim_{n → ∞} {1 \over n} ∑_{k=0}^{n-1} f(k/n).
# $$
# This suggests an algorithm known as the _left-sided rectangular rule_
# for approximating an integral: choose $n$ large and then
# $$
#   ∫_0^1 f(x) {\rm d}x ≈ {1 \over n} ∑_{k=0}^{n-1} f(k/n).
# $$
# To implement this approximation in code we need to turn the sum into a for-loop.
# Let's take as an example $f(x) = \exp(x)$. We can write:

n = 10000     # the number of terms in the summation
ret = 0.0     # ret will store the result, accumulated one argument at a time.
              # The .0 makes it a "real" rather than an "integer".
              # Understanding the "type" will be important later on.
for k = 0:n-1 # k will be equal to 0,1,…,n-1
    ret = ret + exp(k/n) # add exp(k/n) to the result. Now ret = ∑_{j=0}^k f(j/n).
end           # in Julia for-loops are finished with an end
ret/n         # approximates the true answer exp(1) - exp(0) = ℯ-1 = 1.71828… to 4 digits

# It is convenient to wrap this in a function that takes in `f` and `n` and returns
# the left-sided rectangular rule approximation. We can the above routine into a function as follows:

function leftrectangularrule(f, n) # create a function named "leftrectangularrule" that takes in two arguments
    ret = 0.0
    for k = 0:n-1
        ret = ret + f(k/n) # now `f` is the function we put in
    end           
    ret/n   # the last line of a function is returned
end # like for-loops, functions are finished with an end

leftrectangularrule(exp, 100_000_000) # Use n = 100 million points to get an approximation accurate to 8 digits.
                                      # The underscores in numbers are like commas and are ignored.

# Note it is now easy to approximate other functions. For example, the following code computes the
# integral of $x^2$:

function squared(x)
    x^2 # carets ^ mean "to the power of". This is actually a function that just calls x*x.
end
leftrectangularrule(squared, 10_000) # approximates 1/3 to 3 digits

# It is often inconvenient to name a function, and so we might want to integrate a function like $\cos(x^2)$
# by making a so-called anonymous function:

leftrectangularrule(x -> cos(x^2), 10_000) # No nice formula! But I claim we got 4 digits

# **Exercise 1(a)** Complete the following function `rightrectangularrule(f, n)` That approximates
# an integral using the right-sided rectangular rule:
# $$
#   ∫_0^1 f(x) {\rm d}x ≈ {1 \over n} ∑_{k=1}^n f(k/n).
# $$

using Test # Loads the testing packages

function rightrectangularrule(f, n)
    # TODO: return (1/n) * ∑_{k=1}^n f(k/n) computed using a for-loop
end

# Change `@test_broken` to `@test` to check if your solution is correct
@test_broken rightrectangularrule(exp, 1000) ≈ exp(1) - 1 atol=1E-3 # tests that the approximation is accurate to 3 digits after the decimal point
@test_broken leftrectangularrule(exp, 1000) < exp(1) - 1 < rightrectangularrule(exp, 1000) # These two routines bound the true answer. Why is this?

# **Exercise 1(b)** Write a function `trapeziumrule(f, n)` That approximates
# an integral using the trapezium rule:
# $$
#   ∫_0^1 f(x) {\rm d}x ≈ {1 \over n} \left[ f(0)/2 + ∑_{k=1}^{n-1} f(k/n) + f(1)/2 \right]
# $$
# Is it more or less accurate than the rectangular rules?

## 2. Studying errors in approximations



