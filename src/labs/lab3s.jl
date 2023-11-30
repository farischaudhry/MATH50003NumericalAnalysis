# # MATH50003 (2023–24)
# # Lab 3: II.1 Integers and II.2 Reals

# We will use Julia in these notes to explore what is happening as a computer does integer arithmetic.
# We load an external package
# which implements functions `printbits` (and `printlnbits`)
# to print the bits (and with a newline) of numbers in colour:

using ColorBitstring, Test

# ## 1. Integers

# Every primitive number type is stored as a sequence of bits. 
# The number of _bytes_ (i.e. 8-bits) can be deduced using the `sizeof` function:

sizeof(UInt32) # 4 bytes == 4*8 bits == 32 bits

# The function `typeof` can be used to determine the type of a number.
# By default when we write an integer (e.g. `-123`) it is of type `Int`
# (which on 64-bit machines is equivalent to `Int64`):

typeof(5)

# -----
# **Problem 1.1** Use `sizeof` to determine how many bits your machine uses for the type `Int`.

## SOLUTION

sizeof(Int) # returns 8 bytes == 8*8 bits = 64 bits. Some machines may be 32 bits.

## END

# -----

# There are a few ways to create other types of integers. Conversion
# converts between different types:

UInt8(5) # converts an `Int` to an `UInt8`, displaying the result in hex

# This fails if a number cannot be represented as a specified type: e.g. `UInt8(-5)` and `UInt8(2^8)`.

# (These can also be written as e.g. `convert(UInt8, 5)`.)
# We can also create unsigned integers by specifying their bits
# by writing `0b` followed by a sequence of bits:

0b101 # isa UInt8, the smallest type with at least 3 bits
#
0b10111011101 # isa UInt16, the smallest type with at least 11 bits

# Or in base-16 using hexadecimal format (with digits `0–9a–f` following
# an `0x`), where each digit takes 4 bits to represent (since $2^4 = 16$):

0xabcde # isa UInt32, the smallest type with at least 4*5 = 20 bits

# -----
# **Problem 1.2** Use binary format to create an `Int` corresponding to $(101101)_2$.

## SOLUTION

Int(0b101101) # Without the `Int` it would be a UInt8

## END


# -----

# **Problem 1.3** What happens if you specify more than 64 bits using `0b⋅⋅…⋅⋅`? 
# What if you specify more than 128 bits?

## SOLUTION

typeof(0b111111111111111111111111111111111111111111111111111111111111111111111111111111111111) # creates a UInt128

typeof(0b111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111) # creates a BigInt

## END

# -----

# We can also reinterpret a sequence of bits in a different format:

reinterpret(Int8, 0b11111111) # Create an Int8 with the bits 11111111



# Integer arithmetic follows modular arithmetic. The following examples demonstrate this.

# **Example 2 (arithmetic with  8-bit unsigned integers)** 
# If  arithmetic lies between $0$ and $m = 2^8 = 256$ works as expected. 
# For example,
# $$
# \begin{align*}
# 17 ⊕_{256} 3 = 20 ({\rm mod}\ 256) = 20 \\
# 17 ⊖_{256} 3 = 14 ({\rm mod}\ 256) = 14
# \end{align*}
# $$

# This can be seen in Julia:

x = UInt8(17)  # An 8-bit representation of the number 255, i.e. with bits 00010001
y = UInt8(3)   # An 8-bit representation of the number   1, i.e. with bits 00000011
printbits(x); println(" + "); printbits(y); println(" = ")
printlnbits(x + y) # + is automatically modular arithmetic
printbits(x); println(" - "); printbits(y); println(" = ")
printbits(x - y) # - is automatically modular arithmetic



# **Example 3 (overflow with 8-bit unsigned integers)** If we go beyond the range
# the result "wraps around". For example, with integers we have
# $$
# 255 + 1 = (11111111)_2 + (00000001)_2 = (100000000)_2 = 256
# $$
# However, the result is impossible to store in just 8-bits! 
# So as mentioned instead it treats the integers as elements of ${\mathbb Z}_{256}$:
# $$
# 255 ⊕_{256} 1 = 255 + 1 \ ({\rm mod}\ 256) = (00000000)_2 \ ({\rm mod}\ 256) = 0 \ ({\rm mod}\ 256)
# $$
# On the other hand, if we go below $0$ we wrap around from above:
# $$
# 3 ⊖_{256} 5 = -2 ({\rm mod}\ 256) = 254 = (11111110)_2
# $$

# We can see this in  code:

x = UInt8(255) # An 8-bit representation of the number 255, i.e. with bits 11111111
y = UInt8(1)   # An 8-bit representation of the number   1, i.e. with bits 00000001
printbits(x); println(" + "); printbits(y); println(" = ")
printbits(x + y) # + is automatically modular arithmetic


x = UInt8(3) # An 8-bit representation of the number   3, i.e. with bits 00000011
y = UInt8(5) # An 8-bit representation of the number   5, i.e. with bits 00000101
printbits(x); println(" - "); printbits(y); println(" = ")
printbits(x - y) # + is automatically modular arithmetic


# **Example 4 (multiplication of 8-bit unsigned integers)** 
# Multiplication works similarly: for example,
# $$
# 254 ⊗_{256} 2 = 254 * 2 \ ({\rm mod}\ 256) = 252 \ ({\rm mod}\ 256) = (11111100)_2 \ ({\rm mod}\ 256)
# $$
# We can see this behaviour in code by printing the bits:

x = UInt8(254) # An 8-bit representation of the number 254, i.e. with bits 11111110
y = UInt8(2)   # An 8-bit representation of the number   2, i.e. with bits 00000010
printbits(x); println(" * "); printbits(y); println(" = ")
printbits(x * y)



# -----

# **Problem 1.5** Can you predict what the output of the following will be before hitting return?

UInt8(120) + UInt8(10); # Convert to `Int` to see the number printed in decimal
#
Int8(120) + Int8(10);
#
UInt8(2)^7;
#
Int8(2)^7;
#
Int8(2)^8;
#

## SOLUTION

UInt8(120) + UInt8(10) # returns 0x82 = 8*16+2 = 130

Int8(120) + Int8(10) # returns -126 since mod(-126,2^8) == 130

UInt8(2)^7 # Returns 0x80 = 8*16 = 128

Int8(2)^7 # Retuns -128 since mod(-128,2^8) == 128

Int8(2)^8 # Returns 0 since mod(2^8, 2^8) == 0

## END



# ## Division

# In addition to `+`, `-`, and `*` we have integer division `÷`, which rounds towards zero:

5 ÷ 2 # equivalent to div(5,2)

# Standard division `/` (or `\` for division on the right) creates a floating-point number,
# which will be discussed in the next chapter:

5 / 2 # alternatively 2 \ 5


#  We can also create rational numbers using `//`:

(1//2) + (3//4)

# Rational arithmetic often leads to overflow so it
# is often best to combine `big` with rationals:

big(102324)//132413023 + 23434545//4243061 + 23434545//42430534435



# ## 4. Variable bit representation

# An alternative representation for integers uses a variable number of bits,
#     with the advantage of avoiding overflow but with the disadvantage of a substantial
#     speed penalty. In Julia these are `BigInt`s, which we can create by calling `big` on an
#     integer:

x = typemax(Int64) + big(1) # Too big to be an `Int64`

#     Note in this case addition automatically promotes an `Int64` to a `BigInt`.
#     We can create very large numbers using `BigInt`:

x^100

#     Note the number of bits is not fixed, the larger the number, the more bits required 
#     to represent it, so while overflow is impossible, it is possible to run out of memory if a number is
#     astronomically large: go ahead and try `x^x` (at your own risk).
    





# -----

# ## 2. Reals
# 
# Real numbers interpret a sequence of bits in floating point format. 
# 
# -----
# **Problem 2.1** Use `printbits` to guess the binary representation of $1/5$.

## SOLUTION

printbits(1/5) 
## exponent is 0b01111111100 == 1020 so we have 2^(1020 - 1023) = 2^(-3)
## significand is 1.1001100110011001100110011001100110011001100110011010
## guess: 1/5 == 2^(-3) (1.10011001100…)_2 2^(-3) (∑_{k=0}^∞ (2^(-4k) + 2^(-4k-1)))

## END

# -----

# **Problem 2.2** Create a positive `Float64` whose exponent is $q = 156$ and has significand
# bits
# $$
# b_k = \begin{cases}
#     1 & k\hbox{ is prime} \\
#     0 & \hbox{otherwise}
#     \end{cases}
# $$

## SOLUTION

## significand has 52 bits. we can either do it by hand or create a string:

function isprime(k) # quick-and-dirty test for prime
    if k ≤ 1
        return false
    end
    for j=1:k-1
        if gcd(k, j) ≠ 1
            return false
        end
    end
    return true
end

ret = "1" # leading coefficient

for k = 1:52
    global ret # in scripts we need to let Julia know ret is a global variable
    if isprime(k)
        ret *= "1"
    else
        ret *= "0"
    end
end

sig = 2.0^(-52) * parse(Int, ret; base=2)

2.0^(156 - 1023) * sig

## END

# -----

# **Problem 2.3** Create the smallest positive non-zero sub-normal `Float16` by specifying
# its bits.

## SOLUTION
## sign is + so sign bit is 0, exponent is 00000 and significand is all zeros apart from a 1:
reinterpret(Float16, 0b0000000000000001) # == nextfloat(Float16(0))
## END

# -----

# ## 3. Strings and parsing

# Strings are a convenient way of representing arbitrary strings of digits.
# For example we can convert bits of a number to a string of "1"s and "0"s using the function `bitstring`.

# -----

# **Problem 3.1** Can you predict what the output of the following will be before hitting return?

bitstring(11);  # Semi-colon prohibits output, delete to check your answer
#
bitstring(-11);

## SOLUTION
bitstring(11) # "0000000000000000000000000000000000000000000000000000000000001011"
bitstring(-11) # "1111111111111111111111111111111111111111111111111111111111110101"
## this is because mod(-11, 2^64) == 2^64 - 12 == 0b10000…000 - 0b1100 == 0b111…11 - 0b1011 + 0b1
## END

# -----

# We can `parse` a string of digits in base 2 or 10:

parse(Int8, "11"; base=2), 
parse(Int8, "00001011"; base=2)

# Be careful with "negative" numbers, the following will fail: `parse(Int8, "10001011"; base=2)`

# It treats the string as binary digits, NOT bits. That is, negative numbers
# are represented using the minus sign:

parse(Int8, "-00001011"; base=2)

# -----

# **Problem 3.2** Combine `parse`, `reinterpret`, and `UInt8` to convert the
# above string to a (negative) `Int8` with the specified bits.

## SOLUTION

## The above code creates the bits "11110101". Instead, we first parse the bits:

x = reinterpret(Int8, parse(UInt8, "10001011"; base=2)) # -117
bitstring(x)

## END

# -----

# To concatenate strings we use `*` (multiplication is used because string concatenation
# is non-commutative):

"hi" * "bye"

# The string consisting of the first nine characters can be found using `str[1:9]` where `str` is any string:

str="hibye0123445556"
str[1:9]  # returns "hibye0123"

# The string consisting of the 11th through last character can be found using `str[11:end]`:

str="hibye0123445556"
str[11:end]  # returns "45556"

# -----

# **Problem 3.3** Complete the following function that sets the 10th bit of an `Int32` to `1`,
# and returns an `Int32`, assuming that the input is a positive integer, using `bitstring`,
# `parse` and `*`.

function tenthbitto1(x::Int32)
    ## TODO: change the 10th bit of `x` to 1
end

## unit tests are to help you check your result
## Change to `@test` to see if your test passes
@test_broken tenthbitto1(Int32(100)) ≡ Int32(4194404)

## SOLUTION
function tenthbitto1(x::Int32)
    ## TODO: change the 10th bit of `x` to 1
    ret = bitstring(x)
    parse(Int32, ret[1:9] * "1" * ret[11:end]; base=2)
end

## unit tests are to help you check your result
## Change to `@test` to see if your test passes
@test tenthbitto1(Int32(100)) ≡ Int32(4194404)

## END

# -----


# **Problem 3.4**  Modify the previous function to also work with negative numbers. 

function tenthbitto1(x::Int32)
## TODO: change the 10th bit of `x` to 1
end

@test_broken tenthbitto1(Int32(100)) ≡ Int32(4194404)
@test_broken tenthbitto1(-Int32(100000010)) ≡ Int32(-95805706)

## SOLUTION

function tenthbitto1(x::Int32)
    ## TODO: change the 10th bit of `x` to 1
    ret = bitstring(x)
    x = parse(UInt32, ret[1:9] * "1" * ret[11:end]; base=2)
    reinterpret(Int32, x)
end

@test tenthbitto1(Int32(100)) ≡ Int32(4194404)
@test tenthbitto1(-Int32(100000010)) ≡ Int32(-95805706)

## END

# ### Hexadecimal and binary format

# In Julia unsigned integers are displayed in hexadecimal
# form: that is, in base-16.
# Since there are only 10 standard digits (`0-9`) it uses 6 letters (`a–f`) to represent
# 11–16. For example,

UInt8(250)

# because `f` corresponds to 15 and `a` corresponds to 10, and we have
# $$
# 15 * 16 + 10 = 250.
# $$
# The reason for this is that each hex-digit encodes 4 bits (since 4 bits have $2^4 = 16$ possible
# values) and hence two hex-digits are encode 1 byte, and thus the digits correspond
# exactly with how memory is divided into addresses.
# We can create unsigned integers either by specifying their hex format:

0xfa

# Alternatively, we can specify their digits.
# For example, we know $(f)_{16} = 15 = (1111)_2$ and $(a)_{16} = 10 = (1010)_2$ and hence
# $250 = (fa)_{16} = (11111010)_2$ can be written as

0b11111010



# **Example (converting bits to signed integers)** 
# What 8-bit integer has the bits `01001001`? Because the first bit is 0 we know the result is positive.
# Adding the corresponding decimal places we get:

2^0 + 2^3 + 2^6

# What 8-bit (signed) integer has the bits `11001001`? Because the first bit is `1` we know it's a negative 
# number, hence we need to sum the bits but then subtract `2^p`:

2^0 + 2^3 + 2^6 + 2^7 - 2^8

# We can check the results using `printbits`:

printlnbits(Int8(73)) # Int8 is an 8-bit representation of the signed integer 73
printbits(-Int8(55))



# **Example 8 (overflow)** We can find the largest and smallest instances of a type using `typemax` and `typemin`:

printlnbits(typemax(Int8)) # 2^7-1 = 127
printbits(typemin(Int8)) # -2^7 = -128

# As explained, due to modular arithmetic, when we add `1` to the largest 8-bit integer we get the smallest:

typemax(Int8) + Int8(1) # returns typemin(Int8)

# This behaviour is often not desired and is known as _overflow_, and one must be wary
# of using integers close to their largest value.




# In Julia these correspond to 3 different floating-point types:

# 1.  `Float64` is a type representing double precision ($F_{64}$).
# We can create a `Float64` by including a 
# decimal point when writing the number: 
# `1.0` is a `Float64`. Alternatively, one can use scientific notation: `1e0`. 
# `Float64` is the default format for 
# scientific computing (on the _Floating-Point Unit_, FPU).  
# 2. `Float32` is a type representing single precision ($F_{32}$).  We can create a `Float32` by including a 
# `f0` when writing the number: 
# `1f0` is a `Float32` (this is in fact scientific notation so `1f1 ≡ 10f0`). 
# `Float32` is generally the default format for graphics (on the _Graphics Processing Unit_, GPU), 
# as the difference between 32 bits and 64 bits is indistinguishable to the eye in visualisation,
# and more data can be fit into a GPU's limited memory.
# 3.  `Float16` is a type representing half-precision ($F_{16}$).
# It is important in machine learning where one wants to maximise the amount of data
# and high accuracy is not necessarily helpful. 


# We confirm the simple bit representations:

σ,Q,S = 127,8,23 # Float32
εₘ = 2.0^(-S)
printlnbits(Float32(2.0^(1-σ))) # smallest positive normal Float32
printlnbits(Float32(2.0^(2^Q-2-σ) * (2-εₘ))) # largest normal Float32

# For a given floating-point type, we can find these constants using the following functions:

eps(Float32), floatmin(Float32), floatmax(Float32)

# **Example (creating a sub-normal number)** If we divide the smallest normal number by two, we get a subnormal number: 

mn = floatmin(Float32) # smallest normal Float32
printlnbits(mn)
printbits(mn/2)

# Can you explain the bits?




# ### Special numbers

# The special numbers extend the real line by adding $±∞$ but also a notion of "not-a-number" ${\rm NaN}$.
# Whenever the bits of $q$ of a floating-point number are all 1 then they represent an element of $F^{\rm special}$.
# If all $b_k=0$, then the number represents either $±∞$, called `Inf` and `-Inf` for 64-bit floating-point numbers (or `Inf16`, `Inf32`
# for 16-bit and 32-bit, respectively):

printlnbits(Inf16)
printbits(-Inf16)

# All other special floating-point numbers represent ${\rm NaN}$. One particular representation of ${\rm NaN}$ 
# is denoted by `NaN` for 64-bit floating-point numbers (or `NaN16`, `NaN32` for 16-bit and 32-bit, respectively):

printbits(NaN16)

# These are needed for undefined algebraic operations such as:

0/0

# Essentially it is a CPU's way of indicating an error has occurred.


# **Example (many `NaN`s)** What happens if we change some other $b_k$ to be nonzero?
# We can create bits as a string and see:

i = 0b0111110000010001 # an UInt16
reinterpret(Float16, i)

# Thus, there are more than one `NaN`s on a computer.  