### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ 96bc8c74-f10a-11ea-2f86-9f75ccaa5f89
using Primes

# ╔═╡ 986d2a4a-f10c-11ea-0608-c3567591b9f4
using Statistics

# ╔═╡ bfa53152-f10c-11ea-3e60-c70874f29c08
using Gadfly  # the Julia version of the ggplot2 graphics package

# ╔═╡ e166f74e-f10c-11ea-01ec-8332273a87a1
using Colors  # for color keys and color maps in graphics

# ╔═╡ eb919846-f10c-11ea-1eaa-c59b08e5ea3e
using Cairo   # for graphics conversions, eg SVG --> PNG

# ╔═╡ fa380c72-f10c-11ea-3589-f3f3b35f6ced
using DataFrames # 2D arrays with named columns; hard to use Gadfly without them

# ╔═╡ b315f4d6-f10e-11ea-0fdc-9390980b46f1
md"""## Prime After Prime (2020 version for Pluto.jl)

This notebook explores correlations among the residues modulo $m$ of consecutive prime numbers. The notebook is an accompaniment to the article "[Prime After Prime](http://bit-player.org/2016/prime-after-prime)" published at bit-player.org on May 31, 2016. These correlations are *not* my discovery; they were first fully described and analyzed by Robert J. Lemke Oliver and Kannan Soundararajan in March 2016. For the mathematical background, see their arXiv preprint: [Unexpected Biases in the Distribution of Consecutive Primes](http://arxiv.org/abs/1603.03720). 

The code is written in the [Julia programming language](http://julialang.org/). The original version was written in Julia version 0.4, which is now very much out of date. In August of 2020, inspired by an inquiry from a high school student, I set about updating the code for use with Julia 1.0 or later. This was easier than I expected, partly because I discovered I was able to throw away much of the code. Functions I had written are available in Julia packages. Also, I have only included the essentials for exploring the mathematics and plotting results. I have eliminated various false starts and dead ends that serve mainly to distract.

The 2016 version of the notebook, primeafterprime.ipynb, remains available in this repository.

All my development work was done in a Jupyter notebook, but this version has been imported into the new notebook format of [Pluto.jl](https://github.com/fonsp/Pluto.jl).

"""

# ╔═╡ b4672322-f111-11ea-149d-8be5ae03d1d9
md"""## Basic operating instructions

- Install Pluto.js via the Julia Pkg manager
- cd to the directory where this file is stored
- Do ```import Pluto``` and ```Pluto.run()```
- In your default browser, navigate to ```localhost:1234```
- Open this file
- Evaluate all the cells. (This should happen automatically, I think, but I'm just learning Pluto.js. It could take a few minutes.)

To test that the system is working, go to the last three cells in the file, uncomment the lines there, and evaluate them. The first run may take a while, but you should eventually see a pink-and-blue heatmap like the ones in the bit-player article.

For your own experiments, the protocol is as follows:

Generate a list of consecutive primes. For example, get the first 10,000 primes larger than 500,000:

```Julia
p50000 = generate_primes(500000, 10000)
```

Now examine those primes modulo m and analyze correlations between successive residues mod m. For m = 7:

```Julia
p50000mod7 = analyze_primes_mod_m(p50000, 7)
```

The matrices p50000mod7[:nz\_class\_counts] and p50000mod7[:nz\_pair\_counts] are likely to be the main focus of interest.

To produce a graphic heatmap showing the variations in the pair correlation function, do:

```Julia
hm = heatmap(p50000mod7)
```

(Or do something else with the data that you find interesting!)
"""

# ╔═╡ c6901830-f10a-11ea-0e71-9b095034a67a
# Generate a list of 'n' primes starting with the first prime
# greater than or equal to 'start'. The result is returned as an
# array of integers of the same type as the type of 'start'.

# CAUTION: If 'start' is just below a boundary between integer types,
# this procedure may throw an OverflowError or, worse, return a
# list of primes that wraps around to those at the origin of the 
# number line. The boundaries in question are:
# typemax(Int64) = 9223372036854775807
# typemax(Int128) = 170141183460469231731687303715884105727
# If you are working in that neighborhood, it's best to explicitly
# coerce 'start' to the larger type, i.e. Int128(9223372036854775807) or
# BigInt(170141183460469231731687303715884105727)

function generate_primes(start, n)
    primelist = Array{typeof(start)}(undef, n)
    p = nextprime(start)
    for i = 1:n
        primelist[i] = p
        p = nextprime(p + 1)
    end
    return primelist
end

# ╔═╡ d8b4f6a2-f10a-11ea-0263-bb54f81fc3cd
# Auxiliary to 'analyze_primes_mod_m'

# Given a modulus m, return a vector of all possible
# congruence classes for p mod m, where p is a prime > m.
# When m is prime, the vector will consist of the range
# 1:m-1, but for composite m there may be missing values.
# For example, congruence_classes(10) = [1, 3, 7, 9]

function collect_congruence_classes(m)
    classes = Array{Int64}(undef, 0)
    for i = 1:m
        if gcd(i, m) == 1
            push!(classes, i)
        end
    end
    return classes
end

# ╔═╡ 7587d2c8-f10c-11ea-02fd-ad6dce0ead8a
# Auxiliary to 'analyze_primes_mod_m' 

function nonzero_class_counts(class_counts, classes)
    n = length(classes)
    nz_counts = zeros(Int64, n)
    for i = 1:n
        nz_counts[i] = class_counts[classes[i] + 1]
    end
    return nz_counts
end

# ╔═╡ 82e434c0-f10c-11ea-2d53-f796bcc82a68
# Auxiliary to 'analyze_primes_mod_m'

function nonzero_pair_counts(pair_counts, classes)
    n = length(classes)
    nz_counts = zeros(Int64, n, n)
    for i = 1:n
        for j = 1:n
            nz_counts[i, j] = pair_counts[classes[i] + 1, classes[j] + 1]
        end
    end
    return nz_counts
end

# ╔═╡ 8da414dc-f10c-11ea-2511-119c66bf5293
# Auxiliary to 'analyze_primes_mod_m'

# calculate the mean value of all off-diagonal elements and the
# mean of the on-diagonal elements; return the ratio of off/on

function diagonal_ratio(count_matrix)
    n = size(count_matrix, 1)
    diagonal_sum = 0
    for i = 1:n
        diagonal_sum += count_matrix[i, i]
    end
    off_diagonal_sum = sum(count_matrix) - diagonal_sum
    return (off_diagonal_sum / (n^2 - n)) / (diagonal_sum / n)
end

# ╔═╡ b0dc1758-f10c-11ea-000f-a7b0d08cd063
# given a list of primes and an integer modulus m,
# returns a statistical analysis of the primes mod m, organized
# as a data dictionary

function analyze_primes_mod_m(primelist, m)
    
    # associative structure to hold all data on primes
    d = Dict()
    
    # basic info about the set of primes and the modulus
    d[:first_prime] = primelist[1]
    d[:last_prime] = primelist[end]
    d[:median_prime] = primelist[div(end, 2)]
    d[:prime_count] = length(primelist)    
    d[:m] = m
    
    # set of integers that are possible residues for a prime p mod m
    d[:congruence_classes] = collect_congruence_classes(m)
    
    # log10 of the median will become x coordinate in certain graphs
    d[:magnitude] = Float64(log10(d[:median_prime]))
    d[:digits] = round(Int, (ceil(d[:magnitude])))
    
    # arrays of zeros to hold counts of residue classes and consecutive pairs
    class_counts = zeros(Int, m)
    pair_counts = zeros(Int, (m, m))
    
    # loop through the primes mod m, recording the number in
    # each congruence class and the number of each p_{i}, p_{i+1} 
    # pair of successive primes; adjust indices to translate
    # from 0:m-1 range of modular values to 1:m range of Julia arrays.
    prev = primelist[1] % m
    class_counts[prev + 1] += 1
    for i = 2:length(primelist)
        next = primelist[i] % m
        class_counts[next + 1] += 1
        pair_counts[prev + 1, next + 1] += 1
        prev = next
    end
    
    # store the counts in the dictionary
    d[:class_counts] = class_counts
    d[:pair_counts] = pair_counts
    
    # results above include zero entries for impossible congruence
    # classes; here we separately record a vector and a matrix without
    # the zero entries
    d[:nz_class_counts] = nonzero_class_counts(class_counts, d[:congruence_classes])
    d[:nz_pair_counts] = nonzero_pair_counts(pair_counts, d[:congruence_classes])
    
    # normalize the results to allow comparisons of samples of different sizes
    mc = mean(d[:nz_class_counts])
    mp = mean(d[:nz_pair_counts])
    d[:norm_class_counts] = map(x -> (x - mc) / mc, d[:class_counts])
    d[:norm_pair_counts] = map(x -> (x - mp) / mp, d[:pair_counts])
    d[:norm_nz_class_counts] = map(x -> (x - mc) / mc, d[:nz_class_counts])
    d[:norm_nz_pair_counts] = map(x -> (x - mp) / mp, d[:nz_pair_counts])
    
    # elementary stats on the normalized results
    d[:norm_class_counts_mean] = mean(d[:norm_nz_class_counts])
    d[:norm_class_counts_std] = std(d[:norm_nz_class_counts])
    d[:norm_pair_counts_mean] = mean(d[:norm_nz_pair_counts])
    d[:norm_pair_counts_std] = std(d[:norm_nz_pair_counts])
    d[:diag_ratio] = diagonal_ratio(d[:nz_pair_counts])
    return d
end

# ╔═╡ 0218ae7e-f10d-11ea-3916-c3f4ee27d9b6
Gadfly.set_default_plot_size(400pt, 400pt)

# ╔═╡ 11310d90-f10d-11ea-35a1-4b37f1f49b81
function colorfn(x)
    if x == 0.0
        return RGB(0.5, 0.5, 0.5)
    else
        index = round(Int, min(200, 1 + 200x))
        return Colors.diverging_palette(245.0, 10.0, 200, wcolor=RGB(1,0,0), c=0.7)[index]
    end
end

# ╔═╡ 1f3ee64e-f10d-11ea-3aec-bd0e36046247
function clamp_colors(mat)
    max_value = maximum(mat)
    if max_value > 1
        adjusted_mat = similar(mat)
        for i in eachindex(mat)
            if mat[i] > 0
                adjusted_mat[i] = mat[i] / max_value
            else
                adjusted_mat[i] = mat[i]
            end
        end
        return adjusted_mat
    end
    return mat
end

# ╔═╡ 25a9d4e4-f10d-11ea-1c22-ef9ac3b8fcd4
function heatmap(d::Dict)
    mat = clamp_colors(d[:norm_nz_pair_counts])
    labels = d[:congruence_classes]
    spy(mat,
        Scale.x_discrete(labels = x -> string(labels[x])),
        Scale.y_discrete(labels = y -> string(labels[y])),
        Guide.xlabel("j = second prime mod " * string(d[:m])),
        Guide.ylabel("i = first prime mod " * string(d[:m]), orientation=:vertical),
        Guide.title("normalized counts of consecutive " * string(d[:digits]) * "-digit primes mod " * string(d[:m])),
        Guide.colorkey(""),
        Scale.color_continuous(colormap = colorfn, minvalue = -1.0, maxvalue = 1.0))
end

# ╔═╡ 2fba1106-f10d-11ea-3214-83caf8fe272c
function heatmap_z(d::Dict)
    mat = clamp_colors(d[:norm_pair_counts])
    labels = collect(0:d[:m] - 1)
    spy(mat,
        Scale.x_discrete(labels = x -> string(labels[x])),
        Scale.y_discrete(labels = y -> string(labels[y])),
        Guide.xlabel("j = second prime mod " * string(d[:m])),
        Guide.ylabel("i = first prime mod " * string(d[:m]), orientation=:vertical),
        Guide.title("normalized counts of consecutive " * string(d[:digits]) * "-digit primes mod " * string(d[:m])),
        Guide.colorkey(""),
        Scale.color_continuous(colormap = colorfn, minvalue = -1.0, maxvalue = 1.0))
end

# ╔═╡ 38afa348-f10d-11ea-1619-59dc225fa5a6
# testmod7 = generate_primes(1000000, 100000)

# ╔═╡ 49183112-f10d-11ea-19eb-c3fb0a4852d8
# test7data = analyze_primes_mod_m(testmod7, 7)

# ╔═╡ 5035f328-f10d-11ea-0029-39abf9f73636
# hm7 = heatmap(test7data)

# ╔═╡ Cell order:
# ╟─b315f4d6-f10e-11ea-0fdc-9390980b46f1
# ╟─b4672322-f111-11ea-149d-8be5ae03d1d9
# ╠═96bc8c74-f10a-11ea-2f86-9f75ccaa5f89
# ╠═c6901830-f10a-11ea-0e71-9b095034a67a
# ╠═d8b4f6a2-f10a-11ea-0263-bb54f81fc3cd
# ╠═7587d2c8-f10c-11ea-02fd-ad6dce0ead8a
# ╠═82e434c0-f10c-11ea-2d53-f796bcc82a68
# ╠═8da414dc-f10c-11ea-2511-119c66bf5293
# ╠═986d2a4a-f10c-11ea-0608-c3567591b9f4
# ╠═b0dc1758-f10c-11ea-000f-a7b0d08cd063
# ╠═bfa53152-f10c-11ea-3e60-c70874f29c08
# ╠═e166f74e-f10c-11ea-01ec-8332273a87a1
# ╠═eb919846-f10c-11ea-1eaa-c59b08e5ea3e
# ╠═fa380c72-f10c-11ea-3589-f3f3b35f6ced
# ╠═0218ae7e-f10d-11ea-3916-c3f4ee27d9b6
# ╠═11310d90-f10d-11ea-35a1-4b37f1f49b81
# ╠═1f3ee64e-f10d-11ea-3aec-bd0e36046247
# ╠═25a9d4e4-f10d-11ea-1c22-ef9ac3b8fcd4
# ╠═2fba1106-f10d-11ea-3214-83caf8fe272c
# ╠═38afa348-f10d-11ea-1619-59dc225fa5a6
# ╠═49183112-f10d-11ea-19eb-c3fb0a4852d8
# ╠═5035f328-f10d-11ea-0029-39abf9f73636
