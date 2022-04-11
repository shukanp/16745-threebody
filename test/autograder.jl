using NBInclude
using Test
using Formatting
using Crayons
using Crayons.Box

const INDENT_SIZE = 2

import Test: Test, record, finish
using Test: AbstractTestSet, Result, Pass, Fail, Error
using Test: get_testset_depth, get_testset
mutable struct CustomTestSet <: Test.AbstractTestSet
    description::AbstractString
    foo::Int
    results::Vector
    points_total::Float64
    points_missed::Float64
    depth::Int
    # constructor takes a description string and options keyword arguments
    CustomTestSet(desc; foo=1) = new(desc, foo, [], 0, 0, 0)
end

record(ts::CustomTestSet, child::AbstractTestSet) = push!(ts.results, child)
record(ts::CustomTestSet, res::Result) = push!(ts.results, res)
function record(ts::CustomTestSet, res::Pass) 
    # points,total = getpoints(res.source); 
    # ts.points_passed += points; 
    # ts.points_total += points; 
    push!(ts.results, res)
end
function record(ts::CustomTestSet, res::Union{Fail,Error})
    points = getpoints(res.source); 
    ts.points_missed += points
end
function finish(ts::CustomTestSet)
    # just record if we're not the top-level parent
    depth = get_testset_depth()
    ts.depth = depth
    if depth > 0
        record(get_testset(), ts)
        return ts
    end
    calcscore!(ts)

    ts
end

function calcscore!(ts::CustomTestSet)
    for r in ts.results
        if r isa CustomTestSet
            calcscore!(r)
            # ts.points_total += r.points_total
            ts.points_missed += r.points_missed
        end
    end
end

function getdescriptionwidth(ts::CustomTestSet)
    width = ts.depth * INDENT_SIZE + length(ts.description)
    # if allpassed(ts)
    #     return width
    # end
    for r in ts.results
        if r isa CustomTestSet
            child_width = getdescriptionwidth(r)
            width = max(width, child_width)
        end
    end
    return width
end

function getscorecolor(pass, total)
    score =  pass / total * 100
    if score == 100 
        return GREEN_FG
    elseif score >= 90
        return crayon"118" 
    elseif score >= 80
        return YELLOW_FG 
    elseif score >= 60
        return crayon"208"  # orange 
    else
        return RED_FG
    end
end

allpassed(ts::CustomTestSet) = ts.points_total == ts.points_passed
function printresult(ts::CustomTestSet; width=getdescriptionwidth(ts))
    indent = " "^(ts.depth * INDENT_SIZE)
    fspec1 = "{:<$(width)}"
    fspec2 = "{:<5d}"
    fspec3 = "{:<5}"
    if ts.depth == 0
        println(BOLD, format(fspec1,"\nTest Summary:"), " | ", CYAN_FG(format(fspec3,"Score")), " ", BLUE_FG(format(fspec3, "Total")))
    end
    fail = round(Int,ts.points_missed)
    total = round(Int,ts.points_total)
    pass = total - fail 
    color = getscorecolor(pass, total) 
    println(crayon"!bold", format(fspec1, indent * ts.description), BOLD, " | ", color(format(fspec2,pass)), " ", BLUE_FG(format(fspec2, total)))
    # println(format(fspec, indent * ts.description, ts.points_passed, ts.points_total))
    for r in ts.results
        printresult(r, width=width)
    end
end
printresult(res; width=0) = nothing
printresult(res::Fail; width=0) = nothing 

function getsourceline(src::LineNumberNode)
    filestr = string(src.file)
    if isfile(filestr)
        line = ""
        open(filestr) do file
            for i = 1:src.line-1
                if eof(file)
                    @warn "Reached end of file."
                    return ""
                else
                    readline(file)
                end
            end
            line = readline(file)
        end
        return line
    else
        @warn "'$filestr' Not a valid file."
    end
    return ""
end

function getpoints(s::String)
    idx = findfirst(r"POINTS = ", s)
    if !isnothing(idx)
        textafter = split(s[idx[end]+1:end])
        if length(textafter) > 0 
            return parse(Float32, textafter[1])
        end
    end
    return 0.0f0
end
getpoints(src::LineNumberNode) = getpoints(getsourceline(src))

function gettotalpoints(ts::CustomTestSet, testtext)
    d = ts.description
    lines = filter(testtext) do line
        occursin(d, line)
    end
    if !isempty(lines)
        ts.points_total = getpoints(lines[1])
    end
    for res in ts.results
        gettotalpoints(res, testtext)
    end
end
gettotalpoints(any, testtext) = nothing

function getnotebookfile(notebookname::String; repodir=joinpath(@__DIR__, ".."))
    name,ext = splitext(notebookname)
    if isempty(ext)
        notebookname *= ".ipynb"
    elseif ext != ".ipynb"
        error("Input file must be a Jupyter notebook with extension .ipynb")
    end
    nbfile = joinpath(repodir, "src", notebookname)
    if !isfile(nbfile)
        error("$nbfile does not exist.")
    end
    return nbfile
end

function gradequestion(notebookname::String; verbose=true)
    nbfile = getnotebookfile(notebookname)
    name = splitdir(splitext(nbfile)[1])[2]

    # Convert file to regular Julia file
    #   This is so that we can easily extract lines out of the source code to get 
    #   point values
    outfile = joinpath(dirname(nbfile), "$(name)_tmp.jl")
    nbexport(outfile, nbfile)

    # Run the file in the custom test set
    ts = @testset CustomTestSet "$name" begin
        include(outfile)
    end
    rm(outfile)  # Delete source file

    # Generate a file with just the tests
    testfile = joinpath(dirname(nbfile), "$(name)_testsets.jl")
    nbexport(testfile, nbfile, regex=r"\s*@testset*", markdown=false)
    testtext = readlines(testfile)
    rm(testfile)

    # Add up all the points for the problem
    gettotalpoints(ts, testtext)
    for res in ts.results
        if res isa CustomTestSet
            ts.points_total += res.points_total
        end
    end

    # Print results
    if verbose
        printresult(ts)
    end
    return (ts.points_missed, ts.points_total), ts
end

"""
Check if the test sets of the current repo are equal to the test sets of the 
    solution repo, to make sure they haven't been tampered with.
"""
function checktestsets(notebookname::String, solutiondir::String)
    nbfile = getnotebookfile(notebookname)
    name = splitdir(splitext(nbfile)[1])[2]
    outfile = joinpath(dirname(nbfile), "$(name)_testsets.jl")
    nbexport(outfile, nbfile, regex=r"\s*@testset*", markdown=false)

    solfile = getnotebookfile(notebookname, repodir=solutiondir)
    soltests = joinpath(dirname(nbfile), "solution_testsets.jl")
    nbexport(soltests, solfile, regex=r"\s*@testset*", markdown=false)

    solcontents = readlines(soltests)
    nbcontents = readlines(outfile)
    rm(soltests)
    rm(outfile)
    if (solcontents != nbcontents)
        @warn "Student's testsets do not match the solution test sets for $name. They may have been tampered with."
    end
end
