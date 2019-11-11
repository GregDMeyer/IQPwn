
using ArgParse
include("solver.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-N"
            help = "number of samples to generate"
            arg_type = Int
            default = 4096
        "-o"
            help = "output file path"
            default = "samples.dat"
        "-s"
            help = """print formatted secret vector s instead of generating samples.\noptions are "base64" or "bin" """
        "program"
            help = "path to the X-program file"
            required = true
    end

    return parse_args(s)
end

function printkey(key, fmt)
    if !(fmt in ("base64", "bin"))
        println("""Unrecognized key printing format "$(args["s"])" """)
    end

    println(stderr, """Secret key:""")

    if fmt == "base64"
        println("$(vectob64(key))")
    else
        println("$(vectobin(key))")
    end
end

function main()
    args = parse_commandline()

    println(stderr, """Loading X-program at '$(args["program"])'...""")
    P = readprogram(args["program"])

    println(stderr, "Extracting secret key...")
    s, _, _ = extractkey(P)

    if !isnothing(args["s"])
        printkey(s, args["s"])
    else
        println("Generating samples...")
        samples = gensamples(s, args["N"])
        #samples = papersoln(P, args["N"])
        arraytofile(samples, args["o"])
        println("""Samples written to file '$(args["o"])'""")
    end
end

main()
