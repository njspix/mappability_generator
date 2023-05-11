function reader(infile, out_buff::Channel{Tuple{SubString{String},Int64,Char}})
    chr_name = ""
    pos = 0
    for line in eachline(infile)
        if line[1] == '>'
            chr_name = split(line[2:end], ' ')[1]
            pos = 1
        else
            for c in collect(uppercase(line))
                put!(out_buff, (chr_name, pos, c))
                pos += 1
            end
        end
    end
    close(out_buff)
end

function c2t_put(out_buff::Channel{Tuple{SubString{String},Int64,Char}}; c1, c2)
    c1_chr, c1_pos, c1_base = c1
    c2_chr, c2_pos, c2_base = c2

    if c1_chr == c2_chr
        if c1_base == '~'
            return
        elseif c1_base == 'C' && c2_base != 'G'
            put!(out_buff, (c1_chr, c1_pos, 'T'))
        else
            put!(out_buff, (c1_chr, c1_pos, c1_base))
        end
    else
        # C1 is the last base of its chr
        # Last base is never the C of a CpG
        if c1_base == 'C'
            put!(out_buff, (c1_chr, c1_pos, 'T'))
        else
            put!(out_buff, (c1_chr, c1_pos, c1_base))
        end
    end
end

function g2a_put(out_buff::Channel{Tuple{SubString{String},Int64,Char}}; c1, c2)
    c1_chr, c1_pos, c1_base = c1
    c2_chr, c2_pos, c2_base = c2

    if c1_chr == c2_chr
        if c2_base == '~'
            return
        elseif c2_base == 'G' && c1_base != 'C'
            put!(out_buff, (c2_chr, c2_pos, 'A'))
        else
            put!(out_buff, (c2_chr, c2_pos, c2_base))
        end
    else
        # C2 is the first base of its chr
        # First base is never the G of a CpG
        if c2_base == 'G'
            put!(out_buff, (c2_chr, c2_pos, 'A'))
        else
            put!(out_buff, (c2_chr, c2_pos, c2_base))
        end
    end
end

function converter(in_buff::Channel{Tuple{SubString{String},Int64,Char}}, out_buff::Channel{Tuple{SubString{String},Int64,Char}}, conv_func::Function)
    first_char = take!(in_buff)
    conv_func(out_buff, c1=(first_char[1], first_char[2] - 1, '~'), c2=first_char)
    try
        while true
            next_char = take!(in_buff)
            conv_func(out_buff, c1=first_char, c2=next_char)
            first_char = next_char
        end
    catch e
        if e isa InvalidStateException
            true
        else
            rethrow(e)
        end
    finally
        conv_func(out_buff, c1=first_char, c2=(first_char[1], first_char[2] + 1, '~'))
        close(out_buff)
    end
end

function only_n(lst::Vector{Char})
    for c in lst
        if c != 'N'
            return false
        end
    end
    return true
end

function fastq_format(outfile, pos, a::Vector{Char}, name, qual="I")
    if length(a) > 0 && !only_n(a)
        write(outfile, "$(name):$(pos-length(a)+1)\n$(String(a))\n+\n$(qual^length(a))\n")
        flush(outfile)
    end
end

function readmaker(in_buff::Channel{Tuple{SubString{String},Int64,Char}}, read_length::Int, name_infix::String, outfile::String)
    outfile = open(outfile, "w")
    read = Vector{Char}()
    curr_chr = ""
    pos = 0
    try
        while true
            new_chr, new_pos, new_base = take!(in_buff)
            if new_chr != curr_chr
                fastq_format(outfile, pos, read, "$(curr_chr):$(name_infix)")
                # Clear out the buffer
                read = Vector{Char}()
                push!(read, new_base)
                curr_chr = new_chr
                pos = 1
                println("$(name_infix): Working on contig $(curr_chr)")
            else # we're continuing on the same chr
                if length(read) == read_length
                    fastq_format(outfile, pos, read, "$(curr_chr):$(name_infix)")
                    read = read[2:end]
                    push!(read, new_base)
                    pos += 1
                else
                    push!(read, new_base)
                    pos += 1
                end
            end
            pos % 100000 == 0 && println("$(name_infix): $(round(pos/1e6, digits = 1)) Mb processed (contig $(curr_chr))")
        end
    catch e
        if e isa InvalidStateException
            println("$(name_infix): all done with writing")
        else
            rethrow(e)
        end
    finally
        fastq_format(outfile, pos, read, "$(curr_chr):$(name_infix)")
    end
    close(outfile)
end

function process(infile, read_length, out_pfx, buff_size)
    raw1 = Channel{Tuple{SubString{String},Int64,Char}}(buff_size)
    raw2 = Channel{Tuple{SubString{String},Int64,Char}}(buff_size)
    c2t_buff = Channel{Tuple{SubString{String},Int64,Char}}(buff_size)
    g2a_buff = Channel{Tuple{SubString{String},Int64,Char}}(buff_size)
    errormonitor(@async reader(infile, raw1))
    errormonitor(@async reader(infile, raw2))
    errormonitor(@async converter(raw1, c2t_buff, c2t_put))
    errormonitor(@async converter(raw2, g2a_buff, g2a_put))
    t1 = @task readmaker(c2t_buff, read_length, "c2t", "./$(out_pfx)_c2t.fastq")
    t2 = @task readmaker(g2a_buff, read_length, "g2a", "./$(out_pfx)_g2a.fastq")
    errormonitor(schedule(t1))
    errormonitor(schedule(t2))
    wait(t1)
    wait(t2)
end

using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--read_length", "-l"
    help = "Length of reads to generate"
    arg_type = Int
    default = 100
    "--out_pfx", "-o"
    help = "Prefix for output files"
    arg_type = String
    default = "out"
    "--buff_size", "-b"
    help = "Size of buffer for each channel"
    arg_type = Int
    default = 1000000
    "infile"
    help = "Input fasta file"
    arg_type = String
    required = true
end

a = parse_args(s)
process(a["infile"], a["read_length"], a["out_pfx"], a["buff_size"])