#!/bin/julia
using CairoMakie
using LaTeXStrings

function plot_error_csv(output, plots...)
    fig = Figure(resolution=(1000, 600))
    ax = Axis(fig[1,1], xscale=log2, yscale=log10, xlabel="L", ylabel=L"error", aspect=1, title= L"$L^2$-norm of error")

    for (filename, label) in plots
        ns = Vector{Float64}()
        es = Vector{Float64}()

        for l in readlines(filename)
            if length(l) == 0 || l[1] == "#"
                continue
            end

            a = split(l, ',')
            if length(a) >= 2
                push!(ns, parse(Float64, a[1]))
                push!(es, parse(Float64, a[2]))
            end
        end

        lines!(fig[1,1], ns, es, label=label)
        scatter!(fig[1,1], ns, es)
    end

    Legend(fig[1,2], ax)
    save(output, fig)
end

#plot_error_csv("report/t_100_plot.png", ("report/t_100_re_1e2.csv", L"Re = 10^2, t = 100"), ("report/t_100_re_1e4.csv", L"Re = 10^4, t = 100"))
#plot_error_csv("report/t_10_plot.png",  ("report/t_10_re_1e2.csv", L"Re =10^2, t = 10"), ("report/t_10_re_1e4.csv", L"Re = 10^4, t = 10"))
#plot_error_csv("report/t_100_plot_kbc.png",  ("report/t_100_re_1e12_kbc.csv", L"Re=10^{12},t=100"), ("report/t_100_re_1e5_kbc.csv", L"Re=10^{5},t=100"))
plot_error_csv("report/plots/re_30000_lbgk.png",  ("report/data/re_30000_lbgk.csv", L"Re=30000,t=20"))

