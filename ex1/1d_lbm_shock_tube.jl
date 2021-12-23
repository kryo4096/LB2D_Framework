using GLMakie


# parameters
beta = 0.99 # beta for relaxation
s = false # set to true to plot results instead of showing them live

c = [-1,0,1]
c_s = 1 / √3

W = [1.0 / 6., 2.0 / 3., 1.0 / 6.]

N_nodes = 800
N_steps = 100000

F = zeros(N_nodes + 2, 3)

f_eq(rho, u) = rho * W .* (1 .+ c * u / c_s^2 + (c.^2 .- c_s^2) * u^2 / (2 * c_s^4))

rho = Observable(zeros(N_nodes))
u = Observable(zeros(N_nodes))
P = Observable(zeros(N_nodes))



for i in 1:N_nodes
    rho.val[i] = i > N_nodes / 2 ? 1.5 : 1.0
    u.val[i] = 0.0
    F[i+1,:] = f_eq(rho.val[i], u.val[i])
end

fig = Figure()

ax1 = Axis(fig[1,1], ylabel = "u", title="beta = $beta")
lines!(fig[1,1], 1:N_nodes, u)
autolimits!(ax1)

ax2 = Axis(fig[2,1], ylabel="ρ")
lines!(fig[2,1], 1:N_nodes, rho)
autolimits!(ax2)

display(fig)

u_ = zeros(N_nodes)
rho_ = zeros(N_nodes)

if s
	save("report/plots/shock_plot_initial_conditions.png", fig)
end

for t ∈ 1:N_steps - 1
 
    for k ∈ 2:N_nodes+1
        F[k - 1, 1] = F[k, 1]
    end

    for k in reverse(2:N_nodes+1)
        F[k + 1, 3] = F[k, 3]
    end

    F[N_nodes+1, 1] = F[N_nodes+2, 3]
    F[2, 3] = F[1, 1]

    for k ∈ 2:N_nodes+1
        rho_[k-1] = sum(F[k, :])
        u_[k-1] = sum(F[k, :] .* c) / rho_[k-1]
        F[k, :] += 2 * beta * (f_eq(rho_[k-1], u_[k-1]) - F[k,:])
    end

    if t%10 == 0
        rho[] = rho_
        u[] = u_
        autolimits!(ax1)
        autolimits!(ax2)
        sleep(0.001)


        if s && t%100 == 0
        	save("report/plots/shock_plot_t_$t.png", fig)
        end
    end
end





