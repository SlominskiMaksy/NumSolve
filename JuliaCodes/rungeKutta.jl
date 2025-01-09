module runge_kutta_metodes
export RungeKutta
using Plots
"""
Runge-Kutta Method to solve Differential Equations numerically.

# Arguments
- `x0::Float64`: Initial value of x at t = t_start.
- `t_start::Float64`: Starting time for the integration.
- `t_end::Float64`: End time for the integration.
- `n::Int64`: Number of evaluation points.
- `f::Function`: The function defining the differential equation, dx/dt = f(t, x).

# Returns
- Store a plot of the solutions at discrete points
- `x_k::Vector`: The numerical solution at the evaluation points.

# Example
```julia
f(t, x) = -2 * t * x
x0 = 1.0
t_start = 0.0
t_end = 1.0
n = 10
x_k = RungeKutta(x0, t_start, t_end, n, f)
"""
    function RungeKutta(x0::Float64,t_start::Float64,t_end::Float64,n::Int64,f::Function)
        t_k::Vector = LinRange(t_start,t_end,n)
        h::Float64 = (t_end - t_start)/n
        x_k::Vector = zeros(n)
        
        for i in 2:n
            k1 = f(t_k[i-1], x_k[i-1])
            k2 = f(t_k[i-1] + h / 2, x_k[i-1] + h / 2 * k1)
            k3 = f(t_k[i-1] + h / 2, x_k[i-1] + h / 2 * k2)
            k4 = f(t_k[i-1] + h, x_k[i-1] + h * k3)
            x_k[i] = x_k[i-1] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        end
        
        plot1 = scatter(t_k, x_k, label="Evaluation Points", title="Solution of dx/dt = f(t,x) with Runge-Kutta")
        plot!(plot1, t_k, x_k, label="Solution Curve")
        savefig(plot1, "graphicSolution.png")

        
        return x_k
    end
end