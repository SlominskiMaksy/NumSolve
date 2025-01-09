module numericIntegration
export trapezInt
export simpsonInt
"""
Numerische Integration einer Funktion f(x) über die Trapezmethode.

# Argumente
- `f::Function`: Die zu integrierende Funktion.
- `a::Float64`: Die untere Grenze des Integrals.
- `b::Float64`: Die obere Grenze des Integrals.
- `n::Int`: Die Anzahl der Intervalle (optional, Standardwert ist 100).

# Rückgabewert
- Das approximierte Integral der Funktion f zwischen `a` und `b`.

# Beispiel
```julia
f(x) = x^2
result = trapez(f, 0.0, 1.0)
println(result)  # Erwartetes Ergebnis: ca. 0.3333
```
"""
function trapezInt(f::Function,a::Float64,b::Float64,n::Int64)
    auswertungsPunkte = LinRange(a,b,n)
    integralWert = 0
    for i in 1:n - 1
        integralWert = integralWert + (auswertungsPunkte[i + 1]-auswertungsPunkte[i])*(1/2 * f(auswertungsPunkte[i]) + 1/2 * f(auswertungsPunkte[i+1]))
    end
    return integralWert
end

"""
Numerische Integration einer Funktion f(x) über die Simpsonmethode.

# Argumente
- `f::Function`: Die zu integrierende Funktion.
- `a::Float64`: Die untere Grenze des Integrals.
- `b::Float64`: Die obere Grenze des Integrals.
- `n::Int`: Die Anzahl der Intervalle (optional, Standardwert ist 100).

# Rückgabewert
- Das approximierte Integral der Funktion f zwischen `a` und `b`.

# Beispiel
```julia
f(x) = x^2
result = trapez(f, 0.0, 1.0)
println(result)  # Erwartetes Ergebnis: ca. 0.3333
```
"""
function simpsonInt(f::Function,a::Float64,b::Float64,n::Int64)
    auswertungsPunkte = LinRange(a,b,n)
    integralWert = 0
    for i in 1:n - 1
        integralWert = integralWert + (auswertungsPunkte[i + 1]-auswertungsPunkte[i])*(1/6 * f(auswertungsPunkte[i]) + 4/6 * f(1/2 * (auswertungsPunkte[i] + auswertungsPunkte[i+1])) + 1/6 * f(auswertungsPunkte[i + 1]))
    end
    return integralWert
end
end