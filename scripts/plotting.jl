using CairoMakie
CairoMakie.activate!(; type="svg", visible=true)
using DataFrames
import JSON3
import JSONTables
using Printf

#------------------------------Load Data------------------------------

df = open("/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/data/data.jsonl","r") do f
    DataFrame(JSON3.read(f; jsonlines=true))
end

gdf = groupby(df, Not(:vs, :ks, :date));
ugdf = subset(gdf, :date => x-> x.==maximum(x); ungroup=false);
df = combine(ugdf, All());

#------------------------------Growth Rate----------------------------

#------------------------------baroclinic instabiliy------------------
df_bar = filter(row -> row.Ri ≈ 100 && row.a ≈ 12.5, df);

f = Figure(size = (1000, 300));
ax = Axis(f[1,1];
          xlabel = "wavenumber / m⁻¹",
          ylabel = "growth rate",
          aspect = 2.5,
          );
for row in eachrow(df_sym)
    N = row[:N]
    f₀ = row[:f0]
    Ri = row[:Ri]
    M2 = √(N^2 * f₀^2 / Ri)
    if row[:type] == "true"
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%s", row[:type]))
    else
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%6.3f km, %s, %s-grid", row[:a], row[:type], row[:grid]))
    end
end
f[1,2] = Legend(f, ax; merge=true, valign=:top);
colsize!(f.layout, 1, Aspect(1, 2.5));
display(f)

#------------------------------symmetric instability------------------
df_sym = filter(row -> row.Ri ≈ 0.5 && row.thetaU ≈ 0.0, df);

f = Figure(size = (1000, 500));
ax = Axis(f[1,1];
          xlabel = "wavenumber / m⁻¹",
          ylabel = "growth rate",
          aspect = 1.5,
          limits = (nothing, 1.0e-3, nothing, nothing),
          );
for row in eachrow(df_sym)
    N = row[:N]
    f₀ = row[:f0]
    Ri = row[:Ri]
    M2 = √(N^2 * f₀^2 / Ri)
    if row[:type] == "true"
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%s", row[:type]))
    else
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%6.3f km, %s, %s-grid", row[:a], row[:type], row[:grid]), linewidth=1)
    end
end
f[1,2] = Legend(f, ax; merge=true, valign=:top);
colsize!(f.layout, 1, Aspect(1, 1.5));
display(f)
save("plots/test.svg",f)

#------------------------------Galilean Invariance--------------------
grid = "tri-C";
type = "ICON-mimetic";
grid = "hex-C";
type = "MPAS";
grid = "tri-B";
type = "direct-advective";
df_gal = filter(row -> row.Ri ≈ 0.5 && row.a ≈ 3.125, df)

f = Figure(size = (1000, 500));
ax = Axis(f[1,1];
          xlabel = "wavenumber / m⁻¹",
          ylabel = "growth rate",
          aspect = 1.5,
          limits = (nothing, 1.0e-3, nothing, nothing),
          );
for row in eachrow(df_gal)
    N = row[:N]
    f₀ = row[:f0]
    Ri = row[:Ri]
    M2 = √(N^2 * f₀^2 / Ri)
    if row[:type] == "true"
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%s", row[:type]))
    else
        lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("u₀=%7.3f ms⁻¹,θU/π=%5.3f,%s", row[:u0], row[:thetaU]/π, row[:grid]), linewidth=1)
    end
end
# add true instability
row = first(filter(row -> row.Ri ≈ 0.5 && row.a < 1e-13, df))
N = row[:N]
f₀ = row[:f0]
Ri = row[:Ri]
M2 = √(N^2 * f₀^2 / Ri)
lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%s", row[:type]))

f[1,2] = Legend(f, ax; merge=true, valign=:top);
colsize!(f.layout, 1, Aspect(1, 1.5));
display(f)
save("plots/gal-tri-B-direct-advective.svg",f)

#------------------------------Mesh alignment-------------------------
grid = "tri-B";
type = "direct-advective";
df_ali = filter(row -> row.Ri ≈ 0.5 && row.a ≈ 3.125 && row.type == type && row.grid == grid, df)
f = Figure(size = (1000, 500));
ax = Axis(f[1,1];
          xlabel = "wavenumber / m⁻¹",
          ylabel = "growth rate",
          aspect = 1.5,
          limits = (nothing, 1.0e-3, nothing, nothing),
          );
for row in eachrow(df_ali)
    N = row[:N]
    f₀ = row[:f0]
    Ri = row[:Ri]
    M2 = √(N^2 * f₀^2 / Ri)
    lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("θU/π=%5.3f", row[:thetaU]/π), linewidth=1)
end
# add true instability
row = first(filter(row -> row.Ri ≈ 0.5 && row.a < 1e-13, df))
N = row[:N]
f₀ = row[:f0]
Ri = row[:Ri]
M2 = √(N^2 * f₀^2 / Ri)
lines!(ax, row[:ks], row[:vs]*N/abs(M2), label=@sprintf("%s", row[:type]))

f[1,2] = Legend(f, ax; merge=true, valign=:top);
colsize!(f.layout, 1, Aspect(1, 1.5));
display(f)
save("plots/ali-tri-B-direct-advective-sym.svg",f)

#------------------------------Real and imaginary part----------------
row = df_bar[2,:]
imax = argmax(row[:vs]) # [df[2,:ks] .> 1e-4]
Nz = 32
u⃗ = reshape(row[:vecs][imax,1:4*Nz], Nz, 2, 2);
u = u⃗[:,1,:];
v = u⃗[:,2,:];
b = row[:vecs][imax,4*Nz+1:end-1];
pu = plot()
plot!(pu, legend=:outerbottom, framestyle=:box)
plot!(pu, real.(u[:,1]), range(start=0., stop=4000., length=Nz), label="real part")
plot!(pu, imag.(u[:,1]), range(start=0., stop=4000., length=Nz), label="imaginary part")

pv = plot()
plot!(pv, legend=:outerbottom, framestyle=:box)
plot!(pv, real.(v[:,1]), range(start=0., stop=4000., length=Nz), label="real part")
plot!(pv, imag.(v[:,1]), range(start=0., stop=4000., length=Nz), label="imaginary part")

pb = plot()
plot!(pb, legend=:outerbottom, framestyle=:box)
plot!(pb, real.(b), range(start=0., stop=4000., length=Nz), label="real part")
plot!(pb, imag.(b), range(start=0., stop=4000., length=Nz), label="imaginary part")

p = plot(pu, pv, pb, layout=(1,3), size=(1000, 600))

display(p)

