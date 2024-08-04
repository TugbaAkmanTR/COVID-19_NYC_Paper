#
# Source:
#
#@misc{cite1,
#  author = "Alexander Demin",
#  date = "2024-05-27",
#  howpublished = "personal communication"
#}

#@Online{cite2,
# author = {Pogudin, Gleb},
# year = {2024},
# title = {Speeding things up via linear first integrals},
# url = {https://github.com/SciML/StructuralIdentifiability.jl/issues/63},
# urldate = {2020-07-01}
#}

using StructuralIdentifiability

#Transformed model

ode = @ODEmodel(
    x2'(t) = ((betae*x2(t) + betae*x3(t) + betas*x4(t)) * (C1 - x2(t) - int_y1(t) / rho)/(C1 - int_y1(t) / rho + x3(t)+x4(t)+x6(t))) - (1/14)*x2(t),
    x3'(t) = (1-rho)*(1/14)*x2(t) - (1/7)*x3(t),
    x4'(t) = rho*(1/14)*x2(t) - thetaa*(4199 + C2 - phi * int_y3(t) - x5(t) + (1/2.9))*x4(t),
    x5'(t) = thetaa*(4199 + C2 - phi * int_y3(t) - x5(t))*x4(t) - (gammah + nu)*x5(t),
    x6'(t) = (1/7)*x3(t) + (1/2.9)*x4(t) + gammah*x5(t),
    int_y1'(t) = (1/14) * rho * x2(t),
    int_y3'(t) = x5(t),
    y1(t) = (1/14)*rho*x2(t),
    y2(t) = nu,
    y3(t) = x5(t),
    y4(t) = int_y1(t),
    y5(t) = int_y3(t),
)

println(assess_identifiability(ode))

assess_identifiability(ode)


