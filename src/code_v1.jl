##        REPLICATION CODE FOR BRUNNERMEIER & SANNIKOV (2014)

using(DifferentialEquations)

        ## Setting up functions

## Define the investment function

# Fabio: Complete Replication of investment.m except minor details:
# Fabio: not sure if theta, Phi and iota need to be defined inside or outside function
# Fabio: not sure if we have to return(Phi_ioto_temp), but let's keep it for now

function investment(q)
    theta = 10
    Phi = (q - 1)/theta                 # Phi is standard investment technology with adjustment cost
    iota = Phi + (theta*Phi^2)/2
    Phi_iota_temp = [Phi, iota]
    return(Phi_iota_temp)
end

## Define the fnct function

# Fabio: Talk to Tancrede about f, fp and F. biggest problem right now: f is nowhere defined. don't get it how it still works in matlab
# Fabio: The commented line of codes are the MATLAB ORIGINALS
dyn = zeros(8)

function fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)
    # [Phi iota] = investment(f(3));
    Phi = investment(f[3])[1]   # changed this to Phi (used to be Phi_fnct)
    iota = investment(f[3])[2]  # changed this to iota (used to be iota_fnct)
    Phi_iota = [Phi, iota]   # I am sure there is a better way to code this. also not sure if we can't take Phi_iota_temp...

    # psi_L = eta; psi_R = min(f(3)/f(4) + eta, 1);
    psi_L = eta
    psi_R = min(f[3]/f[4] + eta, 1)

    for n = 1:50
        global psi = (psi_L + psi_R)/2
        amplification = 1 - f[4]/f[3]*(psi - eta)

        # VOLATILITY COMPUTATION
        global sigma_eta_eta = sigma*(psi - eta)/amplification  # sigma_eta times eta
        global sigma_q = sigma_eta_eta*f[4]/f[3]
        global sigma_theta = sigma_eta_eta*f[2]/f[1]
        global risk_premium = - sigma_theta*(sigma + sigma_q)

        global household_premium = (a_ - a)/f[3] + delta - delta_ + risk_premium

        if household_premium > 0 # households want to hold more
            psi_R = psi
        else
            psi_L = psi
        end
    end

    mu_q = r - (a - iota)/f[3] - Phi + delta - sigma*sigma_q + risk_premium
    mu_eta_eta = -(psi - eta)*(sigma + sigma_q)*(sigma + sigma_q + sigma_theta) + eta*(a - iota)/f[3] + eta*(1 - psi)*(delta_ - delta)
    qpp = 2*(mu_q*f[3] - f[4]*mu_eta_eta)/sigma_eta_eta^2
    thetapp = 2*((rho - r)*f[1] - f[2]*mu_eta_eta)/sigma_eta_eta^2

    local fp = [f[2], thetapp, f[4], qpp]

    # remember: f = [theta(eta), theta'(eta), q(eta), q'(eta)]

    leverage = psi/eta
    rk = r + risk_premium
    r_k = r + household_premium

    global dyn = [psi, sigma_eta_eta, sigma_q, mu_eta_eta, mu_q, iota, leverage, rk, r_k]

    return(fp, dyn)

end

## Setting Parameters

a = 0.11
a_= 0.05
rho = 0.06
r = 0.05
sigma = 0.025
delta = 0.03
delta_ = 0.08

## Some preparation : Determine First-Best Q, and the worst Q (also q(0))

# Determine first-best q, q(eta_star) that is the price of capital when experts own the equilibrium share of capital
global qmax
global QL = 0
global QR = 10

for iter = 1:50
    qmax = (QL + QR)/2
    Phi = investment(qmax)[1]
    iota = investment(qmax)[2]
    value = (a - iota)/(r + delta - Phi)
    if iota > a # price is too large
        global QR = qmax
    elseif value > qmax # price still too large
        global QL = qmax
    elseif  r + delta < Phi
        println("First-best price is infinite")
        global QR = qmax
    else
        global QR = qmax
    end
end

# Determine the worst q, q_ = q(0) that is the price of capital when experts own no capital
global q_
global QL = 0
global QR = 10

for iter = 1:50
    global q_ = (QL + QR)/2
    Phi = investment(q_)[1]
    iota = investment(q_)[2]
    value = (a_ - iota)/(r + delta_ - Phi)
    if iota > a_ # price is too large
        global QR = q_
    elseif value < q_ # price still too large
        global QR = q_
    else
        global QL = q_
    end
end

## Main part of the code

# Things to do: Look at fnct.m and see how it relates to the paper at 16 point (i)

etaspan = (0.0, 1.0)
F0 = [1, -1e+10, q_, 0]  # [theta(0), theta'(0), q(0), q'(0)]

#options = odeset('RelTol',1e-08,'AbsTol',1e-10, 'events','evntfcn');

# Defining stopping function

# f = [theta(eta), theta'(eta), q(eta), q'(eta)]
# condition (returns true or false), relates to events
function condition(f, eta, integrator)
    return integrator.u[3] > qmax || integrator.u[2] == 0 || integrator.u[4] == 0
end

function affect!(integrator)
    println("terminating")
    terminate!(integrator)
end

cb = ContinuousCallback(condition, affect!)

# u p t
odefun(f, dyn, eta) = fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)[1]
parameters(f, dyn, eta) = fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)[2]

#sols = []
#time = zeros(1)
global QL = 0
global QR = 1e+15
for iter = 1:50
    F0[4] = (QL + QR)/2
    prob = ODEProblem(odefun, F0, etaspan, parameters)
    global sol = solve(prob, Tsit5(), RelTol = 1e-08, abstol = 1e-10)  # sol = solve(prob, Tsit5(), callback = cb, RelTol = 1e-08, abstol = 1e-10)
    #push!(sols, sol)
    #global t = sol.u
    #push!(time,t)
    if sol.u[end][4] == 0 # if q'(eta) has reached zero, we
        global QL = F0[4]  # increase q'(0)
    else        # if q(eta) reached qmax or theta'(0) reached 0
        global QR = F0[4]
    end
end

#    if IE == 3, % if q'(eta) has reached zero, we
#        QL = F0(4);  % increase q'(0)
#    else        % if q(eta) reached qmax or theta'(0) reached 0
#        QR = F0(4);  % we reduce q'(0)
#    end

using Plots
#plot(sols[1])
#plot!(sols[2])
plot(sol, xaxis = "η")
#plot(sols[end], xaxis = "η")


## Computing other variables

#N = length(sols)
N = length(sol.t)
dynout = zeros(N, 9)

for n = 1:N
    fp = fnct(sol.t[n], sol.u[n], r, rho, a, a_, delta, delta_, sigma)[1]
    dynout[n,:] = fnct(sol.t[n], sol.u[n], r, rho, a, a_, delta, delta_, sigma)[2]
end

# We normalize the graphs
fout = hcat(sol.u...)'
normalization = fout[N, 1]
fout[:,1] = fout[:,1]./normalization
fout[:,2] = fout[:,2]./normalization

# Figure 1
p1 = plot(sol.t, fout[:,3], xaxis = "η", yaxis = "q")
p2 = plot(sol.t, fout[:,1], xaxis = "η", yaxis = "theta")
p3 = plot(sol.t[1:N-1], dynout[1:N-1,1], xaxis = "η", yaxis = "psi")
p4 = plot(sol.t[1:N-1], dynout[1:N-1,4], xaxis = "η", yaxis = "η mu^η")
p5 = plot(sol.t[1:N-1], dynout[1:N-1,2], xaxis = "η", yaxis = "η sigma^η")
p6 = plot(sol.t[1:N-1], dynout[1:N-1,3], xaxis = "η", yaxis = "sigma^q")
p7 = plot(sol.t[1:N-1], dynout[1:N-1,6], xaxis = "η", yaxis = "iota")
p8 = plot(sol.t[1:N-1], dynout[1:N-1,7], xaxis = "η", yaxis = "expert leverage")
p9 = plot(sol.t[1:N-1], dynout[1:N-1,8], xaxis = "η", yaxis = "returns")
p10 = plot!(p9, sol.t[1:N-1], dynout[1:N-1,9])
plot(p1, p2, p3, p4, p5, p6, p7, p8, p10, layout = (3, 3), legend = false)

# Figure 2
expert = sol.t.*fout[:,3].*fout[:,1]
household = (1 .- sol.t).*fout[:,3]
plot(expert,household, xaxis = "expert utility", yaxis = "household utility")
