##        REPLICATION CODE FOR BRUNNERMEIER & SANNIKOV (2014)


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

## Define the event function 

# Fabio: Following chunk doesn't work. probably because F isn't defined anywhere...
# Fabio: In my opinion capital F is the derivative of f. but then it's the same thing as fp?
# Fabio: Anyway F is for sure a 4x1 vector just like f or fp is
function eventfcn(eta, F)
    global qmax # not sure
    value = [(qmax - F[3]), F[2], F[4]]
    isterminal = [1, 1, 1]
    direction = [0, 0, 0]
    return (value, isterminal, direction)
end



## Define the fnct function

# Fabio: Talk to Tancrede about f, fp and F. biggest problem right now: f is nowhere defined. don't get it how it still works in matlab
# Fabio: The commented line of codes are the MATLAB ORIGINALS
function fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)
    # [Phi iota] = investment(f(3));
    Phi_fnct = investment(f[3])[1]
    iota_fnct = investment(f[3])[2]
    Phi_iota_fnct = [Phi_fnct, iota_fnct]   # I am sure there is a better way to code this. also not sure if we can't take Phi_iota_temp...

    # psi_L = eta; psi_R = min(f(3)/f(4) + eta, 1);
    psi_L = eta 
    psi_R = min(f[3]/f[4] + eta, 1) 

    for n = 1:50
        psi = (psi_L + psi_R)/2
        amplification = 1 - f[4]/f[3]*(psi - eta)
            
        # VOLATILITY COMPUTATION
        sigma_eta_eta = sigma*(psi - eta)/amplification  # sigma_eta *times* eta
        sigma_q = sigma_eta_eta*f[4]/f[3] 
        sigma_theta = sigma_eta_eta*f[2]/f[1]
        risk_premium = - sigma_theta*(sigma + sigma_q)
            
        household_premium = (a_ - a)/f[3] + delta - delta_ + risk_premium
                   
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

    fp = [f[2], thetapp, f[4], qpp]

    leverage = psi/eta
    rk = r + risk_premium
    r_k = r + household_premium

    dyn = [psi, sigma_eta_eta, sigma_q, mu_eta_eta, mu_q, iota, leverage, rk, r_k]

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
QL = 0
QR = 10 

for iter = 1:50
    qmax = (QL + QR)/2
    Phi = investment(qmax)[1]
    iota = investment(qmax)[2]
    value = (a - iota)/(r + delta - Phi)
    if iota > a # price is too large
        QR = qmax 
    elseif value > qmax # price still too large
        QL = qmax
    elseif  r + delta < Phi 
        println("First-best price is infinite")
        QR = qmax 
    else
        QR = qmax
    end
end

# Determine the worst q, q_ = q(0) that is the price of capital when experts own no capital 
QL = 0 
QR = 10 

for iter = 1:50 
    q_ = (QL + QR)/2
    Phi = investment(q_)[1]
    iota = investment(q_)[2]
    value = (a_ - iota)/(r + delta_ - Phi)
    if iota > a_ # price is too large
        QR = q_
    elseif value < q_ # price still too large
        QR = q_
    else 
        QL = q_
    end
end


        ## Main part of the code 


# Things to do: Look at fnct.m and see how it relates to the paper at 16 point (i)

etaspan = (0.0, 1.0)
F0 = [1, -1e+10, q_, 0]  # [theta(0), theta'(0), q(0), q'(0)]  

#options = odeset('RelTol',1e-08,'AbsTol',1e-10, 'events','evntfcn');


