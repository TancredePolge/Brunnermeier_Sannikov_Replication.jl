##        REPLICATION CODE FOR BRUNNERMEIER & SANNIKOV (2014)


        ## Setting up functions

## Define the investment function

function investment(q)
    theta = 10
    Phi = (q - 1)/theta                 # Phi is standard investment technology with adjustment cost
    iota = Phi + (theta*Phi^2)/2
    Phi_iota_temp = [Phi, iota]
    return(Phi_iota_temp)
end 

## Define the event function 


## Define the fnct function





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