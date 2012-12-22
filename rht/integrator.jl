#load("Winston")
#using Winston
load("Plot")


function leapfrogintegrator(gradU, p, q, epsilon, L)
    # L number of leapfrog steps
    # epsilon leapfrog stepsize

    # make a half step at the beginning
    p -= epsilon * gradU(q) / 2
    q += epsilon * p

    for i=1:L-1
        p += -epsilon * gradU(q)
        q += epsilon * p
    end


    # make a half step for momentum at the end
    p -= epsilon * gradU(q) / 2

    return p, q
end

x = linspace( 0, 3pi, 100 )
f, a = leapfrogintegrator(u -> -u, x, x*3, .1, 100)

p = FramedPlot()
setattr(p, "title", "title!")

setattr(p, "xlabel", L"\Sigma x^2_i")
setattr(p, "ylabel", L"\Theta_i")

add(p, Curve(x, f, "color", "red") )
x11(p)

#file(p, "example1.eps")
#file(p, "example1.png")
