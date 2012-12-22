# The program below is extracted from the following reference:

# Lepage, G. Peter. 1998. Lattice QCD for novices. In <em>Strong Interactions
# at Low and Intermediate Energies: Proceedings of the 13th Annual HUGS
# at CEBAF, Jefferson Laboratory</em> (edited by  J. L. Goity), pp. 49–90.
# River Edge, N.J.: World Scientific.
# Preprint available at arxiv.org/abs/hep-lat/0506036

# To understand what the code does and how to interpret the results, please
# see the Lepage article.

# The only changes made in the published program are minor modifications
# needed to make it compatible with more-recent versions of Python. In
# particular, the older Numeric package has been replaced with the modern
# NumPy package.

tic()

function S(j,x)         # harm. osc. S
    if j == 1
        jm = N
    else
        jm = (j-1) % N  # previous site
    end
    if j == N
        jp = 1
    else
        jp = (j+1)  # next site
    end
    return a*x[j]^2/2 + x[j]*(x[j]-x[jp]-x[jm])/a
end


function update(x,repeat)
    eps = .1
    for i = 1:repeat
        for j = 1:N
            old_x = x[j]
            old_Sj = S(j,x)
            x[j] += eps*(2*rand() - 1)    # update x[j]
            dS = S(j,x) - old_Sj                # change in action
            if (dS > 0) & (exp(-dS) < rand())
                x[j] = old_x                    # restore old value
            end
        end
    end
end


function compute_G(x)
    return [sum([x[j] * x[(j + n) % N] for j = 1:N-1])/N for n = 1:N-1]
end

function MCaverage(x,G)
    update(x, 5 * N_cor)   # thermalize x
    for alpha in 1:N_cf     # loop on random paths
        update(x, N_cor)
        G[alpha] = compute_G(x)
    end
    for n in 1:N            # compute MC averages
        avg_G = sum([G[alpha][n] for alpha = 1:N_cf]) / N_cf
        #print("G(%d) = %g" % (n,avg_G) )
    end
    return x, avg_G
end

function bootstrap(G)
    N_cf = len(G)
    # choose random config from G ensemble
    G_bootstrap = [G[randi(N_cf)] for i = 1:N_cf]
    #G_bootstrap = array(G)[randint(N_cf, size=N_cf)]
    return G_bootstrap
end

function bin(G,binsize)
    G_binned = []                       # binned ensemble
    for i = len(G):binsize   # loop on bins
        G_avg = sum([G[i + j] for j in 1:binsize])  # loop on bin elements
        append_any(G_binned,G_avg/binsize)  # keep bin avg
    end
    return G_binned
end


# set parameters:

N = 20
N_cor = 50
N_cf = 100
a = 0.5
eps = 1.4

# create arrays:
x = zeros(N)
G = zeros((N_cf,N))
# do the simulation:

x, avg_G = MCaverage(x,G)
println(x)

exit()
