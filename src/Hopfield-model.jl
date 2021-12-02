using Distributions
using SpecialFunctions
using JLD2
using FileIO
using Distributed
using LinearAlgebra




function start_pattern(l,up) # initiate 1 pattern
    sign.(rand(l).-up)
end

function get_patterns(l,n,up) # initiate n patterns
    xis=zeros(l,n)
    for ii in 1:n
        xis[:,ii]=start_pattern(l,up)
    end
    return xis
end







function updateactivestate(states)  #Get the state that is presented to the network
    return rand(1:states)
end








function build_matrix_xi_x_runin(l,n,x,γ,runin)   # initiate J_{ij} and n patterns  % γ = (1-λ)
        xi= get_patterns(l,n,0.5)
            xi,J = let_run_in_Hop_uniform(xi,γ,x,runin,n)
    return J,xi
end


function let_run_in_Hop_uniform(xi,γ,x,runs,states)  # initiate J_{ij}
    activestate=updateactivestate(states)


    J = start_J(xi[:,activestate],γ)

    for ii in 1:runs
        xi[:]=new_XI_rand(xi[:], x)
        J=new_J(J,xi[:,activestate],γ)
        activestate=updateactivestate(states)

    end
    return xi,J
end



function start_J(xi,γ)   # Start initial J_ij
    l =length(xi)
    J=zeros(l,l)
    for jj in 1:l
        for ii in 1:jj-1
            J[ii,jj] = (1-γ)* xi[ii]*xi[jj]
        end
        for ii in jj+1:l
            J[ii,jj] = (1-γ)* xi[ii]*xi[jj]
        end

    end
    return J
end






function new_XI_rand(old_XI, x) # random mutations of patterns
    l=length(old_XI)
    for ii in 1:l
        if (rand() < x)
            old_XI[ii] *=-1
        end
    end
    return old_XI
end





function new_J(old_J::Array{Float64,2},xi::Array{Float64,1},γ::Float64)  # update interaction with stadard Hebbian learning # γ = (1-λ)
    l=length(xi)



    for jj in 1:l
        for ii in 1:jj-1
            old_J[ii,jj] = γ* old_J[ii,jj] + (1-γ)* xi[ii]*xi[jj]
        end
        for ii in jj+1:l
            old_J[ii,jj] = γ* old_J[ii,jj] + (1-γ)* xi[ii]*xi[jj]
        end
    end
    return old_J
end

function energie(J,xi)  #energy of a pettern xi in network J
    l=length(xi)

    return (transpose(xi)*J*xi)/l/2
end










function localfiels(xi,J)
    return sum(sign.(xi.*(J*xi)))
end








function Support_entropy(xi,evecs,energies)   ## Support of the patterns in the eigenspace of the network given in units of entropy
    n=length(energies)
    l=length(xi[:,1])
    entropy_vec=Vector{Float64}()
    entropy_vec_normed=Vector{Float64}()
    overlaps=((transpose(xi)*evecs)./(sqrt(l))).^2

    for p in 1:n
        local_entropy=Vector{Float64}()
        local_entropy_normed=Vector{Float64}()
        norm=sum(overlaps[:,p])
        for ii in 1:n
            if (overlaps[ii,p]>0.000001)
                append!(local_entropy,-overlaps[ii,p]* log2(overlaps[ii,p]))
                append!(local_entropy_normed,-overlaps[ii,p]/norm* log2(overlaps[ii,p]/norm))
            end
        end
        append!(entropy_vec,sum(local_entropy))
        append!(entropy_vec_normed,sum(local_entropy_normed))
    end




    entropy_patts=Vector{Float64}()


    overlepsp=transpose(overlaps)




    for p in 1:n
        local_entropy=Vector{Float64}()
        for ii in 1:n
            if (overlepsp[ii,p]>0.000001)
                append!(local_entropy,-overlepsp[ii,p]* log2(overlepsp[ii,p]))
            end
        end
        append!(entropy_patts,sum(local_entropy))
    end


    return entropy_vec,entropy_vec[1],entropy_vec[n],entropy_vec_normed,entropy_vec_normed[1],entropy_vec_normed[n],entropy_patts,entropy_patts[argmax(energies)],entropy_patts[argmin(energies)]
end






function Support_partisipation(xi,evecs,energies)  ## Support of the patterns in the eigenspace of the network given in units of partisipation
    n=length(energies)
    l=length(xi[:,1])
    entropy_vec=Vector{Float64}()

    overlaps=((transpose(xi)*evecs)./(sqrt(l))).^2

    for p in 1:n
        norm=sum(overlaps[:,p].^4)
        partis=sum(overlaps[:,p].^2)^2
        append!(entropy_vec,partis/norm)
        end




    entropy_patts=Vector{Float64}()


    overlepsp=((transpose(evecs)*xi)./(sqrt(l))).^2




    for p in 1:n
        norm=sum(overlepsp[:,p].^4)
        partis=sum(overlepsp[:,p].^2)^2
        append!(entropy_patts,partis/norm)
    end


    return mean(entropy_vec),std(entropy_vec),entropy_vec[1],entropy_vec[n],mean(entropy_patts),std(entropy_patts),entropy_patts[argmax(energies)],entropy_patts[argmin(energies)]
end
