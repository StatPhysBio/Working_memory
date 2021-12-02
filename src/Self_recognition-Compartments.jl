

include("Hopfield-model.jl")
include("Analysis_tools.jl")







function equiliabration(J,xi,β,steps)  # equilibration of a single pattern in a given network J
    for ii in 1:steps


        switchup=rand(1:length(xi)) # choose a given spinflip

        if (exp(-β*(2*(xi[ switchup])*(transpose(J[:,switchup])*xi)/length(xi)))>rand()) # test the energy diffrence from this spinflip

            xi[switchup]*=-1
        end
    end
        return xi
    end




    function Vec_equiliabration_compartment(J,xi,β,βeq,steps)   # equilibration of a set of patterns for multiple compartments
         xi_final=copy(xi)
        for ii in 1:length(xi_final[1,:])

            xi_final[:,ii]=equiliabration(J[:,:,givecompartment(J,xi[:,ii],β)],xi[:,ii],βeq,steps)
        end
        return xi_final
    end






    function Vec_equiliabration(J,xi,β,steps)  # equilibration of a set of patterns for a single network
         xi_final=copy(xi)
        for ii in 1:length(xi_final[1,:])
            xi_final[:,ii]=equiliabration(J,xi[:,ii],β,steps)
        end
        return xi_final
    end









function energies_to_compartments(J,xi)  # energies for the C compartments (compartments)
    es=Vector{Float64}()
    for ii in 1:length(J[1,1,:])
        append!(es,energie(J[:,:,ii],xi))
    end
    return es
end


function givecompartment(J,xi,β) # choose a compartment to update/ find stored state
    probs=exp.(β .*energies_to_compartments(J,xi))
    d = Categorical(probs./sum(probs))
    return rand(d)
end

function givecompartment_stats(J,xi,β) # choose a compartment to update/ find stored state and return the energy
    es=energies_to_compartments(J,xi)
    probs=exp.(β .* es)
    d = Categorical(probs./sum(probs))
    pos=rand(d)
    return pos,es[pos]
end


function update_compartments(J,xi,γ,β) # update a system with many compartments
    compartment=givecompartment(J,xi,β)    # find the compartment
    J[:,:,compartment]=new_J(J[:,:,compartment],xi,γ) # update the compartment
    return J
end

function update_compartments_stats(J,xi,γ,β) # update a system with many compartments and give the energy

    compartment,e_take=givecompartment_stats(J,xi,β)

    J[:,:,compartment]=new_J(J[:,:,compartment],xi,γ)
    return J,compartment,e_take
end


function start_compartment(xi,γ,compartments) # initiate system with multiple compartmetns only intitial state J(0)
    l=length(xi[:,1])
    J=zeros(l,l,compartments)
    for ii in 1:compartments
        J[:,:,ii]=start_J(xi[:,ii],1- 0.95*compartments/length(xi[1,:]) )
    end
    for ii in compartments:length(xi[1,:])
        b=convert(Int64,((ii+(compartments-1)) % compartments) +1 )
        J[:,:,b]=new_J(J[:,:,b],xi[:,ii],1- 0.95*compartments/length(xi[1,:]) )
    end
    return J
end





function build_compartment(l,n,x,γ,β,runin,compartments) # give states and system with multiple compartments after inital steps
        if (compartments > n)
            println(" to many compartments the maximum is the number of patters here $(n) you have selected $(compartments) ")
            return -1
        end

        xi= get_patterns(l,n,0.5)
            xi,J = let_run_in_compartment(xi,γ,x,β,runin,n,compartments)
    return J,xi
end


function let_run_in_compartment(xi,γ,x,β,runin,n,compartments) # initian of system and patterns and making it independed of starting conditions

        J= start_compartment(xi,γ,compartments)

    for ii in 1:runin
        activestate=updateactivestate(n)
        xi[:]=new_XI_rand(xi[:], x)
        J=update_compartments(J,xi[:,activestate],γ,β)
    end
    return xi,J
end




function get_false_compartment_withe_compartment(l,n,x,γ,β,runin,compartments,steps,reps)

    states=Vector{Int16}()
    compartment=Vector{Int16}()
    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()


    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)

            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)
                append!(states,activestate)
                append!(compartment,buc )
                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end


    end
    return states,compartment,energies_take ,energies_rand
end





function get_false_compartment_stats_compartment(l,n,x,γ,β,runin,compartments,steps,reps)

    states=Vector{Int16}()
    compartment=Vector{Int16}()
    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()


    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)

            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)
                append!(states,activestate)
                append!(compartment,buc )
                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end


    end
    return sum(sign.(abs.((convert.(Int64,(states.+(compartments-1)) .% compartments .+1 )) .-compartment)))/length(compartment),mean(1 ./ (1 .+exp.(β .*(energies_take.-energies_rand))))
end



function get_false_compartment_stats_compartment_alongx(l,n,xs,γ,β,runin,compartments,steps,reps)
    d1=get_false_compartment_stats_compartment.(l,n,xs,γ,β,runin,compartments,steps,reps)
    e1s=zeros(length(xs))
e2s=zeros(length(xs))
for ii in 1 :(length(xs))
    e1s[ii] = d1[ii][1]
    e2s[ii] = d1[ii][2]
end
return e1s,e2s
end

function get_false_compartment_stats_compartment_alonggamma(l,n,x,γs,β,runin,compartments,steps,reps)
    d1=get_false_compartment_stats_compartment.(l,n,x,γs,β,runin,compartments,steps,reps)
    e1s=zeros(length(γs))
e2s=zeros(length(γs))
for ii in 1 :(length(γs))
    e1s[ii] = d1[ii][1]
    e2s[ii] = d1[ii][2]
end
return e1s,e2s
end




function get_statistics_withe_compartment(l,n,x,γ,β,runin,compartments,steps,reps,eqsteps,βeq)

    states=Vector{Int16}()
    compartment=Vector{Int16}()
    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)
            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)
                append!(states,activestate)
                append!(compartment,buc )
                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end

            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return states,compartment,energies_take ,energies_rand,overlaps
end

function get_statistics_compartment(l,n,x,γ,β,runin,compartments,steps,reps,eqsteps,βeq)

    states=Vector{Int16}()
    compartment=Vector{Int16}()
    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)
            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)
                append!(states,activestate)
                append!(compartment,buc )
                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end

            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return states,compartment,sum(sign.(energies_take .-maximum(energies_rand)).-1)./(length(energies_take)*(-2)),sum(sign.(energies_rand .-minimum(energies_take)).+1)/(2*length(energies_rand)),overlaps
end





function get_statistics_compartment_no_use_estats(l,n,x,γ,β,runin,compartments,steps,reps,eqsteps,βeq)


    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)
            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)

                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end

            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return sum(sign.(energies_take .-maximum(energies_rand)).-1)./(length(energies_take)*(-2)),sum(sign.(energies_rand .-minimum(energies_take)).+1)/(2*length(energies_rand)),mean(energies_take),std(energies_take),lower_std(energies_take),upper_std(energies_take),overlaps
end



function get_statistics_compartment_no_use(l,n,x,γ,β,runin,compartments,steps,reps,eqsteps,βeq)


    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_compartment(l,n,x,γ,β,runin,compartments)
            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J,buc,e_take=update_compartments_stats(J,xi[:,activestate],γ,β)

                append!(energies_take,e_take)
                append!(energies_rand,energie(J[:,:,buc],start_pattern(l,0.5)))
            end

            finalXI_vec=Vec_equiliabration_compartment(J,xi,β,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return sum(sign.(energies_take .-maximum(energies_rand)).-1)./(length(energies_take)*(-2)),sum(sign.(energies_rand .-minimum(energies_take)).+1)/(2*length(energies_rand)),overlaps
end


function get_statistics_One(l,n,x,γ,runin,steps,reps,eqsteps,βeq)


    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_matrix_xi_x_runin(l,n,x,γ,runin)
            finalXI_vec=Vec_equiliabration(J,xi,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                append!(energies_take,energie(J,xi[:,activestate]))
                append!(energies_rand,energie(J,start_pattern(l,0.5)))
                J=new_J(J,xi[:,activestate],γ)


            end

            finalXI_vec=Vec_equiliabration(J,xi,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return sum(sign.(energies_take .-maximum(energies_rand)).-1)./(length(energies_take)*(-2)),sum(sign.(energies_rand .-minimum(energies_take)).+1)/(2*length(energies_rand)),overlaps
end


function get_statistics_One_estats(l,n,x,γ,runin,steps,reps,eqsteps,βeq)


    energies_take=Vector{Float64}()
    energies_rand=Vector{Float64}()
    overlaps=Vector{Int16}()

    for run in 1:reps
            J,xi=build_matrix_xi_x_runin(l,n,x,γ,runin)
            finalXI_vec=Vec_equiliabration(J,xi,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                append!(energies_take,energie(J,xi[:,activestate]))
                append!(energies_rand,energie(J,start_pattern(l,0.5)))
                J=new_J(J,xi[:,activestate],γ)


            end

            finalXI_vec=Vec_equiliabration(J,xi,βeq,eqsteps)
            data=abs.(transpose(xi)*(finalXI_vec))
            for ii in 1:length(data[:,1])
                append!(overlaps, data[ii,ii])
            end
    end
    return sum(sign.(energies_take .-maximum(energies_rand)).-1)./(length(energies_take)*(-2)),sum(sign.(energies_rand .-minimum(energies_take)).+1)/(2*length(energies_rand)),mean(energies_take),std(energies_take),lower_std(energies_take),upper_std(energies_take),overlaps
end



function get_mean_overlaps_One(l,n,x,γ,runin,steps,reps,eqsteps,βeq)


     dumi1,dumi2,dumi3,dumi4,dumi5,dumi6,overlaps=get_statistics_One_estats(l,n,x,γ,runin,steps,reps,eqsteps,βeq)
     return mean(overlaps)/l
end



function get_statistics_One_allstats(l,n,x,γ,runin,steps,reps,eqsteps,βeq)   # produce all statistics for one compartment at a given point



    dist_start=Vector{Int16}()
    dist_next=Vector{Int16}()
    open_paths=Vector{Int16}()

    dist_start_rand=Vector{Int16}()
    dist_next_rand=Vector{Int16}()
    open_paths_rand=Vector{Int16}()

    partision_mean_vecs=Vector{Float64}()
    partision_std_vecs=Vector{Float64}()
    partision_first_vecs=Vector{Float64}()
    partision_last_vecs=Vector{Float64}()

    partision_mean_patts=Vector{Float64}()
    partision_std_patts=Vector{Float64}()
    partision_first_patts=Vector{Float64}()
    partision_last_patts=Vector{Float64}()

    best_overlaps=Vector{Array{Float64,1}}()
    mean_current_direction=Vector{Array{Float64,1}}()
    best_next_direciton=Vector{Array{Float64,1}}()


    for run in 1:reps
            J,xi=build_matrix_xi_x_runin(l,n,x,γ,runin)

            xi[:]=new_XI_rand(xi[:], x)
            finalXI_vec=Vec_equiliabration(J,copy(xi),βeq,eqsteps)

            F=eigen(J)
                 evecs=zeros(l,n)
                 evlues=zeros(n)
                 for jj in 0:n-1
                    evecs[:,jj+1] = F.vectors[:,l-jj]
                    evlues[jj+1]=F.values[l-jj]
                  end

            dist=abs.(transpose(finalXI_vec)*xi)
            for ii in 1:length(dist[:,1])
                append!(dist_start, convert(Int16,dist[ii,ii]))
                 dist[ii,ii]=0
                 append!(dist_next,convert(Int16,maximum(abs.(dist[ii,:]))))
                 append!(open_paths,convert(Int16,sum(-sign.(xi[ :,ii].*(transpose(J[:,:])*xi[:,ii])).+1)/2))
                 d1,d2,d3=produce_overlap_stats(xi[:,ii],evecs,evlues)
                 append!(best_overlaps,[d1])
                 append!(mean_current_direction,[d2])
                 append!(best_next_direciton,[d3])
            end
            energies=Vector{Float64}()
            for ii in 1:n
               append!(energies,energie(J,xi[:,ii]))
            end

            p1,p2,p3,p4,p5,p6,p7,p8=Support_partisipation(xi,evecs,energies)

            append!(partision_mean_vecs,p1)
            append!(partision_std_vecs,p2)
            append!(partision_first_vecs,p3)
            append!(partision_last_vecs,p4)

            append!(partision_mean_patts,p5)
            append!(partision_std_patts,p6)
            append!(partision_first_patts,p7)
            append!(partision_last_patts,p8)




            xi_rand= get_patterns(l,n,0.5)
            finalXI_rand=Vec_equiliabration(J,copy(xi_rand),βeq,eqsteps)

            dist_rand=abs.(transpose(finalXI_rand)*xi_rand)

            dist_rand_pats=abs.(transpose(finalXI_rand)*xi)
            for ii in 1:length(dist[:,1])
                append!(dist_start_rand, convert(Int16,dist_rand[ii,ii]))

                append!(dist_next_rand,convert(Int16,maximum(abs.(dist_rand_pats[ii,:]))))
                append!(open_paths_rand,convert(Int16,sum(-sign.(xi_rand[ :,ii].*(transpose(J[:,:])*xi_rand[:,ii])).+1)/2))

            end
            for ii in 1:steps

                activestate=updateactivestate(n)
                xi[:]=new_XI_rand(xi[:], x)
                J=new_J(J,xi[:,activestate],γ)


            end
            xi[:]=new_XI_rand(xi[:], x)
            finalXI_vec=Vec_equiliabration(J,xi,βeq,eqsteps)



            F=eigen(J)
                 evecs=zeros(l,n)
                 evlues=zeros(n)
                 for jj in 0:n-1
                    evecs[:,jj+1] = F.vectors[:,l-jj]
                    evlues[jj+1]=F.values[l-jj]
                  end

            dist=abs.(transpose(finalXI_vec)*xi)
            for ii in 1:length(dist[:,1])
                append!(dist_start, convert(Int16,dist[ii,ii]))
                 dist[ii,ii]=0
                 append!(dist_next,convert(Int16,maximum(abs.(dist[ii,:]))))
                 append!(open_paths,convert(Int16,sum(-sign.(xi[ :,ii].*(transpose(J[:,:])*xi[:,ii])).+1)/2))
                 d1,d2,d3=produce_overlap_stats(xi[:,ii],evecs,evlues)
                 append!(best_overlaps,[d1])
                 append!(mean_current_direction,[d2])
                 append!(best_next_direciton,[d3])
            end
            energies=Vector{Float64}()
           for ii in 1:n
               append!(energies,energie(J,xi[:,ii]))
           end
            p1,p2,p3,p4,p5,p6,p7,p8=Support_partisipation(xi,evecs,energies)

            append!(partision_mean_vecs,p1)
            append!(partision_std_vecs,p2)
            append!(partision_first_vecs,p3)
            append!(partision_last_vecs,p4)

            append!(partision_mean_patts,p5)
            append!(partision_std_patts,p6)
            append!(partision_first_patts,p7)
            append!(partision_last_patts,p8)







    end
    return dist_start,dist_next,open_paths,     dist_start_rand,dist_next_rand,open_paths_rand,   mean(partision_mean_vecs),sqrt(mean(partision_std_vecs.^2)),mean(partision_first_vecs),mean(partision_last_vecs),     mean(partision_mean_patts),sqrt(mean(partision_std_patts.^2)),mean(partision_first_patts),mean(partision_last_patts),best_overlaps,mean_current_direction,best_next_direciton


end







function produce_overlap_stats(xi,evecs,evlues)
    n=length(evlues)
    l=length(xi)
    α=n/l
    states=convert.(Int64,range(1,stop=n,length=n))
    best_overlaps=Vector{Float64}()
    best_next_direciton=Vector{Float64}()
    mean_current_direction=Vector{Float64}()
    o_squared_curent=0
    lam_current=0

    overlaps=(transpose(xi)*evecs)./(sqrt(l))

    for ii in 1:min(10,n-4)

        positionnow=states[argmax(abs.(overlaps[states]))]

        lam_current=copy(lam_current)+overlaps[positionnow]^2*evlues[positionnow]
        o_squared_curent=copy(o_squared_curent)+copy(overlaps[positionnow])^2

        append!(best_overlaps,o_squared_curent)
        append!(mean_current_direction,lam_current*(o_squared_curent))
        states=setdiff(states,positionnow)
        append!(best_next_direciton,maximum(abs.(transpose(evlues[states]).*evecs[:,states]))*sqrt(abs.(1-o_squared_curent))/sqrt(α)*(sqrt(n/(n-ii))))
    end
    return best_overlaps,mean_current_direction,best_next_direciton
end















function build_rec_for_gamma_One(l,n,x,runin,steps,βeq,equil_steps,γ_steps,γ_min,reps)

    pointsgamtake=convert.(Int16,range(1,stop=20,length=20))
    pointsgamtake2=convert.(Int16,range(1,stop=10,length=10))
        γ=Vector{Float64}()

        if(γ_min< 1 - 0.02)
            append!(γ,range(0.9*γ_min,stop=γ_min,length=21)[pointsgamtake])
            append!(γ,range(γ_min,stop=0.99999,length=convert(Int,γ_steps-20)))

        else

            append!(γ,range(0.9*γ_min,stop=1-0.02,length=21)[pointsgamtake])

            append!(γ,range(1-0.02,stop=min(0.999,1 - (1-γ_min)/4),length=convert(Int,γ_steps-30)))
            append!(γ,range(min(0.999,1 - (1-γ_min)/4),stop=0.99999,length=convert(Int,11))[pointsgamtake2])
        end

    e_rac=zeros(γ_steps)
    e_rand=zeros(γ_steps)
    runin_e=max.(10*n,min.(floor.(Int64,log(0.00001)./ log.(γ) .+1),runin))
    steps_e=min.(runin_e,steps)
    overlaps=zeros(Int16,n*reps*2,γ_steps)
    for ii in 1:γ_steps
        e_rac[ii],e_rand[ii],overlaps[:,ii]=get_statistics_One(l,n,x,γ[ii],runin_e[ii],steps_e[ii],reps,equil_steps,βeq)
    end
    return e_rac,e_rand,overlaps,γ
end


function build_rec_for_gamma_compartment(l,n,x,runin,steps,β,βeq,equil_steps,γ_steps,γ_min,reps,compartments)
    pointsgamtake=convert.(Int16,range(1,stop=20,length=20))
    pointsgamtake2=convert.(Int16,range(1,stop=10,length=5))
        γ=Vector{Float64}()

        if(γ_min< 1 - 0.02)
            append!(γ,range(0.9*γ_min,stop=γ_min,length=21)[pointsgamtake])
            append!(γ,range(γ_min,stop=0.99999,length=convert(Int,γ_steps-20)))

        else

            append!(γ,range(0.9*γ_min,stop=1-0.02,length=21)[pointsgamtake])

            append!(γ,range(1-0.02,stop=min(0.999,1 - (1-γ_min)/4),length=convert(Int,γ_steps-25)))
            append!(γ,range(min(0.999,1 - (1-γ_min)/4),stop=0.99999,length=convert(Int,6))[pointsgamtake2])
        end


    pattern_used=zeros(Int16,steps*reps,γ_steps)
    compartment_used=zeros(Int16,steps*reps,γ_steps)
    runin_e=max.(10*n,min.(floor.(Int64,2*compartments*(log(0.00001)./ log.(γ) .+1)),runin))
    e_rac=zeros(γ_steps)
    e_rand=zeros(γ_steps)
    overlaps=zeros(Int16,n*reps*2,γ_steps)
    for ii in 1:γ_steps
        pattern_used[:,ii],compartment_used[:,ii],e_rac[ii],e_rand[ii],overlaps[:,ii]=get_statistics_compartment(l,n,x,γ[ii],β,runin_e[ii],compartments,steps,reps,equil_steps,βeq)
    end
    return pattern_used,compartment_used,e_rac,e_rand,overlaps,γ
end

function build_rec_for_gamma_compartment_no_use(l,n,x,runin,steps,β,βeq,equil_steps,γ_steps,γ_min,reps,compartments)

    pointsgamtake=convert.(Int16,range(1,stop=20,length=20))
    pointsgamtake2=convert.(Int16,range(2,stop=6,length=5))
        γ=Vector{Float64}()

        if(γ_min< 1 - 0.02)
            append!(γ,range(0.001,stop=γ_min,length=21)[pointsgamtake])
            append!(γ,range(γ_min,stop=0.99999,length=convert(Int,γ_steps-20)))

        else

            append!(γ,range(0.01,stop=1-0.02,length=21)[pointsgamtake])

            append!(γ,range(1-0.02,stop=min(0.999,1 - (1-γ_min)/4),length=convert(Int,γ_steps-25)))
            append!(γ,range(min(0.999,1 - (1-γ_min)/4),stop=0.99999,length=convert(Int,6))[pointsgamtake2])
        end


    runin_e=max.(10*n,min.(floor.(Int64,2*compartments*(log(0.00001)./ log.(γ) .+1)),runin))
    steps_e=min.(runin_e,steps)
    e_rac=zeros(γ_steps)
    e_rand=zeros(γ_steps)
    overlaps=zeros(Int16,n*reps*2,γ_steps)
    for ii in 1:γ_steps
        e_rac[ii],e_rand[ii],overlaps[:,ii]=get_statistics_compartment_no_use(l,n,x,γ[ii],β,runin_e[ii],compartments,steps_e[ii],reps,equil_steps,βeq)
    end
    return e_rac,e_rand,overlaps,γ
end





function build_rec_for_gamma_x_One(l,n,runin,steps,βeq,equil_steps,γ_steps,x_max,xsteps,reps)
    x=Vector{Float64}()
    append!(x,0.0)
    append!(x,exp.(range(log(0.0001),stop=log(x_max),length=xsteps-1)))
    x=x./n


    γ_min=0.01


        γ_min= max.(1-1.999 /n,1 .-3 .*sqrt.(8 .*x./(n -1)))


    γ_min[1]=0.95  # reset for x=0


    e_rac=zeros(γ_steps,xsteps)
    e_rand=zeros(γ_steps,xsteps)
    overlaps=zeros(Int16,n*reps*2,γ_steps,xsteps)
    γ=zeros(γ_steps,xsteps)
    for ii in 1:xsteps
        e_rac[:,ii],e_rand[:,ii],overlaps[:,:,ii],γ[:,ii]=build_rec_for_gamma_One(l,n,x[ii],runin,steps,βeq,equil_steps,γ_steps,γ_min[ii],reps)
    end

    return e_rac,e_rand,overlaps,γ,x
end


function build_rec_for_gamma_x_compartment(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
    x=Vector{Float64}()
    append!(x,0.0)
    append!(x,exp.(range(log(0.0001),stop=log(x_max),length=xsteps-1)))   # set the smallest number of mutation rates
    x=x./n


    γ_min=0.01 .*ones(xsteps)

    if(compartments<n)
        γ_min= max.(1-1.999 * compartments/n,1 .-3 .*sqrt.(8 .*x./(n/compartments -1)))
    end

    γ_min[1]=0.95  # reset for x=0



    pattern_used=zeros(Int16,steps*reps,γ_steps,xsteps)
    compartment_used=zeros(Int16,steps*reps,γ_steps,xsteps)
    e_rac=zeros(γ_steps,xsteps)
    e_rand=zeros(γ_steps,xsteps)
    overlaps=zeros(Int16,n*reps*2,γ_steps,xsteps)
        γ=zeros(γ_steps,xsteps)
    for ii in 1:xsteps
        e_rac[:,ii],e_rand[:,ii],overlaps[:,:,ii],γ[:,ii]=build_rec_for_gamma_compartment_no_use(l,n,x[ii],runin,steps,β,βeq,equil_steps,γ_steps,γ_min[ii],reps,compartments)
    end

    return e_rac,e_rand,overlaps,γ,x
end








function build_rec_for_opt_gamma_One(l,n,x,γ,runin,steps,βeq,equil_steps,reps)  # give the results and statistics along an (optimal) learning rate
    γ_steps=length(γ)
    e_rac=zeros(γ_steps)
    e_rand=zeros(γ_steps)
    e_mean=zeros(γ_steps)
    e_std_low=zeros(γ_steps)
    e_std_high=zeros(γ_steps)
    overlaps=zeros(Int16,n*reps*2,γ_steps)
    runin_e=max.(10*n,min.(floor.(Int64,log(0.00001)./ log.(γ) .+1),runin))
    steps_e=min.(runin_e,steps)
    for ii in 1:γ_steps
        e_rac[ii],e_rand[ii],e_mean[ii],stdnormal,e_std_low[ii],e_std_high[ii],overlaps[:,ii]=get_statistics_One_estats(l,n,x[ii],γ[ii],runin_e[ii],steps_e[ii],reps,equil_steps,βeq)
    end
    return e_rac,e_rand,e_mean,e_std_low,e_std_high,overlaps
end





function build_rec_for_opt_gamma_compartment(l,n,x,γ,runin,steps,β,βeq,equil_steps,reps,compartments)  # give the results and statistics along an (optimal) learning rate  for multiple compartments
    γ_steps=length(γ)
    pattern_used=zeros(Int16,steps*reps,γ_steps)
    compartment_used=zeros(Int16,steps*reps,γ_steps)

    e_rac=zeros(γ_steps)
    e_rand=zeros(γ_steps)
    overlaps=zeros(Int16,n*reps*2,γ_steps)
    for ii in 1:γ_steps
        pattern_used[:,ii],compartment_used[:,ii],e_rac[ii],e_rand[ii],overlaps[:,ii]=get_statistics_compartment(l,n,x[ii],γ[ii],β,runin,compartments,steps,reps,equil_steps,βeq)
    end
    return pattern_used,compartment_used,e_rac,e_rand,overlaps
end








function distribution_only_allEV(l,n,x,γ,runin,reps)


    evalue_stats=Vector{Array{Float64,1}}()
        for run in 1:reps

                J,xi=build_matrix_xi_x_runin(l,n,x,γ,runin)

                 F=eigen(J)
                 evlues=zeros(l)
                 for jj in 0:l-1
                    evlues[jj+1]=F.values[l-jj]
                  end
                append!(evalue_stats,[evlues])
            end
    evalue_means=zeros(l)

    for ii in 1:l
        evslocal=Vector{Float64}()
        for run in 1:reps
            append!(evslocal,evalue_stats[run][ii])
        end
        evalue_means[ii]=mean(evslocal)
    end


    return evalue_stats,evalue_means
end




function distribution_EV(l,n,x,γ,runin,reps)





    evalue_stats=Vector{Array{Float64,1}}()
    entropies=Vector{Tuple{Array{Float64,1},Float64,Float64,Array{Float64,1},Float64,Float64,Array{Float64,1},Float64,Float64}}()
    for run in 1:reps

                J,xi=build_matrix_xi_x_runin(l,n,x,γ,runin)

                 F=eigen(J)
                 evecs=zeros(l,n)
                 evlues=zeros(n)
                 for jj in 0:n-1
                    evecs[:,jj+1] = F.vectors[:,l-jj]
                    evlues[jj+1]=F.values[l-jj]
                  end
                append!(evalue_stats,[evlues])



                energies=Vector{Float64}()
               for ii in 1:n
                   append!(energies,energie(J,xi[:,ii]))
               end
                local_entropy=Support_entropy(xi,evecs,energies)

                append!(entropies,[local_entropy])

    end


    return evalue_stats,entropies
end


function distribution_EV_for_x(l,n,xsteps,xmin,xmax,runin_max,reps)
    x=exp.(range(log(xmin),stop=log(0.01),length=4))
    if (xsteps>5)
        append!(x,range(0.015,stop=xmax,length=xsteps-4))
    end
    γ=1 .- sqrt.(8 .*x./(n.-1))
    runin=min.(runin_max,floor.(Int64,log(0.00001)./ log.(γ) .+1))
    data=distribution_EV.(l,n,0,γ,runin,reps)
    return data
end
