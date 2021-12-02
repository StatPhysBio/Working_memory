

include("Self_recognition-Compartments.jl")






function simulate_self_recognition(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)


    println("test1")
    filename = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_lowx.jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else
        if !isdir("../Data")
            mkdir("../Data")
        end

        if !isdir("../Data/Workingmemory_compartment")
            mkdir("../Data/Workingmemory_compartment")
        end

        if !isdir("../Data/Workingmemory_compartment/L_$(l)")
            mkdir("../Data/Workingmemory_compartment/L_$(l)")
        end
        if !isdir("../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)")
            mkdir("../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)")
        end
        println("I'm starting the simulation")
        β_scaled=10^(β/10)/(1 /(l/2/n*compartments -sqrt(l/π/2))*(compartments*l)^2 /n/n)
        β_eq_scaled=10^(βeq/10)*n/compartments
        if (compartments==1)
            println("only one Compartment")
            e_rac,e_rand,overlapsload,γ,x = build_rec_for_gamma_x_One(l,n,runin,steps,β_eq_scaled,equil_steps,γ_steps,x_max,xsteps,reps)

            overlaps=zeros(γ_steps,xsteps)
            notrecognized=zeros(γ_steps,xsteps)
            optG=zeros(xsteps)
            min_gamma=zeros(xsteps)
            max_gamma=zeros(xsteps)
            for xx in 1:xsteps
                for gg in 1:γ_steps
                    overlaps[gg,xx] =mean(overlapsload[:,gg,xx])
                    for attempt in 1:length(overlapsload[:,gg,xx])
                            if (overlapsload[attempt,gg,xx]/l <0.8)
                                notrecognized[gg,xx]+=1
                            end
                    end
                end

                errorregime=sqrt(max(1,2*notrecognized[argmax(overlaps[:,xx]),xx]))
                threshEqual= 1/(2*n*reps)*errorregime
                optG[xx],min_gamma[xx],max_gamma[xx],dummi= opt_lambda_analysis(overlaps[:,xx],γ[:,xx],threshEqual)
            end


            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename, "overlaps",overlaps, "gamma",γ,"x",x,"parameters",parameters,"gammaopt",[optG,min_gamma,max_gamma],"notrecognized",notrecognized) ## Here we don't save the energyies to minimize the data size
            #save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",γ,"x",x,"parameters",parameters)

        else
            if (compartments > n)
                println("to many compartments" )
                return -1
            end
            println("you used $(compartments) compartments for $(n) patterns" )
            e_rac,e_rand,overlapsload,γ,x = build_rec_for_gamma_x_compartment(l,n,runin,steps,β_scaled,β_eq_scaled,equil_steps,γ_steps,compartments,x_max,xsteps,reps)

            overlaps=zeros(γ_steps,xsteps)
            notrecognized=zeros(γ_steps,xsteps)
            optG=zeros(xsteps)
            min_gamma=zeros(xsteps)
            max_gamma=zeros(xsteps)
            for xx in 1:xsteps
                for gg in 1:γ_steps
                    overlaps[gg,xx] =mean(overlapsload[:,gg,xx])
                    for attempt in 1:length(overlapsload[:,gg,xx])
                            if (overlapsload[attempt,gg,xx]/l <0.8)
                                notrecognized[gg,xx]+=1
                            end
                    end
                end

                errorregime=sqrt(max(1,2*notrecognized[argmax(overlaps[:,xx]),xx]))
                threshEqual= 1/(2*n*reps)*errorregime
                optG[xx],min_gamma[xx],max_gamma[xx],dummi= opt_lambda_analysis(overlaps[:,xx],γ[:,xx],threshEqual)
            end

            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename, "overlaps",overlaps, "gamma",γ,"x",x,"parameters",parameters,"gammaopt",[optG,min_gamma,max_gamma],"notrecognized",notrecognized) ## Here we don't save the energyies to minimize the data size
            #save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",γ,"x",x,"parameters",parameters)

        end

        println("I'm saving the file $(filename)")
    end
end







function simulate_self_recognition_alongopt(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)


    filename = "../Data/Workingmemory_compartment_opt/L_$(l)/reps_$(reps)/Optimal-compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_lowx.jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_opt")
            mkdir("../Data/Workingmemory_compartment_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_opt/L_$(l)")
            mkdir("../Data/Workingmemory_compartment_opt/L_$(l)")
        end
        if !isdir("../Data/Workingmemory_compartment_opt/L_$(l)/reps_$(reps)")
            mkdir("../Data/Workingmemory_compartment_opt/L_$(l)/reps_$(reps)")
        end
        println("Ich starte mit der siumulation")

        filenameload = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_lowx.jld2"

        β_scaled=10^(β/10)/(1 /(l/2/n*compartments -sqrt(l/π/2))*(compartments*l)^2 /n/n)
        β_eq_scaled=10^(βeq/10)*n/compartments
        if (compartments==1)
            println("only one compartment")
            optg1,dumi1,dumi2=load(filenameload,"gammaopt")

            x=load(filenameload,"x")



            e_rac,e_rand,meane,std_low,std_high,overlaps = build_rec_for_opt_gamma_One(l,n,x,optg1,runin,steps,β_eq_scaled,equil_steps*20,reps)
            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",optg1,"x",x,"meanE",meane,"std_low",std_low,"std_high",std_high,"parameters",parameters)

        else
            if (compartments > n)
                println("to many compartments" )
                return -1
            end
            overlaps=load(filenameload,"overlaps")

            x=load(filenameload,"x")


            optg1,dumi1,dumi2=load(filenameload,"gammaopt")
            pattern_used,compartment_used,e_rac,e_rand,overlaps = build_rec_for_opt_gamma_compartment(l,n,x,optg1,runin,steps,β_scaled,β_eq_scaled,equil_steps*20,reps,compartments)
            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",optg1,"x",x,"pattern_used",pattern_used,"compartment_used",compartment_used,"parameters",parameters)

        end

        println("I'm saving the file $(filename)")
    end
end













function simulate_self_recognition_alongopt_no_evolution(l,n,runin,steps,β,βeq,equil_steps,γ_steps,γ_min,compartments,x_max,xsteps,reps)


    println("test1")
    #filename = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-eqtemp-$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n).jld2"
    filename = "../Data/Workingmemory_compartment_opt/L_$(l)/reps_$(reps)/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_small_no_evolution.jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_opt")
            mkdir("../Data/Workingmemory_compartment_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_opt/L_$(l)")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/L_$(l)")
        end
        if !isdir("../Data/Workingmemory_compartment_pt/L_$(l)/reps_$(reps)")
            mkdir("../Data/Workingmemory_compartment_/L_$(l)/reps_$(reps)")
        end
        println("Ich starte mit der siumulation")

        filenameload = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_lowx.jld2"

        β_scaled=10^(β/10)/(1 /(l/2/n*compartments -sqrt(l/π/2))*(compartments*l)^2 /n/n)
        β_eq_scaled=10^(βeq/10)*n/compartments
        if (compartments==1)
            println("only one compartment")
            optg1,dumi1,dumi2=load(filenameload,"gammaopt")

            x=load(filenameload,"x")


            e_rac,e_rand,meane,std_low,std_high,overlaps = build_rec_for_opt_gamma_One(l,n,0,optg1,runin,steps,β_eq_scaled,equil_steps*20,reps)
            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",optg1,"x",x,"meanE",meane,"std_low",std_low,"std_high",std_high,"parameters",parameters)

        else
            if (compartments > n)
                println("to many compartments" )
                return -1
            end
            overlaps=load(filenameload,"overlaps")

            x=load(filenameload,"x")


            optg1,dumi1,dumi2=load(filenameload,"gammaopt")
            pattern_used,compartment_used,e_rac,e_rand,overlaps = build_rec_for_opt_gamma_compartment(l,n,0,optg1,runin,steps,β_scaled,β_eq_scaled,equil_steps*20,reps,compartments)
            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,compartments,x_max,xsteps,reps)
            save(filename,"e_rac",e_rac, "e_rand",e_rand, "overlaps",overlaps, "gamma",optg1,"x",x,"pattern_used",pattern_used,"compartment_used",compartment_used,"parameters",parameters)

        end

        println("I'm saving the file $(filename)")
    end
end









function simulate_self_recognition_singlex_alongopts_stats(l,n,runin,steps,β,βeq,equil_steps,γ_steps,γ_min,compartments,x_max,xsteps,reps,this_x)


    println("test1")
        filename = "../Data/Workingmemory_compartment_singleX_opt/L_$(l)/reps_$(reps)/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_thisX-$(this_x)_small_stats.jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_singleX_opt")
            mkdir("../Data/Workingmemory_compartment_singleX_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_singleX_opt/L_$(l)")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/L_$(l)")
        end
        if !isdir("../Data/Workingmemory_compartment_singleX_opt/L_$(l)/reps_$(reps)")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/L_$(l)/reps_$(reps)")
        end
        println("Ich starte mit der siumulation")

        filenameload = "../Data/Workingmemory_compartment_singleX/L_$(l)/reps_50/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_thisX-$(this_x)_small.jld2"
        println("Ich lade datei:")
        println(filenameload)
        β_scaled=10^(β/10)/(50)
        β_eq_scaled=10^(βeq/10)*n/compartments

        if (compartments==1)
            println("only one compartment")
            overlaps=load(filenameload,"overlaps")
            γ=load(filenameload,"gamma")
            x=load(filenameload,"x")


            #m,s,opt,optg1,opt2,optg2=give_stats_Overlaps(overlaps./l,γ);
            optg1=γ[argmax(overlaps)]
            runin_e=max.(10*n,min.(floor.(Int64,log(0.00001)./ log.(optg1) .+1),runin))
            steps_e=min.(runin_e,steps)

            dist_start,dist_next,open_paths,     dist_start_rand,dist_next_rand,open_paths_rand,   partision_mean_vecs,partision_std_vecs,partision_first_vecs,partision_last_vecs,     partision_mean_patts,partision_std_patts,partision_first_patts,partision_last_patts,best_overlaps,mean_current_direction,best_next_direciton = get_statistics_One_allstats(l,n,x[this_x],optg1[this_x],runin_e,steps_e,reps,equil_steps*20,β_eq_scaled)
            parameters=(l,n,runin,steps,β,βeq,equil_steps,γ_steps,γ_min,compartments,x_max,xsteps,reps,this_x)
            save(filename,"dist_start",dist_start, "dist_next",dist_next,"open_paths" ,open_paths,   "dist_start_rand",dist_start_rand,"dist_next_rand",dist_next_rand,"open_paths_rand",open_paths_rand,   "partision_mean_vecs",partision_mean_vecs,"partision_std_vecs",partision_std_vecs,"partision_first_vecs",partision_first_vecs,"partision_last_vecs",partision_last_vecs,     "partision_mean_patts",partision_mean_patts,"partision_std_patts",partision_std_patts,"partision_first_patts",partision_first_patts,"partision_last_patts",partision_last_patts,"best_overlaps",best_overlaps,"mean_current_direction",mean_current_direction,"best_next_direciton",best_next_direciton      ,"gamma",optg1,"x",x,"parameters",parameters)

        else
            println("sorry this only works for one compartment at this point try again another time")
        end

        println("I'm saving the file $(filename)")
    end
end




function Make_eigenvalue_eigenvectors_stats_along_x(l,n,runin,steps,equil_steps,β,βeq,compartments,x_max,xsteps,reps,this_x)

    println("test1")
    filename = "../Data/Workingmemory_compartment_singleX_opt/EV_stats/ev_stats_L-$(l)_N-$(n)_xsteps$(xsteps)_reps=$(reps)_thisX_$(this_x).jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_singleX_opt")
            mkdir("../Data/Workingmemory_compartment_singleX_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
        end


            println("Ich starte mit der siumulation")

            filenameload = "../Data/Workingmemory_compartment_singleX/L_$(l)/reps_50/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_thisX-$(this_x)_small.jld2"
            println("Ich lade datei:")
            println(filenameload)
            β_scaled=10^(β/10)/(1 /(l/2/n*compartments -sqrt(l/π/2))*(compartments*l)^2 /n/n)
            β_eq_scaled=10^(βeq/10)*n/compartments

            if (compartments==1)
                println("only one compartment")
                overlaps=load(filenameload,"overlaps")
                γ=load(filenameload,"gamma")
                x=load(filenameload,"x")


                #m,s,opt,optg1,opt2,optg2=give_stats_Overlaps(overlaps./l,γ);
                optg1=γ[argmax(overlaps)]
                runin_e=max.(10*n,min.(floor.(Int64,log(0.00001)./ log.(optg1) .+1),runin))

                data =    distribution_only_allEV(l,n,x[this_x],optg1,runin_e,reps)
                info=" data[1][rep][ii] is the iits eigenvalue of the data with ii in 1:n;  data[2][ii] is the mean of the evs"
                save(filename,"../Data",data,"gamma",optg1, "x",x,"this_x",this_x, "info",info)

            else
                println("sorry this only works for one compartment at this point try again a nother time")
            end

            println("I'm saving the file $(filename)")

    end



    println("Im now doing the non evolving one")
    #filename = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-eqtemp-$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n).jld2"
    filename = "../Data/Workingmemory_compartment_singleX_opt/EV_stats/ev_stats_L-$(l)_N-$(n)_xsteps$(xsteps)_reps=$(reps)_thisX_$(this_x)_noevo.jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_singleX_opt")
            mkdir("../Data/Workingmemory_compartment_singleX_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
        end


            println("Ich starte mit der siumulation")

            filenameload = "../Data/Workingmemory_compartment_singleX/L_$(l)/reps_50/Scaled_temp_compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-elqutemp_$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n)_thisX-$(this_x)_small.jld2"
            println("Ich lade datei:")
            println(filenameload)
            β_scaled=10^(β/10)/(1 /(l/2/n*compartments -sqrt(l/π/2))*(compartments*l)^2 /n/n)
            β_eq_scaled=10^(βeq/10)*n/compartments

            if (compartments==1)
                println("only one compartment")
                overlaps=load(filenameload,"overlaps")
                γ=load(filenameload,"gamma")
                x=load(filenameload,"x")


                #m,s,opt,optg1,opt2,optg2=give_stats_Overlaps(overlaps./l,γ);
                optg1=γ[argmax(overlaps)]
                runin_e=max.(10*n,min.(floor.(Int64,log(0.00001)./ log.(optg1) .+1),runin))

                data =    distribution_only_allEV(l,n,0,optg1,runin_e,reps)
                info=" data[1][rep][ii] is the iits eigenvalue of the data with ii in 1:n;  data[2][ii] is the mean of the evs"
                save(filename,"../Data",data,"gamma",optg1, "x",x,"this_x",this_x, "info",info)

            else
                println("sorry this only works for one compartment at this point try again a nother time")
            end

            println("I'm saving the file $(filename)")

    end

end









function Make_eigenvalue_eigenvectors_stats(l,n,xsteps,xmin,xmax,runin_max,reps)

    println("test1")
    #filename = "../Data/Workingmemory_compartment/L_$(l)/reps_$(reps)/compartments_$(compartments)-SelfRec_Xsteps_$(xsteps)-temp_$(β)-eqtemp-$(βeq)_Xmax_$(x_max)_runin_$(runin)_eqsteps_$(equil_steps)_Patterns_$(n).jld2"
    filename = "../Data/Workingmemory_compartment_singleX_opt/EV_stats/ev_stats_L-$(l)_N-$(n)_xsteps$(xsteps)_reps=$(reps).jld2"

    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("../Data/Workingmemory_compartment_singleX_opt")
            mkdir("../Data/Workingmemory_compartment_singleX_opt")
        end

        if !isdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
            mkdir("../Data/Workingmemory_compartment_singleX_opt/EV_stats")
        end

        println("Ich starte mit der siumulation")


        data =distribution_EV_for_x(l,n,xsteps,xmin,xmax,runin_max,reps)
        parameters=(l,n,xsteps,xmin,xmax,runin_max,reps)
        x=exp.(range(log(xmin),stop=log(0.01),length=4))
        if (xsteps>5)
            append!(x,range(0.015,stop=xmax,length=xsteps-4))
        end
        γ=1 .- sqrt.(8 .*x./(n.-1))
        info="../Data[g] is the data for gamma[g]; data[g][1][rep][ii] is the iits eigenvalue of the data with ii in 1:n;data[g][2][rep][ee] is the entropy stats of the data with ee= 1,2,3 entorpyevec(mean,1 ,last) , 456= same but normed 789 fror eneergies, where the mean,come with mean,std,and lower-uper std "
        save(filename,"../Data",data,"gamma",γ, "parameters",parameters,"info",info)



        println("I'm saving the file $(filename)")
    end
end
