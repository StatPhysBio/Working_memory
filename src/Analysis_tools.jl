using Glob



function get_files(direct,filenames)
    filenames=glob(filenames,direct)
    return filenames
end


function upper_std(data)
    m=mean(data)
    return sqrt(sum(((data.-m).*(0.5 .+ sign.(data.-m)./2)).^2 ./sum(0.5 .+ sign.(data.-m)./2)))
end

function lower_std(data)
    m=mean(data)
    return sqrt(sum(((data.-m).*(0.5 .- sign.(data.-m)./2)).^2 ./sum(0.5 .- sign.(data.-m)./2)))
end


function opt_lambda_analysis(overlaps,γ,threshEqual)
    optg0=argmax(overlaps)
    optper=copy(overlaps[optg0])



    currentp=copy(optper)
    threshold=currentp-threshEqual
    counter=copy(optg0)
    couterinhigh=copy(optg0)
   while ( counter < length(γ) )
       counter+=1
        currentp=copy(overlaps[counter])
        if (currentp >= threshold)
            couterinhigh=copy(counter)
        end

    end

    counter=copy(optg0)
    couterinlow=copy(optg0)
   while ( counter > 1 )
       counter-=1
        currentp=copy(overlaps[counter])
        if (currentp >= threshold)
            couterinlow=copy(counter)
        end

    end

    optG=(γ[couterinlow]+γ[couterinhigh])/2

    optg1=convert(Int,floor((couterinlow+couterinhigh)/2))




    return optG,γ[max(1,min(couterinlow,optg1-1))],γ[min(length(γ),max(couterinhigh,optg1+1))],optg1

end



function normedlog(x) # to give limit 0 * log(0 ) = 0
    if x>0
        return log2(x)
    else
        return 0
    end
end

function entropymatrix(counts) #calculation of the entropy of patterns in the compartments
    entropy=Vector{Float64}()
    norms=Vector{Float64}()
    for ii in 1:length(counts[1,:])
        norm=sum(counts[:,ii]).+0.0000000001  # takeing care of zero norms

        p=counts[:,ii]./norm

        append!(entropy,-sum(p .*normedlog.(p)))
        append!(norms,norm)
    end
    return sum(entropy .*norms ./sum(norms))
end



function entropy_compartments_vs_patterns(compartments,patterns)
    bs=convert(Int64,maximum(compartments))
    ps=convert(Int64,maximum(patterns))
    probs=zeros(bs,bs)
    puses=zeros(ps)
    entropies=zeros(ps)
    for patternanalyse in 1:ps
            Compartment_anlyse=Vector{Int16}()
            for ii in 1:length(patterns)
                if (patterns[ii]==patternanalyse)
                    append!(Compartment_anlyse,compartments[ii])
                    puses[patternanalyse]+=1
                end
            end
            for ii in 1:length(Compartment_anlyse)-1
                probs[Compartment_anlyse[ii],Compartment_anlyse[ii+1]]+=1
            end

           entropies[patternanalyse]=entropymatrix(probs)
    end
    return sum(puses.*entropies./sum(puses))
end


function round_digest(data,diget)
    return round(data*10^diget)/10^diget
end
