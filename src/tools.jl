function improved(b_new, b_old)
    @show b_new
    println("optimized improved:", mean(b_old.times)/mean(b_new.times), " times")
end
