    custom_colors = [
    "#1f77b4",  # niebieski (klasyczny, spokojny)
    "#ff7f0e",  # pomarańczowy (ciepły, kontrast)
    "#2ca02c",  # zielony (czytelny i łagodny)
    "#d62728",  # czerwony (intensywny)
    "#9467bd",  # fioletowy (kontrastowy)
    "#8c564b",  # brązowawy (dojrzały)
    "#e377c2",  # różowy (nietypowy, ale widoczny)
    "#7f7f7f",  # szary (neutralny)
    "#bcbd22",  # oliwkowy
    "#17becf"   # błękitno-morski
    ]

function fig_1(h, c_list, custom_colors)
    N, M = size(c_list[1])
    n_variants = length(c_list)
    layout = (2, 2)

    y_min = minimum([minimum(c) for c in c_list])
    y_max = maximum([maximum(c) for c in c_list])

    plt = plot(layout=layout, size=(1000, 800), legend=false)

    for i in 1:n_variants
        c = c_list[i]
        for j in 1:M
            plot!(plt[i], h.W, c[:, j],
            color = custom_colors[mod1(j, length(custom_colors))],
            label = false,
            linewidth = 2)
        end

        plot!(plt[i],
            xlabel = "Wealth (W)",
            ylabel = "Consumption (c)",
            title = "Panel $i",
            grid = true,
            framestyle = :box,
            xlims = (minimum(h.W), maximum(h.W)),
            ylims = (y_min, y_max)
        )
    end

    display(plt)
    return plt
end

function fig_2(h, φ_list, custom_colors)
    N, M = size(φ_list[1])
    layout = (1, 2)

    plt = plot(layout=layout, size=(1000, 400), legend=false)

    panel_titles = ["Panel 1: Homothetic preferences", "Panel 2: Non-homothetic preferences"]

    for i in 1:2
        φ = φ_list[i]
        for j in 1:M
            plot!(plt[i], h.W, φ[:, j],
                color = custom_colors[mod1(j, length(custom_colors))],
                label = false,
                linewidth = 2)
        end

        plot!(plt[i],
            xlabel = "Wealth (W)",
            ylabel = "Share of risky asset (φ)",
            title = panel_titles[i],
            grid = true,
            framestyle = :box,
            xlims = (0.0, N/2),
            ylims = (0.0, 1.0),
            left_margin = 8mm,
            bottom_margin = 8mm
        )
    end

    display(plt)
    return plt
end

function plot_mpc_and_elasticity(W, mpc, elasticity; z_labels=nothing)
    N, M = size(mpc)
    W_max = W[end]

    if z_labels === nothing
        z_labels = ["z = $j" for j in 1:M]
    end

    plt1 = plot(title="Marginal Propensity to Consume (MPC)",
                xlabel="Wealth", ylabel="MPC", legend=:topright,
                linewidth=2, grid=true,
                xlims = (-20, W_max/2),
                left_margin = 8mm,
                bottom_margin = 8mm)

    for j in 1:M
        plot!(plt1, W, mpc[:, j], 
            color = custom_colors[mod1(j, length(custom_colors))],
            label=z_labels[j])
    end

    plt2 = plot(title="Wealth Elasticity of Consumption",
                xlabel="Wealth", ylabel="Elasticity", legend=:topright,
                linewidth=2, grid=true,
                xlims = (-20, W_max/2),
                left_margin = 8mm,
                bottom_margin = 8mm)

    for j in 1:M
        plot!(plt2, W, elasticity[:, j], 
            color = custom_colors[mod1(j, length(custom_colors))],
            label=z_labels[j])
    end

    plt = plot(plt1, plt2, layout=(1, 2), size=(1000, 400))

    display(plt)
    return plt
end

function plot_mobility_results(results_50, results_150, custom_colors)

    labels_order = [:Median, :Mean, :P90, :P99, :P99_9]
    x_labels = ["Median", "Mean", "P90", "P99", "P99.9"]
    x_pos = 1:1.2:5.8
    bar_width = 0.25
    shifts = [-0.375, -0.125, 0.125, 0.375]

    function extract_series(results)
        p_top01 = [results[l].p_top01 for l in labels_order]
        p_top1  = [results[l].p_top1  for l in labels_order]
        p_top10 = [results[l].p_top10 for l in labels_order]
        p_top50 = [results[l].p_top50 for l in labels_order]
        return p_top01, p_top1, p_top10, p_top50
    end

    p01_50, p1_50, p10_50, p50_50 = extract_series(results_50)
    p01_150, p1_150, p10_150, p50_150 = extract_series(results_150)

    plt = plot(layout = (2,1), size=(700, 600), legend = :topleft,
               titlefontsize=12, guidefontsize=11, tickfontsize=10)

    plot!(plt[1], xlabel = "Initial wealth level", ylabel = "Probability after 50 years",
          xticks = (x_pos, x_labels), title = "Mobility after 50 years", ylim = (0, 1))
    bar!(plt[1], x_pos .+ shifts[1], p01_50, bar_width=bar_width, label="Top 0.1%", color=custom_colors[1])
    bar!(plt[1], x_pos .+ shifts[2], p1_50,  bar_width=bar_width, label="Top 1%",   color=custom_colors[2])
    bar!(plt[1], x_pos .+ shifts[3], p10_50, bar_width=bar_width, label="Top 10%",  color=custom_colors[3])
    bar!(plt[1], x_pos .+ shifts[4], p50_50, bar_width=bar_width, label="Top 50%",  color=custom_colors[4])

    plot!(plt[2], xlabel = "Initial wealth level", ylabel = "Probability after 150 years",
          xticks = (x_pos, x_labels), title = "Mobility after 150 years", ylim = (0, 1))
    bar!(plt[2], x_pos .+ shifts[1], p01_150, bar_width=bar_width, label="Top 0.1%", color=custom_colors[1])
    bar!(plt[2], x_pos .+ shifts[2], p1_150,  bar_width=bar_width, label="Top 1%",   color=custom_colors[2])
    bar!(plt[2], x_pos .+ shifts[3], p10_150, bar_width=bar_width, label="Top 10%",  color=custom_colors[3])
    bar!(plt[2], x_pos .+ shifts[4], p50_150, bar_width=bar_width, label="Top 50%",  color=custom_colors[4])

    display(plt)
    return plt
end





