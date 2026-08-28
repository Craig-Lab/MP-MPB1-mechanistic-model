
# Plot function for all models of tumour growth
function plotall_growth_only(data, sol, error, title, color; seelegend = true)

    time = data[1]; tumour_size = data[2];
    
    trace_data = scatter(
        x = time,
        y = tumour_size,
        mode = "markers",
        name = "Data",
        marker = attr(size = 18, color = "black"),
        error_y = attr(
            type = "data",
            array = error[1],
            arrayminus = error[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    exp_tumour_size = sol[2]; log_tumour_size = sol[3]; gom_tumour_size = sol[4];

    trace_exp = scatter(
        x = sim_time,
        y = exp_tumour_size,
        mode = "lines",
        name = "Exponential",
        line = attr(width = 6, color = color[1])
    )

    trace_log = scatter(
        x = sim_time,
        y = log_tumour_size,
        mode = "lines",
        name = "Logistic",
        line = attr(width = 6, color = color[2])
    )

    trace_gom = scatter(
        x = sim_time,
        y = gom_tumour_size,
        mode = "lines",
        name = "Gompertz",
        line = attr(width = 6, color = color[3])
    )

    if seelegend == true
        yatt = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = attr(
            title = "Time (days)",
            titlefont = attr(size = 28), 
            tickvals = [0, 2, 4, 6, 8, 10], 
            range = [-0.5, 10.5],
            tickfont = attr(size = 28)
        ),
        yaxis = yatt,
        width = 800,
        height = 600,
        legend = attr(
            x=0.15, y=1, 
            xanchor="center", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    fig = plot([trace_exp, trace_log, trace_gom, trace_data], layout);

    display(fig)
    savefig(fig, joinpath(@__DIR__, "Figures", title * ".svg"); width = 800, height = 600)
end

## -------------------------------------

# Plot function for logistic growth with all treatments
function plot_log_treatments(data, sol, error, title, color; seelegend = true)

    time = data[1]; 

    tumour_veh = data[2];
    tumour_cis = data[3];
    tumour_icb = data[4];
    tumour_cis_icb = data[5];

    error_veh = error[1];
    error_cis = error[2];
    error_icb = error[3];
    error_cis_icb = error[4];
    
    veh_data = scatter(
        x = time,
        y = tumour_veh,
        mode = "markers",
        name = "Vehicle Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = error_veh[1],
            arrayminus = error_veh[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    cis_data = scatter(
        x = time,
        y = tumour_cis,
        mode = "markers",
        name = "Cis Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = error_cis[1],
            arrayminus = error_cis[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    icb_data = scatter(
        x = time,
        y = tumour_icb,
        mode = "markers",
        name = "ICB Data",
        showlegend = false,
        marker = attr(size = 18, color = color[3]),
        error_y = attr(
            type = "data",
            array = error_icb[1],
            arrayminus = error_icb[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    cis_icb_data = scatter(
        x = time,
        y = tumour_cis_icb,
        mode = "markers",
        name = "Cis + ICB Data",
        showlegend = false,
        marker = attr(size = 18, color = color[4]),
        error_y = attr(
            type = "data",
            array = error_cis_icb[1],
            arrayminus = error_cis_icb[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour_veh = sol[2]; sim_tumour_cis = sol[3]; 
    sim_tumour_icb = sol[4]; sim_tumour_cis_icb = sol[5];

    trace_veh = scatter(
        x = sim_time,
        y = sim_tumour_veh,
        mode = "lines",
        name = "Vehicle",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    trace_cis = scatter(
        x = sim_time,
        y = sim_tumour_cis,
        mode = "lines",
        name = "Cis",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    trace_icb = scatter(
        x = sim_time,
        y = sim_tumour_icb,
        mode = "lines",
        name = "ICB",
        showlegend = false,
        line = attr(width = 6, color = color[3])
    )

    trace_cis_icb = scatter(
        x = sim_time,
        y = sim_tumour_cis_icb,
        mode = "lines",
        name = "Cis + ICB",
        showlegend = false,
        line = attr(width = 6, color = color[4])
    )

    if seelegend == true

        # Fake treatment legends
        legend_traces = 
            [
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Vehicle", marker = attr(color = vehicle_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin", marker = attr(color = cis_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "ICB", marker = attr(color = icb_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin + ICB", marker = attr(color = cis_icb_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Data",legendrank = 1, marker = attr(color = "black", size = 18)),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Model",legendrank = 2, line = attr(color = "black", width = 6))
            ];

        yatt = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = attr(
            title = "Time (days)",
            titlefont = attr(size = 28), 
            tickvals = [0, 2, 4, 6, 8, 10], 
            range = [-0.5, 10.5],
            tickfont = attr(size = 28)
        ),
        yaxis = yatt,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation = "v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    if seelegend == true
        fig = plot([veh_data, cis_data, icb_data, cis_icb_data, trace_veh, trace_cis, trace_icb, trace_cis_icb, legend_traces...], layout);
    else
        fig = plot([veh_data, cis_data, icb_data, cis_icb_data, trace_veh, trace_cis, trace_icb, trace_cis_icb], layout);
    end

    display(fig)
    savefig(fig, joinpath(@__DIR__, "Figures", title * ".svg"); width = 800, height = 600)
end

## -------------------------------------

# Function plot cisplatin treatment PKPD vs logistic model
function plot_pkpd_cis(data, sol, error, title, color; seelegend = true)

    time = data[1]; 

    tumour_veh = data[2];
    tumour_pkpd = data[3];

    error_veh = error[1];
    error_pkpd = error[2];
    
    veh_data = scatter(
        x = time,
        y = tumour_veh,
        mode = "markers",
        name = "Vehicle Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = error_veh[1],
            arrayminus = error_veh[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    pkpd_data = scatter(
        x = time,
        y = tumour_pkpd,
        mode = "markers",
        name = "Cis Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = error_pkpd[1],
            arrayminus = error_pkpd[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour_veh = sol[2]; sim_tumour_pkpd = sol[3]; 

    trace_veh = scatter(
        x = sim_time,
        y = sim_tumour_veh,
        mode = "lines",
        name = "Vehicle",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    trace_pkpd = scatter(
        x = sim_time,
        y = sim_tumour_pkpd,
        mode = "lines",
        name = "PKPD",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    if seelegend == true

        # Fake treatment legends
        legend_traces = 
            [
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Vehicle (Logistic)", marker = attr(color = vehicle_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin (PKPD)", marker = attr(color = cis_colour, size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Data",legendrank = 1, marker = attr(color = "black", size = 18)),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Model",legendrank = 2, line = attr(color = "black", width = 6))
            ];

        yatt = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = attr(
            title = "Time (days)",
            titlefont = attr(size = 28), 
            tickvals = [0, 2, 4, 6, 8, 10], 
            range = [-0.5, 10.5],
            tickfont = attr(size = 28)
        ),
        yaxis = yatt,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    if seelegend == true
        fig = plot([veh_data, pkpd_data, trace_veh, trace_pkpd, legend_traces...], layout);
    else
        fig = plot([veh_data, pkpd_data, trace_veh, trace_pkpd], layout);
    end

    display(fig)
    savefig(fig, joinpath(@__DIR__, "Figures", title * ".svg"); width = 800, height = 600)
end

## -------------------------------------

# Vehicle immune model plot function
function plot_immune(data, sol, error, title, color; seelegend = true, t2treatment = nothing)

    time = data[1]; 
    tumour = data[2]; T = data[3]; NK = data[4];

    error_tumour = error[1];
    error_T  = error[2];
    error_NK = error[3];
    
    tumour_shift_data = scatter(
        x = time,
        y = tumour,
        mode = "markers",
        name = "Data",
        marker = attr(size = 18, color = color),
        error_y = attr(
            type = "data",
            array = error_tumour[1],
            arrayminus = error_tumour[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    T_shift_data = scatter(
        x = [time[end]],
        y = [T],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color),
        error_y = attr(
            type = "data",
            array = [error_T],
            arrayminus = [error_T],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NK_shift_data = scatter(
        x = [time[end]],
        y = [NK],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color),
        error_y = attr(
            type = "data",
            array = [error_NK],
            arrayminus = [error_NK],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour = sol[2]; sim_T = sol[3]; sim_NK = sol[4]; 

    tumour_shift_trace = scatter(
        x = sim_time,
        y = sim_tumour,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color)
    )

    T_shift_trace = scatter(
        x = sim_time,
        y = sim_T,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color)
    )

    NK_shift_trace = scatter(
        x = sim_time,
        y = sim_NK,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color)
    )

    t2treatment_trace = if t2treatment !== nothing
        scatter(
            x = [t2treatment, t2treatment],
            y = [0, 500000],
            mode = "lines",
            name = "Treatment indication",
            showlegend = false,
            line = attr(width = 6,  dash = "dash", color = "#cbcbcb")
        )
    else
        nothing
    end

    if t2treatment === nothing
        xatt = attr(
            title = "Time (days)",
            titlefont = attr(size = 28), 
            tickvals = [0, 2, 4, 6, 8, 10], 
            range = [-0.5, 10.5],
            tickfont = attr(size = 28)
        )
    else
        xatt = attr(
            title = "Time from tumour initiation (days)",
            titlefont = attr(size = 28), 
            tickvals = [0, 20, 40, 60], 
            range = [-0.5, 60.5],
            tickfont = attr(size = 28)
        )
    end

    # Tumour volume
    if seelegend == true
        yatt_tumour = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt_tumour = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout_tumour = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = xatt,
        yaxis = yatt_tumour,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_tumour = [tumour_shift_data, tumour_shift_trace];
    if t2treatment_trace !== nothing
        push!(traces_tumour, t2treatment_trace)
    end
    fig_tumour = plot(traces_tumour, layout_tumour);
    display(fig_tumour)
    savefig(fig_tumour, joinpath(@__DIR__, "Figures", "Tumour volume - " * title * ".svg"); width = 800, height = 600)

    # T cells
    if seelegend == true
        yatt_T = attr(
            title = "CD8+ T cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 40000, 80000, 120000, 160000], 
            range = [0, 160000],
            ticktext = ["0", "0.4", "0.8", "1.2", "1.6"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="×10<sup>$(Int.(floor.(log10.(160000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_T = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 40000, 80000, 120000, 160000], 
            range = [0, 160000],
            ticktext = ["0", "0.4", "0.8", "1.2", "1.6"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_T = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_T],
        xaxis = xatt,
        yaxis = yatt_T,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_T = [T_shift_data, T_shift_trace];
    if t2treatment_trace !== nothing
        push!(traces_T, t2treatment_trace)
    end
    fig_T = plot(traces_T, layout_T);
    display(fig_T)
    savefig(fig_T, joinpath(@__DIR__, "Figures", "CD8+ T cells - " * title * ".svg"); width = 800, height = 600)

    # NK cells
    if seelegend == true
        yatt_NK = attr(
            title = "NK cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="×10<sup>$(Int.(floor.(log10.(500000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_NK = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_NK = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_NK],
        xaxis = xatt,
        yaxis = yatt_NK,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_NK = [NK_shift_data, NK_shift_trace];
    if t2treatment_trace !== nothing
        push!(traces_NK, t2treatment_trace)
    end
    fig_NK = plot(traces_NK, layout_NK);
    display(fig_NK)
    savefig(fig_NK, joinpath(@__DIR__, "Figures", "NK cells - " * title * ".svg"); width = 800, height = 600)

end

## -------------------------------------

# Immune model + cisplatin treatement plot function
function plot_immune_cis(data, sol, error, title, color; seelegend = true)

    time = data[1]; 
    tumour = data[2]; tumourt = data[3];
    T = data[4];      Tt = data[5];
    NK = data[6];     NKt = data[7];

    error_tumour = error[1]; error_tumourt = error[2];
    error_T  = error[3];     error_Tt = error[4];
    error_NK = error[5];     error_NKt = error[6];
    
    tumour_data = scatter(
        x = time,
        y = tumour,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = error_tumour[1],
            arrayminus = error_tumour[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    tumourt_data = scatter(
        x = time,
        y = tumourt,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = error_tumourt[1],
            arrayminus = error_tumourt[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    T_data = scatter(
        x = [time[end]],
        y = [T],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_T],
            arrayminus = [error_T],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    Tt_data = scatter(
        x = [time[end]],
        y = [Tt],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = [error_Tt],
            arrayminus = [error_Tt],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NK_data = scatter(
        x = [time[end]],
        y = [NK],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_NK],
            arrayminus = [error_NK],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NKt_data = scatter(
        x = [time[end]],
        y = [NKt],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = [error_NKt],
            arrayminus = [error_NKt],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour = sol[2]; sim_tumourt = sol[3]; sim_tumourt_nofit = sol[4]; 
    sim_T = sol[5];      sim_Tt = sol[6];      sim_Tt_nofit = sol[7];
    sim_NK = sol[8];     sim_NKt = sol[9];     sim_NKt_nofit = sol[10];

    tumour_trace = scatter(
        x = sim_time,
        y = sim_tumour,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[1]),
        showlegend = false
    )

    tumourt_trace = scatter(
        x = sim_time,
        y = sim_tumourt,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[2]),
        showlegend = false
    )

    tumourt_nofit_trace = scatter(
        x = sim_time,
        y = sim_tumourt_nofit,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[2], dash = "dot"),
        showlegend = false
    )

    T_trace = scatter(
        x = sim_time,
        y = sim_T,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    Tt_trace = scatter(
        x = sim_time,
        y = sim_Tt,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    Tt_nofit_trace = scatter(
        x = sim_time,
        y = sim_Tt_nofit,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[2], dash = "dot")
    )

    NK_trace = scatter(
        x = sim_time,
        y = sim_NK,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    NKt_trace = scatter(
        x = sim_time,
        y = sim_NKt,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )
    NKt_nofit_trace = scatter(
        x = sim_time,
        y = sim_NKt_nofit,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[2], dash = "dot")
    )

    xatt = attr(
        title = "Time (days)",
        titlefont = attr(size = 28), 
        tickvals = [0, 2, 4, 6, 8, 10], 
        range = [-0.5, 10.5],
        tickfont = attr(size = 28)
    )

    # Tumour volume
    if seelegend == true

        # Fake treatment legends
        legend_traces = 
            [
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Vehicle", marker = attr(color = color[1], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin", marker = attr(color = color[2], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Cisplatin (no refit)", line = attr(width = 6, color = color[2], dash = "dot")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Data", legendrank = 1, marker = attr(color = "black", size = 18)),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Model", legendrank = 2, line = attr(color = "black", width = 6))
            ];

        yatt_tumour = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt_tumour = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout_tumour = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = xatt,
        yaxis = yatt_tumour,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_tumour = [tumour_data, tumourt_data, tumour_trace, tumourt_trace, tumourt_nofit_trace];
    if seelegend == true
        push!(traces_tumour, legend_traces...)
    end
    fig_tumour = plot(traces_tumour, layout_tumour);
    display(fig_tumour)
    savefig(fig_tumour, joinpath(@__DIR__, "Figures", "Tumour volume - " * title * ".svg"); width = 800, height = 600)

    # T cells
    if seelegend == true
        yatt_T = attr(
            title = "CD8+ T cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 60000, 120000, 180000, 240000], 
            range = [0, 240000],
            ticktext = ["0", "0.6", "1.2", "1.8", "2.4"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="×10<sup>$(Int.(floor.(log10.(240000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_T = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 60000, 120000, 180000, 240000], 
            range = [0, 240000],
            ticktext = ["0", "0.6", "1.2", "1.8", "2.4"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_T = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_T],
        xaxis = xatt,
        yaxis = yatt_T,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_T = [T_data, Tt_data, T_trace, Tt_trace, Tt_nofit_trace];
    fig_T = plot(traces_T, layout_T);
    display(fig_T)
    savefig(fig_T, joinpath(@__DIR__, "Figures", "CD8+ T cells - " * title * ".svg"); width = 800, height = 600)

    # NK cells
    if seelegend == true
        yatt_NK = attr(
            title = "NK cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="×10<sup>$(Int.(floor.(log10.(500000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_NK = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_NK = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_NK],
        xaxis = xatt,
        yaxis = yatt_NK,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_NK = [NK_data, NKt_data, NK_trace, NKt_trace, NKt_nofit_trace];
    fig_NK = plot(traces_NK, layout_NK);
    display(fig_NK)
    savefig(fig_NK, joinpath(@__DIR__, "Figures", "NK cells - " * title * ".svg"); width = 800, height = 600)

end

## -------------------------------------

# Immune model + cisplatin treatement plot function
function plot_immune_icb(data, sol, error, title, color; seelegend = true)

    time = data[1]; 
    tumour = data[2]; tumourt = data[3];
    T = data[4];      
    NK = data[5];     

    error_tumour = error[1]; error_tumourt = error[2];
    error_T  = error[3];     
    error_NK = error[4];     
    
    tumour_data = scatter(
        x = time,
        y = tumour,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = error_tumour[1],
            arrayminus = error_tumour[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    tumourt_data = scatter(
        x = time,
        y = tumourt,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = error_tumourt[1],
            arrayminus = error_tumourt[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    T_data = scatter(
        x = [time[end]],
        y = [T],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_T],
            arrayminus = [error_T],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NK_data = scatter(
        x = [time[end]],
        y = [NK],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_NK],
            arrayminus = [error_NK],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour = sol[2]; sim_tumourt = sol[3]; 
    sim_T = sol[4];      sim_Tt = sol[5];      
    sim_NK = sol[6];     sim_NKt = sol[7];     

    tumour_trace = scatter(
        x = sim_time,
        y = sim_tumour,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[1]),
        showlegend = false
    )

    tumourt_trace = scatter(
        x = sim_time,
        y = sim_tumourt,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[2]),
        showlegend = false
    )

    T_trace = scatter(
        x = sim_time,
        y = sim_T,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    Tt_trace = scatter(
        x = sim_time,
        y = sim_Tt,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    NK_trace = scatter(
        x = sim_time,
        y = sim_NK,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    NKt_trace = scatter(
        x = sim_time,
        y = sim_NKt,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )
   
    xatt = attr(
        title = "Time (days)",
        titlefont = attr(size = 28), 
        tickvals = [0, 2, 4, 6, 8, 10], 
        range = [-0.5, 10.5],
        tickfont = attr(size = 28)
    )

    # Tumour volume
    if seelegend == true

        # Fake treatment legends
        legend_traces = 
            [
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Vehicle", marker = attr(color = color[1], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "ICB", marker = attr(color = color[2], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Data", legendrank = 1, marker = attr(color = "black", size = 18)),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Model", legendrank = 2, line = attr(color = "black", width = 6))
            ];

        yatt_tumour = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt_tumour = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout_tumour = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = xatt,
        yaxis = yatt_tumour,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_tumour = [tumour_data, tumourt_data, tumour_trace, tumourt_trace];
    if seelegend == true
        push!(traces_tumour, legend_traces...)
    end
    fig_tumour = plot(traces_tumour, layout_tumour);
    display(fig_tumour)
    savefig(fig_tumour, joinpath(@__DIR__, "Figures", "Tumour volume - " * title * ".svg"); width = 800, height = 600)

    # T cells
    if seelegend == true
        yatt_T = attr(
            title = "CD8+ T cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 100000, 200000, 300000, 400000], 
            range = [0, 400000],
            ticktext = ["0", "1", "2", "3", "4"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="×10<sup>$(Int.(floor.(log10.(400000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_T = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 100000, 200000, 300000, 400000], 
            range = [0, 400000],
            ticktext = ["0", "1", "2", "3", "4"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_T = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_T],
        xaxis = xatt,
        yaxis = yatt_T,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_T = [T_data, T_trace, Tt_trace];
    fig_T = plot(traces_T, layout_T);
    display(fig_T)
    savefig(fig_T, joinpath(@__DIR__, "Figures", "CD8+ T cells - " * title * ".svg"); width = 800, height = 600)

    # NK cells
    if seelegend == true
        yatt_NK = attr(
            title = "NK cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="×10<sup>$(Int.(floor.(log10.(500000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_NK = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_NK = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_NK],
        xaxis = xatt,
        yaxis = yatt_NK,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_NK = [NK_data, NK_trace, NKt_trace];
    fig_NK = plot(traces_NK, layout_NK);
    display(fig_NK)
    savefig(fig_NK, joinpath(@__DIR__, "Figures", "NK cells - " * title * ".svg"); width = 800, height = 600)

end

## -------------------------------------

# Immune model + cisplatin treatement plot function
function plot_immune_all(data, sol, error, title, color; seelegend = true)

    time = data[1]; 
    tumour = data[2]; tumourt = data[3]; tumourt2 = data[4]; tumourt3 = data[5];
    T = data[6];      Tt = data[7];
    NK = data[8];     NKt = data[9];

    error_tumour = error[1]; error_tumourt = error[2]; error_tumourt2 = error[3]; error_tumourt3 = error[4];
    error_T = error[5];      error_Tt = error[6];
    error_NK = error[7];     error_NKt = error[8];  
    
    tumour_data = scatter(
        x = time,
        y = tumour,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = error_tumour[1],
            arrayminus = error_tumour[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    tumourt_data = scatter(
        x = time,
        y = tumourt,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = error_tumourt[1],
            arrayminus = error_tumourt[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    tumourt2_data = scatter(
        x = time,
        y = tumourt2,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[3]),
        error_y = attr(
            type = "data",
            array = error_tumourt2[1],
            arrayminus = error_tumourt2[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    tumourt3_data = scatter(
        x = time,
        y = tumourt3,
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[4]),
        error_y = attr(
            type = "data",
            array = error_tumourt3[1],
            arrayminus = error_tumourt3[2],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    T_data = scatter(
        x = [time[end]],
        y = [T],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_T],
            arrayminus = [error_T],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    Tt_data = scatter(
        x = [time[end]],
        y = [Tt],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = [error_Tt],
            arrayminus = [error_Tt],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NK_data = scatter(
        x = [time[end]],
        y = [NK],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[1]),
        error_y = attr(
            type = "data",
            array = [error_NK],
            arrayminus = [error_NK],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    NKt_data = scatter(
        x = [time[end]],
        y = [NKt],
        mode = "markers",
        name = "Data",
        showlegend = false,
        marker = attr(size = 18, color = color[2]),
        error_y = attr(
            type = "data",
            array = [error_NKt],
            arrayminus = [error_NKt],
            symmetric = false,
            visible = true,
            thickness = 4,
            width = 8
        )
    )

    sim_time = sol[1]; 
    sim_tumour = sol[2]; sim_tumourt = sol[3]; sim_tumourt2 = sol[4]; sim_tumourt3 = sol[5]; 
    sim_T = sol[6];      sim_Tt = sol[7];      sim_Tt2 = sol[8];      sim_Tt3 = sol[9];     
    sim_NK = sol[10];   sim_NKt = sol[11];    sim_NKt2 = sol[12];    sim_NKt3 = sol[13];     

    tumour_trace = scatter(
        x = sim_time,
        y = sim_tumour,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[1]),
        showlegend = false
    )

    tumourt_trace = scatter(
        x = sim_time,
        y = sim_tumourt,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[2]),
        showlegend = false
    )

    tumourt2_trace = scatter(
        x = sim_time,
        y = sim_tumourt2,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[3]),
        showlegend = false
    )

    tumourt3_trace = scatter(
        x = sim_time,
        y = sim_tumourt3,
        mode = "lines",
        name = "Model",
        line = attr(width = 6, color = color[4]),
        showlegend = false
    )

    T_trace = scatter(
        x = sim_time,
        y = sim_T,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    Tt_trace = scatter(
        x = sim_time,
        y = sim_Tt,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    Tt2_trace = scatter(
        x = sim_time,
        y = sim_Tt2,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[3])
    )

    Tt3_trace = scatter(
        x = sim_time,
        y = sim_Tt3,
        mode = "lines",
        name = "CD8+ T cell",
        showlegend = false,
        line = attr(width = 6, color = color[4])
    )

    NK_trace = scatter(
        x = sim_time,
        y = sim_NK,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[1])
    )

    NKt_trace = scatter(
        x = sim_time,
        y = sim_NKt,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[2])
    )

    NKt2_trace = scatter(
        x = sim_time,
        y = sim_NKt2,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[3])
    )

    NKt3_trace = scatter(
        x = sim_time,
        y = sim_NKt3,
        mode = "lines",
        name = "NK cell",
        showlegend = false,
        line = attr(width = 6, color = color[4])
    )
   
    xatt = attr(
        title = "Time (days)",
        titlefont = attr(size = 28), 
        tickvals = [0, 2, 4, 6, 8, 10], 
        range = [-0.5, 10.5],
        tickfont = attr(size = 28)
    )

    # Tumour volume
    if seelegend == true

        # Fake treatment legends
        legend_traces = 
            [
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Vehicle", marker = attr(color = color[1], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin", marker = attr(color = color[2], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "ICB", marker = attr(color = color[3], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Cisplatin + ICB", marker = attr(color = color[4], size = 18, symbol = "square")),
                scatter(x = [NaN], y = [NaN], mode = "markers", name = "Data", legendrank = 1, marker = attr(color = "black", size = 18)),
                scatter(x = [NaN], y = [NaN], mode = "lines", name = "Model", legendrank = 2, line = attr(color = "black", width = 6))
            ];

        yatt_tumour = attr(
            title = "Tumour volume (mm³)", 
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    else
        yatt_tumour = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 400, 800, 1200, 1600], 
            range = [0, 1670],
            tickfont = attr(size = 28)
        )
    end

    layout_tumour = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        xaxis = xatt,
        yaxis = yatt_tumour,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_tumour = [tumour_data, tumourt_data, tumourt2_data, tumourt3_data, tumour_trace, tumourt_trace, tumourt2_trace, tumourt3_trace];
    if seelegend == true
        push!(traces_tumour, legend_traces...)
    end
    fig_tumour = plot(traces_tumour, layout_tumour);
    display(fig_tumour)
    savefig(fig_tumour, joinpath(@__DIR__, "Figures", "Tumour volume - " * title * ".svg"); width = 800, height = 600)

    # T cells
    if seelegend == true
        yatt_T = attr(
            title = "CD8+ T cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="×10<sup>$(Int.(floor.(log10.(500000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_T = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_T = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_T = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_T],
        xaxis = xatt,
        yaxis = yatt_T,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_T = [T_data, Tt_data, T_trace, Tt_trace, Tt2_trace, Tt3_trace];
    fig_T = plot(traces_T, layout_T);
    display(fig_T)
    savefig(fig_T, joinpath(@__DIR__, "Figures", "CD8+ T cells - " * title * ".svg"); width = 800, height = 600)

    # NK cells
    if seelegend == true
        yatt_NK = attr(
            title = "NK cell count", 
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="×10<sup>$(Int.(floor.(log10.(500000))))</sup>",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    else
        yatt_NK = attr(
            title = "", 
            showticklabels = false,
            titlefont = attr(size = 28), 
            tickvals = [0, 125000, 250000, 375000, 500000], 
            range = [0, 500000],
            ticktext = ["0", "1.25", "2.5", "3.75", "5"],
            tickfont = attr(size = 28)
        )
        exp_NK = attr(
            text="",
            x=0, y=1.0,
            xref = "paper", yref = "paper",
            xanchor="left", yanchor="bottom",
            showarrow=false,
            font=attr(size=24)
        )
    end

    layout_NK = Layout(
        template = craig_lab_template,
        title = attr(
            text = title,
            font = attr(size = 32),
            x = 0.5,
            xanchor = "center"
        ),
        annotations = [exp_NK],
        xaxis = xatt,
        yaxis = yatt_NK,
        width = 800,
        height = 600,
        legend = attr(
            x=0.02, y=1, 
            xanchor = "left", 
            orientation="v",
            font = attr(size = 22)
        ),
        margin = attr(l = 130, r = 20, t = 80, b = 80),
        showlegend = seelegend
    )

    traces_NK = [NK_data, NKt_data, NK_trace, NKt_trace, NKt2_trace, NKt3_trace];
    fig_NK = plot(traces_NK, layout_NK);
    display(fig_NK)
    savefig(fig_NK, joinpath(@__DIR__, "Figures", "NK cells - " * title * ".svg"); width = 800, height = 600)

end