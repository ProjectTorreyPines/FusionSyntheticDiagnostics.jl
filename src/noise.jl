import DSP: DSP, AbstractFFTs

mutable struct Noise
    var"pgram"::DSP.Periodograms.Periodogram
    var"time"::Vector{Float64}
    var"noise"::Vector{Float64}
end

function Noise(time::Vector{Float64}, signal::Vector{Float64}, binwidth::Float64)
    fs = 1 / (time[2] - time[1])
    nperseg = round(Int64, fs / binwidth)
    noverlap = div(nperseg, 2)
    pgram = DSP.Periodograms.welch_pgram(signal, nperseg, noverlap; fs)
    return Noise(pgram, time, signal)
end

function Noise(pgram::DSP.Periodograms.Periodogram)
    return Noise(pgram, zeros(Float64, 0), zeros(Float64, 0))
end

function Noise(power_spectrum::Vector{Float64}, freq::AbstractRange)
    pgram = DSP.Periodograms.Periodogram(power_spectrum, freq)
    return Noise(pgram)
end

function generate_noise(n::Noise, t::Vector{Float64})::Vector{Float64}
    signal = zeros(Float64, length(t))
    amp_spec = sqrt.(n.pgram.power * 2)
    for (ii, ff) ∈ enumerate(n.pgram.freq)
        signal .+= amp_spec[ii] .* cos.(2π .* (ff .* t .+ randn()))
    end
    return signal
end

function generate_noise!(n::Noise, t::Vector{Float64})::Noise
    signal = generate_Noise(n, t)
    n.time = t
    n.noise = signal
    return n
end
