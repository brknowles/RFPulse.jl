"""
Creates RF Pulses based on JSON files.
see testing/ for example files
"""
module RFPulse

# requirements
using JSON;
using BlochSim;
using PyPlot; pygui(true);

import Base: getindex, display;

export  AbstractRFPulse,
	HyperbolicSecant,
	ExternalWaveform,
	display, 
	read_pulse, 
	write_pulse, 
	getindex, 
	simulate!,
	plot_pulse,
	plot_simulation;

#
# TODOs:
#	Other Pulses
#	VERSE
#
#

"""
 Abstract RF Pulse:
 Dictionary of RF Properties
 Interface: 
  "RFShape" => Complex RF pulse samples
  "GradShape" => Floating Point Gradient Samples
  "Time" => Vector of time points
"""
abstract type AbstractRFPulse end;

type HyperbolicSecant <: AbstractRFPulse
	properties::Dict{String,Any};
  
 
"""
  Create the HS Pulse
"""
  function HyperbolicSecant(config::Dict{String,Any})
  
    sysconfig   = config["SystemParameters"];
    pulseconfig = config["PulseParameters"];
    rfshape     = config["RFShape"];
  
    # system pars
    Gamma    = sysconfig["Gamma"];
    sys_max_grad = sysconfig["MaxGrad"];
    
    #RF pars
    ns    = pulseconfig["Samples"];
    dur   = pulseconfig["Duration"];
    slice = pulseconfig["SliceThickness"];
    Beta  = pulseconfig["HyperbolicSecantParameters"]["Beta"];
    Mu    = pulseconfig["HyperbolicSecantParameters"]["Mu"];
    #ref_b1= pulseconfig["RefB1"];
    b1_factor = pulseconfig["HyperbolicSecantParameters"]["B1ScaleFactor"];

    #adiabacity
    b1_hs_limit = sqrt(Mu) * Beta / (Gamma*2*pi);  #B1 must be above this limit for adiabacity
    ref_b1= b1_hs_limit * b1_factor;
    
    # create the time vector
    t = linspace(-dur/2.,dur/2,ns);

    b1_hs = ref_b1.*sech.(Beta.*t);         # Amplitude
    df_hs = -1.0.*Mu.*Beta.*tanh.(Beta.*t); # Frequency

    # bandwidth
    bw = Mu*Beta/pi;

    # Gradient Strength
    g_hs = bw/(Gamma*slice);  #mT/m

    # gradient vector
    refgrad  = g_hs;
    grad     = ones(ns); # const. gradient for HS pulses

    # rf vector
    rf_amp   = b1_hs./ref_b1;
    rf_phase = (df_hs .* t);
    rf       = cat(2,rf_amp,rf_phase);

    
    # set some properties
    config["PulseParameters"]["Bandwidth"] = bw;
    config["PulseParameters"]["RefGrad"] = refgrad;
    config["PulseParameters"]["MinSlice"] = 1.0e-3; # default
    config["PulseParameters"]["MaxSlice"] = bw/(Gamma*sys_max_grad);
    config["PulseParameters"]["RefB1"] = ref_b1;
    
    # add to dictionary
    props = config;
    props["RFShape"] = rf;
    props["GradShape"] = grad;
    props["Time"] = Array(linspace(0.0,1.0,ns));

    #finally calculate flip angle
    calculate_flip_angle!(props);

    new(props);
  end

end # HyperbolicSecant

type ExternalWaveform <: AbstractRFPulse
  properties::Dict{String,Any};

  """
  create it
  """
  function ExternalWaveform(config::Dict{String,Any});
    new(config);
  end
end

#getindex for all abstract RF pulses
getindex(p::AbstractRFPulse,idx) = p.properties[idx];

"""
 Reads a JSON config file and returns the requested RF Pulse Type
"""
function read_pulse(filename::String)
  fid = open(filename,"r");
  rawdata = read(fid);
  strdata = String(rawdata);
  config = JSON.parse(strdata);


  # create the RF pulse
  if config["PulseType"] == "HyperbolicSecant"
    pulse = HyperbolicSecant(config);
  else
    error("UnKnown Pulse type");
  end
  return pulse;
end

"""
  Writes the RFPulse type to a JSON file
"""
function write_pulse(pulse::AbstractRFPulse,filename::String)

	fid = open(filename,"w");
	
	JSON.print(fid,pulse.properties);
	
	close(fid);

end

"""
RF Pulse Interface to BlochSim.jl
performs an inplace modification to the RF pulse
"""
function simulate!(pulse::AbstractRFPulse)
  # parse the config dict for the simulation parameters

  ref_b1   = pulse["PulseParameters"]["RefB1"];
  ref_grad = pulse["PulseParameters"]["RefGrad"];
  dur      = pulse["PulseParameters"]["Duration"];

  Gamma    = pulse["SystemParameters"]["Gamma"];
  sim_pars = pulse["SimulationParameters"];
  ns       = sim_pars["PositionSamples"];
  nf       = sim_pars["FrequencySamples"];
  t1       = sim_pars["T1"];
  t2       = sim_pars["T2"];
  frange   = sim_pars["FrequencyRange"];
  prange   = sim_pars["PositionRange"];

  rf_magn  = pulse["RFShape"][:,1] .* ref_b1;
  rf_phase = pulse["RFShape"][:,2];
  rf       = rf_magn .* exp.(im .* rf_phase);
  grad     = pulse["GradShape"] .* ref_grad;
  t        = pulse["Time"] .* dur;
  
  pos_vec  = linspace(-prange/2.,prange/2.,ns);

  #TODO, fix frequency range and think about how to represent a 2D sim in the JSON file
  freq_vec = 0.0; #linspace(-frange/2.,frange/2.,nf);

  (m,m2) = BlochSim.sliceprofile(rf,grad,t,t1,t2,pos_vec,freq_vec,Gamma);

  output = Dict{String,Any}("Mx" => m[1,:],
			    "My" => m[2,:],
			    "Mz" => m[3,:],
			    "X"  => pos_vec,
			    "F"  => freq_vec);
  pulse["SimulationParameters"]["Output"] = output;
end

"""
Plots the RF Pulse
"""
function plot_pulse(pulse::AbstractRFPulse)
  ref_b1   = pulse["PulseParameters"]["RefB1"] * 1000.0; # convert to uT
  ref_grad = pulse["PulseParameters"]["RefGrad"];
  dur      = pulse["PulseParameters"]["Duration"] *1000.0; # convert to ms for viewing

  # get pulse shape from the dictionary
  t    = pulse["Time"];
  rf_magn    = pulse["RFShape"][:,1] .* ref_b1;
  rf_phase   = pulse["RFShape"][:,2];
  rf = rf_magn .* exp.(im .* rf_phase);
  grad = pulse["GradShape"] .* ref_grad;

  subplot(1,2,1);
  plot(t,abs.(rf));
  plot(t,real.(rf));
  plot(t,imag.(rf));
  title("RF Pulse");
  ylabel("B1 [uT]")
  xlabel("Time [ms]");
  subplot(1,2,2);
  plot(t,grad);
  title("gradient");
  ylabel("Grad [mT]")
  xlabel("Time [ms]");
end


"""
Plots the simulation results
"""
function plot_simulation(pulse::AbstractRFPulse)
  # get results from the dictionary
  sim = pulse["SimulationParameters"]["Output"]; 
  st  = pulse["PulseParameters"]["SliceThickness"]; 
  x   = sim["X"];
  mz = sim["Mz"];
  my = sim["My"];
  mx = sim["Mx"];
  mt = abs.(mx + im.*my);

  # ideal slice profile
  ideal_sp = ones(size(x));
  ideal_sp[ find((f)->abs.(f)<st/2.0,x)] = -1.0;

  subplot(1,2,1);
  plot(x,mz);
  plot(x,ideal_sp);
  title("Slice Profile");
  subplot(1,2,2);
  plot(x,mt);
  title("Signal");
end

"""
output pulse parameters to something more human readable
"""
function display(pulse::AbstractRFPulse)
  @printf("RF Pulse: %s \n",pulse["PulseType"]);
  @printf("\nSystem Parameters: \n");
  syspars = pulse["SystemParameters"];
  for (k,v) in syspars
    println(k,": ", v);
  end
  @printf("\nPulseParameters:\n");
  pulsepars = pulse["PulseParameters"];
  for (k,v) in pulsepars
    println(k,": ", v);
  end
end


"""
 Calculates the Ffip angle by integrating the RF Pulse envelope
 args can be RFpulse of config dictionary
"""
function calculate_flip_angle!(config::Union{AbstractRFPulse,Dict{String,Any}})
  Gamma  = config["SystemParameters"]["Gamma"];
  rf_amp = config["RFShape"][:,1] .* config["PulseParameters"]["RefB1"];
  t      = config["Time"] .* config["PulseParameters"]["Duration"];

  fa = 2pi* Gamma * sum(rf_amp .* t);
  config["PulseParameters"]["FlipAngleDeg"] = rad2deg(fa);
end

end #module
