module RFPulse

# requirements
using JSON;


#
# TODOs:
#	Everything
#
#
const default_config = Dict{String,Any}("max_b1"=>20e-6,"gamma"=>42.57e3,"max_grad"=>24.0,"slice_thickness"=>5.0);

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
    max_b1   = sysconfig["MaxB1"];
    gamma    = sysconfig["Gamma"];
    max_grad = sysconfig["MaxGrad"];
    
    #RF pars
    ns    = pulseconfig["Samples"];
    dur   = pulseconfig["Duration"];
    slice = pulseconfig["SliceThickness"];
    beta  = pulseconfig["HyperbolicSecantParameters"]["Beta"];
    mu    = pulseconfig["HyperbolicSecantParameters"]["Mu"];

    # create the time vector
    t = linspace(-dur/2.,dur/2,ns);

    # lots of ugly vectorized julia code
    b1_hs = max_b1.*sech.(beta.*t);         # Amplitude
    df_hs = -1.0.*mu.*beta.*tanh.(beta.*t); # Frequency

    # bandwidth
    bw = mu*beta/pi;

    # Gradient Strength
    g_hs = 2.0*pi*bw/(gamma*slice);  #mT/m

    #adiabacity
    b1_hs_limit = sqrt(mu) * beta / gamma;  #B1 must be above this limit for adiabacity
    
    # gradient vector
    grad = g_hs.*ones(ns);

    # rf vector
    rf = b1_hs .* exp.(-im .* df_hs .* t );
    
    # new dict
    props = config;
    rfshape["Real"] = real(rf);
    rfshape["Imag"] = imag(rf);
    rfshape["Grad"] = grad;
    
    props["RFShape"] = rfshape;

    new(props);
  end

end # HyperbolicSecant

function write(pulse::AbstractRFPulse,filename::String)

	fid = open(filename,"w");
	
	JSON.print(fid,pulse.properties);
	
	close(fid);

end

end #module



