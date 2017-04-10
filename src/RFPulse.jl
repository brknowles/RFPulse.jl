module RFPulse


#
# TODOs:
#	Everything
#
#
const default_config = Dict{String,Any}("max_b1"=>20e-6,"gamma"=>42.57e3,"max_grad"=>24.0,"slice_thickness"=>5.0);

abstract AbstractRFPulse;

type HyperbolicSecant <: AbstractRFPulse
  samples::Integer;
  duration::AbstractFloat; #s
  beta::AbstractFloat; # rad/s, related to cycles
  mu::AbstractFloat; # related to sharpness
  bandwidth::AbstractFloat # Hz
  limit::AbstractFloat; # Approx. Adiabacity limit, mT
  rf_samples::Array{Complex{AbstractFloat},1}; #mT
  grad_samples::Array{AbstractFloat,1}; #mT
  time::Array{AbstractFloat,1}; #secs
  
 
  """
  Create the HS Pulse
  """
  function HyperbolicSecant(ns::Integer,
			    dur::AbstractFloat,
			    beta::AbstractFloat,
			    mu::AbstractFloat)
			    # config::Dict{String,Any})
    config = [];
    # parse function args
    if isempty(config)
      config = default_config;
    end

    # parse config Dictionary
    max_b1   = config["max_b1"];
    gamma    = config["gamma"];
    max_grad = config["max_grad"];
    slice    = config["slice_thickness"];

    # create the time vector
    t = linspace(-dur/2.,dur/2,ns);

    # lots of ugly vectorized julia code
    b1_hs = max_b1.*sech(beta.*t);         # Amplitude
    df_hs = -1.0.*mu.*beta.*tanh(beta.*t); # Frequency

    # bandwidth
    bw = mu*beta/pi;

    # Gradient Strength
    g_hs = 2.0*pi*bw/(gamma*slice);  #mT/m

    #adiabacity
    b1_hs_limit = sqrt(mu) * beta / gamma;  #B1 must be above this limit for adiabacity
    
    # gradient vector
    grad = g_hs.*ones(ns);

    # rf vector
    rf = b1_hs .* exp(-im .* df_hs .* t );

    new(ns,dur,beta,mu,bw,b1_hs_limit,rf,grad,t+(dur/2.));
  end

end # HyperbolicSecant

type FOCI <: AbstractRFPulse
  samples::Integer;
  duration::AbstractFloat; #s
  beta::AbstractFloat; # rad/s, related to cycles
  mu::AbstractFloat; # related to sharpness
  bandwidth::AbstractFloat # Hz
  limit::AbstractFloat; # Approx. Adiabacity limit, mT
  rf_samples::Array{Complex{AbstractFloat},1}; #mT
  grad_samples::Array{AbstractFloat,1}; #mT
  time::Array{AbstractFloat,1}; #secs
  
 
  """
  Create the FOCI pulse
  """
  function FOCI(ns::Integer,
		dur::AbstractFloat,
		beta::AbstractFloat,
		mu::AbstractFloat)
		# config::Dict{String,Any})
    config = [];
    # parse function args
    if isempty(config)
      config = default_config;
    end

    # create a HS pulse 
    HSpulse = HyperbolicSecant(ns,dur,beta,mu);


    new(ns,dur,beta,mu,bw,b1_hs_limit,rf,grad,t+(dur/2.));
  end

end # FOCI

end #module


# not part of anything useful, should be removed
if 0
# RF Pulse setup
SLT = 5.0e-3; # Slice thickness, mm
Gamma = 42.57e3; # Proton Gamma, Hz/mT
num_steps = 101;
T = 10.0e-3; # pulse duration, s
hT = T/2.0;
t_vec = linspace(-hT,hT,num_steps);

#non-vectorised
# From Geoff's Paper (1997)
function foci_c_envelope(t::Real,Beta::Real,amax::Real)
  arg = sech(Beta.*t);
  if arg > (1.0/amax)
    a = 1.0/arg;
  else
    a = amax;
  end
  return a;
end

# verse sin2 time dialation
function verse_sinsq_td(t,u)
  a = 1.9;
  return a - 2.0*(a-1.0)*cos(pi*t/T).^2.0;
end

# verse triangluar time dialation
function verse_triangular_td(t,u)
  amin = 0.7;
  return amin+(abs(t)/T);
end



# Foci-Constants
B1p = 210e-3; #mT (Peak B1)
Gmax = 29.0; #mT/m
Beta = 500.0; # rad/s
Mu = 5.0;

# HS Pulse design
b1_hs(t) = B1p.*sech(Beta.*t);
df_hs(t)= -1.0.*Mu.*Beta.*tanh(Beta.*t);
bw_hs = Mu*Beta/pi;
g_hs = 2*pi*bw_hs/(Gamma*SLT);  #mT/m
b1_hs_limit = sqrt(Mu).*Beta ./ Gamma;  #B1 must be above this limit for adiabacity
G_hs(t) = g_hs.*ones(size(t));

#checks
#if g_hs > Gmax
#  error("Requested RF pulse exceeds hardware limits, g_hs: ", g_hs);
#end

# FOCI envelope
A(t) = foci_c_envelope.(t,Beta,Gmax/g_hs);

# FOCI Gradient
G_foci(t) = A.(t).* G_hs(t);
b1_foci(t)= A.(t).*B1p.*sech(Beta.*t);
df_foci(t)= -1.0.*A.(t).*Mu.*Beta.*tanh(Beta.*t);


VERSE_MODE = "quadratic";
# VERSE Time-dialation
if VERSE_MODE == "sinsq"
  dtau_dt(t,u) = verse_sinsq_td(t,u);
elseif VERSE_MODE == "triangular"
  dtau_dt(t,u) = verse_triangular_td(t,u);
  prob = ODEProblem(dtau_dt,-hT,(-hT,hT));
  s = solve(prob,dt=(T/100.0));

elseif VERSE_MODE == "quadratic"
  Tau(t,t1,t2) = (t1.*t.^5 + t2.*t.^3 + t) / (t1+t2+1);
  Tau(t) = Tau(t/hT,0.1,0.4).*hT;
  sol(x)= ForwardDiff.derivative(Tau,x)
  dtau_dt = sol.(t_vec) ./ maximum(sol.(t_vec));

  #normalise

  #prob = ODEProblem(dtau_dt,-hT,(-hT,hT));
  #s = solve(prob,Euler(),dt=(T/100.0));
  #verse_td = s.u
else
  error("No verse mode selected");
end

#verse-foci
# FOCI Gradient
A_verse(t) = A(t) .* dtau_dt ;
G_versefoci(t) = A_verse(t) .* G_hs.(t);
b1_versefoci(t)= A_verse(t) .* B1p .* sech(Beta.*Tau(t));
df_versefoci(t)= -1.0.*A_verse(t).*Mu.*Beta.*tanh(Beta.*Tau(t));

end
