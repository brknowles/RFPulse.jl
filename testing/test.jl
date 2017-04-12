
using RFPulse;

println("Reading RF Pulse file");
pulse = RFPulse.read("default.json");

println("Simulating...");
RFPulse.simulate!(pulse);


sim = pulse["SimulationParameters"]["Output"];

println("Test finished");
