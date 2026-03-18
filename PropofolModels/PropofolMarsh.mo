within PropofolModels;
model PropofolMarsh
  // Time is implicitly in minutes
  parameter Real weight = 70 "Weight in kg";

  parameter Real k10 = 0.119;
  parameter Real k12 = 0.112;
  parameter Real k13 = 0.0419;
  parameter Real k21 = 0.055;
  parameter Real k31 = 0.0033;
  parameter Real ke0 = 1.21;
  
  parameter Real V1 = 0.228 * weight;
  parameter Real CL_elim = k10 * V1;
  parameter Real CL12 = k12 * V1;
  parameter Real CL13 = k13 * V1;
  parameter Real V2 = CL12 / k21;
  parameter Real V3 = CL13 / k31;

  // Compartments
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment central(V=V1) annotation(
    Placement(transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}})));
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment peripheral1(V=V2) annotation(
    Placement(transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}})));
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment peripheral2(V=V3) annotation(
    Placement(transformation(origin = {80, -10}, extent = {{-10, -10}, {10, 10}})));
  
  // Transfers
  Pharmacolibrary.Pharmacokinetic.TransferFirstOrderNonSym transfer1(CLa=CL12, CLb=CL12) annotation(
    Placement(transformation(origin = {40, 30}, extent = {{-10, -10}, {10, 10}})));
  Pharmacolibrary.Pharmacokinetic.TransferFirstOrderNonSym transfer2(CLa=CL13, CLb=CL13) annotation(
    Placement(transformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}})));
  
  // Elimination
  Pharmacolibrary.Pharmacokinetic.ClearanceDrivenElimination elim(CL=CL_elim) annotation(
    Placement(transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}})));
  
  // Infusion Source
  Pharmacolibrary.Sources.VariableInfusion infusion_source annotation(
    Placement(transformation(origin = {-40, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression infusion_rate(y = if time < 1.0 then 2.0 * weight else if time < 60.0 then 10.0 * weight / 60.0 else 0.0) annotation(
    Placement(transformation(origin = {-80, 40}, extent = {{-10, -10}, {10, 10}})));
  
  // Effect Site
  Real Ce(start=0.0) "Effect site concentration";

equation
  // Effect site calculation
  der(Ce) = ke0 * (central.cport.c - Ce);
  
  // Connections
  connect(infusion_rate.y, infusion_source.massFlow) annotation(
    Line(points = {{-68, 40}, {-50, 40}}, color = {0, 0, 127}));
  connect(infusion_source.cport, central.cport) annotation(
    Line(points = {{-30, 40}, {0, 40}, {0, 10}}, color = {114, 159, 207}));
  connect(elim.cport, central.cport) annotation(
    Line(points = {{-40, 10}, {-20, 10}, {0, 10}}, color = {114, 159, 207}));
  connect(central.cport, transfer1.cport_b) annotation(
    Line(points = {{0, 10}, {20, 10}, {20, 30}, {30, 30}}, color = {114, 159, 207}));
  connect(transfer1.cport_a, peripheral1.cport) annotation(
    Line(points = {{50, 30}, {80, 30}, {80, 40}}, color = {114, 159, 207}));
  connect(central.cport, transfer2.cport_b) annotation(
    Line(points = {{0, 10}, {20, 10}, {20, -10}, {30, -10}}, color = {114, 159, 207}));
  connect(transfer2.cport_a, peripheral2.cport) annotation(
    Line(points = {{50, -10}, {80, -10}, {80, 0}}, color = {114, 159, 207}));
  
  annotation(
    Documentation(info = "<html><head></head><body><h1>Marsh Propofol Model</h1><p>This model implements the standard Marsh PK model for propofol. It is a 3-compartment linearly scaled model relying primarily on patient weight to expand distribution volumes and proportional clearances. It includes an effect-site (ke0) equilibration lag parameter.</p></body></html>"),
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")));
end PropofolMarsh;
