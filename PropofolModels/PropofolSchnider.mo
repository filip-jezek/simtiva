within PropofolModels;
model PropofolSchnider
  parameter Real weight = 70 "Weight in kg";
  parameter Real age = 40 "Age in years";
  parameter Real height = 170 "Height in cm";
  parameter Integer gender = 0 "0=Male";
  parameter Boolean is_adj_bw = false "Use Adj BW";

  parameter Real ibw_male = 50.0 + 0.9 * (height - 152.0);
  parameter Real ibw_female = 45.5 + 0.9 * (height - 152.0);
  parameter Real ibw = if gender == 0 then ibw_male else ibw_female;
  parameter Real adj_bw = ibw + 0.4 * (weight - ibw);

  parameter Real lbm_male = 1.1 * weight - 128.0 * (weight / height)^2;
  parameter Real lbm_female = 1.07 * weight - 148.0 * (weight / height)^2;
  parameter Real lbm = if gender == 0 then lbm_male else lbm_female;

  parameter Real lbm2_male = 1.1 * adj_bw - 128.0 * (adj_bw / height)^2;
  parameter Real lbm2_female = 1.07 * adj_bw - 148.0 * (adj_bw / height)^2;
  parameter Real lbm2 = if gender == 0 then lbm2_male else lbm2_female;

  parameter Real vc = 4.27;
  parameter Real v2 = 18.9 - 0.391 * (age - 53.0);
  parameter Real v3 = 238.0;

  parameter Real cl1 = if not is_adj_bw then 
                          1.89 + 0.0456 * (weight - 77.0) - 0.0681 * (lbm - 59.0) + 0.0264 * (height - 177.0) 
                       else 
                          1.89 + 0.0456 * (adj_bw - 77.0) - 0.0681 * (lbm2 - 59.0) + 0.0264 * (height - 177.0);

  parameter Real cl2 = 1.29 - 0.024 * (age - 53.0);
  parameter Real cl3 = 0.836;
  parameter Real ke0 = 0.456;

  parameter Real V1 = vc;
  parameter Real CL_elim = cl1;
  parameter Real CL12 = cl2;
  parameter Real CL13 = cl3;
  parameter Real V2_param = v2;
  parameter Real V3_param = v3;

  // Compartments
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment central(V=V1) annotation(
    Placement(transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}})));
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment peripheral1(V=V2_param) annotation(
    Placement(transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}})));
  Pharmacolibrary.Pharmacokinetic.NoPerfusedTissueCompartment peripheral2(V=V3_param) annotation(
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
    Documentation(info = "<html><head></head><body><h1>Schnider Propofol Model</h1><p>This model implements the Schnider PK model for propofol, accounting for lean body mass (LBM) using the James equation and adjusted body weight metrics to prevent overdosing in heavier patients. It provides highly accurate predictive capabilities targeting the effect site.</p></body></html>"),
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")));
end PropofolSchnider;
