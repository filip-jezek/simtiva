within PropofolModels;
model PropofolEleveld
  parameter Real weight = 70 "Weight in kg";
  parameter Real age = 40 "Age in years";
  parameter Real height = 170 "Height in cm";
  parameter Integer gender = 0 "0=Male, 1=Female";
  parameter Boolean opioid = true "Opioid co-administration";

  parameter Real toweeks = 52.1429;
  parameter Real pma = age * toweeks + 40;
  
  parameter Real fcentral_w = weight / (weight + 33.6);
  parameter Real fcentral_70 = 70.0 / (70.0 + 33.6);
  parameter Real vc = 6.28 * fcentral_w / fcentral_70;

  parameter Real v2 = 25.5 * (weight / 70.0) * exp((-0.0156) * (age - 35.0));
  
  parameter Real bmi = weight / (height / 100.0)^2;
  parameter Real ffm_male_factor = 0.88 + 0.12 / (1.0 + 1.0 / (age / 13.4)^12.7);
  parameter Real ffm_female_factor = 1.11 - 0.11 / (1.0 + 1.0 / (age / 7.1)^1.1);
  parameter Real ffm = if gender == 0 then ffm_male_factor * ((9270.0 * weight) / (6680.0 + 216.0 * bmi)) else ffm_female_factor * ((9270.0 * weight) / (8780.0 + 244.0 * bmi));
  
  parameter Real ffmref_factor = 0.88 + 0.12 / (1.0 + 1.0 / (35.0 / 13.4)^12.7);
  parameter Real ffmref = ffmref_factor * ((9270.0 * 70.0) / (6680.0 + 216.0 * 24.22145));
  
  parameter Real ext_age = exp((-0.0138) * age);
  parameter Real v3 = if opioid then 273.0 * ffm / ffmref * ext_age else 273.0 * ffm / ffmref;

  parameter Real fq3_age = (age * toweeks + 40.0) / (age * toweeks + 40.0 + 68.3);
  parameter Real fq3_35 = (35.0 * toweeks + 40.0) / (35.0 * toweeks + 40.0 + 68.3);
  
  parameter Real cl2 = 1.75 * (v2 / 25.5)^0.75 * (1.0 + 1.3 * (1.0 - fq3_age));
  parameter Real cl3 = 1.11 * (v3 / 273.0)^0.75 * (fq3_age / fq3_35);

  parameter Real fcl_pma = 1.0 / (1.0 + (42.3 / pma)^9.06);
  parameter Real pma_ref = 35.0 * toweeks + 40.0;
  parameter Real fcl_ref = 1.0 / (1.0 + (42.3 / pma_ref)^9.06);
  parameter Real mat_factor = fcl_pma / fcl_ref;
  parameter Real age_decay = exp((-0.00286) * age);

  parameter Real cl1_base = if gender == 0 then 1.79 else 2.1;
  parameter Real opioid_factor = if opioid then age_decay else 1.0;
  parameter Real cl1 = cl1_base * (weight / 70.0)^0.75 * mat_factor * opioid_factor;

  parameter Real ke0 = 0.146 / (weight / 70.0)^0.25;

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
    Documentation(info = "<html><head></head><body><h1>Eleveld Propofol Model</h1><p>This model implements the Eleveld PK model for propofol, accommodating a wide range of patients including children, adults, obese individuals, and the elderly. It adjusts clearances and volumes dynamically based on covariates including age, weight, height, gender, and opioid co-administration.</p></body></html>"),
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")));
end PropofolEleveld;
