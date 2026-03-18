within ;
package PropofolModels "A collection of Propofol Pharmacokinetic Models"
  extends Modelica.Icons.Package;
  annotation(
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")),
    Documentation(info = "<html><head></head><body>
      <h1>Propofol Models</h1>
      These models are reimplementations of Propofol PK models from Python, utilizing Pharmacolibrary and its standardized components. Models included:
      <ul>
        <li>Marsh</li>
        <li>Schnider</li>
        <li>Paedfusor</li>
        <li>Eleveld</li>
      </ul>
    </body></html>")
  );
end PropofolModels;
