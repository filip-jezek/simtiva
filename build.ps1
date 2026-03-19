$header = @"
package PropofolModels "A collection of Propofol Pharmacokinetic Models"
  extends Modelica.Icons.Package;
  annotation(
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")),
    Documentation(info = "<html><head></head><body><h1>Propofol Models</h1><p>These models are reimplementations of Propofol PK models, utilizing Pharmacolibrary and its standardized components. Contains Marsh, Eleveld, Paedfusor, and Schnider.</p></body></html>")
  );
"@
$footer = "`nend PropofolModels;"

$marsh = (Get-Content "PropofolModels\PropofolMarsh.mo" | Where-Object { $_ -notmatch '^within' -and $_ -notmatch 'uses\(Modelica' }) -join "`n"
$schnider = (Get-Content "PropofolModels\PropofolSchnider.mo" | Where-Object { $_ -notmatch '^within' -and $_ -notmatch 'uses\(Modelica' }) -join "`n"
$paedfusor = (Get-Content "PropofolModels\PropofolPaedfusor.mo" | Where-Object { $_ -notmatch '^within' -and $_ -notmatch 'uses\(Modelica' }) -join "`n"
$eleveld = (Get-Content "PropofolModels\PropofolEleveld.mo" | Where-Object { $_ -notmatch '^within' -and $_ -notmatch 'uses\(Modelica' }) -join "`n"

$content = $header + "`n" + $marsh + "`n" + $schnider + "`n" + $paedfusor + "`n" + $eleveld + $footer
Set-Content -Path "PropofolModels.mo" -Value $content
Write-Output "Successfully built PropofolModels.mo"
