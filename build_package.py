import os

pkg_header = """package PropofolModels "A collection of Propofol Pharmacokinetic Models"
  extends Modelica.Icons.Package;
  annotation(
    uses(Modelica(version = "4.0.0"), Pharmacolibrary(version = "25.09")),
    Documentation(info = "<html><head></head><body><h1>Propofol Models</h1><p>These models are reimplementations of Propofol PK models, utilizing Pharmacolibrary and its standardized components.</p></body></html>")
  );

"""

pkg_footer = """
end PropofolModels;
"""

files = ['PropofolMarsh.mo', 'PropofolSchnider.mo', 'PropofolPaedfusor.mo', 'PropofolEleveld.mo']

def process_model(filepath):
    lines = open(filepath, 'r').read().splitlines()
    if lines[0].startswith('within'):
        lines = lines[1:]
        
    # Remove the `uses` annotation at the end of the model since the package handles it
    out_lines = []
    skip = False
    for line in lines:
        if 'uses(Modelica' in line:
            # Skip this line and remove trailing comma from previous if necessary
            if out_lines and out_lines[-1].strip().endswith(','):
                out_lines[-1] = out_lines[-1].rstrip(',')
            # Also remove the whole annotation block if it was just uses(...)
            # But the models have Documentation() too. We just replace uses(...) with empty string
            continue
        out_lines.append(line)
        
    content = "\\n".join(out_lines)
    # The models currently end with:
    # annotation(
    #   Documentation(...),
    #   uses(...)
    # );
    # We can just do a regex replacement
    import re
    content = re.sub(r',\\s*uses\(Modelica\(version = "4\.0\.0"\),\s*Pharmacolibrary\(version = "25\.09"\)\)', '', content)
    content = re.sub(r'uses\(Modelica\(version = "4\.0\.0"\),\s*Pharmacolibrary\(version = "25\.09"\)\)', '', content)
    return content

with open('PropofolModels.mo', 'w') as f:
    f.write(pkg_header)
    for model_file in files:
        f.write(process_model(os.path.join('PropofolModels', model_file)))
        f.write("\\n\\n")
    f.write(pkg_footer)

print("Generated PropofolModels.mo successfully.")
