import re

with open('/mnt/c/home/git/simtiva/PropofolModels.mo', 'r') as f:
    content = f.read()

models_to_copy = ['PropofolSchnider', 'PropofolPaedfusor', 'PropofolEleveld']
new_content = content

for model_name in models_to_copy:
    if f"model {model_name}_fixedUnits" in new_content:
        continue
        
    pattern = r'(model ' + model_name + r'\b.*?end ' + model_name + r';)'
    match = re.search(pattern, content, re.DOTALL)
    if not match: continue
    orig_text = match.group(1)
    
    new_text = orig_text.replace(f'model {model_name}', f'model {model_name}_fixedUnits')
    new_text = new_text.replace(f'end {model_name}', f'end {model_name}_fixedUnits')
    
    new_text = re.sub(r'central\(V=([^)]+)\)', r'central(V=\1 / 1000.0)', new_text)
    new_text = re.sub(r'peripheral1\(V=([^)]+)\)', r'peripheral1(V=\1 / 1000.0)', new_text)
    new_text = re.sub(r'peripheral2\(V=([^)]+)\)', r'peripheral2(V=\1 / 1000.0)', new_text)
    
    new_text = re.sub(r'transfer1\(CLa=([^,]+), CLb=([^)]+)\)', r'transfer1(CLa=\1 / 60000.0, CLb=\2 / 60000.0)', new_text)
    new_text = re.sub(r'transfer2\(CLa=([^,]+), CLb=([^)]+)\)', r'transfer2(CLa=\1 / 60000.0, CLb=\2 / 60000.0)', new_text)
    
    new_text = re.sub(r'elim\(CL=([^)]+)\)', r'elim(CL=\1 / 60000.0)', new_text)
    
    inf_orig = r'Modelica\.Blocks\.Sources\.RealExpression infusion_rate\(y = if time < 1\.0 then 2\.0 \* weight else if time < 60\.0 then 10\.0 \* weight / 60\.0 else 0\.0\)'
    inf_new = r'Modelica.Blocks.Sources.RealExpression infusion_rate(y = if time < 60.0 then ((2.0 * weight) / 60.0) * 1e-6 else if time < 3660.0 then ((10.0 * weight) / 3600.0) * 1e-6 else 0.0)'
    new_text = re.sub(inf_orig, inf_new, new_text)
    
    der_orig = r'der\(Ce\) = ke0 \* \(central\.cport\.c - Ce\);'
    der_new  = r'der(Ce) = (ke0 / 60.0) * (central.cport.c - Ce);'
    new_text = re.sub(der_orig, der_new, new_text)
    
    new_content = new_content.replace("end PropofolModels;", new_text + "\nend PropofolModels;")

with open('/mnt/c/home/git/simtiva/PropofolModels.mo', 'w') as f:
    f.write(new_content)
    
print("Appended new isolated variants flawlessly.")
