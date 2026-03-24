import subprocess
import json
import numpy as np

result = subprocess.run(["node", "test_js.js"], capture_output=True, text=True)
if result.returncode != 0:
    print("JS Error:", result.stderr)
    exit(1)

try:
    js_data = json.loads(result.stdout)
except Exception as e:
    print("Failed to parse JS output:", result.stdout)
    exit(1)

from models.marsh import MarshModel
from models.schnider import SchniderModel

marsh_py = MarshModel(weight=70, age=40, height=170, gender=0).get_parameters()
schnider_py = SchniderModel(weight=70, age=40, height=170, gender=0).get_parameters()

def compare(name, dict_js, dict_py):
    print(f"\\n--- Comparing {name} Parameters ---")
    keys = ['vc', 'k10', 'k12', 'k13', 'k21', 'k31']
    for k in keys:
        js_val = dict_js[k]
        py_val = dict_py[k]
        match = abs(js_val - py_val) < 1e-6
        print(f"  {k}: JS = {js_val:.5f} | PY = {py_val:.5f} | Match: {'YES' if match else 'NO'}")

compare("Marsh", js_data["marsh"], marsh_py)
compare("Schnider", js_data["schnider"], schnider_py)

marsh_k = marsh_py
k10, k12, k13, k21, k31 = marsh_k['k10'], marsh_k['k12'], marsh_k['k13'], marsh_k['k21'], marsh_k['k31']
a0 = k10 * k21 * k31
a1 = k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
a2 = k10 + k12 + k13 + k21 + k31
py_roots = np.sort(np.roots([1, a2, a1, a0])).real.tolist()

# The JS sort places them descending magnitude (less negative to more negative) because r.sort() ascending on negative numbers means -10, -5, -1, then it reverses them to r[3]=r[2], r[2]=r[1], r[1]=r[0] meaning they become -1, -5, -10. 
# We'll just compare bounding structurally by sorting descending globally over absolute values:
py_roots = sorted(np.abs(py_roots).tolist(), reverse=True)
js_roots = sorted(np.abs(js_data["cube"]).tolist(), reverse=True)

print(f"\\n--- Cube Roots Validation ---")
print(f"  Python analytical np.roots (absolute magnitudes): {py_roots}")
print(f"  JS native 'cube()' roots (absolute magnitudes):   {js_roots}")
print(f"  Match exact parity:         {'YES' if np.allclose(py_roots, js_roots, atol=1e-5) else 'NO'}")
