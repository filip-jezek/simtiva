import http.server
import time
import json
import threading
import subprocess
import os
import sys

port = 8085
result_data = None

class Handler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass
    def do_POST(self):
        global result_data
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length)
        result_data = json.loads(post_data.decode('utf-8'))
        self.send_response(200)
        self.end_headers()
        self.wfile.write(b"OK")
        threading.Thread(target=self.server.shutdown).start()

os.chdir(r"c:\home\git\simtiva")

html_content = """
<!DOCTYPE html>
<html>
<body>
<script>
    let out = {};
    try {
        // Stub global DOM
        window.document.getElementById = function(id) {
            return { style: {display:''}, innerHTML: '', value: 0 };
        };
        // Provide mock global deps
        window.rnd3 = function(x) { return Math.round(x*1000)/1000; };
        window.displayWarning = function(){};
        window.useAdjBW = 0; window.AdjBW = 0; window.mass = 70; window.age = 40; 
        window.height = 170; window.lbm = 60; window.gender = 0; window.paedi_mode = 0; 
        window.opioid = 0; window.active_drug_set_index=0;
        window.drug_sets = [];
        window.fageing = function(){return 1.0;};
        window.fclmaturation = function(){return 1.0;};
        window.fcentral = function(x){return x;};
        window.fq3maturation = function(){return 1.0;};
        window.fffm = function(){return 1.0;};
        
        let script = document.createElement('script');
        script.src = 'pharmacology.js';
        script.onload = function() {
            try {
                readmodel("Marsh", 0);
                readmodel("Schnider", 1);
                out.marsh = {
                    vc: drug_sets[0].vc, k10: drug_sets[0].k10, k12: drug_sets[0].k12,
                    k13: drug_sets[0].k13, k21: drug_sets[0].k21, k31: drug_sets[0].k31
                };
                out.schnider = {
                    vc: drug_sets[1].vc, k10: drug_sets[1].k10, k12: drug_sets[1].k12,
                    k13: drug_sets[1].k13, k21: drug_sets[1].k21, k31: drug_sets[1].k31
                };
                
                // test cube
                window.r = [0,0,0,0];
                cube(1, -6, 11, -6); // roots should be 1, 2, 3
                out.cube_test = [r[1], r[2], r[3]];
                
                fetch('http://127.0.0.1:8085', {method: 'POST', body: JSON.stringify(out)});
            } catch (e) {
                fetch('http://127.0.0.1:8085', {method: 'POST', body: JSON.stringify({error: e.toString(), stack: e.stack})});
            }
        };
        document.body.appendChild(script);
    } catch(err) {
        fetch('http://127.0.0.1:8085', {method: 'POST', body: JSON.stringify({error: err.toString()})});
    }
</script>
</body>
</html>
"""
with open("runner.html", "w", encoding="utf-8") as f:
    f.write(html_content)

server = http.server.HTTPServer(('127.0.0.1', port), Handler)
t = threading.Thread(target=server.serve_forever)
t.start()

print("Launching msedge headless...")
subprocess.Popen(["msedge", "--headless", f"http://127.0.0.1:{port}/runner.html"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

t.join(timeout=10)
try:
    if result_data:
        print("\nSUCCESS! JS BEHAVIOR CAPTURED:")
        print(json.dumps(result_data, indent=2))
        
        # Compare with python models
        from PKPD_Reimplementation.models import MarshModel, SchniderModel
        from PKPD_Reimplementation.core_solvers import run_algebraic_solver
        
        marsh_py = MarshModel(70, 40, 170, 0).get_parameters()
        schnider_py = SchniderModel(70, 40, 170, 0).get_parameters()
        
        print("\nPYTHON PREDICTIONS:")
        print("Marsh:", {k: marsh_py[k] for k in ['vc','k10','k12','k13','k21','k31']})
        print("Schnider:", {k: schnider_py[k] for k in ['vc','k10','k12','k13','k21','k31']})
        
        # cube
        import numpy as np
        # run_algebraic solver internally solves roots for the analytical equation.
        # we can just use np.roots for cube(1, -6, 11, -6)
        py_roots = np.sort(np.roots([1, -6, 11, -6])).real.tolist()
        print(f"Cube Roots [Python np.roots vs JS cube(1,-6,11,-6)] \nPython: {py_roots}\nJS: {result_data.get('cube_test')}")

    else:
        print("Timeout! The headless browser did not respond.")
except Exception as e:
    import traceback
    traceback.print_exc()

server.server_close()
