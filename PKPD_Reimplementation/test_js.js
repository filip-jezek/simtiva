const fs = require('fs');

global.document = {
    getElementById: function(id) {
        return { style: {display:''}, innerHTML: '', value: 0 };
    }
};
global.rnd3 = function(x) { return Math.round(x*1000)/1000; };
global.displayWarning = function(){};
global.useAdjBW = 0; global.AdjBW = 0; global.mass = 70; global.age = 40; 
global.height = 170; global.gender = 0; global.paedi_mode = 0; 
global.lbm = 1.1 * 70 - 128 * Math.pow(70/170, 2); // James LBM correctly mocked
global.opioid = 0; global.active_drug_set_index = 0;
global.drug_sets = [];
global.fageing = function(){return 1.0;};
global.fclmaturation = function(){return 1.0;};
global.fcentral = function(x){return x;};
global.fq3maturation = function(){return 1.0;};
global.fffm = function(){return 1.0;};
global.parseloading = function(){return 1.0;};

const script = fs.readFileSync('../pharmacology.js', 'utf8');
eval(script);

try {
    let out = {};
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
    
    // Test cube root finder
    let m = out.marsh;
    let my_r = new Array(4);
    cube(m.k10, m.k12, m.k21, m.k13, m.k31, my_r);
    out.cube = [my_r[1], my_r[2], my_r[3]];
    console.log(JSON.stringify(out));
} catch (e) {
    console.error(e);
}
