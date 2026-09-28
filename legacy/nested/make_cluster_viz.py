if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Generate a self-contained HTML cluster visualization of disease likelihood.

Clusters = groups of similar diseases (by body system / mechanism).
Brightness/glow = likelihood (bright -> likely, dim -> unlikely).
Bubble size = relative-probability index. Hue = cluster.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "analysis_reports" / "disease_likelihood.json"
OUT = REPO_ROOT / "analysis_reports" / "disease_clusters.html"

diseases = json.loads(DATA.read_text(encoding="utf-8"))

# Hue per group for the "grouped by similarity" encoding.
GROUP_HUES = {
    "Cancers": 350,
    "Cardiometabolic": 25,
    "Neuro / Psychiatric": 275,
    "Autoimmune / Inflammatory": 200,
    "Eye / Sensory": 55,
    "Metabolic / Digestive": 140,
    "Infectious / Immune susceptibility": 315,
    "Other": 220,
}

html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Disease Likelihood Clusters</title>
<style>
  :root { color-scheme: dark; }
  * { box-sizing: border-box; }
  body {
    margin: 0; background: radial-gradient(ellipse at 50% 0%, #0d1326 0%, #05070f 60%, #03040a 100%);
    color: #e8ecf5; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    min-height: 100vh; padding: 28px 20px 60px;
  }
  h1 { font-size: 22px; font-weight: 650; margin: 0 0 4px; letter-spacing: .2px; }
  .sub { color: #8b93a7; font-size: 13px; margin: 0 0 2px; }
  .disclaimer { color: #6a7governments; }
  .disclaimer { color: #7a8296; font-size: 11.5px; margin: 6px 0 18px; font-style: italic; }
  .legend { display: flex; flex-wrap: wrap; gap: 18px; align-items: center; margin: 0 0 18px;
            font-size: 12px; color: #aab2c5; }
  .legend .chip { display:flex; align-items:center; gap:7px; }
  .grad { width: 120px; height: 12px; border-radius: 7px;
          background: linear-gradient(90deg,#2a3042,#3d5a80,#5fb0ff,#8fe0ff,#ffffff); }
  .dot { width:12px;height:12px;border-radius:50%; }
  #wrap { max-width: 1180px; margin: 0 auto; }
  svg { width: 100%; height: auto; display:block; }
  .clabel { fill: #cfd6e6; font-size: 13px; font-weight: 600; }
  .csub { fill: #79melon; }
  .csub { fill: #7b8398; font-size: 10.5px; }
  .region { fill: rgba(255,255,255,0.018); stroke: rgba(255,255,255,0.06); stroke-width:1; }
  .node text { fill: #f2f5fb; font-size: 9.5px; pointer-events:none; text-anchor: middle; }
  .node { cursor: default; }
  #tip {
    position: fixed; pointer-events:none; opacity:0; transition: opacity .12s;
    background: rgba(10,14,26,.97); border: 1px solid rgba(120,150,220,.35);
    border-radius: 10px; padding: 10px 12px; font-size: 12px; max-width: 260px;
    box-shadow: 0 8px 30px rgba(0,0,0,.5); z-index: 10;
  }
  #tip b { color:#bcd4ff; }
  #tip .row { color:#9aa4bb; margin-top:3px; }
  .bar { height:5px;border-radius:3px;background:#232a3d;margin-top:2px;overflow:hidden; }
  .bar > i { display:block;height:100%; }
</style>
</head>
<body>
<div id="wrap">
  <h1>Disease Risk — Similarity Clusters</h1>
  <p class="sub">Cross-checked merge of both analysis runs. Brightness = real-world <b>Probability</b> where the evidence supports it; solid glowing bubbles are cross-checked estimates.</p>
  <p class="disclaimer">Educational / research use only — approximate single-marker estimates, NOT clinical risk or a diagnosis.</p>
  <div class="legend">
    <div class="chip"><span>Dim → Bright =</span><span class="grad"></span><span>lower → higher real-world probability</span></div>
    <div class="chip"><span class="dot" style="background:#5fb0ff;box-shadow:0 0 6px #5fb0ff"></span>solid = cross-checked real-world estimate</div>
    <div class="chip"><span class="dot" style="background:transparent;border:1.5px dashed #7b8398"></span>dashed = low-confidence only (no real-world estimate)</div>
    <div class="chip"><span class="dot" style="background:#e05a7a"></span>color = disease group · size = genetic lean</div>
  </div>
  <svg id="viz" viewBox="0 0 1180 900" xmlns="http://www.w3.org/2000/svg"></svg>
</div>
<div id="tip"></div>
<script>
const DATA = __DATA__;
const HUES = __HUES__;

// group -> nodes
const groups = {};
DATA.forEach(d => { (groups[d.group] = groups[d.group] || []).push(d); });
const groupNames = Object.keys(groups).sort((a,b)=>{
  const ma = Math.max(...groups[a].map(x=>x.likelihood));
  const mb = Math.max(...groups[b].map(x=>x.likelihood));
  return mb - ma;
});

const svg = document.getElementById('viz');
const SVGNS = 'http://www.w3.org/2000/svg';
const W = 1180;
const cols = 3;
const cellW = W / cols;
const cellH = 300;
const tip = document.getElementById('tip');

function el(name, attrs){ const e=document.createElementNS(SVGNS,name); for(const k in attrs) e.setAttribute(k,attrs[k]); return e; }

// glow filter
const defs = el('defs',{});
[['soft',3.2],['strong',5.5]].forEach(([id,dev])=>{
  const f = el('filter',{id:'glow-'+id,x:'-80%',y:'-80%',width:'260%',height:'260%'});
  f.appendChild(el('feGaussianBlur',{stdDeviation:dev,result:'b'}));
  const m=el('feMerge',{}); m.appendChild(el('feMergeNode',{in:'b'})); m.appendChild(el('feMergeNode',{in:'SourceGraphic'}));
  f.appendChild(m); defs.appendChild(f);
});
svg.appendChild(defs);

function radius(d){ return 12 + (d.genetic_lean||0) * 0.34; }

// Brightness 0..1. Real-world nodes scale by probability (sqrt, ref 30%).
// Low-confidence-only nodes stay deliberately dim and are drawn dashed.
const MAXP = 30;
function brightness(d){
  if(d.has_real_world){
    return Math.max(0.14, Math.min(1, Math.sqrt((d.probability_real_world||0)/MAXP)));
  }
  return 0.08 + (d.likelihood_low_conf||0)*0.34;   // capped ~0.30, always dimmer
}

// Deterministic packed layout per cluster (spiral + relaxation).
function layout(nodes, cx, cy, maxR){
  nodes.sort((a,b)=>brightness(b)-brightness(a));
  const placed=[];
  nodes.forEach((n,i)=>{
    const r=radius(n);
    let best=null;
    for(let k=0;k<420;k++){
      const ang=k*2.399963; const rad=6*Math.sqrt(k);
      const x=cx+rad*Math.cos(ang), y=cy+rad*Math.sin(ang);
      let ok=true;
      for(const p of placed){ const dx=x-p.x,dy=y-p.y; if(Math.hypot(dx,dy) < r+p.r+3){ok=false;break;} }
      if(ok){ best={x,y,r}; break; }
      if(!best) best={x,y,r};
    }
    n._x=best.x; n._y=best.y; n._r=r; placed.push({x:best.x,y:best.y,r});
  });
}

groupNames.forEach((gname, gi)=>{
  const col = gi % cols, row = Math.floor(gi/cols);
  const cx = col*cellW + cellW/2;
  const cy = row*cellH + cellH/2 + 34;
  const nodes = groups[gname];
  const hue = HUES[gname] ?? 220;

  // region backdrop
  svg.appendChild(el('rect',{x:col*cellW+14,y:row*cellH+14,width:cellW-28,height:cellH-20,rx:16,class:'region'}));
  const lbl=el('text',{x:col*cellW+26,y:row*cellH+40,class:'clabel'}); lbl.textContent=gname; svg.appendChild(lbl);
  const sub=el('text',{x:col*cellW+26,y:row*cellH+56,class:'csub'});
  const rw=nodes.filter(n=>n.has_real_world);
  const avg=rw.length? (rw.reduce((s,n)=>s+n.probability_real_world,0)/rw.length).toFixed(1)+'% avg prob'
                     : 'no real-world estimate';
  sub.textContent = nodes.length+' diseases · '+avg; svg.appendChild(sub);

  layout(nodes, cx, cy+6, cellW/2-30);

  nodes.forEach(n=>{
    const B = brightness(n);                    // 0..1 brightness driver
    const rw = n.has_real_world;
    const light = 24 + B*66;                     // lightness %
    const sat = rw ? (34 + B*58) : 16;           // low-conf nodes desaturated
    const op = rw ? (0.34 + B*0.66) : 0.30;
    const fill = `hsl(${hue}, ${sat}%, ${light}%)`;
    const g = el('g',{class:'node'});
    // outer glow only for cross-checked real-world nodes
    if(rw && B>0.30){
      g.appendChild(el('circle',{cx:n._x,cy:n._y,r:n._r, fill:fill, opacity:(B*0.9).toFixed(2),
                                 filter: B>0.55?'url(#glow-strong)':'url(#glow-soft)'}));
    }
    const core = el('circle',{cx:n._x,cy:n._y,r:n._r, fill: rw?fill:'rgba(255,255,255,0.03)',
                              opacity:op.toFixed(2),
                              stroke:`hsl(${hue},${rw?sat:30}%,${Math.min(92,light+16)}%)`,
                              'stroke-width': rw?(0.6+B*1.8).toFixed(2):1.3,
                              'stroke-opacity': rw?op.toFixed(2):0.7});
    if(!rw){ core.setAttribute('stroke-dasharray','3 3'); }
    g.appendChild(core);
    // label if bubble is big enough
    if(n._r>19){
      const words = n.trait.replace(/\\s*\\(.*?\\)/,'').trim();
      const t=el('text',{x:n._x,y:n._y+3}); t.textContent = words.length>16?words.slice(0,15)+'…':words;
      t.setAttribute('opacity', (rw && B>0.3)?0.96:0.55); g.appendChild(t);
    }
    g.addEventListener('mousemove', ev=>{
      tip.style.opacity=1; tip.style.left=(ev.clientX+14)+'px'; tip.style.top=(ev.clientY+14)+'px';
      const ev_s = n.evidence_sources.length? n.evidence_sources.join(', ') : 'no cached external evidence';
      let probLine;
      if(rw){
        probLine = `<div class="row" style="color:#bcd4ff">Real-world probability: <b>${n.probability_real_world}%</b></div>`+
          `<div class="bar"><i style="width:${Math.min(100,n.probability_real_world*3).toFixed(0)}%;background:${fill}"></i></div>`+
          `<div class="row">baseline ${(n.baseline_prevalence*100).toFixed(1)}% × genotype ×${n.genotype_multiplier}</div>`+
          `<div class="row">databanks: ${n.databanks.join('+')}</div>`;
      } else {
        probLine = `<div class="row" style="color:#c7a06a">No real-world estimate</div>`+
          `<div class="row" style="font-size:11px">${n.probability_note}</div>`;
      }
      tip.innerHTML = `<b>${n.trait}</b><div class="row">Group: ${n.group}</div>`+
        probLine+
        `<div class="row" style="margin-top:5px">Likelihood (low-conf): ${(n.likelihood_low_conf*100).toFixed(1)}% · Genetic lean: ${n.genetic_lean}%</div>`+
        `<div class="row">Evidence: ${ev_s}</div>`;
    });
    g.addEventListener('mouseleave', ()=> tip.style.opacity=0);
    svg.appendChild(g);
  });
});

// size viewBox to rows used
const rowsUsed = Math.ceil(groupNames.length/cols);
svg.setAttribute('viewBox', `0 0 ${W} ${rowsUsed*cellH+60}`);
</script>
</body>
</html>
"""

html = html.replace("__DATA__", json.dumps(diseases)).replace("__HUES__", json.dumps(GROUP_HUES))
# fix two intentional-looking typos in CSS placeholders
html = html.replace("#6a7governments", "#7a8296").replace("#79melon", "#7b8398")
OUT.write_text(html, encoding="utf-8")
print("Wrote", OUT, "|", len(diseases), "diseases in", len(GROUP_HUES), "groups")
