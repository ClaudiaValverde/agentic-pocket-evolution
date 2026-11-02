// 3D viewers (3Dmol.js, vendored): protein page (pockets) and case page (co-folding, ancestral vs derived).
(function () {
  if (typeof $3Dmol === "undefined") return;

  function plddtColor(atom) {
    const b = atom.b;
    if (b > 90) return "#0053D6";
    if (b > 70) return "#65CBF3";
    if (b > 50) return "#FFDB13";
    return "#FF7D45";
  }

  // One viewer = one structure + optional highlighted residues {chain, resi[]} + optional ligand emphasis.
  function PEViewer(el, o) {
    o = o || {};
    this.el = el;
    this.v = $3Dmol.createViewer(el, { backgroundColor: "white" });
    this.msg = o.msg || null;
    this.plddt = o.plddt || null;   // checkbox
    this.surf = o.surf || null;     // checkbox
    this.current = null;            // {resi:[...], chain, center:[x,y,z]}
    this.near = !!o.nearLigand;     // show residues within 4.5 A of any HETATM ligand
    this.onload = o.onload || null;
    if (this.plddt) this.plddt.addEventListener("change", () => this.apply());
    if (this.surf) this.surf.addEventListener("change", () => this.apply());
  }
  PEViewer.prototype.say = function (t) { if (this.msg) this.msg.textContent = t || ""; };
  PEViewer.prototype.apply = function () {
    const v = this.v;
    const cartoon = this.plddt && this.plddt.checked ? { colorfunc: plddtColor } : { color: "#c8c8c8" };
    v.setStyle({}, { cartoon: cartoon });
    v.setStyle({ hetflag: true }, { stick: { colorscheme: "greenCarbon", radius: 0.22 } });
    v.removeAllShapes();
    v.removeAllSurfaces();
    if (this.near) {
      v.addStyle({ within: { distance: 4.5, sel: { hetflag: true } }, not: { hetflag: true } },
                 { stick: { colorscheme: "cyanCarbon", radius: 0.12 } });
    }
    const c = this.current;
    if (c) {
      const sel = { chain: c.chain || "A", resi: c.resi };
      v.addStyle(sel, { stick: { colorscheme: "orangeCarbon", radius: 0.2 } });
      if (c.center) {
        v.addSphere({ center: { x: c.center[0], y: c.center[1], z: c.center[2] }, radius: 1.2, color: "magenta", opacity: 0.8 });
      }
      if (this.surf && this.surf.checked) {
        v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.6, color: "orange" }, sel);
      }
    }
    v.render();
  };
  PEViewer.prototype.highlight = function (resi, chain, center, zoom) {
    this.current = resi && resi.length ? { resi: resi, chain: chain, center: center || null } : null;
    this.apply();
    if (this.current && zoom !== false) this.v.zoomTo({ chain: chain || "A", resi: resi });
    this.v.render();
  };
  PEViewer.prototype.load = function (url, fmt) {
    if (!url) return Promise.resolve();
    this.say("loading structure…");
    return fetch(url).then(r => { if (!r.ok) throw new Error(r.status + " " + r.statusText); return r.text(); })
      .then(txt => {
        this.v.clear();
        this.v.addModel(txt, fmt || "cif");
        this.current = null;
        this.apply();
        this.v.zoomTo();
        this.v.render();
        this.say("");
        if (this.onload) this.onload();
      })
      .catch(err => this.say("could not load structure: " + err));
  };
  PEViewer.prototype.focusLigand = function () {
    this.v.zoomTo({ hetflag: true });
    this.v.render();
  };

  // ---------------- protein page
  const el = document.getElementById("viewer");
  if (el && !document.getElementById("cofold-data")) {
    const pv = new PEViewer(el, { msg: document.getElementById("vmsg"), plddt: document.getElementById("plddt"),
                                  surf: document.getElementById("surf") });
    const rows = Array.from(document.querySelectorAll("tr.pocket"));
    function select(tr) {
      rows.forEach(x => x.classList.remove("sel"));
      tr.classList.add("sel");
      const resi = tr.dataset.residues.split(/\s+/).filter(Boolean).map(x => parseInt(x.split(":")[1], 10));
      const c = (tr.dataset.center || "").split(",").map(Number);
      pv.highlight(resi, "A", c.length === 3 && c.every(isFinite) ? c : null);
    }
    rows.forEach(tr => tr.addEventListener("click", () => select(tr)));
    const want = new URLSearchParams(location.search).get("pocket");
    pv.onload = function () {
      const tr = rows.find(r => r.dataset.pid === want) || rows[0];
      if (tr) select(tr);
    };
    const reset = document.getElementById("reset");
    if (reset) reset.addEventListener("click", () => {
      rows.forEach(x => x.classList.remove("sel"));
      pv.highlight(null); pv.v.zoomTo(); pv.v.render();
    });
    pv.load(el.dataset.url);
  }

  // ---------------- case page: two panels (ancestral | derived), each with a model selector
  const dataEl = document.getElementById("cofold-data");
  if (dataEl) {
    const models = JSON.parse(dataEl.textContent);
    document.querySelectorAll(".cofpanel").forEach(panel => {
      const role = panel.dataset.role;
      const sel = panel.querySelector("select");
      const info = panel.querySelector(".cofinfo");
      const pv = new PEViewer(panel.querySelector(".viewer"), {
        msg: panel.querySelector(".vmsg"), plddt: panel.querySelector(".plddt"), nearLigand: true });
      function show(m) {
        if (!m) { info.textContent = "no model"; return; }
        const mt = Object.entries(m.metrics || {}).map(([k, x]) => k + " = " + x).join(" · ");
        info.textContent = m.model_id + " · " + (m.accession || "") + " · ligand " + (m.ligand || "?") +
          " · seed " + (m.seed == null ? "?" : m.seed) + (mt ? " · " + mt : "");
        pv.onload = function () {
          if (m.residues && m.residues.length) pv.highlight(m.residues, "A", null, false);
          if (panel.querySelector(".lig").checked) pv.focusLigand();
        };
        pv.load(m.url);
      }
      const idx = models.findIndex(m => m.role === role);
      models.forEach((m, i) => {
        const o = document.createElement("option");
        o.value = i; o.textContent = m.role + " · " + m.model_id + (m.ligand ? " · " + m.ligand : "");
        sel.appendChild(o);
      });
      sel.addEventListener("change", () => show(models[+sel.value]));
      panel.querySelector(".lig").addEventListener("change", () => { if (panel.querySelector(".lig").checked) pv.focusLigand(); else pv.v.zoomTo(); pv.v.render(); });
      if (idx >= 0) { sel.value = idx; show(models[idx]); }
      else { sel.value = ""; info.textContent = "no " + role + " model recorded in case.yaml"; }
    });
  }
})();
