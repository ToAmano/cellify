import os
import tempfile
import webbrowser

from pymatgen.core import Structure

VIEWER_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>cellify 3D Structure Viewer</title>
  <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600&display=swap" rel="stylesheet">
  <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3dmol/2.4.2/3Dmol-min.js"></script>
  <style>
    * {
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }
    body {
      font-family: 'Outfit', sans-serif;
      background-color: #121214;
      color: #e1e1e6;
      overflow: hidden;
      height: 100vh;
      width: 100vw;
    }
    #gdis {
      width: 100%;
      height: 100vh;
      position: absolute;
      top: 0;
      left: 0;
      z-index: 1;
    }
    .panel {
      position: absolute;
      top: 20px;
      left: 20px;
      z-index: 10;
      width: 320px;
      max-height: calc(100vh - 40px);
      background: rgba(22, 22, 26, 0.75);
      backdrop-filter: blur(12px);
      -webkit-backdrop-filter: blur(12px);
      border: 1px solid rgba(255, 255, 255, 0.08);
      border-radius: 16px;
      padding: 24px;
      display: flex;
      flex-direction: column;
      gap: 20px;
      overflow-y: auto;
      box-shadow: 0 10px 30px rgba(0, 0, 0, 0.5);
      transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }
    .panel:hover {
      border-color: rgba(255, 255, 255, 0.15);
      background: rgba(22, 22, 26, 0.85);
    }
    h1 {
      font-size: 1.25rem;
      font-weight: 600;
      letter-spacing: -0.02em;
      background: linear-gradient(135deg, #fff 0%, #a5a5b0 100%);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
    }
    .section {
      border-top: 1px solid rgba(255, 255, 255, 0.08);
      padding-top: 16px;
    }
    .section-title {
      font-size: 0.85rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: #8c8c99;
      margin-bottom: 12px;
    }
    .meta-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 10px;
      font-size: 0.9rem;
    }
    .meta-item {
      display: flex;
      flex-direction: column;
    }
    .meta-label {
      color: #8c8c99;
      font-size: 0.75rem;
    }
    .meta-value {
      font-weight: 400;
      color: #ffffff;
    }
    .inspector-card {
      background: rgba(255, 255, 255, 0.03);
      border: 1px solid rgba(255, 255, 255, 0.05);
      border-radius: 10px;
      padding: 12px;
      font-size: 0.9rem;
      line-height: 1.4;
    }
    .control-row {
      display: flex;
      align-items: center;
      justify-content: space-between;
      font-size: 0.9rem;
      margin-bottom: 10px;
    }
    select, button {
      background: rgba(255, 255, 255, 0.05);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 8px;
      color: #fff;
      padding: 8px 12px;
      font-family: inherit;
      font-size: 0.85rem;
      cursor: pointer;
      outline: none;
      transition: all 0.2s;
    }
    select:hover, button:hover {
      background: rgba(255, 255, 255, 0.1);
      border-color: rgba(255, 255, 255, 0.2);
    }
    button {
      font-weight: 600;
      width: 100%;
      text-align: center;
    }
    .neighbor-list {
      margin-top: 8px;
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .neighbor-item {
      display: flex;
      justify-content: space-between;
      font-size: 0.8rem;
      color: #a5a5b0;
    }
    .toggle-switch {
      position: relative;
      display: inline-block;
      width: 38px;
      height: 20px;
    }
    .toggle-switch input {
      opacity: 0;
      width: 0;
      height: 0;
    }
    .slider {
      position: absolute;
      cursor: pointer;
      top: 0;
      left: 0;
      right: 0;
      bottom: 0;
      background-color: rgba(255, 255, 255, 0.1);
      transition: .3s;
      border-radius: 20px;
    }
    .slider:before {
      position: absolute;
      content: "";
      height: 14px;
      width: 14px;
      left: 3px;
      bottom: 3px;
      background-color: white;
      transition: .3s;
      border-radius: 50%;
    }
    input:checked + .slider {
      background-color: #4f46e5;
    }
    input:checked + .slider:before {
      transform: translateX(18px);
    }
    .watermark {
      position: absolute;
      bottom: 20px;
      right: 20px;
      z-index: 10;
      pointer-events: none;
      opacity: 0.4;
      font-size: 0.9rem;
      font-weight: 300;
      letter-spacing: 0.1em;
      text-transform: uppercase;
    }
  </style>
</head>
<body>
  <div id="gdis"></div>
  <div class="watermark">cellify WebGL viewer</div>

  <div class="panel">
    <div>
      <h1>cellify structure</h1>
      <div style="font-size: 0.8rem; color: #8c8c99; margin-top: 4px;">{{FORMULA}} • {{NUM_ATOMS}} atoms</div>
    </div>

    <!-- Metadata Section -->
    <div class="section">
      <div class="section-title">Structure Info</div>
      <div class="meta-grid">
        <div class="meta-item">
          <span class="meta-label">Volume</span>
          <span class="meta-value">{{VOLUME}} Å³</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Density</span>
          <span class="meta-value">{{DENSITY}} g/cm³</span>
        </div>
      </div>
      <div class="meta-grid" style="margin-top: 10px;">
        <div class="meta-item">
          <span class="meta-label">Lattice (a, b, c)</span>
          <span class="meta-value">{{LATTICE_ABC}} Å</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Angles (α, β, γ)</span>
          <span class="meta-value">{{LATTICE_ANGLES}}°</span>
        </div>
      </div>
    </div>

    <!-- Inspector Section -->
    <div class="section">
      <div class="section-title">Atom Inspector</div>
      <div id="inspector-content" class="inspector-card">
        Click an atom in the 3D viewer to inspect its details and nearest neighbors.
      </div>
    </div>

    <!-- Controls Section -->
    <div class="section">
      <div class="section-title">Visual Controls</div>

      <div class="control-row">
        <span>Unit Cell Grid</span>
        <label class="toggle-switch">
          <input type="checkbox" id="toggle-cell" checked>
          <span class="slider"></span>
        </label>
      </div>

      <div class="control-row">
        <span>Orthographic Camera</span>
        <label class="toggle-switch">
          <input type="checkbox" id="toggle-camera">
          <span class="slider"></span>
        </label>
      </div>

      <div class="control-row">
        <span>Display Style</span>
        <select id="select-style">
          <option value="ball-stick">Ball & Stick</option>
          <option value="spacefill">Spacefill</option>
          <option value="stick">Sticks</option>
        </select>
      </div>
    </div>

    <div class="section">
      <button id="btn-screenshot">Capture Screenshot</button>
    </div>
  </div>

  <script>
    let viewer = null;
    let model = null;
    let selectedRep = null;
    let unitCellRep = null;

    const cifData = `{{CIF_DATA}}`;

    $(document).ready(function() {
      // 1. Initialize 3Dmol viewer
      viewer = $3Dmol.createViewer($("#gdis"), {
        backgroundColor: "#121214"
      });

      // 2. Load model
      model = viewer.addModel(cifData, "cif");

      // 3. Set default style
      applyStyle("ball-stick");

      // 4. Add unit cell box
      if (model.getAtoms().length > 0) {
        unitCellRep = viewer.addUnitCell(model, {
          box: { color: "#ffffff", opacity: 0.7, lineThickness: 1.0 }
        });
      }

      viewer.zoomTo();
      viewer.render();

      // 5. Setup interaction
      viewer.setClickable({}, true, function(atom, viewer, event, container) {
        if (!atom) return;
        inspectAtom(atom);
      });

      // 6. Handle UI controls
      $("#toggle-cell").change(function() {
        if (this.checked) {
          if (unitCellRep === null && model.getAtoms().length > 0) {
            unitCellRep = viewer.addUnitCell(model, {
              box: { color: "#ffffff", opacity: 0.7, lineThickness: 1.0 }
            });
          }
        } else {
          if (unitCellRep !== null) {
            viewer.removeUnitCell(model);
            unitCellRep = null;
          }
        }
        viewer.render();
      });

      $("#toggle-camera").change(function() {
        if (this.checked) {
          viewer.setCameraParameters({ orthographic: true });
        } else {
          viewer.setCameraParameters({ orthographic: false });
        }
        viewer.render();
      });

      $("#select-style").change(function() {
        applyStyle(this.value);
      });

      $("#btn-screenshot").click(function() {
        let imgData = viewer.png();
        let link = document.createElement("a");
        link.download = "{{FORMULA}}_structure.png";
        link.href = "data:image/png;base64," + imgData;
        link.click();
      });
    });

    function applyStyle(styleName) {
      if (styleName === "ball-stick") {
        viewer.setStyle({}, {
          sphere: { scale: 0.25, colorscheme: "Jmol" },
          stick: { radius: 0.1, colorscheme: "Jmol" }
        });
      } else if (styleName === "spacefill") {
        viewer.setStyle({}, {
          sphere: { scale: 0.8, colorscheme: "Jmol" }
        });
      } else if (styleName === "stick") {
        viewer.setStyle({}, {
          stick: { radius: 0.15, colorscheme: "Jmol" }
        });
      }
      if (selectedRep) {
        viewer.removeRepresentation(selectedRep);
        selectedRep = null;
      }
      viewer.render();
    }

    function inspectAtom(atom) {
      if (selectedRep) {
        viewer.removeRepresentation(selectedRep);
      }
      selectedRep = viewer.addRepresentation("sphere", {
        sel: { index: atom.index },
        color: "#f59e0b",
        opacity: 0.6,
        scale: 1.3
      });
      viewer.render();

      const atoms = model.getAtoms();
      let neighbors = [];
      for (let i = 0; i < atoms.length; i++) {
        let a = atoms[i];
        if (a.index === atom.index) continue;
        let dx = a.x - atom.x;
        let dy = a.y - atom.y;
        let dz = a.z - atom.z;
        let d = Math.sqrt(dx*dx + dy*dy + dz*dz);
        neighbors.push({ index: a.index, elem: a.elem, dist: d });
      }
      neighbors.sort((x, y) => x.dist - y.dist);
      const topNeighbors = neighbors.slice(0, 5);

      let html = `
        <div style="font-weight: 600; font-size: 1rem; margin-bottom: 8px; display: flex; align-items: center; gap: 8px;">
          <span style="display:inline-block; width:12px; height:12px; border-radius:50%; background:#f59e0b;"></span>
          Atom #${atom.index} (${atom.elem})
        </div>
        <div style="color: #a5a5b0; font-size: 0.85rem; margin-bottom: 12px;">
          X: ${atom.x.toFixed(4)} Å<br>
          Y: ${atom.y.toFixed(4)} Å<br>
          Z: ${atom.z.toFixed(4)} Å
        </div>
        <div style="font-weight: 600; font-size: 0.8rem; text-transform: uppercase;
                    color: #8c8c99; border-top: 1px solid rgba(255,255,255,0.05); padding-top: 8px;">
          Nearest Neighbors
        </div>
        <div class="neighbor-list">
      `;

      if (topNeighbors.length === 0) {
        html += '<div style="font-size:0.8rem; color:#8c8c99;">No neighbors found</div>';
      } else {
        topNeighbors.forEach(n => {
          html += `
            <div class="neighbor-item">
              <span>Atom #${n.index} (${n.elem})</span>
              <span style="color:#ffffff; font-weight:600;">${n.dist.toFixed(3)} Å</span>
            </div>
          `;
        });
      }
      html += '</div>';

      $("#inspector-content").html(html);
    }
  </script>
</body>
</html>
"""


def open_browser_viewer(structure: Structure) -> None:
    """
    Saves the structure to a temporary HTML file utilizing 3Dmol.js WebGL rendering
    and opens it in the default browser.
    """
    # 1. Convert pymatgen Structure to CIF string
    cif_data = structure.to(fmt="cif")

    # 2. Extract structure metadata
    formula = structure.composition.reduced_formula
    num_atoms = len(structure)
    volume = f"{structure.volume:.2f}"
    density = f"{structure.density:.3f}"
    lattice_abc = ", ".join(f"{x:.3f}" for x in structure.lattice.abc)
    lattice_angles = ", ".join(f"{x:.2f}" for x in structure.lattice.angles)

    # 3. Fill the template
    html_content = VIEWER_TEMPLATE
    html_content = html_content.replace("{{CIF_DATA}}", cif_data.replace("`", "\\`"))
    html_content = html_content.replace("{{FORMULA}}", formula)
    html_content = html_content.replace("{{NUM_ATOMS}}", str(num_atoms))
    html_content = html_content.replace("{{VOLUME}}", volume)
    html_content = html_content.replace("{{DENSITY}}", density)
    html_content = html_content.replace("{{LATTICE_ABC}}", lattice_abc)
    html_content = html_content.replace("{{LATTICE_ANGLES}}", lattice_angles)

    # 4. Write to a temporary file
    # We do not delete the file automatically upon close, so the browser can read it asynchronously.
    # The temporary directory will naturally be cleaned up by the operating system.
    with tempfile.NamedTemporaryFile(
        suffix="_cellify_viewer.html", delete=False, mode="w", encoding="utf-8"
    ) as temp_file:
        temp_file.write(html_content)
        temp_file_path = temp_file.name

    # 5. Open in browser
    print(f"Structure WebGL viewer generated at: file://{temp_file_path}")
    webbrowser.open("file://" + os.path.abspath(temp_file_path))
