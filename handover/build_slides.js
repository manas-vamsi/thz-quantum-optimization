// Progress report deck: THz QUBO array design.
const pptxgen = require("pptxgenjs");
const path = require("path");

const REPO = "D:/thz-quantum-optimization";
const FIG = (p) => path.join(REPO, p);
const OUT = path.join(REPO, "handover", "Progress_Report_Vishwas.pptx");
const URL = "github.com/manas-vamsi/thz-quantum-optimization";

// palette: night sky + warm signal
const NAVY = "0B1D3A", INK = "1B2A41", AMBER = "F2A541", TEAL = "2A9D8F",
      RED = "C8553D", LIGHT = "F4F6FA", MUTED = "5B6B82", WHITE = "FFFFFF",
      CARD = "E8EDF5";

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9";                 // 10 x 5.625 in
pres.author = "Manasa Vamsi";
pres.title = "THz QUBO Array Design - Progress Report";
pres.theme = { headFontFace: "Cambria", bodyFontFace: "Calibri" };

// ---------- layouts ----------
pres.defineSlideMaster({
  title: "DARK",
  background: { color: NAVY },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 1.6, w: 8.8, h: 1.4,
        fontFace: "Cambria", fontSize: 36, bold: true, color: WHITE, valign: "top" },
        text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: WHITE },
  margin: [0.5, 0.5, 0.5, 0.5],
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 9.0, h: 0.75,
        fontFace: "Cambria", fontSize: 28, bold: true, color: NAVY, valign: "middle", margin: 0 },
        text: "" } },
    { text: { text: URL, options: { x: 0.5, y: 5.25, w: 6, h: 0.25, fontSize: 9, color: MUTED,
        fontFace: "Calibri", margin: 0 } } },
  ],
  slideNumber: { x: 9.1, y: 5.22, w: 0.4, h: 0.3, fontSize: 9, color: MUTED, fontFace: "Calibri" },
});

const T = (s, text, o) => s.addText(text, Object.assign({ isTextBox: true, fontFace: "Calibri",
  color: INK, fontSize: 14, margin: 0 }, o));

function card(s, x, y, w, h, fill, name) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, fill: { color: fill || CARD },
    line: { color: fill || CARD }, rectRadius: 0.08, objectName: name });
}
function stat(s, x, y, w, big, label, col) {
  T(s, big, { x, y, w, h: 0.8, fontFace: "Cambria", fontSize: 40, bold: true, color: col || NAVY });
  T(s, label, { x, y: y + 0.8, w, h: 0.55, fontSize: 12, color: MUTED });
}
function badge(s, x, y, n, col) {
  s.addShape(pres.shapes.OVAL, { x, y, w: 0.42, h: 0.42, fill: { color: col }, line: { color: col } });
  T(s, String(n), { x, y, w: 0.42, h: 0.42, align: "center", valign: "middle",
    fontSize: 14, bold: true, color: WHITE });
}
function table(s, rows, o) {
  const head = rows[0].map((t) => ({ text: t, options: { bold: true, color: WHITE, fill: { color: NAVY } } }));
  const body = rows.slice(1).map((r, i) => r.map((t) => {
    const bold = typeof t === "string" && t.startsWith("**");
    return { text: bold ? t.replace(/\*\*/g, "") : String(t),
      options: { bold, fill: { color: i % 2 ? WHITE : LIGHT }, color: INK } };
  }));
  s.addTable([head, ...body], Object.assign({ fontFace: "Calibri", fontSize: 11,
    border: { type: "solid", pt: 0.5, color: "D5DCE6" }, valign: "middle", margin: 0.06 }, o));
}
function content(title, section) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  return s;
}

// ================= 1. Title =================
pres.addSection({ title: "Introduction" });
let s = pres.addSlide({ masterName: "DARK", sectionTitle: "Introduction" });
s.addText("Quantum Optimization of Terahertz Interferometric Architectures", { placeholder: "title" });
T(s, "Progress report: formulation, solver benchmark, and three findings", {
  x: 0.6, y: 3.05, w: 8.8, h: 0.5, fontSize: 18, color: AMBER, italic: true });
T(s, "Manasa Vamsi", { x: 0.6, y: 3.9, w: 6, h: 0.35, fontSize: 16, color: WHITE, bold: true });
T(s, "For: Vishwas, Space Applications Centre, ISRO", { x: 0.6, y: 4.25, w: 8, h: 0.3, fontSize: 13, color: "AFC0D8" });
T(s, URL, { x: 0.6, y: 4.6, w: 8, h: 0.3, fontSize: 12, color: "AFC0D8" });
s.addNotes("Thank you for the time. I will cover three things: what you asked for and where it is, what I built around it so it could be tested, and three results that change the plan. Two of them are not what we expected. Everything is in the repository and every number comes from a script you can re-run.");

// ================= 2. What was asked =================
s = content("What was asked", "Introduction");
card(s, 0.5, 1.3, 4.2, 3.6, NAVY, "request");
T(s, "THE REQUEST", { x: 0.8, y: 1.55, w: 3.6, h: 0.3, fontSize: 12, bold: true, color: AMBER, charSpacing: 2 });
T(s, "The mathematical formulation of the objective function for QUBO application",
  { x: 0.8, y: 1.95, w: 3.6, h: 1.6, fontSize: 20, fontFace: "Cambria", color: WHITE });
T(s, "QUBO = Quadratic Unconstrained Binary Optimization - the form a quantum annealer accepts",
  { x: 0.8, y: 3.75, w: 3.6, h: 0.9, fontSize: 11, color: "AFC0D8", italic: true });
T(s, "DELIVERED", { x: 5.1, y: 1.35, w: 4.4, h: 0.3, fontSize: 12, bold: true, color: TEAL, charSpacing: 2 });
const deliv = [["Full derivation", "docs/objective_function.md"],
  ["Implementation", "src/thz_opt/qubo/"], ["Write-up", "paper/manuscript.pdf, Section 2"]];
deliv.forEach(([a, b], i) => {
  const y = 1.8 + i * 0.85;
  badge(s, 5.1, y, i + 1, TEAL);
  T(s, a, { x: 5.7, y: y - 0.02, w: 3.8, h: 0.3, fontSize: 15, bold: true });
  T(s, b, { x: 5.7, y: y + 0.28, w: 3.8, h: 0.3, fontSize: 12, color: MUTED, fontFace: "Courier New" });
});
T(s, "Doing it properly turned up two things about the onboarding note.",
  { x: 5.1, y: 4.45, w: 4.4, h: 0.45, fontSize: 13, italic: true, color: RED });
s.addNotes("The brief was to formulate the objective so it can be written as a QUBO and handed to a quantum annealer. That is done: derivation, implementation and write-up. Doing it carefully turned up two things about the onboarding note, on the next two slides.");

// ================= 3. Section 4 claim proven =================
pres.addSection({ title: "The formulation" });
s = content("Your Section 4 claim is correct - and now proven", "The formulation");
T(s, "Claim in the note", { x: 0.5, y: 1.3, w: 4.3, h: 0.3, fontSize: 12, bold: true, color: MUTED });
T(s, "Minimising UV auto-correlation minimises PSF sidelobes",
  { x: 0.5, y: 1.6, w: 4.3, h: 0.9, fontSize: 20, fontFace: "Cambria", color: NAVY });
T(s, [
  { text: "Follows from Parseval's theorem", options: { bullet: true, breakLine: true } },
  { text: "Dirty-beam energy = sum of squared gridded sampling counts", options: { bullet: true, breakLine: true } },
  { text: "So 'fewer sidelobes' and 'less UV overlap' are the same number", options: { bullet: true } },
], { x: 0.5, y: 2.65, w: 4.3, h: 1.6, fontSize: 14, paraSpaceAfter: 6 });
card(s, 5.3, 1.3, 4.2, 3.5, LIGHT, "proof");
T(s, "Numerical check", { x: 5.6, y: 1.55, w: 3.6, h: 0.3, fontSize: 12, bold: true, color: MUTED });
T(s, "1.000000000000", { x: 5.4, y: 2.15, w: 4.0, h: 0.9, fontSize: 34, bold: true, fontFace: "Cambria", color: TEAL, align: "center" });
T(s, "ratio of the two quantities on a real layout", { x: 5.6, y: 3.05, w: 3.6, h: 0.35, fontSize: 12, color: MUTED, align: "center" });
T(s, "Why it matters: this is what licenses the whole objective - it is the thing itself, not a proxy.",
  { x: 5.6, y: 3.7, w: 3.6, h: 0.8, fontSize: 12, italic: true, color: INK });
s.addNotes("The note says minimising UV overlap minimises sidelobes. That is true but not obvious. Through Parseval's theorem the dirty-beam energy equals the sum of squared sample counts per UV cell, so they are literally the same quantity. Checked numerically: ratio one to twelve decimal places. Source: src/thz_opt/qubo/baseline_qubo.py, sidelobe_energy_terms.");

// ================= 4. Q_ij cannot be written =================
s = content("Your Q_ij recipe cannot be written as stated", "The formulation");
card(s, 0.5, 1.3, 4.3, 1.75, "FBEAE6", "problem");
T(s, "THE PROBLEM", { x: 0.75, y: 1.45, w: 3.8, h: 0.3, fontSize: 12, bold: true, color: RED, charSpacing: 2 });
T(s, "Baseline overlap needs 4 antenna indices (i, j, k, l). Q_ij has only 2 slots. Written directly the objective is quartic - not a QUBO.",
  { x: 0.75, y: 1.8, w: 3.85, h: 1.15, fontSize: 13 });
card(s, 0.5, 3.25, 4.3, 1.75, "E3F2EF", "fix");
T(s, "THE FIX (core contribution)", { x: 0.75, y: 3.4, w: 3.8, h: 0.3, fontSize: 12, bold: true, color: TEAL, charSpacing: 2 });
T(s, [
  { text: "One activation variable y_k per antenna pair", options: { bullet: true, breakLine: true } },
  { text: "Enforce y_k = x_i * x_j by Rosenberg quadratization", options: { bullet: true, breakLine: true } },
  { text: "Overlap becomes y_k * y_l - quadratic again", options: { bullet: true } },
], { x: 0.75, y: 3.75, w: 3.85, h: 1.2, fontSize: 13, paraSpaceAfter: 3 });
// index diagram
const lbl = ["i", "j", "k", "l"];
lbl.forEach((c, k) => {
  const x = 5.4 + k * 1.0;
  s.addShape(pres.shapes.OVAL, { x, y: 1.5, w: 0.7, h: 0.7, fill: { color: k < 2 ? AMBER : TEAL }, line: { color: WHITE } });
  T(s, c, { x, y: 1.5, w: 0.7, h: 0.7, align: "center", valign: "middle", fontSize: 20, bold: true, color: WHITE, fontFace: "Cambria" });
});
T(s, "baseline 1", { x: 5.4, y: 2.3, w: 1.7, h: 0.3, fontSize: 11, color: MUTED, align: "center" });
T(s, "baseline 2", { x: 7.4, y: 2.3, w: 1.7, h: 0.3, fontSize: 11, color: MUTED, align: "center" });
T(s, "4 indices  >  Q_ij's 2 slots", { x: 5.3, y: 2.75, w: 4.2, h: 0.4, fontSize: 15, bold: true, color: RED, align: "center" });
T(s, "M + M(M-1)/2", { x: 5.4, y: 3.45, w: 2.0, h: 0.6, fontSize: 20, bold: true, color: NAVY, fontFace: "Cambria" });
T(s, "variables after the fix", { x: 5.4, y: 4.1, w: 2.0, h: 0.5, fontSize: 12, color: MUTED });
T(s, "~15,000", { x: 7.5, y: 3.45, w: 2.0, h: 0.6, fontSize: 28, bold: true, color: RED, fontFace: "Cambria" });
T(s, "variables for all 174 ALMA pads", { x: 7.5, y: 4.1, w: 2.0, h: 0.5, fontSize: 12, color: MUTED });
s.addNotes("This is the core contribution. Whether two baselines overlap depends on four antennas; Q_ij can reference only two, so the objective is fourth order, not a QUBO. The fix: a variable per antenna pair saying whether that baseline is active, forced to equal the product of the two pad variables by a Rosenberg penalty. The penalty weights are derived, not guessed, and a test shows the constraint breaks below them. The cost is more variables: about fifteen thousand for the full ALMA field.");

// ================= 5. Real data =================
pres.addSection({ title: "What was built" });
s = content("Real data, checked against published values", "What was built");
table(s, [
  ["Input", "Source", "Check"],
  ["174 ALMA antenna pads", "CASA configuration files", "16,195 m longest baseline (published ~16 km)"],
  ["27 VLA pads (control)", "CASA configuration files", "36,623 m (published ~36 km)"],
  ["Atmospheric phase", "ALMA Memo 624 - 17,000 observations", "Memo's factors 1.51/1.59/1.22 reproduced as 1.516/1.595/1.223"],
  ["Water vapour", "Cortes 2020 - 20 years at Chajnantor", "Monthly percentiles"],
  ["Terrain (Ladakh)", "SRTM / ASTER elevation", "Hanle 4,491 m (pub. 4,500); ALMA 5,036 m (pub. 5,059)"],
  ["Science source", "HL Tau, ALMA 2015", "Gap radii 13.2 / 32.3 / 64.2 au"],
], { x: 0.5, y: 1.3, w: 9.0, colW: [2.1, 2.9, 4.0], rowH: 0.48 });
T(s, "The VLA row is the load-bearing check: geocentric coordinates had to be rotated into a local frame, and recovering 36 km proves the transform.",
  { x: 0.5, y: 4.75, w: 9.0, h: 0.45, fontSize: 12, italic: true, color: MUTED });
s.addNotes("Nothing in the results is invented. Each input is checked against something published. The VLA check proves the coordinate transformation, not just the file. The atmosphere memo publishes its own scaling factors and our implementation reproduces them to three decimals. Source: data/external/SOURCES.md.");

// ================= 6. What I built =================
s = content("What I built around the formulation", "What was built");
const built = [
  ["Optimisers", "greedy, local search, annealing"], ["Exact solver", "MILP with proven dual bound"],
  ["QUBO solver", "dwave-samplers (D-Wave interface)"], ["Science objective", "Fisher information on a real source"],
  ["Imaging", "sky > visibilities > CLEAN > compare"], ["Instrument model", "primary beam, smearing, mJy"],
  ["Site design", "terrain-constrained Ladakh pads"], ["CLI tool", "pad coordinates in seconds"],
];
built.forEach(([a, b], i) => {
  const col = i % 4, row = Math.floor(i / 4);
  const x = 0.5 + col * 2.28, y = 1.3 + row * 1.2;
  card(s, x, y, 2.1, 1.05, LIGHT, "built" + i);
  T(s, a, { x: x + 0.15, y: y + 0.12, w: 1.8, h: 0.35, fontSize: 14, bold: true, color: NAVY });
  T(s, b, { x: x + 0.15, y: y + 0.5, w: 1.8, h: 0.5, fontSize: 11, color: MUTED });
});
[["272", "automated tests"], ["17", "experiments"], ["56", "modules"], ["47", "commits"]].forEach(([n, l], i) => {
  stat(s, 0.5 + i * 2.28, 3.85, 2.1, n, l, i === 0 ? TEAL : NAVY);
});
s.addNotes("The code could score a layout but not search for one, so the first job was optimisers, and an exact solver that returns a proven optimum with a certificate. That gives every other method something to be measured against. The rest followed from asking how we would know if each step is right. 272 tests pass.");

// ================= 7. Finding 1a =================
pres.addSection({ title: "Findings" });
s = content("Finding 1: the QUBO formulation is correct", "Findings");
T(s, "Is the mathematics right? Checked exhaustively on a 16-pad instance of real ALMA geometry.",
  { x: 0.5, y: 1.25, w: 9, h: 0.4, fontSize: 14, color: MUTED });
card(s, 0.5, 1.85, 2.9, 2.6, LIGHT, "s1");
stat(s, 0.75, 2.2, 2.5, "2,324", "feasible selections enumerated - every one", NAVY);
card(s, 3.55, 1.85, 2.9, 2.6, "E3F2EF", "s2");
stat(s, 3.8, 2.2, 2.5, "0 cells", "formulation gap: QUBO optimum = true optimum", TEAL);
card(s, 6.6, 1.85, 2.9, 2.6, LIGHT, "s3");
stat(s, 6.85, 2.2, 2.5, "0.929", "correlation, QUBO score vs true coverage", NAVY);
T(s, "Whatever happens next is about the solver, not the mathematics.",
  { x: 0.5, y: 4.65, w: 9, h: 0.4, fontSize: 15, bold: true, color: NAVY });
s.addNotes("Two separate questions: is the formulation right, and does the solver work? For a small instance you can check exhaustively: all 2,324 legal ways to choose 6 of 16 real pads. What the QUBO thinks is best is what is actually best. Formulation gap zero. Source: experiments/12_qubo_solver.py.");

// ================= 8. Finding 1b chart =================
s = content("Finding 1: but the QUBO route loses", "Findings");
s.addChart(pres.charts.BAR, [
  { name: "QUBO annealing", labels: ["16 pads, pick 6", "24 pads, pick 8", "40 pads, pick 10"], values: [90.8, 92.0, 85.9] },
  { name: "Direct classical search", labels: ["16 pads, pick 6", "24 pads, pick 8", "40 pads, pick 10"], values: [97.7, 98.8, 98.0] },
], { x: 0.5, y: 1.2, w: 5.6, h: 3.9, barDir: "col", barGrouping: "clustered",
  chartColors: [AMBER, TEAL], showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 10,
  dataLabelFormatCode: "0.0", valAxisMinVal: 80, valAxisMaxVal: 100, valAxisTitle: "% of proven optimum",
  showValAxisTitle: true, valAxisTitleFontSize: 10, showLegend: true, legendPos: "b", legendFontSize: 10,
  catAxisLabelFontSize: 10, valAxisLabelFontSize: 9, catAxisLabelColor: MUTED, valAxisLabelColor: MUTED,
  valGridLine: { color: "E1E6EE", size: 0.5 }, catGridLine: { style: "none" },
  showTitle: true, title: "Accuracy against the certified optimum", titleFontSize: 12, titleColor: NAVY,
  catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt", legendFontFace: "+mn-lt", titleFontFace: "+mn-lt" });
card(s, 6.4, 1.2, 3.1, 1.6, NAVY, "speed");
T(s, "35-72x", { x: 6.6, y: 1.3, w: 2.8, h: 0.8, fontSize: 38, bold: true, color: AMBER, fontFace: "Cambria" });
T(s, "faster: classical takes <0.4 s, annealing 3-35 s", { x: 6.6, y: 2.1, w: 2.8, h: 0.6, fontSize: 12, color: WHITE });
T(s, [
  { text: "Every annealing sample was feasible - not a tuning failure", options: { bullet: true, breakLine: true } },
  { text: "Gap widens with problem size", options: { bullet: true, breakLine: true } },
  { text: "40 pads already = 820 variables, 33,689 couplings", options: { bullet: true } },
], { x: 6.4, y: 3.0, w: 3.1, h: 2.0, fontSize: 12, paraSpaceAfter: 5 });
s.addNotes("The uncomfortable result. The QUBO ran through dwave-samplers, the same interface a real D-Wave machine uses. Annealing reaches 86 to 92 percent of the proven optimum; plain greedy plus swaps reaches 98 percent in under half a second. Every annealing sample was feasible, so the formulation works; converting and annealing is simply worse here. A quantum result must beat a 0.35 second classical method already within 2 percent of optimal.");

// ================= 9. Finding 2 =================
s = content("Finding 2: the golden ratio is not the point", "Findings");
T(s, "Named layouts change two things at once - how radius grows and how angle advances. So I crossed 5 radial laws x 5 angular laws, everything else fixed.",
  { x: 0.5, y: 1.2, w: 4.6, h: 0.95, fontSize: 13, color: MUTED });
s.addChart(pres.charts.DOUGHNUT, [{ name: "Share of variance", labels: ["Radial law  84.6%", "Angular law  11.2%", "Interaction  3.5%", "Noise  0.7%"], values: [84.6, 11.2, 3.5, 0.7] }],
  { x: 0.4, y: 2.15, w: 4.7, h: 3.0, holeSize: 55, chartColors: [NAVY, AMBER, TEAL, "C9D2DE"],
    showPercent: false, showValue: false, showLabel: false, showLegend: true,
    legendPos: "r", legendFontSize: 11, showTitle: true, title: "Variance in UV coverage", titleFontSize: 12,
    titleColor: NAVY, legendFontFace: "+mn-lt", titleFontFace: "+mn-lt", dataLabelFontFace: "+mn-lt" });
card(s, 5.4, 1.2, 4.1, 1.75, LIGHT, "r1");
stat(s, 5.65, 1.35, 3.6, "2.7x", "uniform-in-area vs Gaussian radial profile", NAVY);
card(s, 5.4, 3.1, 4.1, 1.75, "FDF1DF", "r2");
stat(s, 5.65, 3.25, 3.6, "2.5%", "golden angle vs a degenerate 5-bearing angle", AMBER);
s.addNotes("A golden spiral changes radius and angle together, so a ranking cannot say which matters. Separated with a factorial design and variance decomposition: radial distribution explains 85 percent, angle 11. The golden angle beats a deliberately bad 5-bearing angle by 2.5 percent. Caveat: sidelobes behave differently, with a large interaction. Source: experiments/16_factorial_layout.py.");

// ================= 10. Finding 3 =================
s = content("Finding 3: Hanle versus the onboarding note", "Findings");
T(s, "The note: Hanle is 'a comparable high-altitude, low-PWV site' to ALMA. Measured against the note's own 1 mm threshold:",
  { x: 0.5, y: 1.2, w: 9, h: 0.5, fontSize: 13, color: MUTED });
s.addChart(pres.charts.BAR, [{ name: "% time PWV < 1 mm",
  labels: ["ALMA, Chajnantor", "Ladakh site A", "Ladakh site B", "Merak", "Hanle"], values: [50, 23, 19, 8, 5] }],
  { x: 0.5, y: 1.75, w: 5.6, h: 3.35, barDir: "bar", chartColors: [NAVY, TEAL, TEAL, AMBER, RED],
    invertedColors: [NAVY], showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0\"%\"",
    dataLabelFontSize: 10, showLegend: false, catAxisLabelFontSize: 11, valAxisHidden: true,
    valGridLine: { style: "none" }, catGridLine: { style: "none" }, catAxisLabelColor: INK,
    catAxisLabelFontFace: "+mn-lt", dataLabelFontFace: "+mn-lt", varyColors: true });
card(s, 6.4, 1.75, 3.1, 3.35, LIGHT, "frame");
T(s, "NOT a rejection of Hanle", { x: 6.6, y: 1.9, w: 2.8, h: 0.35, fontSize: 13, bold: true, color: TEAL });
T(s, [
  { text: "Only candidate with road, power, observatory", options: { bullet: true, breakLine: true } },
  { text: "Flattest ground: 15 m relief vs 57 m at site A", options: { bullet: true, breakLine: true } },
  { text: "Argue it as infrastructure vs dryness - not parity with ALMA", options: { bullet: true } },
], { x: 6.6, y: 2.3, w: 2.8, h: 2.0, fontSize: 12, paraSpaceAfter: 5 });
T(s, "Raghunath et al. 2026, arXiv:2604.13487 - 184 months of ERA5", { x: 6.6, y: 4.5, w: 2.8, h: 0.5, fontSize: 9, color: MUTED, italic: true });
s.addNotes("Careful framing: this is about framing, not the site. Against the note's own threshold Hanle qualifies about 5 percent of the time and Chajnantor about 50, an order of magnitude. Hanle still has road, power and the flattest ground. Argue it as infrastructure against dryness with numbers stated; if the paper claims parity a reviewer will check.");

// ================= 11. Ladakh design =================
pres.addSection({ title: "Design and validation" });
s = content("A real array design for Ladakh", "Design and validation");
s.addImage({ path: FIG("ladakh_array_design/outputs/fig_hanle_extended_terrain.png"), x: 0.5, y: 1.2, w: 5.2, h: 5.2 / 2.315, objectName: "terrain" });
T(s, "Hanle terrain (left) and buildable ground with chosen pads (right)", { x: 0.5, y: 3.5, w: 5.2, h: 0.3, fontSize: 10, color: MUTED, italic: true });
table(s, [
  ["", "Compact", "Intermediate", "Extended"],
  ["Baselines", "32-232 m", "64-878 m", "338-2,788 m"],
  ["Resolves to", "1.16\"", "0.31\"", "**0.096\""],
  ["Sees up to", "**5.01\"", "2.50\"", "0.48\""],
  ["Relief", "11 m", "9 m", "15 m"],
], { x: 5.9, y: 1.2, w: 3.6, colW: [0.9, 0.85, 0.95, 0.9], rowH: 0.36, fontSize: 10 });
T(s, [
  { text: "Continuous coverage 0.096\" to 5.01\" (52x)", options: { bullet: true, breakLine: true } },
  { text: "Pad lat / long / elevation, CASA format", options: { bullet: true, breakLine: true } },
  { text: "Centre moved 2,049 m off the summit onto the plain", options: { bullet: true } },
], { x: 5.9, y: 3.15, w: 3.6, h: 1.3, fontSize: 12, paraSpaceAfter: 4 });
T(s, "Preliminary geometry - not a construction plan (no land rights, roads, power).",
  { x: 0.5, y: 4.6, w: 9, h: 0.4, fontSize: 12, italic: true, color: RED });
s.addNotes("Real elevation data, filtered for buildable slope and confined to one landform. Outputs are latitude, longitude and elevation for every pad in ALMA's own file format. Three configurations because sixteen antennas over three kilometres cannot also give 32 metre baselines; I checked they overlap so nothing falls between them. Not a construction plan.");

// ================= 12. Fisher objective =================
s = content("Replacing an invented score with a derived one", "Design and validation");
card(s, 0.5, 1.25, 4.3, 1.4, "FBEAE6", "prob");
T(s, "Problem", { x: 0.7, y: 1.35, w: 3.9, h: 0.3, fontSize: 13, bold: true, color: RED });
T(s, "The score was 'occupied UV cells' - chosen for convenience. Optimising a score you invented is circular.", { x: 0.7, y: 1.7, w: 3.9, h: 0.9, fontSize: 12 });
const chain = ["HL Tau disc model", "Visibility", "Fisher information", "Error bar (Cramer-Rao)"];
chain.forEach((c, i) => {
  const y = 2.85 + i * 0.55;
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.5, y, w: 4.3, h: 0.42, fill: { color: i === 3 ? TEAL : LIGHT }, line: { color: i === 3 ? TEAL : LIGHT }, rectRadius: 0.06 });
  T(s, (i + 1) + ".  " + c, { x: 0.7, y, w: 3.9, h: 0.42, valign: "middle", fontSize: 13, bold: i === 3, color: i === 3 ? WHITE : INK });
});
table(s, [
  ["Parameter", "Error-bar gain", "UV cells lost"],
  ["Gap depth", "10.4%", "43%"],
  ["**Gap radius", "**16.0%", "9%"],
  ["Gap width", "10.9%", "45%"],
], { x: 5.2, y: 1.25, w: 4.3, colW: [1.6, 1.4, 1.3], rowH: 0.42, fontSize: 12 });
T(s, "Arrays now compare by the error bar on a physical quantity. Bonus: this objective needs no Rosenberg trick at all.",
  { x: 5.2, y: 3.15, w: 4.3, h: 1.0, fontSize: 13, italic: true, color: NAVY });
s.addNotes("Nobody builds a telescope to touch UV cells. The replacement comes from estimation theory: for a real source, HL Tau, Fisher information says exactly how much each UV sample is worth for measuring a ring's depth, radius or width. Optimising for the measurement beats cell count by 10 to 16 percent on the error bar, at a cost in coverage: the estimation versus imaging trade-off.");

// ================= 13. Imaging validation =================
s = content("Does UV coverage actually predict image quality?", "Design and validation");
T(s, "Built end-to-end imaging: sky > visibilities > real coverage > noise > CLEAN > compare. Transform verified to 2e-16.",
  { x: 0.5, y: 1.2, w: 9, h: 0.45, fontSize: 13, color: MUTED });
table(s, [
  ["Sky model", "Cells predict STRUCTURE", "Cells predict FLUX"],
  ["Point source", "**+1.000", "+1.000"],
  ["Close double", "**+1.000", "+0.600"],
  ["Smooth Gaussian", "**+1.000", "**0.000"],
  ["Disc with gap", "**+1.000", "+0.300"],
], { x: 0.5, y: 1.8, w: 5.2, colW: [1.8, 1.8, 1.6], rowH: 0.45, fontSize: 12 });
card(s, 6.0, 1.8, 3.5, 1.4, "E3F2EF", "ok");
T(s, "Structure: validated", { x: 6.2, y: 1.9, w: 3.1, h: 0.35, fontSize: 14, bold: true, color: TEAL });
T(s, "Perfect rank correlation against every sky model.", { x: 6.2, y: 2.3, w: 3.1, h: 0.8, fontSize: 12 });
card(s, 6.0, 3.35, 3.5, 1.4, "FBEAE6", "bad");
T(s, "Flux: not predicted", { x: 6.2, y: 3.45, w: 3.1, h: 0.35, fontSize: 14, bold: true, color: RED });
T(s, "Zero correlation on extended sources.", { x: 6.2, y: 3.85, w: 3.1, h: 0.8, fontSize: 12 });
s.addNotes("Everything so far optimises a stand-in for image quality. I built an end-to-end imaging pipeline to test it. The transform returns the input sky to one part in ten to the sixteen. For structure, occupied cells predict perfectly. For flux on extended sources they predict nothing. Source: experiments/14_imaging_validation.py.");

// ================= 14. Flux failure =================
s = content("The failure a cell count cannot see", "Design and validation");
card(s, 0.5, 1.3, 4.3, 2.4, LIGHT, "a");
T(s, "Searched, maximum cells", { x: 0.75, y: 1.45, w: 3.8, h: 0.35, fontSize: 14, bold: true, color: NAVY });
T(s, "5,204 UV cells", { x: 0.75, y: 1.8, w: 3.8, h: 0.3, fontSize: 12, color: MUTED });
T(s, "42%", { x: 0.75, y: 2.15, w: 3.8, h: 1.0, fontSize: 60, bold: true, color: RED, fontFace: "Cambria" });
T(s, "of the source's flux recovered", { x: 0.75, y: 3.2, w: 3.8, h: 0.35, fontSize: 12, color: MUTED });
card(s, 5.2, 1.3, 4.3, 2.4, LIGHT, "b");
T(s, "Golden spiral, on real pads", { x: 5.45, y: 1.45, w: 3.8, h: 0.35, fontSize: 14, bold: true, color: NAVY });
T(s, "3,350 UV cells", { x: 5.45, y: 1.8, w: 3.8, h: 0.3, fontSize: 12, color: MUTED });
T(s, "91%", { x: 5.45, y: 2.15, w: 3.8, h: 1.0, fontSize: 60, bold: true, color: TEAL, fontFace: "Cambria" });
T(s, "of the source's flux recovered", { x: 5.45, y: 3.2, w: 3.8, h: 0.35, fontSize: 12, color: MUTED });
T(s, "55% more cells, less than half the source. Maximising cells pushes antennas outward and resolves extended emission out.",
  { x: 0.5, y: 3.9, w: 9, h: 0.5, fontSize: 13 });
T(s, "This qualifies our own headline: the +10.1% from search is a gain in structure, not in flux.",
  { x: 0.5, y: 4.45, w: 9, h: 0.45, fontSize: 13, bold: true, color: NAVY });
s.addNotes("The cell-maximising array has 55 percent more cells and recovers 42 percent of a smooth source where the spiral recovers 91. More cells live at large UV radius, so maximising cells drives antennas outward, and long baselines resolve out extended emission. The paper now says the 10 percent search gain is structural, not flux, because a reviewer would find it otherwise.");

// ================= 15. Instrumental =================
s = content("Limits the coverage metrics cannot see", "Design and validation");
table(s, [
  ["Setup", "Bandwidth", "Time", "Primary beam", "Binds"],
  ["Spectral line", "108\"", "427\"", "9.7\"", "primary beam"],
  ["Narrow continuum", "5.4\"", "142\"", "9.7\"", "**bandwidth"],
  ["Wide continuum", "**0.67\"", "28\"", "9.7\"", "**bandwidth"],
], { x: 0.5, y: 1.25, w: 9, colW: [2.2, 1.6, 1.4, 1.8, 2.0], rowH: 0.42, fontSize: 12 });
s.addImage({ path: FIG("figures/fig31_field_of_view_limits.png"), x: 0.5, y: 3.05, w: 4.8, h: 4.8 / 2.401, objectName: "fov" });
card(s, 5.6, 3.05, 3.9, 2.0, NAVY, "tension");
T(s, "The tension", { x: 5.8, y: 3.15, w: 3.5, h: 0.35, fontSize: 14, bold: true, color: AMBER });
T(s, "Smearing scales with the synthesised beam. Extend the array for resolution and the usable field shrinks in proportion.",
  { x: 5.8, y: 3.5, w: 3.5, h: 1.45, fontSize: 12, color: WHITE });
s.addNotes("Three effects were missing: primary beam, bandwidth and time smearing. Which one limits the usable field depends on the observation: the dish for spectral line, bandwidth for continuum, collapsing to two-thirds of an arcsecond for wide continuum. Both smearing limits scale with the synthesised beam, so extending the array shrinks the usable field. Sensitivity is now in physical units: about 30 microjansky per hour at Hanle in good winter conditions.");

// ================= 16. CLI =================
s = content("A tool anyone can use", "Design and validation");
card(s, 0.5, 1.25, 9.0, 1.05, INK, "code");
T(s, "python design_array.py --site alma  --n 20 --dish 12\npython design_array.py --site hanle --n 16 --dish 8 --max-baseline 3000",
  { x: 0.7, y: 1.35, w: 8.6, h: 0.85, fontFace: "Courier New", fontSize: 13, color: AMBER });
const outs = [["pads.csv", "latitude, longitude, elevation"], ["baselines.csv", "every pairwise distance"],
  ["array.cfg", "CASA format, loads in simobserve"], ["run.yaml", "every parameter, exact rerun"]];
outs.forEach(([a, b], i) => {
  const x = 0.5 + i * 2.28;
  card(s, x, 2.55, 2.1, 1.15, LIGHT, "out" + i);
  T(s, a, { x: x + 0.15, y: 2.65, w: 1.8, h: 0.35, fontSize: 13, bold: true, color: NAVY, fontFace: "Courier New" });
  T(s, b, { x: x + 0.15, y: 3.05, w: 1.8, h: 0.55, fontSize: 11, color: MUTED });
});
T(s, "No penalty-weight input, deliberately.", { x: 0.5, y: 3.95, w: 9, h: 0.35, fontSize: 14, bold: true, color: RED });
T(s, "The classical search never breaks the antenna count, so nothing needs penalising. QUBO penalties are derived and shown with --show-penalties: a hand-typed value below them makes the formulation wrong silently.",
  { x: 0.5, y: 4.3, w: 9, h: 0.75, fontSize: 12 });
s.addNotes("Everything is wrapped in a command-line tool: give it a site, antenna count and dish size, get pad coordinates and every pairwise distance. There is no penalty input because the classical search never changes the count; on the QUBO route the code derives the smallest safe penalty and shows it, rather than inviting a wrong value.");

// ================= 17. Status =================
pres.addSection({ title: "Next steps" });
s = content("Where it stands", "Next steps");
[["13", "manuscript pages"], ["17", "figures"], ["25", "references, 0 broken"], ["272", "tests passing"]].forEach(([n, l], i) => {
  card(s, 0.5 + i * 2.28, 1.3, 2.1, 1.5, LIGHT, "st" + i);
  stat(s, 0.7 + i * 2.28, 1.45, 1.8, n, l, i === 3 ? TEAL : NAVY);
});
card(s, 0.5, 3.05, 4.3, 1.9, "E3F2EF", "ready");
T(s, "READY", { x: 0.75, y: 3.2, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: TEAL, charSpacing: 2 });
T(s, "The science. Every figure and table is generated by a script; paper/ is self-contained and checksummed.", { x: 0.75, y: 3.55, w: 3.8, h: 1.3, fontSize: 13 });
card(s, 5.2, 3.05, 4.3, 1.9, "FBEAE6", "notready");
T(s, "NOT READY - 7 fields", { x: 5.45, y: 3.2, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: RED, charSpacing: 2 });
T(s, "Affiliation, email, coauthors, funding, acknowledgements, venue, DOI.\npaper/build_pdf.py --check-metadata", { x: 5.45, y: 3.55, w: 3.8, h: 1.3, fontSize: 13 });
s.addNotes("The paper is thirteen pages, every number generated not typed, and the paper folder is self-contained with checksummed code. What blocks submission is not science but seven metadata fields I cannot fill in. The build refuses to call itself ready while any is blank.");

// ================= 18. Asks =================
s = pres.addSlide({ masterName: "DARK", sectionTitle: "Next steps" });
T(s, "What I need from you", { x: 0.6, y: 0.4, w: 8.8, h: 0.7, fontSize: 32, bold: true, color: WHITE, fontFace: "Cambria" });
const asks = [
  ["Coauthorship", "Do you want to be on the paper? Anyone else at SAC?"],
  ["Seven metadata fields", "Affiliation, email, funding, acknowledgements, interests, venue, DOI"],
  ["Which story do we tell?", "'Formulated, solved, quantum route does not currently pay' is defensible. 'Quantum optimization of THz arrays' is not."],
  ["What should this array observe?", "The one input I could not derive. Everything downstream follows from it."],
];
asks.forEach(([a, b], i) => {
  const y = 1.35 + i * 1.0;
  badge(s, 0.6, y + 0.05, i + 1, AMBER);
  T(s, a, { x: 1.25, y, w: 8.2, h: 0.35, fontSize: 17, bold: true, color: WHITE });
  T(s, b, { x: 1.25, y: y + 0.38, w: 8.2, h: 0.55, fontSize: 12, color: "AFC0D8" });
});
s.addNotes("Four things. Authorship: your entry is prepared but nobody goes on a paper without reading it. The metadata I cannot supply. The real decision: the evidence supports an honest paper that says the quantum route does not currently pay, not one titled around quantum succeeding. And the single most valuable input: what should this array observe? One sentence, measure X in sources of type Y, unblocks more than any further optimisation.");

// ================= Backup: glossary =================
pres.addSection({ title: "Backup" });
s = content("Backup: glossary", "Backup");
table(s, [
  ["Term", "Plain meaning", "Our value"],
  ["Lambda (wavelength)", "Size of the radio wave", "1.3 mm at 230 GHz; 1 mm at 300 GHz"],
  ["Baseline", "Distance between two antennas", "32-2,788 m (Ladakh); up to 16,195 m (ALMA)"],
  ["N / M", "Antennas placed / candidate spots", "e.g. 16 of 544"],
  ["Resolution", "Finest detail = lambda / longest baseline (smaller better)", "0.096\""],
  ["Largest scale", "Biggest structure = 0.6 lambda / shortest baseline", "5.01\""],
  ["Primary beam", "Field of view of one dish ~ 1.13 lambda / D", "38\" (8 m at 230 GHz)"],
  ["Sidelobe", "False ghost near a real source (lower better)", "0.08-0.10"],
  ["Shadowing", "Minimum separation = 1.5 x dish", "12 m for 8 m dishes"],
], { x: 0.5, y: 1.2, w: 9, colW: [1.9, 4.2, 2.9], rowH: 0.4, fontSize: 11 });
s.addNotes("Reference if asked about terms. One arcsecond is 1/3600 of a degree, roughly a one-rupee coin seen from 4 km.");

// ================= Backup: not claimed =================
s = content("Backup: what is deliberately not claimed", "Backup");
const nc = [
  ["No quantum hardware", "QUBO ran on a classical annealer behind the hardware interface"],
  ["No real visibility data", "Imaging is simulated - correctly, but simulated"],
  ["Synthetic pads in Sections 4.1-4.6", "By design, so N and extent stay fixed across comparisons"],
  ["No Ladakh phase measurement", "Needs an interferometer on site for years; ALMA values transferred and labelled a model"],
  ["Terrain only for Ladakh", "No land rights, roads, power or geotechnics"],
];
nc.forEach(([a, b], i) => {
  const y = 1.25 + i * 0.75;
  badge(s, 0.5, y + 0.05, i + 1, MUTED);
  T(s, a, { x: 1.1, y, w: 8.4, h: 0.32, fontSize: 14, bold: true, color: NAVY });
  T(s, b, { x: 1.1, y: y + 0.32, w: 8.4, h: 0.32, fontSize: 12, color: MUTED });
});
s.addNotes("Have this ready if the quantum angle is pushed. Each is stated in the manuscript.");

pres.writeFile({ fileName: OUT }).then((f) => console.log("wrote " + f));
