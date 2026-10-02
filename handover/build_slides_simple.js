// Short, professional, plain-English progress deck.
const pptxgen = require("pptxgenjs");
const path = require("path");

const REPO = "D:/thz-quantum-optimization";
const OUT = path.join(REPO, "handover", "Progress_Report_Final.pptx");
const URL = "github.com/manas-vamsi/thz-quantum-optimization";

const NAVY = "0B1D3A", INK = "1B2A41", AMBER = "F2A541", TEAL = "2A9D8F",
      RED = "C8553D", LIGHT = "F4F6FA", MUTED = "5B6B82", WHITE = "FFFFFF",
      GOOD_BG = "E3F2EF", BAD_BG = "FBEAE6";

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9";
pres.author = "Manasa Vamsi";
pres.title = "THz Array Design with QUBO - Progress Report";
pres.theme = { headFontFace: "Cambria", bodyFontFace: "Calibri" };

pres.defineSlideMaster({
  title: "DARK", background: { color: NAVY },
  objects: [{ placeholder: { options: { name: "title", type: "title", x: 0.6, y: 1.4, w: 8.8, h: 1.5,
    fontFace: "Cambria", fontSize: 36, bold: true, color: WHITE, valign: "top" }, text: "" } }],
});
pres.defineSlideMaster({
  title: "CONTENT", background: { color: WHITE },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 9.0, h: 0.75,
      fontFace: "Cambria", fontSize: 28, bold: true, color: NAVY, valign: "middle", margin: 0 }, text: "" } },
    { text: { text: URL, options: { x: 0.5, y: 5.25, w: 6, h: 0.25, fontSize: 9, color: MUTED, fontFace: "Calibri", margin: 0 } } },
  ],
  slideNumber: { x: 9.1, y: 5.22, w: 0.4, h: 0.3, fontSize: 9, color: MUTED, fontFace: "Calibri" },
});

const T = (s, text, o) => s.addText(text, Object.assign({ isTextBox: true, fontFace: "Calibri",
  color: INK, fontSize: 15, margin: 0 }, o));
const card = (s, x, y, w, h, fill) => s.addShape(pres.shapes.ROUNDED_RECTANGLE,
  { x, y, w, h, fill: { color: fill }, line: { color: fill }, rectRadius: 0.08 });
const bullets = (items) => items.map((t, i) => ({ text: t,
  options: { bullet: true, breakLine: i < items.length - 1 } }));
function slide(title, section) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  return s;
}
function circle(s, x, y, d, label, fill) {
  s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill }, line: { color: fill } });
  T(s, label, { x, y, w: d, h: d, align: "center", valign: "middle", fontSize: Math.round(d * 30), bold: true, color: WHITE });
}
// verdict panel: strength (green) on top, limitation (red) below
function verdict(s, x, y, w, good, bad) {
  card(s, x, y, w, 1.85, GOOD_BG);
  circle(s, x + 0.2, y + 0.2, 0.42, "+", TEAL);
  T(s, "STRENGTH", { x: x + 0.75, y: y + 0.25, w: w - 0.9, h: 0.35, fontSize: 14, bold: true, color: TEAL });
  T(s, bullets(good), { x: x + 0.25, y: y + 0.75, w: w - 0.45, h: 1.0, fontSize: 13, paraSpaceAfter: 4 });
  card(s, x, y + 2.0, w, 1.85, BAD_BG);
  circle(s, x + 0.2, y + 2.2, 0.42, "-", RED);
  T(s, "LIMITATION", { x: x + 0.75, y: y + 2.25, w: w - 0.9, h: 0.35, fontSize: 14, bold: true, color: RED });
  T(s, bullets(bad), { x: x + 0.25, y: y + 2.75, w: w - 0.45, h: 1.0, fontSize: 13, paraSpaceAfter: 4 });
}
const chartText = { catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
  dataLabelFontFace: "+mn-lt", legendFontFace: "+mn-lt", titleFontFace: "+mn-lt" };
function hbar(s, title, labels, values, colors, fmt) {
  s.addChart(pres.charts.BAR, [{ name: "v", labels, values }], Object.assign({
    x: 0.5, y: 1.2, w: 5.0, h: 3.9, barDir: "bar", chartColors: colors, varyColors: true,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: fmt, dataLabelFontSize: 13,
    showLegend: false, catAxisLabelFontSize: 12, valAxisHidden: true, valGridLine: { style: "none" },
    catGridLine: { style: "none" }, catAxisLabelColor: INK, showTitle: true, title,
    titleFontSize: 13, titleColor: NAVY }, chartText));
}

// 1 Title ---------------------------------------------------------
pres.addSection({ title: "Introduction" });
let s = pres.addSlide({ masterName: "DARK", sectionTitle: "Introduction" });
s.addText("QUBO-Based Design of a Terahertz Telescope Array", { placeholder: "title" });
T(s, "Progress Report", { x: 0.6, y: 3.0, w: 8.8, h: 0.5, fontSize: 22, color: AMBER, italic: true });
T(s, "Manasa Vamsi", { x: 0.6, y: 3.9, w: 6, h: 0.4, fontSize: 18, bold: true, color: WHITE });
T(s, URL, { x: 0.6, y: 4.35, w: 8, h: 0.35, fontSize: 13, color: "AFC0D8" });
s.addNotes("Agenda: objective, method, results, and the strength and limitation of each result.");

// 2 Objective -----------------------------------------------------
s = slide("Project Objective", "Introduction");
card(s, 0.5, 1.3, 4.3, 3.6, NAVY);
T(s, "OBJECTIVE", { x: 0.8, y: 1.55, w: 3.7, h: 0.35, fontSize: 14, bold: true, color: AMBER });
T(s, "Write the antenna-placement problem as a QUBO, so it can be solved on a quantum computer",
  { x: 0.8, y: 2.0, w: 3.7, h: 2.2, fontSize: 20, fontFace: "Cambria", color: WHITE });
T(s, "DELIVERABLES", { x: 5.2, y: 1.35, w: 4.3, h: 0.35, fontSize: 14, bold: true, color: TEAL });
T(s, bullets([
  "QUBO formula - completed",
  "One error in the formula - corrected",
  "Working code - 272 tests pass",
  "Tested on real telescope data",
  "First array design for Ladakh",
]), { x: 5.2, y: 1.8, w: 4.3, h: 3.0, fontSize: 17, paraSpaceAfter: 10 });
s.addNotes("QUBO means every choice is 0 or 1. Here: 1 = antenna on this pad, 0 = empty. Quantum annealers only accept this form.");

// 3 Terms ---------------------------------------------------------
pres.addSection({ title: "Background" });
s = slide("Key Terms and Formulas", "Background");
const rows = [
  ["Term", "Meaning", "Formula", "Our value"],
  ["Wavelength (lambda)", "Size of the radio wave", "lambda = c / f", "1.3 mm at 230 GHz"],
  ["Baseline", "Distance between 2 antennas", "N(N-1)/2 baselines", "16 antennas = 120"],
  ["Resolution", "Smallest detail seen", "lambda / longest baseline", "0.096 arcsec"],
  ["Largest scale", "Biggest object seen", "0.6 lambda / shortest baseline", "5 arcsec"],
  ["UV coverage", "Share of picture sampled", "(u,v) = baseline / lambda", "about 50%"],
  ["QUBO", "Form a quantum computer accepts", "E = sum Q_ij x_i x_j, x = 0 or 1", "-"],
];
const head = rows[0].map((t) => ({ text: t, options: { bold: true, color: WHITE, fill: { color: NAVY } } }));
const body = rows.slice(1).map((r, i) => r.map((t, j) => ({ text: t,
  options: { fill: { color: i % 2 ? WHITE : LIGHT }, color: INK, bold: j === 0,
    fontFace: j === 2 ? "Courier New" : "Calibri" } })));
s.addTable([head, ...body], { x: 0.5, y: 1.2, w: 9.0, colW: [1.8, 2.4, 2.9, 1.9], rowH: 0.42,
  fontSize: 12, fontFace: "Calibri", border: { type: "solid", pt: 0.5, color: "D5DCE6" }, valign: "middle", margin: 0.06 });
T(s, "Long baselines give sharp detail.  Short baselines show big objects.  A good array needs both.",
  { x: 0.5, y: 4.65, w: 9, h: 0.4, fontSize: 15, bold: true, color: TEAL });
s.addNotes("Worked example: lambda = 3e8 / 230e9 = 1.3 mm. Resolution = 1.3 mm / 2788 m = 0.096 arcsec. 1 arcsec = 1/3600 degree. Note: lambda is also the name of the penalty weights in the QUBO - ask which one is meant.");

// 4 Parseval ------------------------------------------------------
s = slide("Parseval's Theorem and How We Used It", "Background");
card(s, 0.5, 1.25, 4.3, 1.75, LIGHT);
T(s, "WHAT IT SAYS", { x: 0.75, y: 1.4, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: NAVY });
T(s, "Total energy is the same in space and in frequency.", { x: 0.75, y: 1.75, w: 3.8, h: 0.5, fontSize: 15 });
T(s, "sum|image|^2 = sum|UV|^2 / N", { x: 0.75, y: 2.35, w: 3.8, h: 0.4, fontSize: 13, fontFace: "Courier New", color: TEAL, bold: true });
card(s, 0.5, 3.15, 4.3, 1.9, LIGHT);
T(s, "HOW WE USED IT", { x: 0.75, y: 3.3, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: NAVY });
T(s, bullets([
  "Image and UV plane are a Fourier pair",
  "So ghost-image energy = UV overlap count",
  "Reduce UV overlap = reduce ghost images",
]), { x: 0.75, y: 3.65, w: 3.8, h: 1.3, fontSize: 13, paraSpaceAfter: 4 });
// flow diagram on the right
const flow = [["UV overlap", "sum of n^2"], ["Parseval", "equal"], ["Ghost images", "sidelobe energy"]];
flow.forEach(([a, b], i) => {
  const y = 1.25 + i * 1.0;
  card(s, 5.2, y, 4.3, 0.8, i === 1 ? AMBER : NAVY);
  T(s, a, { x: 5.4, y: y + 0.08, w: 3.9, h: 0.35, fontSize: 16, bold: true, color: WHITE });
  T(s, b, { x: 5.4, y: y + 0.43, w: 3.9, h: 0.3, fontSize: 12, color: i === 1 ? WHITE : "AFC0D8" });
});
card(s, 5.2, 4.25, 4.3, 0.8, GOOD_BG);
T(s, "Checked by computer: ratio = 1.000000000000", { x: 5.4, y: 4.25, w: 3.9, h: 0.8, valign: "middle", fontSize: 14, bold: true, color: TEAL });
s.addNotes("Why it matters: ghost images are expensive to compute and not quadratic. UV overlap is simple counting and quadratic, so it fits a QUBO. Parseval proves the two are the same, so we can optimise the easy one. This also proves the guide's Section 4 claim is correct.");

// 5 Formulation fix -----------------------------------------------
s = slide("Correction to the QUBO Formula", "Background");
["i", "j", "k", "l"].forEach((c, k) => circle(s, 0.7 + k * 1.0, 1.35, 0.7, c, k < 2 ? AMBER : TEAL));
T(s, "baseline 1", { x: 0.7, y: 2.15, w: 1.7, h: 0.3, fontSize: 11, color: MUTED, align: "center" });
T(s, "baseline 2", { x: 2.7, y: 2.15, w: 1.7, h: 0.3, fontSize: 11, color: MUTED, align: "center" });
card(s, 0.5, 2.6, 4.3, 1.1, BAD_BG);
T(s, "PROBLEM", { x: 0.75, y: 2.7, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: RED });
T(s, "Overlap needs 4 antennas. Q_ij holds only 2. It does not fit.", { x: 0.75, y: 3.05, w: 3.8, h: 0.6, fontSize: 14 });
card(s, 0.5, 3.85, 4.3, 1.2, GOOD_BG);
T(s, "FIX (Rosenberg method)", { x: 0.75, y: 3.95, w: 3.8, h: 0.3, fontSize: 13, bold: true, color: TEAL });
T(s, "Add one 0/1 switch per antenna pair:  y = x_i * x_j.  Now it fits.", { x: 0.75, y: 4.3, w: 3.8, h: 0.7, fontSize: 14 });
verdict(s, 5.2, 1.2, 4.3,
  ["Formula is now correct", "Works with quantum software"],
  ["Many more variables", "About 15,000 for 174 pads"]);
s.addNotes("Formula: E(x) = sum Q_ii x_i + sum Q_ij x_i x_j. Full objective H = H_select + H_shadow + H_uv + H_phase. H_select = lambda_s (sum x - N)^2 forces N antennas. Rosenberg penalty = lambda_R (x_i x_j - 2 x_i y - 2 x_j y + 3y). Penalty weights are calculated, not guessed.");

// 6 Data ----------------------------------------------------------
pres.addSection({ title: "Method" });
s = slide("Data Sources (All Real and Verified)", "Method");
const data = [
  ["Antenna pads", "174 real ALMA pads, Chile"],
  ["Atmosphere", "17,000 real ALMA observations"],
  ["Water vapour", "20 years of measurements"],
  ["Ladakh land", "Real height maps (SRTM)"],
  ["Test star", "HL Tau, a real young star"],
];
data.forEach(([a, b], i) => {
  const y = 1.3 + i * 0.66;
  circle(s, 0.5, y, 0.5, String(i + 1), TEAL);
  T(s, a, { x: 1.2, y: y + 0.07, w: 2.3, h: 0.4, fontSize: 16, bold: true, color: NAVY });
  T(s, b, { x: 3.5, y: y + 0.07, w: 3.3, h: 0.4, fontSize: 15 });
});
card(s, 7.0, 1.3, 2.5, 3.2, NAVY);
T(s, "VERIFIED", { x: 7.2, y: 1.45, w: 2.1, h: 0.35, fontSize: 14, bold: true, color: AMBER });
T(s, "ALMA size\nOurs: 16,195 m\nPublished: ~16 km\n\nHanle height\nOurs: 4,491 m\nPublished: 4,500 m",
  { x: 7.2, y: 1.85, w: 2.1, h: 2.5, fontSize: 13, color: WHITE });
s.addNotes("No data is invented. Each input is checked against a published value.");

// 7 Method --------------------------------------------------------
s = slide("Methodology", "Method");
const steps = ["Write QUBO formula", "Build search code", "Find proven best", "Run QUBO solver", "Make test images", "Design for Ladakh"];
steps.forEach((t, i) => {
  const x = 0.5 + i * 1.53;
  card(s, x, 1.35, 1.38, 1.45, LIGHT);
  T(s, String(i + 1), { x, y: 1.45, w: 1.38, h: 0.5, align: "center", fontSize: 24, bold: true, fontFace: "Cambria", color: TEAL });
  T(s, t, { x: x + 0.05, y: 2.0, w: 1.28, h: 0.7, align: "center", fontSize: 13, bold: true, color: NAVY });
});
T(s, bullets([
  "Search code: greedy, swap and annealing methods",
  "Proven best: exact solver (MILP) - the answer key to compare against",
  "QUBO solver: D-Wave software on a normal computer (not quantum hardware)",
  "Test images: check if our score gives good pictures",
]), { x: 0.5, y: 3.1, w: 9, h: 1.9, fontSize: 15, paraSpaceAfter: 8 });
s.addNotes("Each step has automated tests - 272 in total.");

// 8 Result 1 ------------------------------------------------------
pres.addSection({ title: "Results" });
s = slide("Result 1: QUBO Solver vs Classical Search", "Results");
s.addChart(pres.charts.BAR, [
  { name: "QUBO solver", labels: ["16 pads", "24 pads", "40 pads"], values: [90.8, 92.0, 85.9] },
  { name: "Classical search", labels: ["16 pads", "24 pads", "40 pads"], values: [97.7, 98.8, 98.0] },
], Object.assign({ x: 0.5, y: 1.2, w: 5.0, h: 3.9, barDir: "col", barGrouping: "clustered",
  chartColors: [AMBER, TEAL], showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 11,
  dataLabelFormatCode: "0.0", valAxisMinVal: 80, valAxisMaxVal: 100, showLegend: true, legendPos: "b",
  legendFontSize: 11, catAxisLabelFontSize: 12, valAxisLabelFontSize: 10, catAxisLabelColor: MUTED,
  valAxisLabelColor: MUTED, valGridLine: { color: "E1E6EE", size: 0.5 }, catGridLine: { style: "none" },
  showTitle: true, title: "% of the best possible answer", titleFontSize: 13, titleColor: NAVY }, chartText));
verdict(s, 5.8, 1.2, 3.7,
  ["QUBO formula is correct", "Checked all 2,324 options"],
  ["Classical is more accurate", "Classical is 35-72x faster"]);
s.addNotes("Classical: under 0.4 s at 98%. QUBO: 3-35 s at 86-92%. Quantum is not faster on this problem yet.");

// 9 Result 2 ------------------------------------------------------
s = slide("Result 2: What Makes a Good Antenna Layout", "Results");
hbar(s, "Share of effect on picture quality",
  ["Both together", "Angle (golden ratio)", "Distance from centre"], [3.5, 11.2, 84.6],
  [TEAL, AMBER, NAVY], "0.0\"%\"");
verdict(s, 5.8, 1.2, 3.7,
  ["Clear design rule found", "Spread antennas evenly by area"],
  ["Golden ratio helps only ~2.5%", "Rule differs for ghost images"]);
s.addNotes("Tested 5 distance rules x 5 angle rules = 25 layouts, all else fixed. Method: ANOVA (separates which factor causes the difference).");

// 10 Result 3 -----------------------------------------------------
s = slide("Result 3: Site Assessment - Hanle", "Results");
hbar(s, "% of time air is dry enough (PWV < 1 mm)",
  ["Hanle", "Merak", "Ladakh site B", "Ladakh site A", "ALMA (Chile)"], [5, 8, 19, 23, 50],
  [RED, AMBER, TEAL, TEAL, NAVY], "0\"%\"");
verdict(s, 5.8, 1.2, 3.7,
  ["Roads, power, observatory", "Very flat land (15 m)"],
  ["Dry only 5% vs ALMA 50%", "Not equal to ALMA"]);
s.addNotes("Source: Raghunath et al. 2026, 184 months of weather data. PWV = water in the air; water absorbs these radio waves. Suggested wording: good infrastructure, less dry air.");

// 11 Result 4 -----------------------------------------------------
s = slide("Result 4: First Array Design for Ladakh", "Results");
s.addImage({ path: path.join(REPO, "ladakh_array_design/outputs/fig_hanle_extended_terrain.png"),
  x: 0.5, y: 1.25, w: 5.0, h: 5.0 / 2.315 });
T(s, "Left: land height.  Right: buildable land and chosen antenna spots.",
  { x: 0.5, y: 3.5, w: 5.0, h: 0.4, fontSize: 11, color: MUTED, italic: true });
T(s, bullets(["16 antennas, 8 m dishes", "3 layouts: small, medium, large", "Detail from 0.096 to 5 arcsec"]),
  { x: 0.5, y: 3.95, w: 5.0, h: 1.1, fontSize: 14, paraSpaceAfter: 3 });
verdict(s, 5.8, 1.2, 3.7,
  ["Exact location of each antenna", "Built on real land data"],
  ["First design only", "No land, road or power survey"]);
s.addNotes("Three layouts because one layout cannot give both sharp detail and big-scale view. ALMA also moves antennas between layouts.");

// 12 Result 5 -----------------------------------------------------
s = slide("Result 5: Image Quality Check", "Results");
card(s, 0.5, 1.25, 2.4, 2.2, LIGHT);
T(s, "Best-score array", { x: 0.65, y: 1.35, w: 2.1, h: 0.35, fontSize: 13, bold: true, color: NAVY });
T(s, "42%", { x: 0.65, y: 1.8, w: 2.1, h: 0.9, fontSize: 48, bold: true, color: RED, fontFace: "Cambria" });
T(s, "brightness captured", { x: 0.65, y: 2.75, w: 2.1, h: 0.4, fontSize: 12, color: MUTED });
card(s, 3.1, 1.25, 2.4, 2.2, LIGHT);
T(s, "Golden spiral", { x: 3.25, y: 1.35, w: 2.1, h: 0.35, fontSize: 13, bold: true, color: NAVY });
T(s, "91%", { x: 3.25, y: 1.8, w: 2.1, h: 0.9, fontSize: 48, bold: true, color: TEAL, fontFace: "Cambria" });
T(s, "brightness captured", { x: 3.25, y: 2.75, w: 2.1, h: 0.4, fontSize: 12, color: MUTED });
T(s, bullets(["Made real test pictures from 4 sky models", "Compared each picture with the true sky"]),
  { x: 0.5, y: 3.7, w: 5.0, h: 1.3, fontSize: 14, paraSpaceAfter: 5 });
verdict(s, 5.8, 1.2, 3.7,
  ["Our score predicts shape perfectly", "Validates the whole method"],
  ["Misses brightness of big objects", "Needs short baselines too"]);
s.addNotes("Pipeline: sky model, visibilities, noise, CLEAN, compare. Shape correlation = 1.000. Brightness correlation for big blurry sources = 0.");

// 13 Summary ------------------------------------------------------
pres.addSection({ title: "Summary" });
s = slide("Summary", "Summary");
const sum = [
  ["Result", "Strength", "Limitation"],
  ["QUBO formula", "Correct, fixed one error", "Many variables"],
  ["QUBO solver", "Works, all answers valid", "Slower than classical"],
  ["Antenna layout", "Clear rule: spread by area", "Golden ratio adds little"],
  ["Hanle site", "Roads, power, flat land", "Dry only 5% of time"],
  ["Ladakh design", "Exact antenna positions", "No survey yet"],
  ["Image check", "Shape predicted perfectly", "Brightness not predicted"],
];
const sh = sum[0].map((t) => ({ text: t, options: { bold: true, color: WHITE, fill: { color: NAVY } } }));
const sb = sum.slice(1).map((r, i) => r.map((t, j) => ({ text: t, options: {
  bold: j === 0, color: j === 1 ? "1E6B61" : j === 2 ? "A2412E" : INK,
  fill: { color: j === 1 ? GOOD_BG : j === 2 ? BAD_BG : (i % 2 ? WHITE : LIGHT) } } })));
s.addTable([sh, ...sb], { x: 0.5, y: 1.2, w: 9.0, colW: [2.4, 3.3, 3.3], rowH: 0.48,
  fontSize: 14, fontFace: "Calibri", border: { type: "solid", pt: 0.5, color: WHITE }, valign: "middle", margin: 0.08 });
s.addNotes("Paper: 13 pages. 272 tests pass. All code and results are in the repository.");

pres.writeFile({ fileName: OUT }).then((f) => console.log("wrote " + f));
