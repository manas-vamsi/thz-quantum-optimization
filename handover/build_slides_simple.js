// Short, plain-English version of the progress deck.
const pptxgen = require("pptxgenjs");
const path = require("path");

const REPO = "D:/thz-quantum-optimization";
const OUT = path.join(REPO, "handover", "Progress_Report_Simple.pptx");
const URL = "github.com/manas-vamsi/thz-quantum-optimization";

const NAVY = "0B1D3A", INK = "1B2A41", AMBER = "F2A541", TEAL = "2A9D8F",
      RED = "C8553D", LIGHT = "F4F6FA", MUTED = "5B6B82", WHITE = "FFFFFF",
      GREEN_BG = "E3F2EF", RED_BG = "FBEAE6";

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9";
pres.author = "Manasa Vamsi";
pres.title = "THz Array Design - Progress (short)";
pres.theme = { headFontFace: "Cambria", bodyFontFace: "Calibri" };

pres.defineSlideMaster({
  title: "DARK", background: { color: NAVY },
  objects: [{ placeholder: { options: { name: "title", type: "title", x: 0.6, y: 1.5, w: 8.8, h: 1.4,
    fontFace: "Cambria", fontSize: 36, bold: true, color: WHITE, valign: "top" }, text: "" } }],
});
pres.defineSlideMaster({
  title: "CONTENT", background: { color: WHITE },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 9.0, h: 0.75,
      fontFace: "Cambria", fontSize: 30, bold: true, color: NAVY, valign: "middle", margin: 0 }, text: "" } },
    { text: { text: URL, options: { x: 0.5, y: 5.25, w: 6, h: 0.25, fontSize: 9, color: MUTED, fontFace: "Calibri", margin: 0 } } },
  ],
  slideNumber: { x: 9.1, y: 5.22, w: 0.4, h: 0.3, fontSize: 9, color: MUTED, fontFace: "Calibri" },
});

const T = (s, text, o) => s.addText(text, Object.assign({ isTextBox: true, fontFace: "Calibri",
  color: INK, fontSize: 16, margin: 0 }, o));
const card = (s, x, y, w, h, fill) => s.addShape(pres.shapes.ROUNDED_RECTANGLE,
  { x, y, w, h, fill: { color: fill }, line: { color: fill }, rectRadius: 0.08 });
const bullets = (items) => items.map((t, i) => ({ text: t,
  options: { bullet: true, breakLine: i < items.length - 1 } }));
function slide(title, section) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  return s;
}
function big(s, x, y, w, num, label, col) {
  T(s, num, { x, y, w, h: 0.85, fontFace: "Cambria", fontSize: 44, bold: true, color: col });
  T(s, label, { x, y: y + 0.85, w, h: 0.6, fontSize: 14, color: MUTED });
}
const chartText = { catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
  dataLabelFontFace: "+mn-lt", legendFontFace: "+mn-lt", titleFontFace: "+mn-lt" };

// 1 ---------------------------------------------------------------
pres.addSection({ title: "Start" });
let s = pres.addSlide({ masterName: "DARK", sectionTitle: "Start" });
s.addText("Designing a Terahertz Telescope Array with QUBO", { placeholder: "title" });
T(s, "What I did, what I found, what comes next", { x: 0.6, y: 3.0, w: 8.8, h: 0.5, fontSize: 20, color: AMBER, italic: true });
T(s, "Manasa Vamsi", { x: 0.6, y: 3.9, w: 6, h: 0.4, fontSize: 18, bold: true, color: WHITE });
T(s, URL, { x: 0.6, y: 4.35, w: 8, h: 0.35, fontSize: 13, color: "AFC0D8" });
s.addNotes("Today: what you asked, how I did it, three results, and what I need from you. Everything is in the repo.");

// 2 ---------------------------------------------------------------
s = slide("What you asked, and what I did", "Start");
card(s, 0.5, 1.3, 4.3, 3.6, NAVY);
T(s, "YOU ASKED", { x: 0.8, y: 1.55, w: 3.7, h: 0.35, fontSize: 14, bold: true, color: AMBER });
T(s, "Write the maths (objective function) so the problem can run on a quantum computer (QUBO)",
  { x: 0.8, y: 2.0, w: 3.7, h: 2.2, fontSize: 20, fontFace: "Cambria", color: WHITE });
T(s, "I DID", { x: 5.2, y: 1.35, w: 4.3, h: 0.35, fontSize: 14, bold: true, color: TEAL });
T(s, bullets([
  "Wrote the maths - done",
  "Fixed one problem in it",
  "Built the code and tested it",
  "Ran it on real telescope data",
  "Made a real design for Ladakh",
]), { x: 5.2, y: 1.8, w: 4.3, h: 3.0, fontSize: 18, paraSpaceAfter: 10 });
s.addNotes("Main point: the request is done. Plus I tested it and used real data.");

// 3 ---------------------------------------------------------------
pres.addSection({ title: "Basics" });
s = slide("Basic terms and formulas", "Basics");
const rows = [
  ["Term", "Meaning", "Formula", "Our value"],
  ["Lambda (wavelength)", "Size of the radio wave", "lambda = c / f", "1.3 mm at 230 GHz"],
  ["Baseline", "Distance between 2 antennas", "N antennas give N(N-1)/2", "16 antennas = 120"],
  ["Resolution", "Smallest detail we can see", "lambda / longest baseline", "0.096 arcsec"],
  ["Largest scale", "Biggest thing we can see", "0.6 lambda / shortest baseline", "5 arcsec"],
  ["UV coverage", "How much of the picture we sample", "(u, v) = baseline / lambda", "~50% of grid"],
  ["QUBO", "Form a quantum computer accepts", "E = sum Q_ij x_i x_j,  x = 0 or 1", "-"],
];
const head = rows[0].map((t) => ({ text: t, options: { bold: true, color: WHITE, fill: { color: NAVY } } }));
const body = rows.slice(1).map((r, i) => r.map((t, j) => ({ text: t,
  options: { fill: { color: i % 2 ? WHITE : LIGHT }, color: INK, bold: j === 0,
    fontFace: j === 2 ? "Courier New" : "Calibri" } })));
s.addTable([head, ...body], { x: 0.5, y: 1.25, w: 9.0, colW: [1.8, 2.5, 2.8, 1.9], rowH: 0.42,
  fontSize: 12, fontFace: "Calibri", border: { type: "solid", pt: 0.5, color: "D5DCE6" }, valign: "middle", margin: 0.06 });
T(s, "Rule: long baselines = sharp detail.  Short baselines = see big things.  We need both.",
  { x: 0.5, y: 4.65, w: 9, h: 0.4, fontSize: 15, bold: true, color: TEAL });
s.addNotes("Example: lambda = 3e8 / 230e9 = 1.3 mm. Resolution = 1.3 mm / 2788 m = 0.096 arcsec. One arcsecond is 1/3600 of a degree.");

// 4 ---------------------------------------------------------------
s = slide("Your guide: one part right, one part fixed", "Basics");
card(s, 0.5, 1.3, 4.3, 3.6, GREEN_BG);
T(s, "RIGHT", { x: 0.8, y: 1.5, w: 3.7, h: 0.4, fontSize: 18, bold: true, color: TEAL });
T(s, bullets([
  "Less UV overlap = fewer false ghost images (sidelobes)",
  "I proved it with Parseval's theorem",
  "Checked by computer: 1.000000000000",
]), { x: 0.8, y: 2.0, w: 3.7, h: 2.7, fontSize: 16, paraSpaceAfter: 10 });
card(s, 5.2, 1.3, 4.3, 3.6, RED_BG);
T(s, "FIXED", { x: 5.5, y: 1.5, w: 3.7, h: 0.4, fontSize: 18, bold: true, color: RED });
T(s, bullets([
  "Overlap needs 4 antennas (i, j, k, l)",
  "Q_ij has only 2 slots - does not fit",
  "Fix: add one switch per antenna pair",
  "Now it fits QUBO again",
]), { x: 5.5, y: 2.0, w: 3.7, h: 2.7, fontSize: 16, paraSpaceAfter: 10 });
s.addNotes("The fix is called Rosenberg quadratization. It adds y_k = x_i * x_j for each pair. Cost: more variables, about 15,000 for all ALMA pads.");

// 5 ---------------------------------------------------------------
pres.addSection({ title: "Method" });
s = slide("Data I used - all real", "Method");
const data = [
  ["ALMA antenna pads", "174 real pads, Chile"],
  ["Atmosphere", "17,000 real ALMA observations"],
  ["Water vapour", "20 years of measurements"],
  ["Ladakh land", "Real height maps (SRTM)"],
  ["Test star", "HL Tau, a real young star"],
];
data.forEach(([a, b], i) => {
  const y = 1.3 + i * 0.66;
  s.addShape(pres.shapes.OVAL, { x: 0.5, y, w: 0.5, h: 0.5, fill: { color: TEAL }, line: { color: TEAL } });
  T(s, String(i + 1), { x: 0.5, y, w: 0.5, h: 0.5, align: "center", valign: "middle", fontSize: 16, bold: true, color: WHITE });
  T(s, a, { x: 1.2, y: y + 0.05, w: 2.6, h: 0.4, fontSize: 17, bold: true, color: NAVY });
  T(s, b, { x: 3.8, y: y + 0.05, w: 3.0, h: 0.4, fontSize: 16 });
});
card(s, 7.0, 1.3, 2.5, 3.2, NAVY);
T(s, "Checked", { x: 7.2, y: 1.45, w: 2.1, h: 0.35, fontSize: 15, bold: true, color: AMBER });
T(s, "Our ALMA size: 16,195 m\nPublished: ~16 km\n\nOur Hanle height: 4,491 m\nPublished: 4,500 m",
  { x: 7.2, y: 1.9, w: 2.1, h: 2.4, fontSize: 13, color: WHITE });
s.addNotes("Nothing is made up. Every input is checked against a published number.");

// 6 ---------------------------------------------------------------
s = slide("The process, step by step", "Method");
const steps = ["Write the maths", "Build optimiser", "Find true best", "Run QUBO", "Make images", "Design Ladakh"];
steps.forEach((t, i) => {
  const x = 0.5 + i * 1.53;
  card(s, x, 1.4, 1.38, 1.4, i === 3 ? AMBER : LIGHT);
  T(s, String(i + 1), { x, y: 1.5, w: 1.38, h: 0.5, align: "center", fontSize: 24, bold: true, fontFace: "Cambria", color: i === 3 ? WHITE : TEAL });
  T(s, t, { x: x + 0.05, y: 2.05, w: 1.28, h: 0.65, align: "center", fontSize: 13, bold: true, color: i === 3 ? WHITE : NAVY });
});
T(s, bullets([
  "Step 2: search methods - greedy, swap, annealing",
  "Step 3: exact solver (MILP) - gives the proven best answer to compare against",
  "Step 4: QUBO run with D-Wave software (simulated, not real quantum hardware)",
  "Step 5: make real pictures to check if our score means anything",
]), { x: 0.5, y: 3.1, w: 9, h: 1.9, fontSize: 15, paraSpaceAfter: 8 });
s.addNotes("272 automated tests check every step.");

// 7 ---------------------------------------------------------------
pres.addSection({ title: "Results" });
s = slide("Result 1: QUBO maths is right, but slower", "Results");
s.addChart(pres.charts.BAR, [
  { name: "QUBO", labels: ["16 pads", "24 pads", "40 pads"], values: [90.8, 92.0, 85.9] },
  { name: "Simple classical search", labels: ["16 pads", "24 pads", "40 pads"], values: [97.7, 98.8, 98.0] },
], Object.assign({ x: 0.5, y: 1.2, w: 5.4, h: 3.9, barDir: "col", barGrouping: "clustered",
  chartColors: [AMBER, TEAL], showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 11,
  dataLabelFormatCode: "0.0", valAxisMinVal: 80, valAxisMaxVal: 100, showLegend: true, legendPos: "b",
  legendFontSize: 11, catAxisLabelFontSize: 12, valAxisLabelFontSize: 10, catAxisLabelColor: MUTED,
  valAxisLabelColor: MUTED, valGridLine: { color: "E1E6EE", size: 0.5 }, catGridLine: { style: "none" },
  showTitle: true, title: "% of the best possible answer", titleFontSize: 13, titleColor: NAVY }, chartText));
T(s, bullets([
  "QUBO maths is correct (checked all 2,324 options)",
  "But simple search is better",
  "And 35-72x faster",
  "Quantum is not winning here yet",
]), { x: 6.2, y: 1.5, w: 3.3, h: 3.3, fontSize: 16, paraSpaceAfter: 12 });
s.addNotes("Classical takes under 0.4 seconds, QUBO 3-35 seconds.");

// 8 ---------------------------------------------------------------
s = slide("Result 2: golden ratio is not the key", "Results");
s.addChart(pres.charts.BAR, [{ name: "Share", labels: ["Both together", "Angle (golden ratio)", "Distance from centre"], values: [3.5, 11.2, 84.6] }],
  Object.assign({ x: 0.5, y: 1.2, w: 5.4, h: 3.9, barDir: "bar", chartColors: [TEAL, AMBER, NAVY], varyColors: true,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.0\"%\"", dataLabelFontSize: 13,
    showLegend: false, catAxisLabelFontSize: 13, valAxisHidden: true, valGridLine: { style: "none" },
    catGridLine: { style: "none" }, catAxisLabelColor: INK, showTitle: true,
    title: "What decides how good an array is", titleFontSize: 13, titleColor: NAVY }, chartText));
T(s, bullets([
  "How far antennas are from the centre matters most",
  "The golden-ratio angle helps only a little (~2.5%)",
  "Best: spread antennas evenly by area",
]), { x: 6.2, y: 1.5, w: 3.3, h: 3.3, fontSize: 16, paraSpaceAfter: 12 });
s.addNotes("Tested 5 radius rules x 5 angle rules = 25 combinations, everything else fixed.");

// 9 ---------------------------------------------------------------
s = slide("Result 3: Hanle is dry enough only 5% of the time", "Results");
s.addChart(pres.charts.BAR, [{ name: "%", labels: ["Hanle", "Merak", "Ladakh site B", "Ladakh site A", "ALMA (Chile)"], values: [5, 8, 19, 23, 50] }],
  Object.assign({ x: 0.5, y: 1.2, w: 5.4, h: 3.9, barDir: "bar", chartColors: [RED, AMBER, TEAL, TEAL, NAVY], varyColors: true,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0\"%\"", dataLabelFontSize: 13,
    showLegend: false, catAxisLabelFontSize: 13, valAxisHidden: true, valGridLine: { style: "none" },
    catGridLine: { style: "none" }, catAxisLabelColor: INK, showTitle: true,
    title: "% of time the air is dry enough (PWV < 1 mm)", titleFontSize: 13, titleColor: NAVY }, chartText));
T(s, bullets([
  "Guide says Hanle is like ALMA - data says no",
  "But Hanle has roads, power, flat land",
  "Say: good infrastructure, less dry air",
]), { x: 6.2, y: 1.5, w: 3.3, h: 3.3, fontSize: 16, paraSpaceAfter: 12 });
s.addNotes("Source: Raghunath et al. 2026, 184 months of weather data. Hanle is still a good choice for practical reasons.");

// 10 --------------------------------------------------------------
s = slide("A real array design for Ladakh", "Results");
s.addImage({ path: path.join(REPO, "ladakh_array_design/outputs/fig_hanle_extended_terrain.png"),
  x: 0.5, y: 1.2, w: 5.4, h: 5.4 / 2.315 });
T(s, "Left: land height.  Right: buildable land and chosen antenna spots (stars).",
  { x: 0.5, y: 3.6, w: 5.4, h: 0.5, fontSize: 11, color: MUTED, italic: true });
T(s, bullets([
  "16 antennas, 8 m dishes",
  "3 layouts: small, medium, large",
  "See detail 0.096 to 5 arcsec",
  "Exact latitude / longitude for each antenna",
]), { x: 6.2, y: 1.3, w: 3.3, h: 2.8, fontSize: 16, paraSpaceAfter: 10 });
T(s, "First design only - no land, road or power survey yet.", { x: 0.5, y: 4.5, w: 9, h: 0.4, fontSize: 13, italic: true, color: RED });
s.addNotes("Three layouts because one layout cannot give both sharp detail and big-scale view. Files: ladakh_array_design/outputs.");

// 11 --------------------------------------------------------------
s = slide("Result 4: does our score give good pictures?", "Results");
card(s, 0.5, 1.3, 4.3, 3.0, GREEN_BG);
T(s, "SHAPE: YES", { x: 0.8, y: 1.5, w: 3.7, h: 0.45, fontSize: 20, bold: true, color: TEAL });
T(s, "More UV coverage = correct shape in the picture. Checked on 4 test images.",
  { x: 0.8, y: 2.05, w: 3.7, h: 1.3, fontSize: 16 });
card(s, 5.2, 1.3, 4.3, 3.0, RED_BG);
T(s, "BRIGHTNESS: NO", { x: 5.5, y: 1.5, w: 3.7, h: 0.45, fontSize: 20, bold: true, color: RED });
T(s, "Best-score array caught only 42% of a big blurry source. Golden spiral caught 91%.",
  { x: 5.5, y: 2.05, w: 3.7, h: 1.3, fontSize: 16 });
T(s, "Lesson: we need short baselines too, not only long ones.",
  { x: 0.5, y: 4.55, w: 9, h: 0.45, fontSize: 16, bold: true, color: NAVY });
s.addNotes("Pipeline: sky model, visibilities, noise, CLEAN, compare with truth.");

// 12 --------------------------------------------------------------
pres.addSection({ title: "Next" });
s = pres.addSlide({ masterName: "DARK", sectionTitle: "Next" });
T(s, "What I need from you", { x: 0.6, y: 0.4, w: 8.8, h: 0.7, fontSize: 32, bold: true, color: WHITE, fontFace: "Cambria" });
[["Coauthor?", "Do you want to be on the paper?"],
 ["Details for the paper", "Affiliation, email, funding"],
 ["Which story?", "Honest result: quantum is not faster here yet"],
 ["Science target", "What should this telescope look at?"]].forEach(([a, b], i) => {
  const y = 1.35 + i * 0.95;
  s.addShape(pres.shapes.OVAL, { x: 0.6, y, w: 0.5, h: 0.5, fill: { color: AMBER }, line: { color: AMBER } });
  T(s, String(i + 1), { x: 0.6, y, w: 0.5, h: 0.5, align: "center", valign: "middle", fontSize: 16, bold: true, color: WHITE });
  T(s, a, { x: 1.35, y: y - 0.02, w: 8, h: 0.4, fontSize: 20, bold: true, color: WHITE });
  T(s, b, { x: 1.35, y: y + 0.38, w: 8, h: 0.4, fontSize: 15, color: "AFC0D8" });
});
T(s, "Paper: 13 pages ready  |  272 tests pass  |  Code: " + URL, { x: 0.6, y: 5.0, w: 8.8, h: 0.35, fontSize: 12, color: "AFC0D8" });
s.addNotes("The science target is the most useful input: everything else follows from it.");

pres.writeFile({ fileName: OUT }).then((f) => console.log("wrote " + f));
